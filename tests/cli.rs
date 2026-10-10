//! Tests of the command line interface.
#![cfg(feature = "cli")]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use yaml_rust2::{Yaml, YamlLoader};

/// Path to an example input.
fn example(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(name)
}

/// Run the binary with the given arguments and standard input.
fn run(args: &[&str], stdin: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_cygv"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to run the cygv binary");
    child
        .stdin
        .as_mut()
        .expect("failed to open stdin")
        .write_all(stdin.as_bytes())
        .expect("failed to write to stdin");
    child.wait_with_output().expect("failed to read the output")
}

/// Check that the given documents contain results for the expected invariants.
fn check_results(docs: &[Yaml], expected: &[&str]) {
    assert_eq!(docs.len(), expected.len());
    for (doc, key) in docs.iter().zip(expected.iter()) {
        assert_eq!(doc["invariants"].as_str(), Some(*key));
        let results = doc["results"].as_vec().expect("missing results");
        assert!(!results.is_empty());
        assert!(!doc["grading_vector"]
            .as_vec()
            .expect("missing grading vector")
            .is_empty());
        for entry in results {
            assert!(!entry["curve_class"]
                .as_vec()
                .expect("missing class")
                .is_empty());
            assert!(entry["degree"].as_i64().is_some_and(|d| d > 0));
            assert!(!entry[*key].is_badvalue());
        }
    }
}

#[test]
fn test_stdin_to_stdout() {
    let input = std::fs::read_to_string(example("threefold.yaml")).unwrap();
    let output = run(&["--threads", "2"], &input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let docs = YamlLoader::load_from_str(&String::from_utf8(output.stdout).unwrap()).unwrap();
    check_results(&docs, &["gv", "gw"]);

    // The GW invariants are written as fractions, and every GV invariant of this CY
    // that is not filtered out is a nonzero integer.
    let gv = &docs[0]["results"][0]["gv"];
    assert!(gv.as_i64().is_some_and(|c| c != 0), "{gv:?}");
    assert!(docs[1]["results"][0]["gw"].as_str().is_some());
}

#[test]
fn test_file_to_file() {
    let out_path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("fourfold.yaml");
    let output = run(
        &[
            "--file",
            example("fourfold.yaml").to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
        ],
        "",
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());

    let docs = YamlLoader::load_from_str(&std::fs::read_to_string(&out_path).unwrap()).unwrap();
    check_results(&docs, &["gv"]);
    assert_eq!(docs[0]["name"].as_str(), Some("fourfold"));
    assert_eq!(docs[0]["is_threefold"].as_bool(), Some(false));
    // The invariants of CYs that are not threefolds carry a reference surface index.
    let first = &docs[0]["results"][0];
    assert_eq!(first["surface_index"].as_i64(), Some(4));
    // These are the invariants of the complete intersection; read as a hypersurface, the
    // same data gives 79 instead, starting with a class that has none here.
    let class: Vec<i64> = first["curve_class"]
        .as_vec()
        .unwrap()
        .iter()
        .map(|c| c.as_i64().unwrap())
        .collect();
    assert_eq!(class, [0, 0, 0, 0, -2, 1]);
    assert_eq!(first["gv"].as_i64(), Some(-64));
    assert_eq!(docs[0]["results"].as_vec().unwrap().len(), 56);
}

#[test]
fn test_invalid_input() {
    let output = run(&[], "generators: [[1, 0], [0, 1]]\n");
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("grading_vector"), "{stderr}");
    assert!(output.stdout.is_empty());
}

#[test]
fn test_missing_file() {
    let output = run(&["--file", "does_not_exist.yaml"], "");
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("does_not_exist.yaml"), "{stderr}");
}

#[test]
fn test_methods() {
    // Both methods give the same invariants of the quintic, written the same way.
    let quintic = "\
generators: [[1]]
grading_vector: [1]
q: [[1, 1, 1, 1, 1]]
intnums: [[0, 0, 0, 5]]
max_deg: 6
";
    let input = format!("---\n{quintic}---\n{quintic}invariants: gw\n");
    let outputs: Vec<String> = ["reference", "modular", "auto"]
        .iter()
        .map(|method| {
            let output = run(&["--method", method], &input);
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8(output.stdout).unwrap()
        })
        .collect();
    assert!(outputs[0].contains("{curve_class: [1], degree: 1, gv: 2875}"));
    assert!(outputs[0].contains("{curve_class: [2], degree: 2, gw: '4876875/8'}"));
    assert_eq!(outputs[1], outputs[0]);
    assert_eq!(outputs[2], outputs[0]);

    // The modular method needs max_deg, whether it is asked for in the input or with
    // the flag, which overrides the input.
    let input = quintic.replace("max_deg: 6", "min_points: 4");
    for (args, data) in [
        (&[][..], format!("{input}method: modular\n")),
        (
            &["--method", "modular"][..],
            format!("{input}method: auto\n"),
        ),
    ] {
        let output = run(args, &data);
        assert!(!output.status.success());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("min_points"), "{stderr}");
    }
    let output = run(&["--method", "fast"], quintic);
    assert!(!output.status.success());
}

#[test]
fn test_cone_data() {
    // The example with cone data and no generators gives the same invariants as the
    // reference method with the generators.
    let input = std::fs::read_to_string(example("cones.yaml")).unwrap();
    let with_cones = run(&[], &input);
    assert!(
        with_cones.status.success(),
        "{}",
        String::from_utf8_lossy(&with_cones.stderr)
    );
    let head = &input[..input.find("# The inequalities").unwrap()];
    let reference = run(
        &["--method", "reference"],
        &format!("{head}generators: [[1, 0], [0, 1]]\n"),
    );
    assert!(reference.status.success());
    assert_eq!(with_cones.stdout, reference.stdout);
    let docs = YamlLoader::load_from_str(&String::from_utf8(with_cones.stdout).unwrap()).unwrap();
    check_results(&docs, &["gv"]);

    // Without generators the reference method cannot run.
    let output = run(&["--method", "reference"], &input);
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("generators"), "{stderr}");
}

#[test]
fn test_vex_cone_data() {
    // The vex example, with the cone of its stratum, gives the same invariants as the
    // reference method with the generators.
    let input = std::fs::read_to_string(example("vex.yaml")).unwrap();
    let with_cones = run(&[], &input);
    assert!(
        with_cones.status.success(),
        "{}",
        String::from_utf8_lossy(&with_cones.stderr)
    );
    let head = &input[..input.find("# The inequalities").unwrap()];
    let reference = run(
        &["--method", "reference"],
        &format!("{head}generators: [[0, 1], [1, -1]]\n"),
    );
    assert!(reference.status.success());
    assert_eq!(with_cones.stdout, reference.stdout);
    let out = String::from_utf8(with_cones.stdout).unwrap();
    assert!(
        out.contains("{curve_class: [1, -1], degree: 1, gv: 1}"),
        "{out}"
    );

    // Without the cone of the stratum, the cone data is incomplete.
    let without = input.replace(
        "  - {divisors: [1, 2, 4], hilbert_basis: [[0, 1], [1, -1]]}\n",
        "",
    );
    assert_ne!(without, input);
    let output = run(&[], &without);
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("stratum"), "{stderr}");
}
