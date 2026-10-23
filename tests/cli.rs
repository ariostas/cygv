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
    // A binary that rejects its arguments exits without reading its input, which may close
    // the pipe before the input is written.
    match child
        .stdin
        .as_mut()
        .expect("failed to open stdin")
        .write_all(stdin.as_bytes())
    {
        Err(e) if e.kind() != std::io::ErrorKind::BrokenPipe => {
            panic!("failed to write to stdin: {e}")
        }
        _ => {}
    }
    child.wait_with_output().expect("failed to read the output")
}

/// The output without its report lines, which say how the results were computed.
fn without_report(output: &[u8]) -> String {
    String::from_utf8_lossy(output)
        .lines()
        .filter(|l| !l.starts_with("report:"))
        .map(|l| format!("{l}\n"))
        .collect()
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
    assert_eq!(docs[0]["dimension"].as_i64(), Some(4));
    assert_eq!(docs[0]["report"]["method"].as_str(), Some("reference"));
    // The invariants of CYs that are not threefolds carry a reference surface index.
    let first = &docs[0]["results"][0];
    assert_eq!(first["surface"].as_i64(), Some(4));
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
    assert!(stderr.contains("\"q\""), "{stderr}");
    assert!(output.stdout.is_empty());

    // Every document is checked before any is computed.
    let quintic = "q: [[1, 1, 1, 1, 1]]\nintnums: [[0, 0, 0, 5]]\ngenerators: [[1]]\n";
    let output = run(&[], &format!("---\n{quintic}max_deg: 3\n---\n{quintic}"));
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("max_deg"), "{stderr}");
    assert!(output.stdout.is_empty());
}

#[test]
fn test_verbose() {
    // The progress goes to stderr, and the results are those of a quiet run.
    let quintic = "q: [[1, 1, 1, 1, 1]]\nintnums: [[0, 0, 0, 5]]\ngenerators: [[1]]\nmax_deg: 4\n";
    let quiet = run(&[], quintic);
    assert!(quiet.status.success());
    assert!(quiet.stderr.is_empty());
    let output = run(&["-v"], quintic);
    assert!(output.status.success());
    assert_eq!(output.stdout, quiet.stdout);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("pass over primes"), "{stderr}");
    assert!(!stderr.contains("degree 2:"), "{stderr}");
    let output = run(&["-vv", "--method", "modular"], quintic);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("degree 2:"), "{stderr}");
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
            let report = format!("report: {{method: {}", method.replace("auto", "modular"));
            assert!(String::from_utf8_lossy(&output.stdout).contains(&report));
            without_report(&output.stdout)
        })
        .collect();
    assert!(outputs[0].contains("{curve_class: [1], degree: 1, gv: 2875}"));
    assert!(outputs[0].contains("{curve_class: [2], degree: 2, gw: '4876875/8'}"));
    assert_eq!(outputs[1], outputs[0]);
    assert_eq!(outputs[2], outputs[0]);

    // Certified, as the flag or in the input, the results are the same.
    let certified = run(&["--certify"], &input);
    assert!(certified.status.success());
    let out = String::from_utf8_lossy(&certified.stdout);
    assert!(out.contains("certified: true"), "{out}");
    assert_eq!(without_report(&certified.stdout), outputs[0]);
    let in_input = input.replace("max_deg: 6\n", "max_deg: 6\ncertify: true\n");
    let certified = run(&[], &in_input);
    assert!(certified.status.success());
    assert_eq!(without_report(&certified.stdout), outputs[0]);
    let uncertified = run(&["--no-certify"], &in_input);
    let out = String::from_utf8_lossy(&uncertified.stdout);
    assert!(out.contains("certified: false"), "{out}");

    // So are those of the low-memory mode, as the flag or in the input.
    let low = run(&["--low-memory"], &input);
    assert!(low.status.success());
    assert_eq!(without_report(&low.stdout), outputs[0]);
    let in_input = input.replace("max_deg: 6\n", "max_deg: 6\nlow_memory: true\n");
    let low = run(&[], &in_input);
    assert!(low.status.success());
    assert_eq!(without_report(&low.stdout), outputs[0]);

    // The maximum degree and the invariants can be set for every document.
    let output = run(&["--max-deg", "2", "--invariants", "gw"], &input);
    assert!(output.status.success());
    let out = String::from_utf8(output.stdout).unwrap();
    assert!(!out.contains("gv:") && !out.contains("degree: 3"), "{out}");
    assert!(
        out.contains("{curve_class: [2], degree: 2, gw: '4876875/8'}"),
        "{out}"
    );
    let output = run(
        &["--max-deg", "2"],
        &quintic.replace("max_deg: 6", "min_points: 4"),
    );
    assert!(!output.status.success());

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
    assert_eq!(
        without_report(&with_cones.stdout),
        without_report(&reference.stdout)
    );
    let docs = YamlLoader::load_from_str(&String::from_utf8(with_cones.stdout).unwrap()).unwrap();
    check_results(&docs, &["gv"]);

    // The reference method takes the Mori cone from its inequalities, and has no use for the
    // Hilbert bases.
    let output = run(&["--method", "reference"], &input);
    assert!(output.status.success());
    assert_eq!(output.stdout, reference.stdout);
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
    assert_eq!(
        without_report(&with_cones.stdout),
        without_report(&reference.stdout)
    );
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

#[test]
fn test_fan() {
    // The vex example given by its fan gives the same output as with its strata and cone
    // data, since the grading vector derived from the fan is the same.
    let with_fan = run(&["--file", example("fan.yaml").to_str().unwrap()], "");
    assert!(
        with_fan.status.success(),
        "{}",
        String::from_utf8_lossy(&with_fan.stderr)
    );
    let with_cones = run(&["--file", example("vex.yaml").to_str().unwrap()], "");
    assert_eq!(
        without_report(&with_fan.stdout),
        without_report(&with_cones.stdout)
    );
    let out = String::from_utf8_lossy(&with_fan.stdout);
    assert!(out.contains("grading_vector_derived: true"), "{out}");

    // A vex phase cannot be certified.
    let certified = run(
        &["--certify", "--file", example("fan.yaml").to_str().unwrap()],
        "",
    );
    assert!(!certified.status.success());
    let stderr = String::from_utf8(certified.stderr).unwrap();
    assert!(stderr.contains("vex phase"), "{stderr}");
}
