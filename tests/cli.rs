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
    assert!(docs[0]["results"][0]["surface_index"].as_i64().is_some());
}

#[test]
fn test_float_arithmetic() {
    let input = std::fs::read_to_string(example("threefold.yaml")).unwrap();
    let exact = run(&[], &input);
    let approx = run(&[], &format!("{input}prec: 300\n"));
    assert!(
        approx.status.success(),
        "{}",
        String::from_utf8_lossy(&approx.stderr)
    );

    let exact = YamlLoader::load_from_str(&String::from_utf8(exact.stdout).unwrap()).unwrap();
    let approx = YamlLoader::load_from_str(&String::from_utf8(approx.stdout).unwrap()).unwrap();
    // Only the last document carries the added precision field.
    assert_eq!(exact[0]["results"], approx[0]["results"]);
    assert_eq!(
        exact[1]["results"].as_vec().map(Vec::len),
        approx[1]["results"].as_vec().map(Vec::len)
    );
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
