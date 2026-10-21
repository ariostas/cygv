//! Frozen fixtures: GV invariants that neither pipeline computed (layer 4 of the testing
//! strategy in `PLAN.md` 5.2), so that they also catch a bug the two share. The inputs are CLI
//! documents in `tests/data/refs.yaml` (computed by cygv 0.2.3) and `tests/data/refs_vex.yaml`
//! (vex phases and FRSTs, computed by the C program of PR #94), and the expected results are
//! in the matching `*_gv.yaml`, in the output format of the CLI. The files say where each
//! geometry comes from.
#![cfg(feature = "cli")]

use cygv::io::Input;
use cygv::{compute_gv, Method, Options};
use std::collections::HashMap;
use std::path::Path;
use yaml_rust2::{Yaml, YamlLoader};

/// The fixtures that the reference method also runs, which take it well under a second each
/// in a debug build.
const REFERENCE: &[&str] = &[
    "quintic_D10",
    "h10_72_3998_0_default_K5",
    "h10_92_226_0_default_K5",
    "h10_120_48_0_plike_K7",
    "liam_h2_fan0",
    "liam_h2_fan1",
    "liam_h3_fan0",
    "liam_h3_fan1",
    "liam_h3_fan2",
    "liam_h3_fan3",
    "liam_h3_fan4",
    "liam_h3_fan5",
    "liam_h3_fan6",
    "liam_h3_fan7",
    "liam_h3_fan8",
    "liam_h3_fan9",
];

/// GV invariants by curve class, with the values as strings.
type Table = HashMap<Vec<i32>, String>;

fn read(file: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(file);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn ints(y: &Yaml) -> Vec<i32> {
    y.as_vec()
        .unwrap()
        .iter()
        .map(|x| x.as_i64().unwrap() as i32)
        .collect()
}

/// The grading vector and the invariants of each expected document, by name.
fn expected(text: &str) -> HashMap<String, (Vec<i32>, Table)> {
    YamlLoader::load_from_str(text)
        .unwrap()
        .iter()
        .map(|doc| {
            let table = doc["results"]
                .as_vec()
                .unwrap()
                .iter()
                .map(|e| {
                    // Values beyond an i64 are read as reals that keep their digits.
                    let value = match &e["gv"] {
                        Yaml::Integer(v) => v.to_string(),
                        Yaml::Real(v) => v.clone(),
                        v => panic!("{v:?}"),
                    };
                    (ints(&e["curve_class"]), value)
                })
                .collect();
            let name = doc["name"].as_str().unwrap().to_owned();
            (name, (ints(&doc["grading_vector"]), table))
        })
        .collect()
}

fn check(file: &str) {
    let inputs = Input::load_all(&read(&format!("{file}.yaml"))).unwrap();
    let expected = expected(&read(&format!("{file}_gv.yaml")));
    assert_eq!(inputs.len(), expected.len());
    for input in &inputs {
        let name = input.name.as_deref().unwrap();
        let (grading_vector, table) = &expected[name];
        let mut methods = vec![Method::Auto];
        if REFERENCE.contains(&name) {
            methods.push(Method::Reference);
        }
        for method in methods {
            let options = Options {
                method,
                ..input.options.clone()
            };
            let gv = compute_gv(&input.problem, &options)
                .unwrap_or_else(|e| panic!("{name} ({method:?}): {e}"));
            assert_eq!(&gv.grading_vector, grading_vector, "{name}");
            let computed: Table = gv
                .entries
                .into_iter()
                .map(|e| (e.curve, e.value.to_string()))
                .collect();
            let differ = table.iter().filter(|(c, v)| computed.get(*c) != Some(v));
            assert!(
                computed.len() == table.len() && differ.clone().next().is_none(),
                "{name} ({method:?}): {} computed, {} expected, {} of those differ, e.g. {:?}",
                computed.len(),
                table.len(),
                differ.clone().count(),
                differ.take(3).collect::<Vec<_>>()
            );
        }
    }
}

#[test]
fn test_refs() {
    check("refs");
}

#[test]
fn test_refs_vex() {
    check("refs_vex");
}

#[test]
fn test_reference_names() {
    let names: Vec<String> = ["refs_gv.yaml", "refs_vex_gv.yaml"]
        .iter()
        .flat_map(|f| expected(&read(f)).into_keys())
        .collect();
    for name in REFERENCE {
        assert!(names.iter().any(|n| n == name), "{name}");
    }
}
