//! Frozen fixtures: GV invariants that neither pipeline computed, so that they also catch a
//! bug the two share. The inputs are CLI documents in `tests/data/refs.yaml` (computed by
//! cygv 0.2.3) and `tests/data/refs_vex.yaml` (vex phases and FRSTs, computed by cgv, the C
//! program the modular pipeline reimplements), and the expected results are in the matching
//! `*_gv.yaml`, in the output format of the CLI. `tests/data/refs_030.yaml` adds an n-fold,
//! GW invariants and target points, computed by the released cygv 0.3.0. The files say where
//! each geometry comes from.
#![cfg(feature = "cli")]

use cygv::io::Input;
use cygv::{compute, Computed, Method};
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
    "fourfold_gv",
    "fourfold_gw",
    "p11222_gw",
    "p11222_targets",
];

/// GV or GW invariants by curve class and surface (none for a threefold), with the values as
/// strings.
type Table = HashMap<(Vec<i32>, Option<usize>), String>;

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

/// The grading vector, if given, and the invariants of each expected document, by name.
fn expected(text: &str) -> HashMap<String, (Option<Vec<i32>>, Table)> {
    YamlLoader::load_from_str(text)
        .unwrap()
        .iter()
        .map(|doc| {
            let table = doc["results"]
                .as_vec()
                .unwrap()
                .iter()
                .map(|e| {
                    // Values beyond an i64 are read as reals that keep their digits, and GW
                    // invariants are strings.
                    let value = match (&e["gv"], &e["gw"]) {
                        (Yaml::Integer(v), _) => v.to_string(),
                        (Yaml::Real(v), _) | (_, Yaml::String(v)) => v.clone(),
                        v => panic!("{v:?}"),
                    };
                    let surface = e["surface"].as_i64().map(|s| s as usize);
                    ((ints(&e["curve_class"]), surface), value)
                })
                .collect();
            let name = doc["name"].as_str().unwrap().to_owned();
            let grading_vector =
                (!doc["grading_vector"].is_badvalue()).then(|| ints(&doc["grading_vector"]));
            (name, (grading_vector, table))
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
            let mut options = input.options.clone();
            options.method = method;
            let res = compute(&input.problem, input.kind, &options)
                .unwrap_or_else(|e| panic!("{name} ({method:?}): {e}"));
            let (grading, computed): (_, Table) = match res {
                Computed::Gv(gv) => (
                    gv.grading_vector,
                    gv.entries
                        .into_iter()
                        .map(|e| ((e.curve, e.surface), e.value.to_string()))
                        .collect(),
                ),
                Computed::Gw(gw) => (
                    gw.grading_vector,
                    gw.entries
                        .into_iter()
                        .map(|e| ((e.curve, e.surface), e.value.to_string()))
                        .collect(),
                ),
            };
            if let Some(g) = grading_vector {
                assert_eq!(&grading, g, "{name}");
            }
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
fn test_refs_030() {
    check("refs_030");
}

#[test]
fn test_reference_names() {
    let names: Vec<String> = ["refs_gv.yaml", "refs_vex_gv.yaml", "refs_030_gv.yaml"]
        .iter()
        .flat_map(|f| expected(&read(f)).into_keys())
        .collect();
    for name in REFERENCE {
        assert!(names.iter().any(|n| n == name), "{name}");
    }
}
