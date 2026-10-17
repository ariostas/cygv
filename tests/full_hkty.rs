use cygv::malachite::{Integer, Rational};
use cygv::{
    compute_gv, compute_gw, Error, Geometry, Invariants, Method, MoriCone, Options, Problem,
    Selection,
};
use nalgebra::{dmatrix, DMatrix};
use std::collections::HashMap;

/// A problem on a hypersurface or complete intersection with the given generators.
fn problem(
    q: &DMatrix<i32>,
    nefpart: &[Vec<usize>],
    intnums: &HashMap<(usize, usize, usize), i32>,
    generators: &DMatrix<i32>,
    grading_vector: &[i32],
    selection: Selection,
) -> Problem {
    Problem {
        geometry: Geometry {
            q: q.clone(),
            intnums: intnums.clone(),
            nefpart: nefpart.to_vec(),
            ..Geometry::default()
        },
        cone: MoriCone::Generators(generators.clone()),
        selection,
        grading_vector: Some(grading_vector.to_vec()),
    }
}

fn method(method: Method) -> Options {
    Options {
        method,
        ..Options::default()
    }
}

/// The curve classes, surfaces and values of the invariants.
fn values<T: Clone>(inv: &Invariants<T>) -> Vec<(Vec<i32>, Option<usize>, T)> {
    inv.entries
        .iter()
        .map(|e| (e.curve.clone(), e.surface, e.value.clone()))
        .collect()
}

#[test]
fn test_threefold() {
    let generators = dmatrix![
         0, 1;
        -1, 2;
    ];
    let q = dmatrix![
        1,  0;
        1,  0;
        1, -1;
        0,  1;
        1,  1;
        2, -1;
    ];
    let intnums = HashMap::from([
        ((0, 0, 0), 2),
        ((0, 0, 1), 1),
        ((0, 1, 1), -1),
        ((1, 1, 1), -5),
    ]);
    let problem = problem(
        &q,
        &[],
        &intnums,
        &generators,
        &[3, -1],
        Selection::MinPoints(100),
    );

    // For now, these are just smoke tests
    let gv = compute_gv(&problem, &Options::default()).unwrap();
    assert_eq!(gv.dimension, 3);
    assert_eq!(gv.report.method, Method::Reference);
    assert!(gv.entries.iter().all(|e| e.surface.is_none()));
    compute_gw(&problem, &Options::default()).unwrap();
}

#[test]
fn test_fourfold() {
    let generators = dmatrix![
        0,1,0,0,0,0;
        0,0,0,1,-2,0;
        0,-4,1,0,0,0;
        -2,0,1,0,1,0;
        3,0,0,1,0,-2;
        0,0,0,-2,2,1;
    ];
    let grading_vector = vec![73, -30, 18, -17, -11, -21];

    let q = dmatrix![
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        8,-4,2,-2,-2,-3;
        12,-6,3,-3,-2,-4;
        4,-3,0,-2,-1,-2;
        0,0,0,1,0,0;
        0,0,0,0,1,0;
        0,0,0,0,0,1;
    ];
    let nefpart = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7, 8, 9, 10, 11]];

    let intnums = HashMap::from([
        ((0, 0, 2), 1),
        ((0, 1, 1), -8),
        ((0, 1, 3), 4),
        ((0, 1, 5), 8),
        ((0, 2, 2), -4),
        ((0, 3, 3), -8),
        ((0, 4, 4), -16),
        ((0, 4, 5), 8),
        ((0, 5, 5), -16),
        ((1, 0, 1), -8),
        ((1, 0, 3), 4),
        ((1, 0, 5), 8),
        ((1, 1, 3), 16),
        ((1, 1, 5), -32),
        ((1, 3, 3), -16),
        ((1, 5, 5), 64),
        ((2, 0, 0), 1),
        ((2, 0, 2), -4),
        ((2, 2, 2), 16),
        ((3, 0, 1), 4),
        ((3, 0, 3), -8),
        ((3, 1, 1), 16),
        ((3, 1, 3), -16),
        ((4, 0, 4), -16),
        ((4, 0, 5), 8),
        ((4, 4, 4), -128),
        ((4, 4, 5), 32),
        ((5, 0, 1), 8),
        ((5, 0, 4), 8),
        ((5, 0, 5), -16),
        ((5, 1, 1), -32),
        ((5, 1, 5), 64),
        ((5, 4, 4), 32),
        ((5, 5, 5), -128),
        ((6, 0, 0), -8),
        ((6, 0, 3), 16),
        ((6, 0, 5), -32),
        ((6, 1, 1), -128),
        ((6, 1, 3), 64),
        ((6, 1, 5), 128),
        ((6, 3, 3), -64),
        ((6, 5, 5), -256),
        ((7, 0, 0), 4),
        ((7, 0, 1), 16),
        ((7, 0, 3), -16),
        ((7, 1, 1), 64),
        ((7, 1, 3), -64),
        ((7, 3, 3), 64),
        ((8, 0, 0), 8),
        ((8, 0, 1), -32),
        ((8, 0, 5), 64),
        ((8, 1, 1), 128),
        ((8, 1, 5), -256),
        ((8, 5, 5), 512),
        ((9, 0, 0), -4),
        ((9, 0, 2), 16),
        ((9, 2, 2), -64),
        ((10, 0, 0), -8),
        ((10, 0, 1), -16),
        ((10, 1, 1), -64),
        ((10, 1, 3), 64),
        ((10, 3, 3), -128),
        ((11, 0, 0), -16),
        ((11, 0, 4), -128),
        ((11, 0, 5), 32),
        ((11, 4, 4), -768),
        ((11, 4, 5), 128),
        ((12, 0, 0), 8),
        ((12, 0, 4), 32),
        ((12, 4, 4), 128),
        ((13, 0, 0), -16),
        ((13, 0, 1), 64),
        ((13, 0, 5), -128),
        ((13, 1, 1), -256),
        ((13, 1, 5), 512),
        ((13, 5, 5), -1024),
    ]);

    // For now, these are just smoke tests
    let min_points = problem(
        &q,
        &nefpart,
        &intnums,
        &generators,
        &grading_vector,
        Selection::MinPoints(100),
    );
    let gv = compute_gv(&min_points, &Options::default()).unwrap();
    assert_eq!(gv.dimension, 4);
    assert!(gv.entries.iter().all(|e| e.surface.is_some()));
    compute_gw(&min_points, &Options::default()).unwrap();

    // Both methods give the same invariants of this complete intersection.
    let max_deg = problem(
        &q,
        &nefpart,
        &intnums,
        &generators,
        &grading_vector,
        Selection::MaxDeg(8),
    );
    let gv = |m| values(&compute_gv(&max_deg, &method(m)).unwrap());
    let gw = |m| values(&compute_gw(&max_deg, &method(m)).unwrap());
    let expected = gv(Method::Reference);
    assert!(!expected.is_empty());
    assert_eq!(gv(Method::Modular), expected);
    assert_eq!(gv(Method::Auto), expected);
    let expected = gw(Method::Reference);
    assert!(expected.len() > gv(Method::Reference).len());
    assert_eq!(gw(Method::Modular), expected);
}

/// The quintic, through the entry points with each method.
#[test]
fn test_methods() {
    let q = DMatrix::from_element(5, 1, 1);
    let intnums = HashMap::from([((0, 0, 0), 5)]);
    let quintic = |selection| problem(&q, &[], &intnums, &dmatrix![1], &[1], selection);
    let gv = |m, selection| compute_gv(&quintic(selection), &method(m)).map(|gv| values(&gv));
    let gw = |m| values(&compute_gw(&quintic(Selection::MaxDeg(6)), &method(m)).unwrap());

    let expected = gv(Method::Reference, Selection::MaxDeg(6)).unwrap();
    assert_eq!(expected.len(), 6);
    assert_eq!(expected[0], (vec![1], None, Integer::from(2875)));
    assert_eq!(gv(Method::Modular, Selection::MaxDeg(6)).unwrap(), expected);
    assert_eq!(gv(Method::Auto, Selection::MaxDeg(6)).unwrap(), expected);
    let expected = gw(Method::Reference);
    // GV_2 + GV_1 / 8
    assert_eq!(expected[1].2, Rational::from_signeds(609250 * 8 + 2875, 8));
    assert_eq!(gw(Method::Modular), expected);

    // The report says which method ran.
    let auto = compute_gv(&quintic(Selection::MaxDeg(3)), &Options::default()).unwrap();
    assert_eq!(auto.report.method, Method::Modular);
    assert!(!auto.report.certified);
    let certified = Options {
        certify: true,
        ..Options::default()
    };
    assert!(
        compute_gv(&quintic(Selection::MaxDeg(3)), &certified)
            .unwrap()
            .report
            .certified
    );

    // The modular method only supports max_deg and target_points; auto runs the reference
    // one for the others.
    assert!(matches!(
        gv(Method::Modular, Selection::MinPoints(4)),
        Err(Error::Unsupported(r)) if r.contains("min_points")
    ));
    let expected = gv(Method::Reference, Selection::MinPoints(4)).unwrap();
    assert!(!expected.is_empty());
    assert_eq!(gv(Method::Auto, Selection::MinPoints(4)).unwrap(), expected);
    let reference = compute_gv(&quintic(Selection::MinPoints(4)), &Options::default()).unwrap();
    assert_eq!(reference.report.method, Method::Reference);
    assert!(reference.report.certified);
}

/// Invalid nef partitions, intersection numbers and curve classes are errors with each
/// method, rather than panics of the reference, and the curve classes get zero.
#[test]
fn test_invalid_input() {
    let q = DMatrix::from_element(5, 1, 1);
    let quintic = HashMap::from([((0, 0, 0), 5)]);
    let gv = |selection, nefpart: &[Vec<usize>], intnums: &HashMap<_, _>, m| {
        let problem = problem(&q, nefpart, intnums, &dmatrix![1], &[1], selection);
        compute_gv(&problem, &method(m)).map(|gv| values(&gv))
    };
    let field = |res: Result<_, Error>| match res {
        Err(Error::Input { field, reason }) => (field, reason),
        res => panic!("{res:?}"),
    };
    for m in [Method::Auto, Method::Modular, Method::Reference] {
        let run = |nefpart: &[Vec<usize>], intnums: &HashMap<_, _>| {
            gv(Selection::MaxDeg(3), nefpart, intnums, m)
        };
        for (nefpart, reason) in [
            (vec![vec![0, 1, 2, 3]], "include every row"),
            (vec![vec![0, 1, 2, 3, 4, 4]], "disjoint"),
            (vec![vec![0, 1, 2, 3, 7]], "indices of rows"),
            (vec![vec![0, 1, 2, 3, 4], vec![]], "not be empty"),
        ] {
            let (f, r) = field(run(&nefpart, &quintic));
            assert_eq!(f, "nefpart");
            assert!(r.contains(reason), "{m:?}: {r}");
        }
        let (f, r) = field(run(&[vec![0, 1], vec![2, 3, 4]], &quintic));
        assert_eq!(f, "q");
        assert!(r.contains("but it is 2"), "{r}");
        // A partition into one part is the hypersurface.
        assert_eq!(run(&[vec![4, 0, 1, 2, 3]], &quintic), run(&[], &quintic));
        for intnums in [
            HashMap::from([((0, 0, 0), 5), ((0, 0, 1), 1)]),
            HashMap::from([((0, 0, 0), 0)]),
        ] {
            assert_eq!(field(run(&[], &intnums)).0, "intnums", "{m:?}: {intnums:?}");
        }
    }

    // The curve classes are used as given, with zero added, each once.
    let expected = gv(Selection::MaxDeg(2), &[], &quintic, Method::Reference).unwrap();
    assert_eq!(expected.len(), 2);
    for classes in [dmatrix![1, 2], dmatrix![2, 0, 1], dmatrix![1, 2, 1]] {
        let res = gv(
            Selection::CurveClasses(classes),
            &[],
            &quintic,
            Method::Auto,
        );
        assert_eq!(res.unwrap(), expected);
    }
    assert_eq!(
        field(gv(
            Selection::CurveClasses(dmatrix![1, -1]),
            &[],
            &quintic,
            Method::Auto
        ))
        .0,
        "curve_classes"
    );
}
