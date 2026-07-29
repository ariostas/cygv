//! Wall-clock benchmarks of the full HKTY pipeline.
//!
//! Run with `cargo bench --bench bench`. A filter can be passed to run a subset,
//! e.g. `cargo bench --bench bench -- threefold/deg10`.
//!
//! Two environment variables tune what is run:
//!
//! - `CYGV_BENCH_HEAVY=1` adds the high-degree cases, which take a few minutes
//!   each to sample but are the ones dominated by bignum arithmetic.
//! - `CYGV_BENCH_THREADS=n` pins the number of worker threads, which makes
//!   results less noisy and comparable across machines.

mod common;

use common::Scenario;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use std::time::Duration;

fn benchmark_hkty(c: &mut Criterion) {
    // Scenarios of the same group are emitted contiguously, so a running fold
    // is enough to collect them.
    let mut groups: Vec<(String, Vec<Scenario>)> = Vec::new();
    for scenario in common::scenarios() {
        match groups.last_mut() {
            Some((name, cases)) if *name == scenario.group() => cases.push(scenario),
            _ => groups.push((scenario.group(), vec![scenario])),
        }
    }

    for (name, cases) in groups {
        let model = common::model(cases[0].model);
        let mut group = c.benchmark_group(&name);
        // A run allocates its own number pools, so there is little state for a
        // long warm-up to settle; spend the budget on the samples instead.
        group.warm_up_time(Duration::from_secs(1));

        for scenario in cases {
            group.sample_size(scenario.sample_size);
            group.measurement_time(scenario.measurement_time());
            group.bench_function(scenario.variant.name(), |b| {
                b.iter(|| black_box(scenario.run(&model)))
            });
        }

        group.finish();
    }
}

criterion_group!(benches, benchmark_hkty);
criterion_main!(benches);
