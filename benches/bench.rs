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
//!
//! The cases enabled by `CYGV_BENCH_HUGE` are never run here, since a single
//! run of each takes minutes; `cargo bench --bench memory` runs them once.

mod common;

use common::{Scenario, Tier};
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use std::time::Duration;

fn benchmark_hkty(c: &mut Criterion) {
    // The scenarios of a group are not all next to each other (the modular
    // pipeline's come after the reference ones of the same tier), so they are
    // collected by name.
    let mut groups: Vec<(String, Vec<Scenario>)> = Vec::new();
    for scenario in common::scenarios()
        .into_iter()
        .filter(|s| s.tier != Tier::Huge)
    {
        match groups
            .iter_mut()
            .find(|(name, _)| *name == scenario.group())
        {
            Some((_, cases)) => cases.push(scenario),
            None => groups.push((scenario.group(), vec![scenario])),
        }
    }

    for (name, cases) in groups {
        let model = common::model(cases[0].model);
        let mut group = c.benchmark_group(&name);
        // A run builds all of its state from scratch, so there is little for a
        // long warm-up to settle; spend the budget on the samples instead.
        group.warm_up_time(Duration::from_secs(1));

        for scenario in cases {
            group.sample_size(scenario.sample_size);
            group.measurement_time(scenario.measurement_time());
            let name = format!("{}{}", scenario.variant.name(), scenario.pipeline.suffix());
            group.bench_function(name, |b| b.iter(|| black_box(scenario.run(&model))));
        }

        group.finish();
    }
}

criterion_group!(benches, benchmark_hkty);
criterion_main!(benches);
