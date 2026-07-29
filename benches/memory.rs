//! Memory benchmarks of the full HKTY pipeline.
//!
//! Run with `cargo bench --bench memory`. A filter can be passed to run a
//! subset, e.g. `cargo bench --bench memory -- threefold/deg10`, and the same
//! `CYGV_BENCH_HEAVY` and `CYGV_BENCH_THREADS` variables as the wall-clock
//! benchmarks apply. Each scenario is run exactly once, so enabling the heavy
//! cases here is much cheaper than it is there.
//!
//! Almost all of the memory this crate uses is held by GMP/MPFR numbers, which
//! are allocated by C code and therefore never reach Rust's global allocator.
//! To account for them, this benchmark installs its own allocation functions
//! with `mp_set_memory_functions` (MPFR routes its allocations through the same
//! hooks) alongside a global allocator wrapper, so every allocation on either
//! side is counted. GMP hands the size back on free and realloc, so the tracking
//! needs no bookkeeping of its own beyond a few counters.
//!
//! Each scenario runs in its own child process. That keeps one scenario's peak
//! from leaking into the next one, and it lets the child report a peak RSS that
//! covers everything a counter cannot see (thread stacks, allocator overhead,
//! fragmentation).

mod common;

use std::alloc::{GlobalAlloc, Layout, System};
use std::ffi::c_void;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Bytes currently held by live allocations.
static CURRENT: AtomicUsize = AtomicUsize::new(0);
/// High-water mark of `CURRENT`.
static PEAK: AtomicUsize = AtomicUsize::new(0);
/// Bytes handed out over the whole run, counting only the growth of a realloc.
static TOTAL: AtomicUsize = AtomicUsize::new(0);
/// Number of allocation and reallocation calls.
static COUNT: AtomicUsize = AtomicUsize::new(0);

fn on_alloc(size: usize) {
    let current = CURRENT.fetch_add(size, Ordering::Relaxed) + size;
    PEAK.fetch_max(current, Ordering::Relaxed);
    TOTAL.fetch_add(size, Ordering::Relaxed);
    COUNT.fetch_add(1, Ordering::Relaxed);
}

fn on_free(size: usize) {
    CURRENT.fetch_sub(size, Ordering::Relaxed);
}

fn on_realloc(old_size: usize, new_size: usize) {
    if new_size >= old_size {
        on_alloc(new_size - old_size);
    } else {
        COUNT.fetch_add(1, Ordering::Relaxed);
        on_free(old_size - new_size);
    }
}

// ---------------------------------------------------------------------------
// Rust side
// ---------------------------------------------------------------------------

struct TrackingAllocator;

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            on_alloc(layout.size());
        }
        ptr
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if !ptr.is_null() {
            on_alloc(layout.size());
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        on_free(layout.size());
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new_ptr = unsafe { System.realloc(ptr, layout, new_size) };
        if !new_ptr.is_null() {
            on_realloc(layout.size(), new_size);
        }
        new_ptr
    }
}

#[global_allocator]
static ALLOCATOR: TrackingAllocator = TrackingAllocator;

// ---------------------------------------------------------------------------
// GMP/MPFR side
// ---------------------------------------------------------------------------

/// Alignment of the blocks handed to GMP. It only requires `mp_limb_t`
/// alignment, but matching `max_align_t` keeps this safe for the odd struct
/// MPFR allocates through the same hooks.
const GMP_ALIGN: usize = 16;

/// GMP never asks for zero bytes, but the layout has to be valid regardless,
/// and the same adjustment has to be made when the block comes back.
fn gmp_layout(size: usize) -> Layout {
    Layout::from_size_align(size.max(1), GMP_ALIGN).expect("invalid GMP allocation size")
}

// GMP declares its allocation functions as safe `extern "C"` pointers, so these
// cannot be marked `unsafe` even though they are only ever called by GMP, with
// the pointers and sizes it got from here.
extern "C" fn gmp_alloc(size: usize) -> *mut c_void {
    let ptr = unsafe { System.alloc(gmp_layout(size)) };
    if ptr.is_null() {
        std::alloc::handle_alloc_error(gmp_layout(size));
    }
    on_alloc(size);
    ptr.cast()
}

extern "C" fn gmp_realloc(ptr: *mut c_void, old_size: usize, new_size: usize) -> *mut c_void {
    let new_ptr = unsafe { System.realloc(ptr.cast(), gmp_layout(old_size), new_size.max(1)) };
    if new_ptr.is_null() {
        std::alloc::handle_alloc_error(gmp_layout(new_size));
    }
    on_realloc(old_size, new_size);
    new_ptr.cast()
}

extern "C" fn gmp_free(ptr: *mut c_void, size: usize) {
    on_free(size);
    unsafe { System.dealloc(ptr.cast(), gmp_layout(size)) }
}

/// Routes GMP and MPFR allocations through the counters.
///
/// # Safety
///
/// Must be called before any GMP or MPFR object exists, since blocks allocated
/// by the previous functions would be freed by the new ones.
unsafe fn install_gmp_hooks() {
    unsafe {
        gmp_mpfr_sys::gmp::set_memory_functions(Some(gmp_alloc), Some(gmp_realloc), Some(gmp_free))
    }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

/// Peak resident set size of this process, in bytes, if the platform reports it.
fn peak_rss() -> Option<usize> {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        let line = status.lines().find(|l| l.starts_with("VmHWM:"))?;
        let kib: usize = line.split_whitespace().nth(1)?.parse().ok()?;
        Some(kib * 1024)
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

fn human_bytes(bytes: usize) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{value:.1} {}", UNITS[unit])
    }
}

struct Measurement {
    peak: usize,
    total: usize,
    allocs: usize,
    rss: Option<usize>,
    secs: f64,
    results: usize,
}

impl Measurement {
    /// Serialized on a single line so the parent process can read it back.
    fn encode(&self) -> String {
        format!(
            "RESULT peak={} total={} allocs={} rss={} secs={} results={}",
            self.peak,
            self.total,
            self.allocs,
            self.rss.map_or(-1, |r| r as i64),
            self.secs,
            self.results,
        )
    }

    fn decode(line: &str) -> Option<Measurement> {
        let mut fields = std::collections::HashMap::new();
        for field in line.strip_prefix("RESULT ")?.split_whitespace() {
            let (key, value) = field.split_once('=')?;
            fields.insert(key, value);
        }
        let get = |key: &str| fields.get(key).copied();
        let rss = get("rss")?.parse::<i64>().ok()?;
        Some(Measurement {
            peak: get("peak")?.parse().ok()?,
            total: get("total")?.parse().ok()?,
            allocs: get("allocs")?.parse().ok()?,
            rss: (rss >= 0).then_some(rss as usize),
            secs: get("secs")?.parse().ok()?,
            results: get("results")?.parse().ok()?,
        })
    }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Runs a single scenario and reports what it used.
fn measure(id: &str) -> Measurement {
    let scenario = common::scenarios()
        .into_iter()
        .find(|s| s.id() == id)
        .unwrap_or_else(|| panic!("unknown scenario {id:?}"));
    let model = common::model(scenario.model);

    // Anything allocated before this point (the model, the process itself) is
    // not part of what is being measured.
    let baseline = CURRENT.load(Ordering::Relaxed);
    PEAK.store(baseline, Ordering::Relaxed);
    TOTAL.store(0, Ordering::Relaxed);
    COUNT.store(0, Ordering::Relaxed);

    let start = Instant::now();
    let results = scenario.run(&model);
    let secs = start.elapsed().as_secs_f64();

    Measurement {
        peak: PEAK.load(Ordering::Relaxed).saturating_sub(baseline),
        total: TOTAL.load(Ordering::Relaxed),
        allocs: COUNT.load(Ordering::Relaxed),
        rss: peak_rss(),
        secs,
        results,
    }
}

fn run_child(id: &str) {
    println!("{}", measure(id).encode());
}

fn run_parent(filters: &[String]) {
    let scenarios = common::filtered_scenarios(filters);
    if scenarios.is_empty() {
        eprintln!("no scenarios matched {filters:?}");
        std::process::exit(1);
    }

    let exe = std::env::current_exe().expect("cannot locate the benchmark binary");

    println!(
        "{:<32} {:>12} {:>12} {:>12} {:>12} {:>10} {:>10}",
        "scenario", "peak alloc", "total alloc", "allocs", "peak RSS", "time (s)", "results"
    );

    for scenario in scenarios {
        let id = scenario.id();
        let output = Command::new(&exe)
            .arg("--run")
            .arg(&id)
            .output()
            .expect("failed to spawn the benchmark child process");

        if !output.status.success() {
            eprintln!("{id}: child process failed ({})", output.status);
            eprint!("{}", String::from_utf8_lossy(&output.stderr));
            continue;
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        let Some(measurement) = stdout.lines().find_map(Measurement::decode) else {
            eprintln!("{id}: child process produced no measurement");
            continue;
        };

        println!(
            "{:<32} {:>12} {:>12} {:>12} {:>12} {:>10.3} {:>10}",
            id,
            human_bytes(measurement.peak),
            human_bytes(measurement.total),
            measurement.allocs,
            measurement.rss.map_or("n/a".to_string(), human_bytes),
            measurement.secs,
            measurement.results,
        );
    }
}

fn main() {
    // Nothing may touch GMP or MPFR before this point.
    unsafe { install_gmp_hooks() };

    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.iter().position(|a| a == "--run") {
        Some(pos) => {
            let id = args.get(pos + 1).expect("--run needs a scenario id");
            run_child(id);
        }
        // Cargo passes flags such as `--bench` that do not apply here; anything
        // else is treated as a filter on the scenario ids.
        None => {
            let filters: Vec<String> = args.into_iter().filter(|a| !a.starts_with('-')).collect();
            run_parent(&filters);
        }
    }
}
