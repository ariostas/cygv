//! Interrupting a computation from another thread, as the Python bindings do on ctrl+c.
//!
//! The handle of a computation ([`Options::interrupt`](crate::Options::interrupt)) is made
//! the current one of the thread that runs it and of the threads of the pools it builds
//! ([`pool`]), and the pipelines poll it at their checkpoints ([`check`]): between the
//! degrees and passes of the modular pipeline, and between the degrees of the stages of the
//! reference one. So a computation stops within about the time of one degree.

use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// A handle that interrupts the computations it is given to: once [`Interrupt::interrupt`]
/// is called, from any thread, they stop at their next checkpoint and return
/// [`Error::Interrupted`](crate::Error::Interrupted). Clones share the same flag.
#[derive(Debug, Clone, Default)]
pub struct Interrupt(Arc<AtomicBool>);

impl Interrupt {
    /// A handle that has not been interrupted.
    pub fn new() -> Self {
        Self::default()
    }

    /// Interrupts the computations of this handle.
    pub fn interrupt(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    /// Whether [`Interrupt::interrupt`] was called.
    pub fn is_interrupted(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
}

/// Two handles are equal when they share their flag.
impl PartialEq for Interrupt {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for Interrupt {}

thread_local! {
    /// The handle of the computation that the thread runs.
    static CURRENT: RefCell<Option<Interrupt>> = const { RefCell::new(None) };
}

/// Runs `f` with `interrupt` as the current thread's handle, restoring the previous one
/// after.
pub(crate) fn with<T>(interrupt: Option<&Interrupt>, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Interrupt>);
    impl Drop for Restore {
        fn drop(&mut self) {
            CURRENT.with(|c| *c.borrow_mut() = self.0.take());
        }
    }
    let _restore = Restore(CURRENT.with(|c| c.replace(interrupt.cloned())));
    f()
}

/// A builder of a thread pool whose threads take the current thread's handle.
pub(crate) fn pool(threads: usize) -> rayon::ThreadPoolBuilder {
    let current = CURRENT.with(|c| c.borrow().clone());
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .start_handler(move |_| CURRENT.with(|c| *c.borrow_mut() = current.clone()))
}

/// Whether the computation that the current thread runs was interrupted.
pub(crate) fn interrupted() -> bool {
    CURRENT.with(|c| c.borrow().as_ref().is_some_and(Interrupt::is_interrupted))
}

/// `Err(e)` if the computation that the current thread runs was interrupted.
pub(crate) fn check<E>(e: E) -> Result<(), E> {
    match interrupted() {
        true => Err(e),
        false => Ok(()),
    }
}
