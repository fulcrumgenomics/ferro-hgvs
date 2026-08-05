//! The chunked, order-preserving streaming engine behind the batch APIs (#975).
//!
//! One engine, two public surfaces: [`crate::batch::BatchProcessor`]'s streaming
//! methods and [`crate::parallel`]'s streaming functions. It lives here rather
//! than beside either of them because `parallel` is behind the `parallel`
//! feature while `batch` is not, so `batch` cannot borrow from `parallel` — and a
//! second copy of a chunk-and-drain loop is a second thing to drift.
//!
//! The shape is the CLI's `run_batch`/`flush_chunk` pair (#687), which #975 asks
//! to reuse: read a bounded chunk, map it order-preservingly, drain it, repeat.
//! Peak memory is therefore a function of the chunk size rather than of the
//! stream's length, which is the whole point — the `Vec`-based APIs materialize
//! both their input and their output by contract and no internal chunking can fix
//! that.

use rayon::prelude::*;

/// Build the rayon pool a chunked stream should run on.
///
/// `Some` is a dedicated pool; `None` means "not a dedicated pool", which the
/// `serial` flag then disambiguates between running inline and using the global
/// pool. Built once per stream rather than once per chunk: a
/// `ThreadPoolBuilder::build()` per chunk would spawn and join `num_threads` OS
/// threads every few thousand items.
fn dedicated_pool(num_threads: usize) -> Option<rayon::ThreadPool> {
    if num_threads <= 1 {
        return None;
    }
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Ok(pool) => Some(pool),
        Err(err) => {
            log::warn!(
                "failed to build a dedicated {num_threads}-thread pool ({err}); \
                 falling back to the global rayon pool"
            );
            None
        }
    }
}

/// Chunked, order-preserving streaming map.
///
/// Generic over the mapped output because the two surfaces yield different
/// types — `parallel` yields `Result`, the batch processor yields `ItemResult` —
/// and that is the only thing that differed between the two copies this replaces.
pub(crate) struct StreamingMap<I, S, F, T> {
    input: I,
    map: F,
    chunk_size: usize,
    /// Reused input buffer: one allocation per stream, not per chunk.
    inputs: Vec<S>,
    /// The current chunk's results, drained front-first so order is preserved.
    pending: std::vec::IntoIter<T>,
    pool: Option<rayon::ThreadPool>,
    /// `num_threads == 1`: map inline rather than fanning out at all.
    serial: bool,
}

impl<I, S, F, T> StreamingMap<I, S, F, T>
where
    I: Iterator<Item = S>,
{
    /// `num_threads`: `1` is serial, `0` is the global rayon pool, `n` is a
    /// dedicated `n`-thread pool (falling back to the global pool with a warning
    /// if one cannot be built) — the same contract as `BatchProcessor`'s
    /// `Vec`-based methods.
    pub(crate) fn new<J>(variants: J, num_threads: usize, chunk_size: usize, map: F) -> Self
    where
        J: IntoIterator<Item = S, IntoIter = I>,
    {
        StreamingMap {
            input: variants.into_iter(),
            map,
            chunk_size,
            inputs: Vec::new(),
            pending: Vec::new().into_iter(),
            pool: dedicated_pool(num_threads),
            serial: num_threads == 1,
        }
    }
}

impl<I, S, F, T> Iterator for StreamingMap<I, S, F, T>
where
    I: Iterator<Item = S>,
    S: AsRef<str> + Send + Sync,
    F: Fn(&S) -> T + Sync + Send,
    T: Send,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            if let Some(next) = self.pending.next() {
                return Some(next);
            }
            // The chunk is drained; read the next one. `take` moves the buffer
            // out so `self.input` can be borrowed mutably to refill it.
            let mut inputs = std::mem::take(&mut self.inputs);
            inputs.clear();
            inputs.extend(self.input.by_ref().take(self.chunk_size));
            if inputs.is_empty() {
                self.inputs = inputs;
                return None;
            }
            let results: Vec<T> = if self.serial {
                inputs.iter().map(&self.map).collect()
            } else {
                match self.pool.as_ref() {
                    Some(pool) => pool.install(|| inputs.par_iter().map(&self.map).collect()),
                    None => inputs.par_iter().map(&self.map).collect(),
                }
            };
            self.inputs = inputs;
            self.pending = results.into_iter();
        }
    }
}
