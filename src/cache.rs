//! A module for caching dynamically-allocated numerical variables.
//!
//! This module is useful to avoid frequently allocating and deallocating
//! arbitrary-precision rational and multi-precision floating-point numbers.

/// A structure to cache numbers.
///
/// The structure stores a clonable template number so that we can easily
/// implement creating new numbers when the cache is empty without having to
/// worry about the specifics of how they should be initialized.
#[derive(Clone, Debug)]
pub struct NumCache<T> {
    template_num: T,
    nums: Vec<T>,
    max_nums: usize,
}

impl<T: Clone> NumCache<T> {
    /// Create a new number cache with a given maximum size.
    pub fn new(template_num: T, max_nums: usize) -> Self {
        let nums = Vec::with_capacity(max_nums);
        Self {
            template_num,
            nums,
            max_nums,
        }
    }

    /// Pop a number from the cache. If the cache is empty the template number
    /// is cloned.
    pub fn pop(&mut self) -> T {
        if let Some(c) = self.nums.pop() {
            c
        } else {
            self.template_num.clone()
        }
    }

    /// Push a number to the cache. If the cache is full the number is dropped.
    pub fn push(&mut self, c: T) {
        if self.nums.len() < self.max_nums {
            self.nums.push(c);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Rational;

    #[test]
    fn test_cache() {
        let zero = Rational::new();
        let one = Rational::from(1);
        let two = Rational::from(2);

        let mut cache = NumCache::new(zero.clone(), 2);

        cache.push(one.clone());
        cache.push(two.clone());
        cache.push(two.clone());
        assert_eq!(cache.pop(), two);
        assert_eq!(cache.pop(), one);
        assert_eq!(cache.pop(), zero);
    }
}
