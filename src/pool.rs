//! Object pools of dynamically-allocated numerical variables.

/// A structure to pool numerical variables.
///
/// The structure stores a clonable template number so that we can easily
/// implement creating new numbers when the pool is empty without having to
/// worry about the specifics of how they should be initialized.
#[derive(Clone, Debug)]
pub struct NumberPool<T> {
    template_num: T,
    nums: Vec<T>,
    max_nums: usize,
}

impl<T: Clone> NumberPool<T> {
    /// Create a new number pool with a given maximum size.
    pub fn new(template_num: T, max_nums: usize) -> Self {
        let nums = Vec::with_capacity(max_nums);
        Self {
            template_num,
            nums,
            max_nums,
        }
    }

    /// Pop a number from the pool. If the pool is empty the template number
    /// is cloned.
    pub fn pop(&mut self) -> T {
        if let Some(c) = self.nums.pop() {
            c
        } else {
            self.template_num.clone()
        }
    }

    /// Push a number to the pool. If the pool is full the number is dropped.
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
    fn test_pool() {
        let zero = Rational::new();
        let one = Rational::from(1);
        let two = Rational::from(2);

        let mut pool = NumberPool::new(zero.clone(), 2);

        pool.push(one.clone());
        pool.push(two.clone());
        pool.push(two.clone());
        assert_eq!(pool.pop(), two);
        assert_eq!(pool.pop(), one);
        assert_eq!(pool.pop(), zero);
    }
}
