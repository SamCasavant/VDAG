use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Node<T, U> {
    pub val: T,
    pub next: U,
}

impl<T, U> Node<T, U> {
    pub fn new(val: T, next: U) -> Self {
        Self { val, next }
    }
}
