use std::{
    iter::Enumerate,
    ops::{Index, IndexMut},
    slice::ChunksExact,
};

use serde::{Deserialize, Serialize};

use crate::node::Node;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Layer<T, U>(pub Vec<Node<T, U>>);

impl<T, U> Index<usize> for Layer<T, U> {
    type Output = Node<T, U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T, U> IndexMut<usize> for Layer<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T, U> IntoIterator for Layer<T, U> {
    type Item = Node<T, U>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T, U> Layer<T, U> {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn len(&self) -> usize {
        // Returns length of content
        self.0.len()
    }

    pub fn cubes_iter(&self) -> Enumerate<ChunksExact<'_, Node<T, U>>> {
        // Returns an iterator over 2x2x2 cubes of a layer
        self.0.chunks_exact(8).into_iter().enumerate()
    }

    pub fn free(&mut self) {
        // Deletes content of layer and frees memory
        self.0.clear();
        self.0.shrink_to_fit();
    }

    pub fn shrink_to_fit(&mut self) {
        // Reduces capacity of content to current length
        self.0.shrink_to_fit();
    }

    pub fn push(&mut self, element: Node<T, U>) {
        // Pushes an element to underlying vector
        self.0.push(element);
    }

    pub fn capacity(&self) -> usize {
        self.0.capacity()
    }

    pub unsafe fn get_unchecked(&self, index: usize) -> &Node<T, U> {
        // Gets the value at index without bounds checking
        self.0.get_unchecked(index)
    }
}
