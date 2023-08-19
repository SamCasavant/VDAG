use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Coordinate {
    pub x: i16,
    pub y: i16,
    pub z: i16,
}

#[derive(Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Axis {
    X,
    Y,
    Z,
}
