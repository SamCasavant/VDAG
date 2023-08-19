use crate::{node::Node, utility::Coordinate};

#[derive(PartialEq)]
enum Surface {
    // Represents a surface on a cube
    // Two dimensions represent which plane the surface is on
    // Pos/Neg represents whether the surface is greater or less along the unmentioned dimension
    XZPos,
    XYPos,
    YZPos,
    XZNeg,
    XYNeg,
    YZNeg,
}

#[derive(PartialEq)]
enum AdjacentSurface {
    // Represents a surface that touches another
    // Names are arbitrary
    Up,
    Right,
    Down,
    Left,
}

/*
   Surface::XZPos
   AdjacentSurface::Down
        _ _ _ _ _ _
  +y   /          /
   ^  /  XZPOS   / |  +z
   | /          /  |  ^
   |/_ _  _ _ _/   | /
   |           | R |/
   |           |   /
   |   DOWN    |  /
   |           | /
   |_ _ _ _ _ _|/_ _ _> +x


*/

// u64 data maps:
static ROTATE_90_DEG_CW_X: [u8; 64] = [
    36, 37, 32, 33, 38, 39, 34, 35, 44, 45, 40, 41, 46, 47, 42, 43, 4, 5, 0, 1, 6, 7, 2, 3, 12, 13,
    8, 9, 14, 15, 10, 11, 52, 53, 48, 49, 54, 55, 50, 51, 60, 61, 56, 57, 62, 63, 58, 59, 20, 21,
    16, 17, 22, 23, 18, 19, 28, 29, 24, 25, 30, 31, 26, 27,
];

static ROTATE_90_DEG_CW_Y: [u8; 64] = [
    9, 13, 11, 15, 8, 12, 10, 14, 41, 45, 43, 47, 40, 44, 42, 46, 25, 29, 27, 31, 24, 28, 26, 30,
    57, 61, 59, 63, 56, 60, 58, 62, 1, 5, 3, 7, 0, 4, 2, 6, 33, 37, 35, 39, 32, 36, 34, 38, 17, 21,
    19, 23, 16, 20, 18, 22, 49, 53, 51, 55, 48, 52, 50, 54,
];

static ROTATE_90_DEG_CW_Z: [u8; 64] = [
    9, 11, 8, 10, 13, 15, 12, 14, 25, 27, 24, 26, 29, 31, 28, 30, 1, 3, 0, 2, 5, 7, 4, 6, 17, 19,
    16, 18, 21, 23, 20, 22, 41, 43, 40, 42, 45, 47, 44, 46, 57, 59, 56, 58, 61, 63, 60, 62, 33, 35,
    32, 34, 37, 39, 36, 38, 49, 51, 48, 50, 53, 55, 52, 54,
];

// [Node; 8] maps:

static XZNEG_NODES_TOP: [u8; 8] = [6, 7, 4, 5, 2, 3, 0, 1];
static XYPOS_NODES_TOP: [u8; 8] = [2, 3, 6, 7, 0, 1, 4, 5];
static XYNEG_NODES_TOP: [u8; 8] = [4, 5, 0, 1, 6, 7, 2, 3];
static YZPOS_NODES_TOP: [u8; 8] = [2, 0, 3, 1, 6, 4, 7, 5];
static YZNEG_NODES_TOP: [u8; 8] = [1, 3, 0, 2, 5, 7, 4, 6];

static RIGHT_NODES_FRONT: [u8; 8] = [1, 5, 3, 7, 0, 4, 2, 6];
static UP_NODES_FRONT: [u8; 8] = [5, 4, 7, 6, 1, 0, 3, 2];
static LEFT_NODES_FRONT: [u8; 8] = [4, 0, 6, 2, 5, 1, 7, 3];
// XZ-R * XY+R = YZ-D
// XZ-R * XY-R = YZ+D
// XZ-R * XZ+R = XZ-U
// XZ-R * XZ-R = XZ+D
// XZ-R * YZ+R = XY-U
// XZ-R * YZ-R = XY+U

pub struct Transform {
    top: Surface,
    front: AdjacentSurface,
    reflected: bool,
}

impl Transform {
    // Represents a transformation on a cube
    // Presently represents only orientation
    pub fn from_u8(mut transform: u8) -> Self {
        // u8 = [2 Bits Unused] [3 Bits Top] [2 Bits Front]  [1 Bit Reflection]
        let reflected = (transform % 2) != 0;
        transform /= 2;
        let front = match transform % 4 {
            0 => AdjacentSurface::Up,
            1 => AdjacentSurface::Right,
            2 => AdjacentSurface::Down,
            3 => AdjacentSurface::Left,
            _ => unreachable!(),
        };
        transform /= 4;
        let top = match transform % 8 {
            0 => Surface::XZPos,
            1 => Surface::XYPos,
            2 => Surface::YZPos,
            3 => Surface::XZNeg,
            4 => Surface::XYNeg,
            5 => Surface::YZNeg,
            _ => panic!("Undefined behavior"),
        };
        return Self {
            reflected,
            front,
            top,
        };
    }

    pub fn to_u8(&self) -> u8 {
        let top = match self.top {
            Surface::XZPos => 0,
            Surface::XYPos => 1,
            Surface::YZPos => 2,
            Surface::XZNeg => 3,
            Surface::XYNeg => 4,
            Surface::YZNeg => 5,
        };

        let front = match self.front {
            AdjacentSurface::Up => 0,
            AdjacentSurface::Right => 1,
            AdjacentSurface::Down => 2,
            AdjacentSurface::Left => 3,
        };

        let reflected = self.reflected as u8;

        reflected | front << 1 | top << 3
    }

    fn apply_data(&self, data: u64) -> u64 {
        // Apply transformation to a u64 of data
        let mut transformed = data;
        if (self.reflected) {
            // Reflect across all three axes (maybe one axis would be okay, this makes more sense to me abstractly)
            //  Reflect X
            //      There are 8 regions of 8 bits.
            //      In even regions, indices are 7 or 9 below their x reflection, toggling at each index.
            //      Odd regions are the other half of the reflection, so only even regions need to be considered.
            for region in 0..4 {
                for sub_index in 0..8 {
                    let index = region * 16 + sub_index;
                    transformed = if sub_index % 2 == 0 {
                        Self::swap_bits(transformed, index, index + 9)
                    } else {
                        Self::swap_bits(transformed, index, index + 7)
                    }
                }
            }
            //  Reflect Y
            //      There are 4 regions of 16 bits.
            //      In even regions, indices are 18 or 14 below their y reflection, toggling at each pair of indices.
            for region in 0..2 {
                for sub_index in 0..16 {
                    let index = region * 32 + sub_index;
                    transformed = if (sub_index / 2) % 2 == 0 {
                        Self::swap_bits(transformed, index, index + 18)
                    } else {
                        Self::swap_bits(transformed, index, index + 14)
                    }
                }
            }
            //  Reflect Z
            //      Divide the first 32 bits into groups of 4.
            //      Even groups are 36 below their z reflection, odds are 28 below their z reflection.
            //      The second 32 are the other half of the reflection, and need no consideration.
            for index in 0..32 {
                transformed = if (index / 4) % 2 == 0 {
                    Self::swap_bits(transformed, index, index + 36)
                } else {
                    Self::swap_bits(transformed, index, index + 28)
                }
            }
        }

        // Rotations (oh boy (why am I doing this morton encoded what is wrong with me))
        // Initially, the top is in the XZPos plane
        if (self.top == Surface::XYNeg) {
            // Rotating 90 degrees around X changes the top to the XYNeg plane
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
        } else if (self.top == Surface::XZNeg) {
            // Rotate 180 degrees around the X axis to XZNeg
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
        } else if (self.top == Surface::XYPos) {
            // Rotate 270 degrees around X axis to XYPos
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_X);
        } else if (self.top == Surface::YZNeg) {
            // Rotate 90 degrees around Z axis to ZYNeg
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Z);
        } else if (self.top == Surface::YZPos) {
            // Rotate 270 degrees around Z axis to ZYPos
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Z);
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Z);
            transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Z);
        }

        // Rotate around Y axis until the front is correct
        if (self.front == AdjacentSurface::Down) {
            return transformed;
        }
        transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Y);
        if (self.front == AdjacentSurface::Left) {
            return transformed;
        }
        transformed = Self::map_bits(transformed, ROTATE_90_DEG_CW_Y);
        if (self.front == AdjacentSurface::Up) {
            return transformed;
        }
        return Self::map_bits(transformed, ROTATE_90_DEG_CW_Y);
    }

    pub fn transform_coords(&self, coords: Coordinate) -> Coordinate {
        let mut x = coords.x;
        let mut y = coords.y;
        let mut z = coords.z;

        if self.reflected {
            x = -x;
            y = -y;
            z = -z;
        }

        match self.top {
            Surface::XYPos => {
                let temp = z;
                z = -y;
                y = temp;
            }
            Surface::YZPos => {
                let temp = x;
                x = -y;
                y = temp;
            }
            Surface::XZNeg => {
                z = -z;
                y = -y;
            }
            Surface::XYNeg => {
                let temp = y;
                y = -z;
                z = temp;
            }
            Surface::YZNeg => {
                let temp = y;
                y = -x;
                x = temp;
            }
            Surface::XZPos => {}
        }

        match self.front {
            AdjacentSurface::Down => return Coordinate { x, y, z },
            AdjacentSurface::Right => return Coordinate { x: z, y, z: -x },
            AdjacentSurface::Up => return Coordinate { x: -x, y, z: -z },
            AdjacentSurface::Left => return Coordinate { x: -z, y, z: x },
        }
    }

    pub fn get_symmetries_data(data: u64) -> Vec<Transform> {
        // Brute force search for all transforms that are no-ops for given data
        let mut symmetries = Vec::new();
        for transform in 0b000000..0b101111 {
            if Transform::from_u8(transform).apply_data(data) == data {
                symmetries.push(Transform::from_u8(transform));
            }
        }
        return symmetries;
    }

    fn apply_nodes(&self, nodes: [Node<u8, u32>; 8]) -> [Node<u8, u32>; 8] {
        // Children are ordered like this:
        // (x, y, z) =>
        // [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        let mut transformed = nodes;
        if self.reflected {
            // Exchange every pair of values
            for index in [0, 2, 4, 6] {
                transformed.swap(index, index + 1);
            }
            // Exchange pairs of (index, index+2)
            for index in [0, 1, 4, 5] {
                transformed.swap(index, index + 2);
            }
            // Exchange pairs of (index, index+4)
            for index in [0, 1, 2, 3] {
                transformed.swap(index, index + 4);
            }
        }
        transformed = match self.top {
            Surface::XYNeg => Self::map_elements(transformed, XYNEG_NODES_TOP),
            Surface::XYPos => Self::map_elements(transformed, XYPOS_NODES_TOP),
            Surface::XZNeg => Self::map_elements(transformed, XZNEG_NODES_TOP),
            Surface::XZPos => transformed,
            Surface::YZNeg => Self::map_elements(transformed, YZNEG_NODES_TOP),
            Surface::YZPos => Self::map_elements(transformed, YZPOS_NODES_TOP),
        };

        transformed = match self.front {
            AdjacentSurface::Down => transformed,
            AdjacentSurface::Left => Self::map_elements(transformed, LEFT_NODES_FRONT),
            AdjacentSurface::Right => Self::map_elements(transformed, RIGHT_NODES_FRONT),
            AdjacentSurface::Up => Self::map_elements(transformed, UP_NODES_FRONT),
        };

        // TODO: Don't update transforms for elements on which it has no effect
        // Update the transforms for all of the stored nodes
        for node in transformed {}
        return transformed;
    }

    const fn swap_bits(mut data: u64, first: u8, second: u8) -> u64 {
        // Swaps the bit at {first} position with the bit at {second} in {data}
        let first_bit = (data >> first) & 1;
        let second_bit = (data >> second) & 1;
        if first_bit != second_bit {
            let mut x = first_bit ^ second_bit;
            x = (x << first) | (x << second);
            data ^= x;
        }
        data
    }

    pub fn map_bits(data: u64, map: [u8; 64]) -> u64 {
        let mut mapped = 0;
        for index in map {
            if (Self::get_bit_at_index(data, index)) {
                mapped += 1;
            }
            mapped <<= 1;
        }
        return mapped;
    }
    const fn get_bit_at_index(data: u64, index: u8) -> bool {
        return (data & (1 << index - 1)) != 0;
    }

    pub fn map_elements(elements: [Node<u8, u32>; 8], map: [u8; 8]) -> [Node<u8, u32>; 8] {
        let mut mapped = [Node::new(0, u32::MAX); 8];
        for index in 0..8 {
            mapped[index] = elements[map[index] as usize]
        }
        return mapped;
    }

//    pub fn transform(&self, other: &Self) -> Self {
        // Combines two transform operations
//        let reflected = self.reflected ^ other.reflected;
//        let top;
//        let front;
//        if self.top == Surface::XZPos {
//            top = other.top;
//        } else if (other.top == Surface::XZPos) {
//            top = self.top;
//        } else {
//        }
//    }
}
