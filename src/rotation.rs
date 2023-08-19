use crate::utility::Coordinate;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Rotation {
    // Represents a 90 degree rotation in up to three dimensions
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl Rotation {
    pub const fn to_u8(self) -> u8 {
        self.x as u8 | ((self.y as u8) << 1) | ((self.z as u8) << 2)
    }
    pub const fn from_u8(rotation: u8) -> Self {
        Self {
            x: (rotation % 2) != 0,
            y: ((rotation / 2) % 2) != 0,
            z: ((rotation / 4) % 2) != 0,
        }
    }
    pub const fn rotate_coords(self, coords: Coordinate) -> Coordinate {
        // TODO: Cleanup
        let mut x = coords.x;
        let mut y = coords.y;
        let mut z = coords.z;

        if self.x {
            let temp = z;
            z = y;
            y = -temp;
        }

        if self.y {
            let temp = x;
            x = -z;
            z = temp;
        }

        if self.z {
            let temp = y;
            y = -x;
            x = temp;
        }
        Coordinate { x, y, z }
    }
    pub const fn flip_x(rotation: u8) -> u8 {
        if rotation % 2 == 1 {
            rotation - 1
        } else {
            rotation + 1
        }
    }
    pub const fn flip_y(rotation: u8) -> u8 {
        if (rotation / 2) % 2 == 1 {
            rotation - 2
        } else {
            rotation + 2
        }
    }
    pub const fn flip_z(rotation: u8) -> u8 {
        if (rotation / 4) % 2 == 1 {
            rotation - 4
        } else {
            rotation + 4
        }
    }
}
