use crate::utility::Coordinate;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Reflection {
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl Reflection {
    pub const fn to_u8(self) -> u8 {
        self.x as u8 | ((self.y as u8) << 1) | ((self.z as u8) << 2)
    }
    pub const fn from_u8(reflection: u8) -> Self {
        Self {
            x: (reflection % 2) != 0,
            y: ((reflection / 2) % 2) != 0,
            z: ((reflection / 4) % 2) != 0,
        }
    }
    pub const fn reflect_coords(self, coords: Coordinate) -> Coordinate {
        let x = if self.x {
            (-(coords.x as i32) - 1) as i16
        } else {
            coords.x
        };
        let y = if self.y {
            (-(coords.y as i32) - 1) as i16
        } else {
            coords.y
        };
        let z = if self.z {
            (-(coords.z as i32) - 1) as i16
        } else {
            coords.z
        };
        Coordinate { x, y, z }
    }
    pub const fn flip_x(reflection: u8) -> u8 {
        if reflection % 2 == 1 {
            reflection - 1
        } else {
            reflection + 1
        }
    }
    pub const fn flip_y(reflection: u8) -> u8 {
        if (reflection / 2) % 2 == 1 {
            reflection - 2
        } else {
            reflection + 2
        }
    }
    pub const fn flip_z(reflection: u8) -> u8 {
        if (reflection / 4) % 2 == 1 {
            reflection - 4
        } else {
            reflection + 4
        }
    }
}
