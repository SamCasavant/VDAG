#![feature(test)]
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize)]
struct Voxel {
    x: String,
    y: String,
    z: String,
}

#[derive(Serialize, Deserialize)]
struct Dimension {
    width: String,
    height: String,
    depth: String,
}

#[derive(Serialize, Deserialize)]
struct VoxData {
    dimension: Vec<Dimension>,
    voxels: Vec<Voxel>,
}
fn main() {
    let mut world = UChunkDag::new();

    // let vox_data_str = fs::read_to_string("assets/torus.json").unwrap();
    // let vox_data: VoxData = serde_json::from_str(&vox_data_str).unwrap();
    // let mut count = 0;
    // for voxel in vox_data.voxels.iter() {
    //     let pos_coords = Coordinate {
    //         x: voxel.x.parse::<i16>().unwrap(),
    //         y: voxel.y.parse::<i16>().unwrap(),
    //         z: voxel.z.parse::<i16>().unwrap(),
    //     };
    //     let neg_coords = Coordinate {
    //         x: -voxel.x.parse::<i16>().unwrap(),
    //         y: -voxel.y.parse::<i16>().unwrap(),
    //         z: -voxel.z.parse::<i16>().unwrap(),
    //     };
    //     world.set_pos(pos_coords);
    //     world.set_pos(neg_coords);
    //     count += 2;
    // }
    // println!("{} voxels", count);
    world.set_pos(Coordinate { x: 0, y: 0, z: 0 });

    let mut compressed = ChunkDAG::new(world);
    for n in 1..=8 {
        compressed.compress(n);
    }
    for x in 0..1000 {
        compressed.get_pos(Coordinate { x: 0, y: 0, z: 0 });
        compressed.get_pos(Coordinate { x: 1, y: 1, z: 1 });
        compressed.get_pos(Coordinate {
            x: 100,
            y: 100,
            z: 100,
        });
    }
    // for voxel in vox_data.voxels.iter() {
    //     let pos_coords = Coordinate {
    //         x: voxel.x.parse::<i16>().unwrap(),
    //         y: voxel.y.parse::<i16>().unwrap(),
    //         z: voxel.z.parse::<i16>().unwrap(),
    //     };
    //     let neg_coords = Coordinate {
    //         x: -voxel.x.parse::<i16>().unwrap(),
    //         y: -voxel.y.parse::<i16>().unwrap(),
    //         z: -voxel.z.parse::<i16>().unwrap(),
    //     };
    //     assert!(
    //         compressed.get_pos(pos_coords),
    //         "Failure at {:?}",
    //         pos_coords
    //     );
    //     assert!(
    //         compressed.get_pos(neg_coords),
    //         "Failure at {:?}",
    //         neg_coords
    //     );
    // }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Coordinate {
    x: i16,
    y: i16,
    z: i16,
}

pub struct UChunkDag {
    /// Uncompressed Chunk DAG
    /// This is a wrapper around an array, which provides write access to voxel data based on coordinates.
    /// Indices across each dimension in a 2x2x2 block proceed as follows:
    /// [(-x, -y, -z), (+x, -y, -z), (-x, +y, -z), (+x, +y, -z),
    ///  (-x, -y, +z), (+x, -y, +z), (-x, +y, +z), (+x, +y, +z)]
    /// At higher levels of abstraction, the ordering is the same, so, e.g. a 2x2x2 block whose coordinates
    /// are all less than another 2x2x2 block across x and y but greater across z, will appear later in the array.
    ///
    /// The structure can be imagined as a 3d cube recursively divided into 8 octant (3d quadrant) subregions.
    /// At each layer of abstraction, 8x as many voxels are represented, and each layer shares the index ordering above.
    ///
    /// By convention, the bottom layer contains u64s which are structured as above. Each bit in these ints represents a single voxel being present or absent.
    /// There are 8^9 u64s, which use ~1 GB before compression. Beyond this point, I would use all my RAM.
    /// Because the u64s represent 64 voxels, this struct can represent 8^11 voxels, a 2048x2048x2048 cube.
    data: Box<[u64; 134217728]>,
}
impl UChunkDag {
    pub fn new() -> Self {
        let data = Box::new([0; 134_217_728]);
        Self { data }
    }
    pub fn set_pos(&mut self, coords: Coordinate) {
        //Sets a voxel to 1
        let (index, subindex) = Self::coords_to_indices(coords);
        let updater = 1u64 << (63 - subindex);
        self.data[index] |= updater;
    }
    pub fn set_all(&mut self, coords: Coordinate) {
        //Sets a 4x4x4 chunk to all 1s
        let (index, _subindex) = Self::coords_to_indices(coords);
        let updater = u64::MAX;
        self.data[index] = updater;
    }

    fn coords_to_indices(coords: Coordinate) -> (usize, u8) {
        // TODO: Clean this up
        let (x, y, z) = (coords.x + 1024, coords.y + 1024, coords.z + 1024);

        let xi = (x / 4) as u32;
        let mut x_index = xi;
        for n in 1..=10 {
            x_index += 2u32.pow(3 * (n - 1)) * 6 * (xi / 2u32.pow(n));
        }

        let yi = (y / 4) as u32;
        let mut y_index = 2 * yi;
        for n in 1..10 {
            y_index += 2 * (2u32.pow(3 * (n - 1)) * 6 * (yi / 2u32.pow(n)));
        }

        let zi = (z / 4) as u32;
        let mut z_index = 4 * zi;
        for n in 1..10 {
            z_index += 4 * (2u32.pow(3 * (n - 1)) * 6 * (zi / 2u32.pow(n)));
        }
        let index = (x_index + y_index + z_index) as usize;

        let sub_i_x = x % 4;
        let sub_i_y = y % 4;
        let sub_i_z = z % 4;
        let sub_index_x = sub_i_x + 6 * (sub_i_x / 2);
        let sub_index_y = 2 * sub_i_y + 12 * (sub_i_y / 2);
        let sub_index_z = 4 * sub_i_z + 24 * (sub_i_z / 2);
        let sub_index = (sub_index_x + sub_index_y + sub_index_z) as u8;

        (index, sub_index)
    }
}

#[derive(Debug, PartialEq)]
struct Reflection {
    x: bool,
    y: bool,
    z: bool,
}
impl Reflection {
    fn to_u8(&self) -> u8 {
        self.x as u8 | ((self.y as u8) << 1) | ((self.z as u8) << 2)
    }
    fn from_u8(reflection: u8) -> Self {
        Self {
            x: (reflection % 2) != 0,
            y: ((reflection / 2) % 2) != 0,
            z: ((reflection / 4) % 2) != 0,
        }
    }
    fn reflect_coords(&self, coords: Coordinate) -> Coordinate {
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
    fn flip_x(reflection: u8) -> u8 {
        if reflection % 2 == 1 {
            reflection - 1
        } else {
            reflection + 1
        }
    }
    fn flip_y(reflection: u8) -> u8 {
        if (reflection / 2) % 2 == 1 {
            reflection - 2
        } else {
            reflection + 2
        }
    }
    fn flip_z(reflection: u8) -> u8 {
        if (reflection / 4) % 2 == 1 {
            reflection - 4
        } else {
            reflection + 4
        }
    }
}

#[derive(Debug)]
pub struct ChunkDAG {
    //TODO: Most of these don't need to be u32, but consistency made the code easier to generalize. Generics?
    //TODO: Consider cost/benefits of shortest-common-substring compression (cost: index of higher layer cannot be divided by 8)
    //In progress: Symmetrical compression
    layer9: Vec<(u8, u32)>,
    layer8: Vec<(u8, u32)>,
    layer7: Vec<(u8, u32)>,
    layer6: Vec<(u8, u32)>,
    layer5: Vec<(u8, u32)>,
    layer4: Vec<(u8, u32)>,
    layer3: Vec<(u8, u32)>,
    layer2: Vec<(u8, u32)>,
    //TODO: This is the same size as (u32, u32) and only 3 bits are used in u8, so I can fit 29 more bits/node.
    //       16 for material information, 13 for additional transforms?
    //       Or with u16 indices, there are 13 remaining bits.
    layer1: Vec<(u8, u32)>,
    data: Vec<u64>,
}
impl ChunkDAG {
    fn new(chunk: UChunkDag) -> Self {
        let leaves = Box::new(
            chunk
                .data
                .iter()
                .enumerate()
                .collect::<Vec<(usize, &u64)>>(),
        );
        let mut reflected_leaves = Vec::new();
        for leaf in leaves.iter() {
            let (reflected, reflection) = Self::least_symmetrical_data(*leaf.1);
            reflected_leaves.push((leaf.0, reflected, reflection.to_u8()));
        }

        // Sort in reverse order so we can pop elements without shifting
        reflected_leaves.sort_by(|a, b| {
            if a.1 == b.1 {
                b.0.cmp(&a.0)
            } else {
                b.1.cmp(&a.1)
            }
        });
        let mut layer1_contents = vec![(0, 0); 134_217_728];
        // Data without duplicates
        let mut data = Vec::new();

        println!("Compressing...");
        while !reflected_leaves.is_empty() {
            let active_leaf = reflected_leaves.pop().unwrap();
            let active_index = match active_leaf.1 {
                0 => u32::MAX,
                val => {
                    data.push(val);
                    (data.len() - 1) as u32
                }
            };
            layer1_contents[active_leaf.0] = (active_leaf.2, active_index);
            while !reflected_leaves.is_empty()
                && reflected_leaves[reflected_leaves.len() - 1].1 == active_leaf.1
            {
                let removed = reflected_leaves.pop().unwrap();
                layer1_contents[removed.0] = (removed.2, active_index);
            }
        }
        Self {
            layer9: Vec::new(),
            layer8: Vec::new(),
            layer7: Vec::new(),
            layer6: Vec::new(),
            layer5: Vec::new(),
            layer4: Vec::new(),
            layer3: Vec::new(),
            layer2: Vec::new(),
            layer1: layer1_contents,
            data,
        }
    }
    fn least_symmetrical_data(data: u64) -> (u64, Reflection) {
        // Takes a 4x4x4 block of voxel data and returns the reflection with the lowest value
        // 'least' is chosen arbitrarily, but guarantees a consistent output for a similar data
        let mut min = data;
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        for x in [false, true] {
            for y in [false, true] {
                for z in [false, true] {
                    let test = Self::reflect_data(data, Reflection { x, y, z });
                    if test < min {
                        min = test;
                        reflection = Reflection { x, y, z };
                    }
                }
            }
        }
        (min, reflection)
    }
    fn reflect_data(data: u64, reflection: Reflection) -> u64 {
        // Reflects a 4x4x4 block across a list of axes
        // TODO: Simplify with math! This can probably be cleaner AND faster

        let mut reflected = data;
        if reflection.x {
            // There are 8 regions of 8 bits.
            // In even regions, indices are 7 or 9 below their x reflection, toggling at each index.
            // Odd regions are the other half of the reflection, so only even regions need to be considered.
            for region in 0..4 {
                for sub_index in 0..8 {
                    let index = region * 16 + sub_index;
                    reflected = if sub_index % 2 == 0 {
                        Self::swap_bits(reflected, index, index + 9)
                    } else {
                        Self::swap_bits(reflected, index, index + 7)
                    }
                }
            }
        }
        if reflection.y {
            // There are 4 regions of 16 bits.
            // In even regions, indices are 18 or 14 below their y reflection, toggling at each pair of indices.
            for region in 0..2 {
                for sub_index in 0..16 {
                    let index = region * 32 + sub_index;
                    reflected = if (sub_index / 2) % 2 == 0 {
                        Self::swap_bits(reflected, index, index + 18)
                    } else {
                        Self::swap_bits(reflected, index, index + 14)
                    }
                }
            }
        }
        if reflection.z {
            // Divide the first 32 bits into groups of 4.
            // Even groups are 36 below their z reflection, odds are 28 below their z reflection.
            // The second 32 are the other half of the reflection, and need no consideration.
            for index in 0..32 {
                reflected = if (index / 4) % 2 == 0 {
                    Self::swap_bits(reflected, index, index + 36)
                } else {
                    Self::swap_bits(reflected, index, index + 28)
                }
            }
        }
        reflected
    }
    fn swap_bits(mut data: u64, first: u8, second: u8) -> u64 {
        let first_bit = (data >> first) & 1;
        let second_bit = (data >> second) & 1;
        if first_bit != second_bit {
            let mut x = first_bit ^ second_bit;
            x = (x << first) | (x << second);
            data ^= x;
        }
        data
    }

    fn compress(&mut self, layer: u8) {
        // Compresses layer n, producing layer n+1
        let (lower, upper) = match layer {
            // These two variables have to be assigned in the same statement to
            // make rust accept that I'm not borrowing the same variable twice,
            // even though layer is immutable.
            // THANKS RUST (this code is a little cleaner though, so actually thanks, Rust)
            1 => (&mut self.layer1, &mut self.layer2),
            2 => (&mut self.layer2, &mut self.layer3),
            3 => (&mut self.layer3, &mut self.layer4),
            4 => (&mut self.layer4, &mut self.layer5),
            5 => (&mut self.layer5, &mut self.layer6),
            6 => (&mut self.layer6, &mut self.layer7),
            7 => (&mut self.layer7, &mut self.layer8),
            8 => (&mut self.layer8, &mut self.layer9),
            n => panic!("Can't compress layer {}", n),
        };
        let lower_chunks: Vec<&[(u8, u32)]> = lower.chunks(8).collect();

        let mut lower_data = Vec::new();
        for (i, chunk) in lower_chunks.iter().enumerate() {
            lower_data.push((i, *chunk));
        }

        let mut reflected_lower_data = Vec::new();
        for &children in lower_data.iter() {
            let (reflected, reflection) =
                Self::least_symmetrical_nodes(children.1.try_into().unwrap());
            reflected_lower_data.push((children.0, reflected, reflection.to_u8()));
        }

        // Sort in reverse order so we can pop elements without shifting
        reflected_lower_data.sort_by(|a, b| {
            if a.1 == b.1 {
                if a.2 == b.2 {
                    b.0.cmp(&a.0)
                } else {
                    b.2.cmp(&a.2)
                }
            } else {
                b.1.cmp(&a.1)
            }
        });

        let upper_length = match layer {
            1 => 16_777_216,
            2 => 2_097_152,
            3 => 262_144,
            4 => 32_768,
            5 => 4_096,
            6 => 512,
            7 => 64,
            8 => 8,
            _ => panic!("Can't compress given layer"),
        };

        let mut upper_contents = vec![(0, 0); upper_length];
        // Data without duplicates
        let mut uniq_lower_chunks: Vec<[(u8, u32); 8]> = Vec::new();

        println!("Compressing...");
        while !reflected_lower_data.is_empty() {
            let active_node = reflected_lower_data.pop().unwrap();
            let empty = [(0, u32::MAX); 8];
            let active_index = match active_node.1 {
                null if null == empty => u32::MAX,
                val => {
                    uniq_lower_chunks.push(val);
                    (uniq_lower_chunks.len() - 1) as u32
                }
            };
            upper_contents[active_node.0] = (active_node.2, active_index);
            while !reflected_lower_data.is_empty()
                && (
                    reflected_lower_data[reflected_lower_data.len() - 1].1,
                    reflected_lower_data[reflected_lower_data.len() - 1].2,
                ) == (active_node.1, active_node.2)
            {
                let removed = reflected_lower_data.pop().unwrap();
                upper_contents[removed.0] = (removed.2, active_index);
            }
        }
        let mut lower_contents = Vec::new();
        for chunk in uniq_lower_chunks.iter() {
            for element in chunk.iter() {
                lower_contents.push(*element);
            }
        }
        *lower = lower_contents;
        *upper = upper_contents;
    }
    fn least_symmetrical_nodes(nodes: [(u8, u32); 8]) -> ([(u8, u32); 8], Reflection) {
        // Takes a nodes children and returns the lowest (according to however Rust orders slices) similar children by reflection
        let mut min = nodes;
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        for x in [false, true] {
            for y in [false, true] {
                for z in [false, true] {
                    let test = Self::reflect_nodes(nodes, Reflection { x, y, z });
                    if test < min {
                        min = test;
                        reflection = Reflection { x, y, z };
                    }
                }
            }
        }

        (min, reflection)
    }
    fn reflect_nodes(nodes: [(u8, u32); 8], reflection: Reflection) -> [(u8, u32); 8] {
        // Children are ordered like this:
        // (x, y, z) =>
        // [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        let mut reflected = nodes;
        if reflection.x {
            // Exchange every pair of values and flip the x reflection bit
            // TODO: This accounts for one invariant case (empty), make it account for all of them
            reflected = [
                (
                    if reflected[1].1 != u32::MAX {
                        Reflection::flip_x(reflected[1].0)
                    } else {
                        reflected[1].0
                    },
                    reflected[1].1,
                ),
                (
                    if reflected[0].1 != u32::MAX {
                        Reflection::flip_x(reflected[0].0)
                    } else {
                        reflected[0].0
                    },
                    reflected[0].1,
                ),
                (
                    if reflected[3].1 != u32::MAX {
                        Reflection::flip_x(reflected[3].0)
                    } else {
                        reflected[3].0
                    },
                    reflected[3].1,
                ),
                (
                    if reflected[2].1 != u32::MAX {
                        Reflection::flip_x(reflected[2].0)
                    } else {
                        reflected[2].0
                    },
                    reflected[2].1,
                ),
                (
                    if reflected[5].1 != u32::MAX {
                        Reflection::flip_x(reflected[5].0)
                    } else {
                        reflected[5].0
                    },
                    reflected[5].1,
                ),
                (
                    if reflected[4].1 != u32::MAX {
                        Reflection::flip_x(reflected[4].0)
                    } else {
                        reflected[4].0
                    },
                    reflected[4].1,
                ),
                (
                    if reflected[7].1 != u32::MAX {
                        Reflection::flip_x(reflected[7].0)
                    } else {
                        reflected[7].0
                    },
                    reflected[7].1,
                ),
                (
                    if reflected[6].1 != u32::MAX {
                        Reflection::flip_x(reflected[6].0)
                    } else {
                        reflected[6].0
                    },
                    reflected[6].1,
                ),
            ];
        }
        if reflection.y {
            // Exchange pairs of (index, index+2) and flip the y reflection bit
            reflected = [
                (
                    if reflected[2].1 != u32::MAX {
                        Reflection::flip_y(reflected[2].0)
                    } else {
                        reflected[2].0
                    },
                    reflected[2].1,
                ),
                (
                    if reflected[3].1 != u32::MAX {
                        Reflection::flip_y(reflected[3].0)
                    } else {
                        reflected[3].0
                    },
                    reflected[3].1,
                ),
                (
                    if reflected[0].1 != u32::MAX {
                        Reflection::flip_y(reflected[0].0)
                    } else {
                        reflected[0].0
                    },
                    reflected[0].1,
                ),
                (
                    if reflected[1].1 != u32::MAX {
                        Reflection::flip_y(reflected[1].0)
                    } else {
                        reflected[1].0
                    },
                    reflected[1].1,
                ),
                (
                    if reflected[6].1 != u32::MAX {
                        Reflection::flip_y(reflected[6].0)
                    } else {
                        reflected[6].0
                    },
                    reflected[6].1,
                ),
                (
                    if reflected[7].1 != u32::MAX {
                        Reflection::flip_y(reflected[7].0)
                    } else {
                        reflected[7].0
                    },
                    reflected[7].1,
                ),
                (
                    if reflected[4].1 != u32::MAX {
                        Reflection::flip_y(reflected[4].0)
                    } else {
                        reflected[4].0
                    },
                    reflected[4].1,
                ),
                (
                    if reflected[5].1 != u32::MAX {
                        Reflection::flip_y(reflected[5].0)
                    } else {
                        reflected[5].0
                    },
                    reflected[5].1,
                ),
            ];
        }
        if reflection.z {
            // Exchange pairs of (index, index+4) and flip the z reflection bit
            reflected = [
                (
                    if reflected[4].1 != u32::MAX {
                        Reflection::flip_z(reflected[4].0)
                    } else {
                        reflected[4].0
                    },
                    reflected[4].1,
                ),
                (
                    if reflected[5].1 != u32::MAX {
                        Reflection::flip_z(reflected[5].0)
                    } else {
                        reflected[5].0
                    },
                    reflected[5].1,
                ),
                (
                    if reflected[6].1 != u32::MAX {
                        Reflection::flip_z(reflected[6].0)
                    } else {
                        reflected[6].0
                    },
                    reflected[6].1,
                ),
                (
                    if reflected[7].1 != u32::MAX {
                        Reflection::flip_z(reflected[7].0)
                    } else {
                        reflected[7].0
                    },
                    reflected[7].1,
                ),
                (
                    if reflected[0].1 != u32::MAX {
                        Reflection::flip_z(reflected[0].0)
                    } else {
                        reflected[0].0
                    },
                    reflected[0].1,
                ),
                (
                    if reflected[1].1 != u32::MAX {
                        Reflection::flip_z(reflected[1].0)
                    } else {
                        reflected[1].0
                    },
                    reflected[1].1,
                ),
                (
                    if reflected[2].1 != u32::MAX {
                        Reflection::flip_z(reflected[2].0)
                    } else {
                        reflected[2].0
                    },
                    reflected[2].1,
                ),
                (
                    if reflected[3].1 != u32::MAX {
                        Reflection::flip_z(reflected[3].0)
                    } else {
                        reflected[3].0
                    },
                    reflected[3].1,
                ),
            ];
        }
        reflected
    }
    fn get_pos(&self, mut coords: Coordinate) -> bool {
        if coords.y >= 1024 {
            return false;
        } else if coords.y < -1024 {
            return true;
        } else if (coords.x < -1024 || coords.x >= 1024) || (coords.z < -1024 || coords.z >= 1024) {
            return coords.y < 0;
        }
        let mut range = 2i16.pow(10);
        let mut offset;
        let mut reflection = 0;
        let mut next_index = 0;

        // Safe because these layers are constructed based on the bounds of the lower layer, so bounds checking is redundant
        // This buys a 16% speed boost, which will be very valuable in raycasting

        for layer in (1..=9).rev() {
            range /= 2;
            let next_offset_coords = Self::coords_to_offset(coords, range);
            offset = next_offset_coords.0;
            coords = next_offset_coords.1;
            unsafe {
                let next_data = match layer {
                    9 => self.layer9.get_unchecked(offset),
                    8 => self.layer8.get_unchecked(8 * next_index + offset),
                    7 => self.layer7.get_unchecked(8 * next_index + offset),
                    6 => self.layer6.get_unchecked(8 * next_index + offset),
                    5 => self.layer5.get_unchecked(8 * next_index + offset),
                    4 => self.layer4.get_unchecked(8 * next_index + offset),
                    3 => self.layer3.get_unchecked(8 * next_index + offset),
                    2 => self.layer2.get_unchecked(8 * next_index + offset),
                    1 => self.layer1.get_unchecked(8 * next_index + offset),
                    _ => panic!("Trying to index into nonexistant layer."),
                };

                reflection = next_data.0;
                next_index = next_data.1 as usize;
                if next_index == u32::MAX as usize {
                    return false;
                }
                if reflection != 0 {
                    coords = Reflection::from_u8(reflection).reflect_coords(coords);
                }
            }
        }

        // Data:
        let data = self.data[next_index];
        range /= 2;

        let (offset, coords) = Self::coords_to_offset(coords, range);
        range /= 2;
        let (sub_offset, _coords) = Self::coords_to_offset(coords, range);

        data & (1 << (63 - (offset * 8 + sub_offset))) != 0
    }
    fn coords_to_offset(coords: Coordinate, range: i16) -> (usize, Coordinate) {
        // Indices in layer N point to 8 elements in layer N-1.
        // This function returns a number 0-7 that corresponds to the given coordinates
        // and the coordinates updated to their relative value within layer N-1.
        // Range is the length of a side of the volume represented by layerN.
        let (mut x, mut y, mut z) = (coords.x, coords.y, coords.z);
        let mut offset = 0;
        if x >= 0 {
            offset += 1;
            x -= range;
        } else {
            x += range;
        };
        if y >= 0 {
            offset += 2;
            y -= range;
        } else {
            y += range;
        };
        if z >= 0 {
            offset += 4;
            z -= range;
        } else {
            z += range;
        };
        (offset, Coordinate { x, y, z })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torus() {
        let mut world = UChunkDag::new();
        let vox_data_str = fs::read_to_string("assets/torus.json").unwrap();
        let vox_data: Box<VoxData> = Box::new(serde_json::from_str(&vox_data_str).unwrap());
        for voxel in vox_data.voxels.iter() {
            let pos_coords = Coordinate {
                x: voxel.x.parse::<i16>().unwrap(),
                y: voxel.y.parse::<i16>().unwrap(),
                z: voxel.z.parse::<i16>().unwrap(),
            };
            let neg_coords = Coordinate {
                x: -voxel.x.parse::<i16>().unwrap(),
                y: -voxel.y.parse::<i16>().unwrap(),
                z: -voxel.z.parse::<i16>().unwrap(),
            };
            world.set_pos(pos_coords);
            world.set_pos(neg_coords);
        }
        let mut compressed = ChunkDAG::new(world);
        for n in 1..=8 {
            compressed.compress(n);
        }
        for voxel in vox_data.voxels.iter() {
            let pos_coords = Coordinate {
                x: voxel.x.parse::<i16>().unwrap(),
                y: voxel.y.parse::<i16>().unwrap(),
                z: voxel.z.parse::<i16>().unwrap(),
            };
            let neg_coords = Coordinate {
                x: -voxel.x.parse::<i16>().unwrap(),
                y: -voxel.y.parse::<i16>().unwrap(),
                z: -voxel.z.parse::<i16>().unwrap(),
            };
            assert!(
                compressed.get_pos(pos_coords),
                "Failure at {:?}",
                pos_coords
            );
            assert!(
                compressed.get_pos(neg_coords),
                "Failure at {:?}",
                neg_coords
            );
        }
    }
    #[test]
    fn reflect() {
        assert!(
            Reflection {
                x: true,
                y: true,
                z: true
            }
            .to_u8()
                == 7
        );
        assert!(
            Reflection::from_u8(5)
                == Reflection {
                    x: true,
                    y: false,
                    z: true
                }
        );
        assert!(
            Reflection {
                x: true,
                y: true,
                z: false
            }
            .reflect_coords(Coordinate {
                x: i16::MAX,
                y: i16::MIN,
                z: 0
            }) == Coordinate {
                x: i16::MIN,
                y: i16::MAX,
                z: 0
            }
        );
        assert!(
            ChunkDAG::least_symmetrical_data(0)
                == (
                    0,
                    Reflection {
                        x: false,
                        y: false,
                        z: false
                    }
                )
        );
        assert!(
            ChunkDAG::least_symmetrical_data(u64::MAX)
                == (
                    u64::MAX,
                    Reflection {
                        x: false,
                        y: false,
                        z: false
                    }
                )
        );
    }
    extern crate test;
    use test::Bencher;
    #[bench]
    fn torus_read_bench(b: &mut Bencher) {
        let mut world = UChunkDag::new();
        let vox_data_str = fs::read_to_string("assets/torus.json").unwrap();
        let vox_data: Box<VoxData> = Box::new(serde_json::from_str(&vox_data_str).unwrap());
        for voxel in vox_data.voxels.iter() {
            let pos_coords = Coordinate {
                x: voxel.x.parse::<i16>().unwrap(),
                y: voxel.y.parse::<i16>().unwrap(),
                z: voxel.z.parse::<i16>().unwrap(),
            };
            let neg_coords = Coordinate {
                x: -voxel.x.parse::<i16>().unwrap(),
                y: -voxel.y.parse::<i16>().unwrap(),
                z: -voxel.z.parse::<i16>().unwrap(),
            };
            world.set_pos(pos_coords);
            world.set_pos(neg_coords);
        }
        let mut compressed = ChunkDAG::new(world);
        for n in 1..=8 {
            compressed.compress(n);
        }

        b.iter(|| {
            for x in -128..128 {
                for y in -128..128 {
                    for z in -128..128 {
                        compressed.get_pos(Coordinate {
                            x: 8 * x,
                            y: 8 * y,
                            z: 8 * z,
                        });
                    }
                }
            }
        });
    }
}

// 10.234
// 10.257 ns
// 10.279

// Coordinate { x: 519, y: 144, z: 0 }
// Coordinate { x: 1, y: -2, z: -2 }
