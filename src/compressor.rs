use crate::constructor::UChunkDAG;
use crate::layer::Layer;
use crate::node::Node;
use crate::reflection::Reflection;
use crate::utility::Coordinate;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct ChunkDAG {
    //TODO: Most of these don't need to be u32, but consistency made the code easier to generalize. Generics?
    //TODO: This is the same size as (u32, u32) and only 3 bits are used in u8, so I can fit 29 more bits/node.
    //       16 for material information, 13 for additional transforms?
    //       Or with u16 indices, there are 13 remaining bits.
    pub layers: Vec<Layer<u8, u32>>,
    pub data: Vec<u64>,
}
impl ChunkDAG {
    pub fn new(mut uncompressed: UChunkDAG) -> Self {
        let mut leaves = Vec::with_capacity(134_217_728);
        while let Some(leaf) = uncompressed.data.pop() {
            let (reflected, reflection) = Self::least_symmetrical_data(leaf);
            leaves.push((
                uncompressed.data.len() as u32,
                reflected,
                reflection.to_u8(),
            ));
        }

        // Free memory formerly used by uncompressed data
        uncompressed.data.shrink_to_fit();

        // Sort in reverse order so we can pop elements without shifting
        leaves.sort_unstable_by(|a, b| {
            if a.1 == b.1 {
                b.0.cmp(&a.0)
            } else {
                b.1.cmp(&a.1)
            }
        });
        let mut layer1 = Layer(vec![Node::new(0, u32::MAX); 134_217_728]);
        // Data without duplicates
        let mut data = Vec::with_capacity(128);

        while let Some(active_leaf) = leaves.pop() {
            let active_index = match active_leaf.1 {
                0 => u32::MAX,
                val => {
                    data.push(val);
                    (data.len() - 1) as u32
                }
            };
            if active_index != u32::MAX {
                layer1[active_leaf.0 as usize] = Node::new(active_leaf.2, active_index);
            }
            while !leaves.is_empty() && leaves[leaves.len() - 1].1 == active_leaf.1 {
                let removed = leaves.pop().unwrap();
                if active_index != u32::MAX {
                    layer1[removed.0 as usize] = Node::new(removed.2, active_index);
                }
            }
        }
        // Data will not grow until a merge, so free extra allocated space
        data.shrink_to_fit();

        let layers = vec![layer1];

        let mut result = Self { layers, data };

        // Compress remaining layers
        for _layer in 1..=8 {
            result.compress();
        }
        result
    }
    fn compress(&mut self) {
        // Compresses the highest layer in the DAG
        let layer = self.layers.len();
        // Identify symmetrical regions in layer[n - 1].
        // These inform which transformations can be ignored when checking for symmetries during
        //      the construction of layer[n + 1]
        // let mut active_symmetries = HashSet::new();
        // if layer == 1 {
        //     for (index, elem) in self.data.iter().enumerate() {
        //         let symmetries = Self::get_symmetries_data(*elem);
        //         if symmetries.x {
        //             active_symmetries.insert((index as u32, Axis::X));
        //         }
        //         if symmetries.y {
        //             active_symmetries.insert((index as u32, Axis::Y));
        //         }
        //         if symmetries.z {
        //             active_symmetries.insert((index as u32, Axis::Z));
        //         }
        //     }
        // } else {
        //     for (index, chunk) in self.layers[layer - 2].cubes_iter() {
        //         let symmetries =
        //             Self::get_symmetries_nodes(chunk.try_into().unwrap(), &active_symmetries);
        //         if symmetries.x {
        //             active_symmetries.insert((index as u32, Axis::X));
        //         }
        //         if symmetries.y {
        //             active_symmetries.insert((index as u32, Axis::Y));
        //         }
        //         if symmetries.z {
        //             active_symmetries.insert((index as u32, Axis::Z));
        //         }
        //     }
        // }
        // active_symmetries.insert((u32::MAX, Axis::X));
        // active_symmetries.insert((u32::MAX, Axis::Y));
        // active_symmetries.insert((u32::MAX, Axis::Z));

        let lower = self.layers.get_mut(layer - 1).unwrap();

        if lower.len() % 64 != 0 {
            panic!("Attempted to compress incomplete layer.");
        }
        let mut lower_data = Vec::with_capacity(lower.len() / 8);
        for (i, chunk) in lower.cubes_iter() {
            let (reflected, reflection) = Self::least_symmetrical_nodes(chunk.try_into().unwrap()); //, &active_symmetries);
            lower_data.push((i, reflected, reflection.to_u8()));
        }
        let upper_length = lower.len() / 8;
        lower.free();

        // Sort in reverse order so we can pop elements without shifting
        lower_data.sort_unstable_by(|a, b| {
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

        let mut upper_contents = Layer(vec![Node::new(0, 0); upper_length]);
        let mut index = 0_u32;
        println!("Compressing...");
        while let Some(active_node) = lower_data.pop() {
            let empty = [Node::new(0, u32::MAX); 8];
            let active_index = match active_node.1 {
                null if null == empty => u32::MAX,
                val => {
                    for element in val {
                        lower.push(element);
                    }
                    index += 1;
                    index - 1
                }
            };
            upper_contents[active_node.0] = Node::new(active_node.2, active_index);
            while !lower_data.is_empty()
                && (
                    lower_data[lower_data.len() - 1].1,
                    lower_data[lower_data.len() - 1].2,
                ) == (active_node.1, active_node.2)
            {
                // Delete duplicate elements from lower_data and shift upper_context indices
                let removed = lower_data.pop().unwrap();
                upper_contents[removed.0] = Node::new(removed.2, active_index);
            }
        }
        lower.shrink_to_fit();

        self.layers.push(upper_contents);
    }
    fn get_symmetries_data(data: u64) -> Reflection {
        // Takes a 4x4x4 block of voxel data and returns the axes across which it is symmetrical
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        if Self::reflect_data(
            data,
            Reflection {
                x: true,
                y: false,
                z: false,
            },
        ) == data
        {
            reflection.x = true;
        }
        if Self::reflect_data(
            data,
            Reflection {
                x: false,
                y: true,
                z: false,
            },
        ) == data
        {
            reflection.y = true;
        }
        if Self::reflect_data(
            data,
            Reflection {
                x: false,
                y: false,
                z: true,
            },
        ) == data
        {
            reflection.z = true;
        }
        reflection
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
        for axis in 0b001..=0b111 {
            let test = Self::reflect_data(data, Reflection::from_u8(axis));
            if test < min {
                min = test;
                reflection = Reflection::from_u8(axis);
            }
        }
        (min, reflection)
    }
    fn reflect_data(data: u64, reflection: Reflection) -> u64 {
        // Reflects a 4x4x4 block across a list of axes
        let mut reflected = data;
        if reflection.x {
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
    fn get_symmetries_nodes(
        nodes: [Node<u8, u32>; 8],
        // active_symmetries: &HashSet<(u32, Axis)>,
    ) -> Reflection {
        // Takes a 4x4x4 block of voxel data and returns the axes along which a reflection is identical
        // Note to self: In my mind, each axis can be treated independently, but I don't have hard proof of this
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        if Self::reflect_nodes(
            nodes,
            Reflection {
                x: true,
                y: false,
                z: false,
            }, // active_symmetries,
        ) == nodes
        {
            reflection.x = true;
        }
        if Self::reflect_nodes(
            nodes,
            Reflection {
                x: false,
                y: true,
                z: false,
            },
            //active_symmetries,
        ) == nodes
        {
            reflection.y = true;
        }
        if Self::reflect_nodes(
            nodes,
            Reflection {
                x: false,
                y: false,
                z: true,
            },
            //active_symmetries,
        ) == nodes
        {
            reflection.z = true;
        }
        reflection
    }
    fn least_symmetrical_nodes(
        nodes: [Node<u8, u32>; 8],
        // active_symmetries: &HashSet<(u32, Axis)>,
    ) -> ([Node<u8, u32>; 8], Reflection) {
        // Takes a node's children and returns the lowest (according to however Rust orders slices) similar children by reflection
        let mut min = nodes;
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        for axis in 0b001..=0b111 {
            let test = Self::reflect_nodes(nodes, Reflection::from_u8(axis)); //, active_symmetries);
            if test < min {
                min = test;
                reflection = Reflection::from_u8(axis);
            }
        }

        (min, reflection)
    }
    fn reflect_nodes(
        nodes: [Node<u8, u32>; 8],
        reflection: Reflection,
        // active_symmetries: &HashSet<(u32, Axis)>,
    ) -> [Node<u8, u32>; 8] {
        // Children are ordered like this:
        // (x, y, z) =>
        // [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        let mut reflected = nodes;
        if reflection.x {
            // Exchange every pair of values and flip the x reflection bit
            for index in [0, 2, 4, 6] {
                reflected.swap(index, index + 1);
                // if !active_symmetries.contains(&(reflected[index].next, Axis::X)) {
                reflected[index].val = Reflection::flip_x(reflected[index].val);
                // };
                // if !active_symmetries.contains(&(reflected[index + 1].next, Axis::X)) {
                reflected[index + 1].val = Reflection::flip_x(reflected[index + 1].val);
                // };
            }
        }
        if reflection.y {
            // Exchange pairs of (index, index+2) and flip the y reflection bit
            for index in [0, 1, 4, 5] {
                reflected.swap(index, index + 2);
                // if !active_symmetries.contains(&(reflected[index].next, Axis::Y)) {
                reflected[index].val = Reflection::flip_y(reflected[index].val);
                // };
                // if !active_symmetries.contains(&(reflected[index + 2].next, Axis::Y)) {
                reflected[index + 2].val = Reflection::flip_y(reflected[index + 2].val);
                // };
            }
        }
        if reflection.z {
            // Exchange pairs of (index, index+4) and flip the z reflection bit
            for index in [0, 1, 2, 3] {
                reflected.swap(index, index + 4);
                // if !active_symmetries.contains(&(reflected[index].next, Axis::Z)) {
                reflected[index].val = Reflection::flip_z(reflected[index].val);
                // };
                // if !active_symmetries.contains(&(reflected[index + 4].next, Axis::Z)) {
                reflected[index + 4].val = Reflection::flip_z(reflected[index + 4].val);
                // };
            }
        }
        reflected
    }
    pub fn get_pos(&self, mut coords: Coordinate) -> bool {
        if coords.y >= 1024 {
            return false;
        } else if coords.y < -1024 {
            return true;
        } else if (coords.x < -1024 || coords.x >= 1024) || (coords.z < -1024 || coords.z >= 1024) {
            return coords.y < 0;
        }
        let mut range = 2i16.pow(10);
        let mut offset;
        let mut reflection;
        let mut next_index = 0;

        for layer in (1..=9).rev() {
            range /= 2;
            let next_offset_coords = Self::coords_to_offset(coords, range);
            offset = next_offset_coords.0;
            coords = next_offset_coords.1;
            // SAFETY: because these layers are constructed based on the bounds of the lower layer, bounds checking is redundant
            // This buys a 16% speed boost, which will be very valuable in raycasting
            unsafe {
                let next_data = self.layers[layer - 1].get_unchecked(8 * next_index + offset);

                reflection = next_data.val;
                next_index = next_data.next as usize;
                if next_index == u32::MAX as usize {
                    return false;
                }
                if reflection != 0 {
                    coords = Reflection::from_u8(reflection).reflect_coords(coords);
                }
            }
        }

        let data = self.data[next_index];
        range /= 2;

        (offset, coords) = Self::coords_to_offset(coords, range);
        range /= 2;
        let (sub_offset, _coords) = Self::coords_to_offset(coords, range);

        data & (1 << (63 - (offset * 8 + sub_offset))) != 0
    }
    const fn coords_to_offset(coords: Coordinate, range: i16) -> (usize, Coordinate) {
        // Indices in layer N point to 8 elements in layer N-1.
        // This function returns an index into those elements that corresponds to the given coordinates
        // and the coordinates updated to their relative value within layer N-1.
        // Range is the length of a side of the volume represented by layerN
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
    pub fn size(&self) -> usize {
        println!("Data length: {}", self.data.len());
        self.data.capacity() * 8
            + self.layers[0].capacity() * 8
            + self.layers[1].capacity() * 8
            + self.layers[2].capacity() * 8
            + self.layers[3].capacity() * 8
            + self.layers[4].capacity() * 8
            + self.layers[5].capacity() * 8
            + self.layers[6].capacity() * 8
            + self.layers[7].capacity() * 8
            + self.layers[8].capacity() * 8
    }
}
