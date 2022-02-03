#![feature(test)]
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, fs};

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
    // world.set_pos(Coordinate { x: 1, y: 1, z: 2 });
    // world.set_pos(Coordinate { x: 1, y: 2, z: 1 });
    // world.set_pos(Coordinate { x: 1, y: 2, z: 2 });
    // world.set_pos(Coordinate { x: 2, y: 1, z: 1 });
    // world.set_pos(Coordinate { x: 2, y: 1, z: 2 });
    // world.set_pos(Coordinate { x: 2, y: 2, z: 1 });
    // world.set_pos(Coordinate { x: 2, y: 2, z: 2 });

    let vox_data_str = fs::read_to_string("assets/torus.json").unwrap();
    let vox_data: VoxData = serde_json::from_str(&vox_data_str).unwrap();
    let mut count = 0;
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
        let higher_coords = Coordinate {
            x: voxel.x.parse::<i16>().unwrap(),
            y: voxel.y.parse::<i16>().unwrap() * 2,
            z: voxel.z.parse::<i16>().unwrap(),
        };
        let lower_coords = Coordinate {
            x: voxel.x.parse::<i16>().unwrap(),
            y: -voxel.y.parse::<i16>().unwrap() * 2,
            z: voxel.z.parse::<i16>().unwrap(),
        };
        world.set_pos(lower_coords);
        world.set_pos(higher_coords);
        world.set_pos(pos_coords);
        world.set_pos(neg_coords);
        count += 4;
    }
    println!("{} voxels", count);

    let compressed = ChunkDAG::new(world);
    println!("Layer4: {:?}", compressed.layer4);

    //println!("{:?}", compressed);
    // for each in [
    //     Coordinate { x: 1, y: 1, z: 2 },
    //     Coordinate { x: 1, y: 2, z: 1 },
    //     Coordinate { x: 1, y: 2, z: 2 },
    //     Coordinate { x: 2, y: 1, z: 1 },
    //     Coordinate { x: 2, y: 1, z: 2 },
    //     Coordinate { x: 2, y: 2, z: 1 },
    //     Coordinate { x: 2, y: 2, z: 2 },
    // ] {
    //     println!("{}", compressed.get_pos(each));
    // }
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
        let higher_coords = Coordinate {
            x: voxel.x.parse::<i16>().unwrap(),
            y: voxel.y.parse::<i16>().unwrap() * 2,
            z: voxel.z.parse::<i16>().unwrap(),
        };
        let lower_coords = Coordinate {
            x: voxel.x.parse::<i16>().unwrap(),
            y: -voxel.y.parse::<i16>().unwrap() * 2,
            z: voxel.z.parse::<i16>().unwrap(),
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
        assert!(
            compressed.get_pos(higher_coords),
            "Failure at {:?}",
            pos_coords
        );
        assert!(
            compressed.get_pos(lower_coords),
            "Failure at {:?}",
            neg_coords
        );
    }
    for x in -128..128 {
        for y in -128..128 {
            for z in -128..128 {
                compressed.get_pos(Coordinate { x, y, z });
            }
        }
    }
    println!("(Theoretical) memory usage: {}", compressed.size());
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
    data: Vec<u64>, //Box<[u64; 134217728]>,
}
impl UChunkDag {
    pub fn new() -> Self {
        let data = vec![0; 134_217_728]; //Box::new([0; 134_217_728]);
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
        let (x, y, z) = (
            (coords.x + 1024) as u32,
            (coords.y + 1024) as u32,
            (coords.z + 1024) as u32,
        );

        // Split into first 9 bits / last 2 bits
        let (mut x_index, mut x_sub_index) = (x / 4, x % 4);
        let (mut y_index, mut y_sub_index) = (y / 4, y % 4);
        let (mut z_index, mut z_sub_index) = (z / 4, z % 4);

        // Morton(ish) encoding math
        for n in 0..=9 {
            x_index += 8u32.pow(n) * 6 * (x / 2u32.pow(n + 3));
            y_index += 8u32.pow(n) * 6 * (y / 2u32.pow(n + 3));
            z_index += 8u32.pow(n) * 6 * (z / 2u32.pow(n + 3));
        }
        x_sub_index = x_sub_index + 6 * (x_sub_index / 2);
        y_sub_index = y_sub_index + 6 * (y_sub_index / 2);
        z_sub_index = z_sub_index + 6 * (z_sub_index / 2);

        // Each increment in y advances the index by 2
        // Each increment in z advances the index by 4

        let index = (x_index + y_index * 2 + z_index * 4) as usize;

        let sub_index = (x_sub_index + y_sub_index * 2 + z_sub_index * 4) as u8;

        (index, sub_index)
    }
}
#[derive(Debug, PartialEq, Eq, Hash)]
enum Axis {
    X,
    Y,
    Z,
}
#[derive(Debug, PartialEq, Clone, Copy)]
struct Reflection {
    x: bool,
    y: bool,
    z: bool,
}
impl Reflection {
    const fn to_u8(self) -> u8 {
        self.x as u8 | ((self.y as u8) << 1) | ((self.z as u8) << 2)
    }
    const fn from_u8(reflection: u8) -> Self {
        Self {
            x: (reflection % 2) != 0,
            y: ((reflection / 2) % 2) != 0,
            z: ((reflection / 4) % 2) != 0,
        }
    }
    const fn reflect_coords(self, coords: Coordinate) -> Coordinate {
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
    const fn flip_x(reflection: u8) -> u8 {
        if reflection % 2 == 1 {
            reflection - 1
        } else {
            reflection + 1
        }
    }
    const fn flip_y(reflection: u8) -> u8 {
        if (reflection / 2) % 2 == 1 {
            reflection - 2
        } else {
            reflection + 2
        }
    }
    const fn flip_z(reflection: u8) -> u8 {
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
    active_symmetries: HashSet<(u32, Axis)>,
}
impl ChunkDAG {
    fn new(mut chunk: UChunkDag) -> Self {
        let mut leaves = Vec::with_capacity(134_217_728);
        while let Some(leaf) = chunk.data.pop() {
            let (reflected, reflection) = Self::least_symmetrical_data(leaf);
            leaves.push((chunk.data.len() as u32, reflected, reflection.to_u8()));
        }
        // Free memory formerly used by uncompressed data
        chunk.data.shrink_to_fit();

        // Sort in reverse order so we can pop elements without shifting
        leaves.sort_unstable_by(|a, b| {
            if a.1 == b.1 {
                b.0.cmp(&a.0)
            } else {
                b.1.cmp(&a.1)
            }
        });
        let mut layer1 = vec![(0, u32::MAX); 134_217_728];
        // Data without duplicates
        let mut data = Vec::with_capacity(128);

        println!("Compressing...");
        while let Some(active_leaf) = leaves.pop() {
            let active_index = match active_leaf.1 {
                0 => u32::MAX,
                val => {
                    data.push(val);
                    (data.len() - 1) as u32
                }
            };
            if active_index != u32::MAX {
                layer1[active_leaf.0 as usize] = (active_leaf.2, active_index);
            }
            while !leaves.is_empty() && leaves[leaves.len() - 1].1 == active_leaf.1 {
                let removed = leaves.pop().unwrap();
                if active_index != u32::MAX {
                    layer1[removed.0 as usize] = (removed.2, active_index);
                }
            }
        }
        // Data will not grow until a merge, so free extra allocated space
        data.shrink_to_fit();
        // Produce a set representing which of the indices into data are self-symmetrical
        let mut active_symmetries = HashSet::new();
        for (index, elem) in data.iter().enumerate() {
            let symmetries = Self::get_symmetries_data(*elem);
            if symmetries.x {
                active_symmetries.insert((index as u32, Axis::X));
            }
            if symmetries.y {
                active_symmetries.insert((index as u32, Axis::Y));
            }
            if symmetries.z {
                active_symmetries.insert((index as u32, Axis::Z));
            }
        }
        active_symmetries.insert((u32::MAX, Axis::X));
        active_symmetries.insert((u32::MAX, Axis::Y));
        active_symmetries.insert((u32::MAX, Axis::Z));
        println!("Precompression data symmetries: {:?}", active_symmetries);

        let mut result = Self {
            layer9: Vec::new(),
            layer8: Vec::new(),
            layer7: Vec::new(),
            layer6: Vec::new(),
            layer5: Vec::new(),
            layer4: Vec::new(),
            layer3: Vec::new(),
            layer2: Vec::new(),
            layer1,
            data,
            active_symmetries,
        };

        for layer in 1..=8 {
            result.compress(layer);
        }
        result.active_symmetries = HashSet::new();
        result
    }
    fn get_symmetries_data(data: u64) -> Reflection {
        // Takes a 4x4x4 block of voxel data and returns the axes along which a reflection is identical
        // Note to self: In my mind, each axis can be treated independently, but I don't have hard proof of this
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
    const fn swap_bits(mut data: u64, first: u8, second: u8) -> u64 {
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

        let mut lower_data = Vec::with_capacity(lower.len() / 8);
        for (i, chunk) in lower.chunks_exact(8).into_iter().enumerate() {
            let (reflected, reflection) =
                Self::least_symmetrical_nodes(chunk.try_into().unwrap(), &self.active_symmetries);
            lower_data.push((i, reflected, reflection.to_u8()));
        }
        // Temporarily free some memory
        lower.clear();
        lower.shrink_to_fit();

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
        let mut index = 0_u32;
        println!("Compressing...");
        // Dedup
        while let Some(active_node) = lower_data.pop() {
            let empty = [(0, u32::MAX); 8];
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
            upper_contents[active_node.0] = (active_node.2, active_index);
            while !lower_data.is_empty()
                && (
                    lower_data[lower_data.len() - 1].1,
                    lower_data[lower_data.len() - 1].2,
                ) == (active_node.1, active_node.2)
            {
                let removed = lower_data.pop().unwrap();
                upper_contents[removed.0] = (removed.2, active_index);
            }
        }
        lower.shrink_to_fit(); //Release extra capacity, vector will no longer grow until a merge

        let mut active_symmetries = HashSet::new();
        for (index, chunk) in lower.chunks_exact(8).into_iter().enumerate() {
            let symmetries =
                Self::get_symmetries_nodes(chunk.try_into().unwrap(), &self.active_symmetries);
            if symmetries.x {
                active_symmetries.insert((index as u32, Axis::X));
            }
            if symmetries.y {
                active_symmetries.insert((index as u32, Axis::Y));
            }
            if symmetries.z {
                active_symmetries.insert((index as u32, Axis::Z));
            }
        }
        println!("Node symmetries: {:?}", active_symmetries);
        self.active_symmetries = active_symmetries;
        *upper = upper_contents;
    }
    fn get_symmetries_nodes(
        nodes: [(u8, u32); 8],
        active_symmetries: &HashSet<(u32, Axis)>,
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
            },
            active_symmetries,
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
            active_symmetries,
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
            active_symmetries,
        ) == nodes
        {
            reflection.z = true;
        }
        reflection
    }
    fn least_symmetrical_nodes(
        nodes: [(u8, u32); 8],
        active_symmetries: &HashSet<(u32, Axis)>,
    ) -> ([(u8, u32); 8], Reflection) {
        // Takes a nodes children and returns the lowest (according to however Rust orders slices) similar children by reflection
        let mut min = nodes;
        let mut reflection = Reflection {
            x: false,
            y: false,
            z: false,
        };
        for axis in 0b001..=0b111 {
            let test = Self::reflect_nodes(nodes, Reflection::from_u8(axis), active_symmetries);
            if test < min {
                min = test;
                reflection = Reflection::from_u8(axis);
            }
        }

        (min, reflection)
    }
    fn reflect_nodes(
        nodes: [(u8, u32); 8],
        reflection: Reflection,
        active_symmetries: &HashSet<(u32, Axis)>,
    ) -> [(u8, u32); 8] {
        // Children are ordered like this:
        // (x, y, z) =>
        // [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        let mut reflected = nodes;
        if reflection.x {
            // Exchange every pair of values and flip the x reflection bit
            for index in [0, 2, 4, 6] {
                reflected.swap(index, index + 1);
                if !active_symmetries.contains(&(reflected[index].1, Axis::X)) {
                    reflected[index].0 = Reflection::flip_x(reflected[index].0);
                };
                if !active_symmetries.contains(&(reflected[index + 1].1, Axis::X)) {
                    reflected[index + 1].0 = Reflection::flip_x(reflected[index + 1].0);
                };
            }
        }
        if reflection.y {
            // Exchange pairs of (index, index+2) and flip the y reflection bit
            for index in [0, 1, 4, 5] {
                reflected.swap(index, index + 2);
                if !active_symmetries.contains(&(reflected[index].1, Axis::Y)) {
                    reflected[index].0 = Reflection::flip_y(reflected[index].0);
                };
                if !active_symmetries.contains(&(reflected[index + 2].1, Axis::Y)) {
                    reflected[index + 2].0 = Reflection::flip_y(reflected[index + 2].0);
                };
            }
        }
        if reflection.z {
            // Exchange pairs of (index, index+4) and flip the z reflection bit
            for index in [0, 1, 2, 3] {
                reflected.swap(index, index + 4);
                if !active_symmetries.contains(&(reflected[index].1, Axis::Z)) {
                    reflected[index].0 = Reflection::flip_z(reflected[index].0);
                };
                if !active_symmetries.contains(&(reflected[index + 4].1, Axis::Z)) {
                    reflected[index + 4].0 = Reflection::flip_z(reflected[index + 4].0);
                };
            }
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
    fn size(&self) -> usize {
        return self.data.capacity() * 8
            + self.layer1.capacity() * 8
            + self.layer2.capacity() * 8
            + self.layer3.capacity() * 8
            + self.layer4.capacity() * 8
            + self.layer5.capacity() * 8
            + self.layer6.capacity() * 8
            + self.layer7.capacity() * 8
            + self.layer8.capacity() * 8
            + self.layer9.capacity() * 8;
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
        let compressed = ChunkDAG::new(world);
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
        let compressed = ChunkDAG::new(world);

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
