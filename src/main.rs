#![feature(test)]
use serde::{Deserialize, Serialize};
use std::fs;
pub mod compressor;
pub mod constructor;
pub mod layer;
pub mod node;
pub mod reflection;
pub mod rotation;
pub mod transform;
pub mod utility;
use utility::Coordinate;

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
    let mut world = constructor::UChunkDAG::new();
    let vox_data_str = fs::read_to_string("assets/torus.json").unwrap();
    let vox_data: VoxData = serde_json::from_str(&vox_data_str).unwrap();
    let mut count = 0;
    for voxel in vox_data.voxels.iter() {
        let (x, y, z) = (
            voxel.x.parse::<i16>().unwrap(),
            voxel.y.parse::<i16>().unwrap(),
            voxel.z.parse::<i16>().unwrap(),
        );
        let pos_coords = Coordinate { x, y, z };
        let neg_coords = Coordinate { x, y: -y, z };
        let higher_coords = Coordinate { x, y: y * 2, z };
        let lower_coords = Coordinate { x, y: -y * 2, z };
        let sideways_coords = Coordinate { x, y: z, z: y };
        let neg_sideways_coords = Coordinate { x, y: z, z: -y };
        world.set_pos(lower_coords);
        world.set_pos(higher_coords);
        world.set_pos(pos_coords);
        world.set_pos(neg_coords);
        world.set_pos(sideways_coords);
        world.set_pos(neg_sideways_coords);
        count += 6;
    }
    println!("{} voxels", count);

    let compressed = compressor::ChunkDAG::new(world);

    println!("(Theoretical) memory usage: {}", compressed.size());
    let output_data = bincode::serialize(&compressed).unwrap();
    fs::write("out/compressed", output_data).expect("Unable to write file");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torus() {
        let mut world = constructor::UChunkDAG::new();
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
        let compressed = compressor::ChunkDAG::new(world);
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
            reflection::Reflection {
                x: true,
                y: true,
                z: true
            }
            .to_u8()
                == 7
        );
        assert!(
            reflection::Reflection::from_u8(5)
                == reflection::Reflection {
                    x: true,
                    y: false,
                    z: true
                }
        );
        assert!(
            reflection::Reflection {
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
    }
    extern crate test;
    use test::Bencher;
    #[bench]
    fn torus_read_bench(b: &mut Bencher) {
        let mut world = constructor::UChunkDAG::new();
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
        let compressed = compressor::ChunkDAG::new(world);

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
