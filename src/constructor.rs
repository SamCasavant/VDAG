use crate::utility::Coordinate;

pub struct UChunkDAG {
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
    pub data: Vec<u64>, //Box<[u64; 134217728]>,
}
impl UChunkDAG {
    pub fn new() -> Self {
        let data = vec![0; 134_217_728];
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
