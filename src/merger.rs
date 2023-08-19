use crate::compressor::ChunkDAG;

pub struct WorldDAG {
    // Combines 32 ChunkDAGS
    layer11: Vec<(u8, u32)>,
    layer10: Vec<(u8, u32)>,
    layer9: Vec<(u8, u32)>,
    layer8: Vec<(u8, u32)>,
    layer7: Vec<(u8, u32)>,
    layer6: Vec<(u8, u32)>,
    layer5: Vec<(u8, u32)>,
    layer4: Vec<(u8, u32)>,
    layer3: Vec<(u8, u32)>,
    layer2: Vec<(u8, u32)>,
    layer1: Vec<(u8, u32)>,
    data: Vec<u64>,
}

impl WorldDAG {
    pub fn new(chunks: [ChunkDAG; 8]) -> Self {
        /// [(-x, -y, -z), (+x, -y, -z), (-x, +y, -z), (+x, +y, -z),
        ///  (-x, -y, +z), (+x, -y, +z), (-x, +y, +z), (+x, +y, +z)]
        let mut data = Vec::new(); //chunks[0].data;
        let mut layer1 = Vec::new(); //chunks[0].layer1;
        let mut layer2 = Vec::new(); //chunks[0].layer2;
        let mut layer3 = Vec::new(); //chunks[0].layer3;
        let mut layer4 = Vec::new(); //chunks[0].layer4;
        let mut layer5 = Vec::new(); //chunks[0].layer5;
        let mut layer6 = Vec::new(); //chunks[0].layer6;
        let mut layer7 = Vec::new(); //chunks[0].layer7;
        let mut layer8 = Vec::new(); //chunks[0].layer8;
        let mut layer9 = Vec::new(); //chunks[0].layer9;

        for chunk in chunks.iter() {
            let data_offset = data.len();
            data.extend(chunk.data.clone());

            let layer1_offset = layer1.len();
            layer1.extend(
                chunk.layers[0]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + data_offset as u32)),
            );

            let layer2_offset = layer2.len();
            layer2.extend(
                chunk.layers[1]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer1_offset as u32)),
            );

            let layer3_offset = layer3.len();
            layer3.extend(
                chunk.layers[2]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer2_offset as u32)),
            );

            let layer4_offset = layer4.len();
            layer4.extend(
                chunk.layers[3]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer3_offset as u32)),
            );

            let layer5_offset = layer5.len();
            layer5.extend(
                chunk.layers[4]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer4_offset as u32)),
            );

            let layer6_offset = layer6.len();
            layer6.extend(
                chunk.layers[5]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer5_offset as u32)),
            );

            let layer7_offset = layer7.len();
            layer7.extend(
                chunk.layers[6]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer6_offset as u32)),
            );

            let layer8_offset = layer8.len();
            layer8.extend(
                chunk.layers[7]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer7_offset as u32)),
            );

            layer9.extend(
                chunk.layers[8]
                    .clone()
                    .into_iter()
                    .map(|(y, x)| (y, x + layer8_offset as u32)),
            );
        }
        Self {
            layer11: Vec::new(),
            layer10: Vec::new(),
            layer9,
            layer8,
            layer7,
            layer6,
            layer5,
            layer4,
            layer3,
            layer2,
            layer1,
            data,
        }
    }
}
