use crate::latest::binary_matrix::BinaryMatrix;
use std::simd::{LaneCount, Simd, SupportedLaneCount};

#[derive(Clone, Debug)]
pub struct BinaryMatrixSimd<const LANES: usize>
where
    LaneCount<LANES>: SupportedLaneCount,
{
    rows: usize,
    columns: Vec<Vec<Simd<u64, LANES>>>,
}

impl<const LANES: usize> BinaryMatrixSimd<LANES>
where
    LaneCount<LANES>: SupportedLaneCount,
{
    pub fn new() -> Box<dyn BinaryMatrix> {
        Box::new(BinaryMatrixSimd {
            rows: 0,
            columns: vec![],
        })
    }

    #[inline(always)]
    const fn simd_bits(&self) -> usize {
        (LANES.trailing_zeros() + u64::BITS.trailing_zeros()) as usize
    }

    #[inline(always)]
    const fn zero_simd(&self) -> Simd<u64, LANES> {
        Simd::from_array([0; LANES])
    }

    #[inline(always)]
    const fn simd_lane_mask(&self) -> usize {
        LANES - 1
    }
    #[inline(always)]
    const fn simd_base_mask(&self) -> usize {
        (u64::BITS - 1) as usize
    }
}

impl<const LANES: usize> BinaryMatrix for BinaryMatrixSimd<LANES>
where
    LaneCount<LANES>: SupportedLaneCount,
{
    fn nrows(&self) -> usize {
        self.rows
    }
    fn ncols(&self) -> usize {
        self.columns.len()
    }
    fn insert_rows(&mut self, new_rows: usize, new_cols: usize) {
        let simd_bits = self.simd_bits();
        let zero_simd = self.zero_simd();
        self.rows += new_rows;
        for c in 0..self.columns.len() {
            self.columns[c].resize((self.rows >> simd_bits) + 1, zero_simd);
        }
        for _ in 0..new_cols {
            let mut col = Vec::with_capacity((self.rows >> simd_bits) + 1);
            col.resize((self.rows >> simd_bits) + 1, zero_simd);
            self.columns.push(col);
        }
    }

    // TODO: bit tilt algorithm?
    fn transpose(&self) -> Box<dyn BinaryMatrix> {
        let mut new = BinaryMatrixSimd::new();
        new.insert_rows(self.ncols(), self.nrows());
        for c in 0..self.columns.len() {
            for r in 0..self.rows {
                new.set(c, r, self.get(r, c));
            }
        }
        new
    }

    fn get(&self, r: usize, c: usize) -> u8 {
        assert!(c < self.columns.len());
        assert!(r < self.rows);
        let x = self.columns[c][r >> self.simd_bits()].as_array()
            [(r >> u64::BITS.trailing_zeros()) & self.simd_lane_mask()];
        let shift = r & self.simd_base_mask();
        ((x >> shift) & 1) as u8
    }

    fn set(&mut self, r: usize, c: usize, val: u8) {
        assert!(c < self.columns.len());
        assert!(r < self.rows);
        let shift = r & self.simd_base_mask();
        let row_idx = r >> self.simd_bits();
        let lane_idx = (r >> u64::BITS.trailing_zeros()) & self.simd_lane_mask();
        if val == 1 {
            self.columns[c][row_idx].as_mut_array()[lane_idx] |= 1 << shift;
        } else {
            self.columns[c][row_idx].as_mut_array()[lane_idx] &= !(1 << shift);
        }
    }

    fn swap_columns(&mut self, c1: usize, c2: usize) -> () {
        self.columns.swap(c1, c2);
    }

    fn xor_col(&mut self, c1: usize, c2: usize) -> () {
        let maxc = self.columns[c1].len();
        unsafe {
            let mut x1 = self.columns[c1].as_mut_ptr();
            let mut x2 = self.columns[c2].as_ptr();
            for _ in 0..maxc {
                //self.columns[c1][i] ^= self.columns[c2][i];
                *x1 ^= *x2;
                x1 = x1.offset(1);
                x2 = x2.offset(1);
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_transpose() {
        let mut mat = BinaryMatrixSimd::<16>::new();
        mat.insert_rows(2, 3);
        mat.transpose();
    }
}
