extern crate nalgebra as na;
use crate::latest::binary_matrix::BinaryMatrix;
use na::DMatrix;

impl BinaryMatrix for DMatrix<u8> {
    fn nrows(&self) -> usize {
        DMatrix::nrows(&self)
    }

    fn ncols(&self) -> usize {
        DMatrix::ncols(&self)
    }

    fn insert_rows(&mut self, new_rows: usize, new_cols: usize) {
        self.resize_mut(self.nrows() + new_rows, self.ncols() + new_cols, 0u8);
    }

    fn transpose(&self) -> Box<dyn BinaryMatrix> {
        let t = DMatrix::transpose(&self);
        Box::new(t)
    }

    fn get(&self, r: usize, c: usize) -> u8 {
        self[(r, c)]
    }

    fn set(&mut self, r: usize, c: usize, val: u8) {
        self[(r, c)] = val;
    }

    fn swap_columns(&mut self, c1: usize, c2: usize) -> () {
        DMatrix::swap_columns(self, c1, c2)
    }

    fn xor_col(&mut self, c1: usize, c2: usize) -> () {
        for r in 0..self.nrows() {
            self[(r, c1)] ^= self[(r, c2)];
        }
    }
}
