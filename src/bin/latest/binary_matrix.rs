use core::ops;
use rand::{rngs::SmallRng, Rng, SeedableRng};

pub trait BinaryMatrix {
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    fn insert_rows(&mut self, new_rows: usize, new_cols: usize);
    fn transpose(&self) -> Box<dyn BinaryMatrix>;
    fn get(&self, r: usize, c: usize) -> u8;
    fn set(&mut self, r: usize, c: usize, val: u8);
    fn left_kernel(&self) -> Option<Vec<BinaryVector>> {
        // compute kernel of matrix over GF(2) using Gaussian elimination
        let mut extend = self.transpose();
        // extend out with an identity matrix of size RxR
        let maxr = extend.nrows();
        let maxc = extend.ncols();
        extend.insert_rows(maxc, 0);
        for c in 0..maxc {
            extend.set(maxr + c, c, 1);
        }

        // pivot on a row for each column
        let mut nextc = 0;
        for r in 0..maxr {
            let mut pivot = -1i32;
            for c in nextc..maxc {
                if extend.get(r, c) != 0 {
                    pivot = c as i32;
                    break;
                }
            }
            if pivot == -1 {
                continue;
            }

            extend.swap_columns(nextc, pivot as usize);

            for c in (nextc as usize + 1)..maxc {
                if extend.get(r, c) == 1 {
                    extend.xor_col(c, nextc);
                }
            }
            nextc += 1;
        }

        let mut kernel = vec![];

        for c in 0..maxc {
            if extend.column_part_all_zero(c, maxr) {
                let v = extend.extract_column_part(c, maxr, extend.nrows() - maxr);
                kernel.push(v);
            }
        }

        return Some(kernel);
    }
    fn left_mul(&self, result_vector: &BinaryVector) -> BinaryVector {
        assert_eq!(result_vector.bits.len(), self.nrows());
        let mut result = Vec::with_capacity(self.ncols());
        for c in 0..self.ncols() {
            let mut bit = 0;
            for r in 0..self.nrows() {
                bit ^= self.get(r, c) * result_vector.bits[r];
            }
            result.push(bit);
        }
        BinaryVector { bits: result }
    }
    fn swap_columns(&mut self, c1: usize, c2: usize) -> ();
    fn xor_col(&mut self, c1: usize, c2: usize) -> ();
    fn column_part_all_zero(&self, c: usize, maxr: usize) -> bool {
        for x in 0..maxr {
            if self.get(x, c) == 1 {
                return false;
            }
        }
        return true;
    }
    fn extract_column_part(&self, c: usize, maxr_1: usize, size: usize) -> BinaryVector {
        let mut bits: Vec<u8> = Vec::with_capacity(size);
        for r in maxr_1..maxr_1 + size {
            bits.push(self.get(r, c));
        }
        BinaryVector { bits }
    }
}

#[derive(Clone, Debug)]
pub struct BinaryMatrix64 {
    rows: usize,
    columns: Vec<Vec<u64>>,
}

#[derive(Clone, Debug)]
pub struct BinaryVector {
    pub bits: Vec<u8>,
}

impl BinaryVector {
    pub fn zeros(n: usize) -> BinaryVector {
        let mut v = Vec::with_capacity(n);
        v.resize(n, 0);
        BinaryVector { bits: v }
    }

    pub fn random(n: usize) -> BinaryVector {
        let mut rng = SmallRng::from_entropy();
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(rng.gen());
        }
        BinaryVector { bits: v }
    }

    pub fn all_zero(&self) -> bool {
        self.bits.iter().all(|x| *x == 0)
    }

    pub fn get(&self, i: usize) -> u8 {
        self.bits[i]
    }

    pub fn norm(&self) -> u8 {
        let mut x = 0;
        for b in self.bits.iter() {
            x ^= b;
        }
        x
    }
}

impl ops::Mul<&BinaryVector> for &BinaryVector {
    type Output = u8;

    fn mul(self, rhs: &BinaryVector) -> Self::Output {
        assert_eq!(self.bits.len(), rhs.bits.len());
        let mut x = 0;
        for i in 0..self.bits.len() {
            x ^= self.bits[i] * rhs.bits[i];
        }
        x
    }
}

impl ops::Mul<u8> for &BinaryVector {
    type Output = BinaryVector;

    fn mul(self, rhs: u8) -> Self::Output {
        if rhs == 0 {
            BinaryVector::zeros(self.bits.len())
        } else {
            self.clone()
        }
    }
}

impl ops::Add<&BinaryVector> for &BinaryVector {
    type Output = BinaryVector;

    fn add(self, rhs: &BinaryVector) -> Self::Output {
        assert_eq!(self.bits.len(), rhs.bits.len());
        let mut new_bits = Vec::with_capacity(self.bits.len());
        for i in 0..self.bits.len() {
            new_bits.push(self.bits[i] ^ rhs.bits[i]);
        }
        BinaryVector { bits: new_bits }
    }
}

impl BinaryMatrix64 {
    pub fn new() -> Box<dyn BinaryMatrix> {
        Box::new(BinaryMatrix64 {
            rows: 0,
            columns: vec![],
        })
    }
}

impl BinaryMatrix for BinaryMatrix64 {
    fn nrows(&self) -> usize {
        self.rows
    }
    fn ncols(&self) -> usize {
        self.columns.len()
    }
    fn insert_rows(&mut self, new_rows: usize, new_cols: usize) {
        self.rows += new_rows;
        for c in 0..self.columns.len() {
            self.columns[c].resize((self.rows >> 6) + 1, 0u64);
        }
        for _ in 0..new_cols {
            let mut col = Vec::with_capacity((self.rows >> 6) + 1);
            col.resize((self.rows >> 6) + 1, 0u64);
            self.columns.push(col);
        }
    }

    // TODO: bit tilt algorithm?
    fn transpose(&self) -> Box<dyn BinaryMatrix> {
        let mut new = BinaryMatrix64::new();
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
        let x = self.columns[c][r >> 6];
        let shift = r & 0x3f;
        ((x >> shift) & 1) as u8
    }

    fn set(&mut self, r: usize, c: usize, val: u8) {
        assert!(c < self.columns.len());
        assert!(r < self.rows);
        let shift = r & 0x3f;
        if val == 1 {
            self.columns[c][r >> 6] |= 1 << shift;
        } else {
            self.columns[c][r >> 6] &= !(1 << shift);
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
        let mut mat = BinaryMatrix64::new();
        mat.insert_rows(2, 3);
        mat.transpose();
    }
}
