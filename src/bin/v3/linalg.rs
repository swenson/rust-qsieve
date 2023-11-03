extern crate nalgebra as na;
use na::DMatrix;
use na::DVector;

pub fn left_kernel_gf2(mat: &DMatrix<u8>) -> Option<Vec<DVector<u8>>> {
    // compute kernel of matrix over GF(2) using Gaussian elimination
    let mut extend = mat.clone().transpose();
    // extend out with an identity matrix of size RxR
    let maxr = extend.nrows();
    let maxc = extend.ncols();
    extend = extend.insert_rows(maxr, maxc, 0);
    for c in 0..maxc {
        extend[(maxr + c, c)] = 1;
    }

    // pivot on a row for each column
    let mut nextc = 0;
    for r in 0..maxr {
        let mut pivot = -1i32;
        for c in nextc..maxc {
            if extend[(r, c)] != 0 {
                pivot = c as i32;
                break;
            }
        }
        if pivot == -1 {
            continue;
        }

        extend.swap_columns(nextc, pivot as usize);

        for c in (nextc as usize + 1)..maxc {
            if extend[(r, c)] == 1 {
                xor_col(&mut extend, c, nextc);
            }
        }
        nextc += 1;
    }

    let mut kernel = vec![];

    for c in 0..maxc {
        if extend.column_part(c, maxr).iter().all(|x| *x == 0) {
            let v = DVector::from(extend.column(c).rows(maxr, extend.nrows() - maxr));
            kernel.push(v);
        }
    }

    return Some(kernel);
}

fn xor_col(mat: &mut DMatrix<u8>, dst: usize, op: usize) {
    for r in 0..mat.nrows() {
        mat[(r, dst)] ^= mat[(r, op)];
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn kernel_test() {
        // 0010
        // 1100
        // 0010
        // left kernel = [1 0 1]

        let mat: DMatrix<u8> =
            DMatrix::from_row_slice(3, 4, &vec![0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0]);
        assert_eq!(
            left_kernel_gf2(&mat)
                .unwrap()
                .get(0)
                .unwrap()
                .iter()
                .map(|x| *x)
                .collect::<Vec<u8>>(),
            vec![1, 0, 1]
        );

        // [1 1 0 1 1 1]
        // [1 1 1 0 1 1]
        // [0 1 1 1 0 0]
        // [1 0 1 0 1 1]
        // [0 0 0 1 0 0]
        let mat: DMatrix<u8> = DMatrix::from_row_slice(
            5,
            6,
            &vec![
                1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1,
                0, 0,
            ],
        );
        assert_eq!(
            left_kernel_gf2(&mat)
                .unwrap()
                .get(0)
                .unwrap()
                .iter()
                .map(|x| *x)
                .collect::<Vec<u8>>(),
            vec![1, 0, 1, 1, 0]
        );
    }
}
