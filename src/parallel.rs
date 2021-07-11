//! This module provides functions which can compute many magnetic fields in parallel using Rayon
//! and ndarray.

use crate::axialobject::AxialSystem;
use ndarray::{Array2, ArrayView1, ArrayView2, ErrorKind, Zip};
use rayon::prelude::*;
use std::convert::TryInto;

/// performs the computation of the magnetic field in parallel for multiple positions
/// simultaneously using [Rayon](https://docs.rs/rayon). Acts on the basic rust vector containing
/// length 3 position array.
///
/// # Arguments:
/// 'axialsystem' instance of AxialSystem
/// *'positions' a vector containing a list of positions to compute
pub fn get_b_vec(axialsystem: &AxialSystem, positions: Vec<[f64; 3]>, tol: f64) -> Vec<[f64; 3]> {
    let fields = positions
        .into_par_iter()
        .map(|i| axialsystem.get_field(i, &tol));
    fields.collect()
}

/// parallel iterates over an nd array containing positions.
///
/// Takes a 2d ndarray object and uses a parallel iterator to compute
/// the magnetic field for all rows in the array. Initalises a new array of the
/// same shape to return the fields.
///
/// # Arguments
/// `axialsystem` the AxialSystem to compute fields for
/// `positions` 2d ndarray array with shape (_,3)
/// `tol` the tolerance to compute fields to
pub fn get_b_ndarray(
    axialsystem: &AxialSystem,
    positions: &ArrayView2<f64>,
    tol: f64,
) -> Result<Array2<f64>, ErrorKind> {
    match positions.dim().1 {
        3 => {}
        _ => return Err(ErrorKind::IncompatibleShape),
    }
    let mut fields = Array2::<f64>::zeros(positions.dim()); //fields array to populate
    Zip::from(fields.rows_mut())
        .and(positions.rows())
        .par_for_each(|mut field, position| {
            let res = _get_field_ndarray(position, axialsystem, &tol);
            field[0] = res[0];
            field[1] = res[1];
            field[2] = res[2];
        });
    Ok(fields)
}
// internal function which computes the field for each row in the ndarray, only
// called from a function which has already checked the length of the slices.
// panics if cant convert ArrayView1 to [f64;3]
fn _get_field_ndarray(position: ArrayView1<f64>, axialsystem: &AxialSystem, tol: &f64) -> [f64; 3] {
    let pos_as_array = match position.as_slice() {
        Some(x) => x,
        None => panic!("failed to convert ndarrayview to slice"),
    };
    let pos_as_array: [f64; 3] = match pos_as_array.try_into() {
        Ok(x) => x,
        Err(_) => panic!("failed to convert slice to [f64;3]"),
    };
    axialsystem.get_field(pos_as_array, tol)
}

#[cfg(test)]
mod paralell_test {
    use super::*;
    #[test]
    fn test_axial() {
        let mut myvec: Vec<[f64; 3]> = Vec::new();
        let z = 0.0;
        let y = 0.0;
        let tol = 1e-12;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z]);
        }
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = get_b_vec(&myaxial, myvec.clone(), tol);
        let mut res_seq: Vec<[f64; 3]> = Vec::new();
        for pos in myvec.into_iter() {
            res_seq.push(myaxial.get_field(pos, &tol));
        }
        assert_eq!(res_para, res_seq);
    }
    #[test]
    fn test_off_axial() {
        let mut myvec: Vec<[f64; 3]> = Vec::new();
        let z = 0.2;
        let y = 0.2;
        let tol = 1e-12;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z]);
        }
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = get_b_vec(&myaxial, myvec.clone(), tol);
        let mut res_seq: Vec<[f64; 3]> = Vec::new();
        for pos in myvec.into_iter() {
            res_seq.push(myaxial.get_field(pos, &tol));
        }
        assert_eq!(res_para, res_seq);
    }

    #[test]
    fn test_ndarray() {
        let mut myvec: Vec<[f64; 3]> = Vec::new();
        let z = 0.2;
        let y = 0.2;
        let tol = 1e-12;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z]);
        }
        let mypositions = ndarray::arr2(&myvec);
        let mypositions = mypositions.view();

        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = match get_b_ndarray(&myaxial, &mypositions, tol) {
            Ok(x) => x,
            Err(_) => panic!(),
        };
        let mut res_seq: Vec<[f64; 3]> = Vec::new();
        for pos in myvec.into_iter() {
            res_seq.push(myaxial.get_field(pos, &tol));
        }
        let res_seq = ndarray::arr2(&res_seq);
        assert_eq!(res_para, res_seq);
    }
    #[test]
    fn test_ndarray_wrong_length() {
        let mut myvec: Vec<[f64; 4]> = Vec::new();
        let z = 0.2;
        let y = 0.2;
        let tol = 1e-12;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z, 0.0]);
        }
        let mypositions = ndarray::arr2(&myvec);
        let mypositions = mypositions.view();

        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = get_b_ndarray(&myaxial, &mypositions, tol);
        assert_eq!(res_para, Err(ndarray::ErrorKind::IncompatibleShape));
    }
}
