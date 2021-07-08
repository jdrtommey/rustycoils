use crate::axialobject::AxialSystem;
use ndarray::Zip;
use ndarray::{Array1, Array2, ArrayView1, ArrayView2};
use rayon::prelude::*;

use std::convert::TryInto;
/// performs the computation of the magnetic field in parallel for multiple positions
/// simultaneously using [Rayon](https://docs.rs/rayon).
///
/// # Arguments:
/// 'axialsystem' instance of AxialSystem
/// *'positions' a vector containing a list of positions to compute
pub fn get_b_parallel(
    axialsystem: &AxialSystem,
    positions: Vec<[f64; 3]>,
    tol: f64,
) -> Vec<[f64; 3]> {
    let fields = positions
        .into_par_iter()
        .map(|i| axialsystem.get_field(i, &tol));
    fields.collect()
}
/// parallel iterates over an nd array containing positions.
/// not optimized
pub fn get_b_ndarray(
    axialsystem: &AxialSystem,
    positions: ArrayView2<f64>,
    tol: f64,
) -> Array2<f64> {
    let mut fields = Array2::<f64>::zeros(positions.dim()); //fields array to populate
    Zip::from(fields.rows_mut())
        .and(positions.rows())
        .for_each(|mut field, position| {
            let res = _get_field_ndarray(&position, axialsystem, &tol);
            field[0] = res[0];
            field[1] = 0.1;
            field[2] = res[2];
        });
    fields
}

fn _get_field_ndarray(
    position: &ArrayView1<f64>,
    axialsystem: &AxialSystem,
    tol: &f64,
) -> Array1<f64> {
    let pos_as_array = match position.as_slice() {
        Some(x) => x,
        None => panic!(),
    };
    let pos_as_array: [f64; 3] = match pos_as_array.try_into() {
        Ok(x) => x,
        Err(_) => panic!(),
    };
    let field = axialsystem.get_field(pos_as_array, tol).to_vec();
    Array1::from(field)
}
// added for benchmarking
pub fn get_b_seq(axialsystem: &AxialSystem, positions: Vec<[f64; 3]>, tol: f64) -> Vec<[f64; 3]> {
    let fields = positions
        .into_iter()
        .map(|i| axialsystem.get_field(i, &tol));
    fields.collect()
}
#[cfg(test)]
mod paralell_test {
    use super::*;
    #[test]
    fn test_axial() {
        let mut myvec: Vec<[f64; 3]> = Vec::new();
        let z = 0.0;
        let y = 0.0;
        let tol = 1e-18;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z]);
        }
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = get_b_parallel(&myaxial, myvec.clone(), tol);
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
        let tol = 1e-18;
        for x in -1000..1000 {
            myvec.push([(x as f64) * 0.01, y, z]);
        }
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 0.0, 1.0);
        let res_para = get_b_parallel(&myaxial, myvec.clone(), tol);
        let mut res_seq: Vec<[f64; 3]> = Vec::new();
        for pos in myvec.into_iter() {
            res_seq.push(myaxial.get_field(pos, &tol));
        }
        assert_eq!(res_para, res_seq);
    }
}
