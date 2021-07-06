use crate::axialobject::AxialSystem;
use rayon::prelude::*;
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
