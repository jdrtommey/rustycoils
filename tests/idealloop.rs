///
/// Test the ideal loop case against the analytic formula for the
/// ideal current loop.
///
#[cfg(test)]
extern crate rgsl;
use rustycoils;
const MU0: f64 = 1.25663706212 * 1e-6;
const PI: f64 = std::f64::consts::PI;
fn ellip_k(k: f64) -> f64 {
    rgsl::elliptic::legendre::complete::ellint_Kcomp(k, rgsl::Mode::PrecDouble)
}
fn ellip_e(k: f64) -> f64 {
    rgsl::elliptic::legendre::complete::ellint_Ecomp(k, rgsl::Mode::PrecDouble)
}
fn err_tol(x: f64, y: f64, tol: &f64) -> bool {
    if (x - y).abs() < *tol {
        true
    } else {
        false
    }
}
fn float_compare(a: f64, b: f64, epsilon: f64) -> bool {
    if a == b {
        return true;
    } else if ((a - b) / b).abs() < epsilon {
        return true;
    }
    false
}
#[test]
fn check_external_ellipticals_correct() {
    assert!(err_tol(ellip_k(0.1), 1.5747455615174, &1e-13));
    assert!(err_tol(ellip_k(0.2), 1.5868678474542, &1e-13));
    assert!(err_tol(ellip_k(-0.8), 1.9953027776647, &1e-13));
    assert!(err_tol(ellip_e(0.1), 1.5668619420217, &1e-13));
    assert!(err_tol(ellip_e(0.2), 1.5549685462425, &1e-13));
    assert!(err_tol(ellip_e(-0.8), 1.2763499431699, &1e-13));
}

fn get_x(z: f64, z0: f64, radius: f64) -> f64 {
    (z - z0) / radius
}
fn get_b0(current: f64, radius: f64) -> f64 {
    (MU0 * current) / (2.0 * radius)
}

fn get_nu(r: f64, radius: f64) -> f64 {
    r / radius
}
fn get_k_squared(nu: f64, x: f64) -> f64 {
    (2.0 * nu) / (1.0 + nu.powi(2) + x.powi(2))
}
fn get_t(k_squared: f64) -> f64 {
    ((2.0 * k_squared) / (1.0 + k_squared)).sqrt()
}

fn get_axial_field(r: f64, z: f64, z0: f64, radius: f64, current: f64) -> f64 {
    let x = get_x(z, z0, radius);
    let b0 = get_b0(current, radius);
    if r == 0.0 {
        return b0 * (1.0 / (1.0 + x.powi(2)).powf(1.5));
    } else {
        let nu = get_nu(r, radius);
        let k_squared = get_k_squared(nu, x);
        let t = get_t(k_squared);
        let overall_term = (b0 * t) / (2.0 * PI * nu.sqrt());
        let second_factor = (k_squared - nu) / (nu * (1.0 - k_squared));
        let first_ellip = ellip_k(t);
        let second_ellip = ellip_e(t);
        overall_term * (first_ellip + second_factor * second_ellip)
    }
}
fn get_radial_field(r: f64, z: f64, z0: f64, radius: f64, current: f64) -> f64 {
    let x = get_x(z, z0, radius);
    let b0 = get_b0(current, radius);
    if r == 0.0 {
        return 0.0;
    } else {
        let nu = get_nu(r, radius);
        let k_squared = get_k_squared(nu, x);
        let t = get_t(k_squared);
        let overall_term = (b0 * t * x) / (2.0 * PI * nu * nu.sqrt());
        let second_factor = 1.0 / (1.0 - k_squared);
        let first_ellip = ellip_k(t);
        let second_ellip = ellip_e(t);
        overall_term * (-first_ellip + second_factor * second_ellip)
    }
}

#[test]
fn compare_unit_loop_on_axis_center() {
    let mut myloop = rustycoils::AxialSystem::default();
    let radius = 1.0;
    let current = 1.0;
    let x0 = 0.0;
    let z = 1.0;
    let r = 0.0;
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field_axial(&z, &r, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    assert_eq!(computed_field[0], ana_axial_field);
    assert!(float_compare(computed_field[0], ana_axial_field, 0.001));
    assert!(float_compare(computed_field[1], ana_radial_field, 0.001));
}
#[test]
fn compare_unit_loop_off_axis_center() {
    let mut myloop = rustycoils::AxialSystem::default();
    let radius = 1.0;
    let current = 1.0;
    let x0 = 0.0;
    let z = 1.0;
    let r = 0.1;
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field_axial(&z, &r, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    assert!(float_compare(computed_field[0], ana_axial_field, 0.001));
    assert!(float_compare(computed_field[1], ana_radial_field, 0.001));
}
#[test]
fn compare_unit_loop_off_axis_center_arb_numbers() {
    let mut myloop = rustycoils::AxialSystem::default();
    let radius = 1.9;
    let current = 3.0;
    let x0 = 0.1;
    let z = 1.3;
    let r = 0.4;
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field_axial(&z, &r, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    assert!(float_compare(computed_field[0], ana_axial_field, 0.001));
    assert!(float_compare(computed_field[1], ana_radial_field, 0.001));
}
#[test]
fn compare_unit_loop_off_axis_center_xyz() {
    let mut myloop = rustycoils::AxialSystem::default();
    let radius = 1.9;
    let current = 3.0;
    let x0 = 0.1;
    let z = 1.3;
    let r = 0.4;
    let pos = [z, r, 0.0];
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field(pos, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    assert!(float_compare(computed_field[0], ana_axial_field, 0.001));
    assert!(float_compare(computed_field[1], ana_radial_field, 0.001));
}
#[test]
fn compare_unit_loop_off_axis_center_xyz_rotate() {
    let mut myloop = rustycoils::AxialSystem::default();
    let _res = myloop.transform_y();
    let radius = 1.9;
    let current = 3.0;
    let x0 = 0.1;
    let z = 1.3;
    let r = 0.4;
    let pos = [r, z, 0.0];
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field(pos, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    assert!(float_compare(computed_field[1], ana_axial_field, 0.001));
    assert!(float_compare(computed_field[0], ana_radial_field, 0.001));
}
//place the radial in both x and y.
#[test]
fn compare_unit_loop_off_axis_center_r_in_both() {
    let mut myloop = rustycoils::AxialSystem::default();
    let _res = myloop.transform_z();
    let radius = 1.9;
    let current = 3.0;
    let x0 = 0.1;
    let z = 1.3;
    let r = 0.4;
    let x = (r) / f64::sqrt(2.0);
    let y = x;
    let pos = [x, y, z];
    let _res = myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_field = myloop.get_field(pos, &1e-14);
    let ana_axial_field = get_axial_field(r, z, x0, radius, current);
    let ana_radial_field = get_radial_field(r, z, x0, radius, current);
    let computed_radial =
        (f64::powi(computed_field[0], 2) + f64::powi(computed_field[1], 2)).sqrt();
    assert!(float_compare(computed_field[2], ana_axial_field, 0.001));
    assert!(float_compare(computed_radial, ana_radial_field, 0.001));
}
