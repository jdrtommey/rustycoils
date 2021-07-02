#[cfg(test)]
extern crate rgsl;
use solenoid;
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
#[test]
fn check_external_ellipticals_correct() {
    let x = rgsl::elliptic::legendre::complete::ellint_Kcomp(0.1, rgsl::Mode::PrecDouble);
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
fn get_w0(radius: f64, b0: f64) -> f64 {
    0.5 * radius.powi(2) * b0
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
        overall_term * (first_ellip + second_factor * second_ellip)
    }
}

#[test]
fn compare_unit_loop_on_axis_center() {
    let mut myloop = solenoid::AxialObject::default();
    let radius = 1.0;
    let current = 1.0;
    let x0 = 0.0;
    let pos = (0.0, 0.0, 0.0);
    myloop.add_loop("loop1".to_string(), radius, x0, current);
    print!("{}", myloop);
    let computed_axial_field = myloop.get_field(&pos, &1e-20).0;
    let ana_axial_field = get_axial_field(0.0, x0, 0.0, radius, current);
    assert!(err_tol(computed_axial_field, ana_axial_field, &1e-14));
}
#[test]
fn compare_unit_loop_off_axis() {
    let mut myloop = solenoid::AxialObject::default();
    let radius = 1.0;
    let current = 1.0;
    let x0 = 0.0;
    let pos = (0.0, 0.01, 0.0);
    myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_axial_field = myloop.get_field(&pos, &1e-20).0;
    let ana_axial_field = get_axial_field(0.01, 0.0, x0, radius, current);
    assert_eq!(computed_axial_field, ana_axial_field);
    assert!(err_tol(computed_axial_field, ana_axial_field, &1e-10));
}
#[test]
fn compare_unit_loop_off_axis_radial() {
    let mut myloop = solenoid::AxialObject::default();
    let radius = 1.0;
    let current = 1.0;
    let x0 = 0.0;
    let r = 0.25;
    let z = 1.0;
    let pos = (z, r, 0.0);
    myloop.add_loop("loop1".to_string(), radius, x0, current);
    let computed_axial_field = myloop.get_field(&pos, &1e-20).1;
    let ana_axial_field = get_radial_field(r, z, x0, radius, current);
    assert_eq!(computed_axial_field, ana_axial_field);
    assert!(err_tol(computed_axial_field, ana_axial_field, &1e-9));
}
