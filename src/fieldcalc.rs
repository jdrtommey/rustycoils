// REFERENCE PAPER = "Off-Axis Expansion Solution of Laplace's Equation: Application to Accurate
// and Rapid Calculation of Coil Magnetic Fields" by Robert H. Jackson.
// each primitive impliments Prititive trait which contains a get_fields(z,r,tol) function
// which returns the axial and radial fields.
// Primitives are described in Fig.1. of reference paper and include
// an ideal current loop, an annular, a thin solenoid and a coil.
pub mod primitives {
    use super::polynomials;
    use std::fmt;
    const MU0: f64 = 1.25663706212e-6;
    //trait which acts on all primitive shapes to sum the analtyic polynomials.
    pub trait Primitive {
        // finds the value of the nth derivative at location x, specific to each primitive
        fn get_nth_derivative(&self, n: &u32, x: &f64) -> f64;

        fn get_max_depth(&self) -> u32;
        // impliments equation 14a in reference paper.
        fn get_axial_field(&self, z: &f64, r: &f64, tol: &f64) -> f64 {
            let mut answer = 0.0;
            let mut diff = 1.0;
            let mut counter: i32 = 0;
            while diff > *tol && counter <= self.get_max_depth() as i32 {
                let old_answer = answer;
                let next_derivative = self.get_nth_derivative(&(2 * counter as u32), &z);
                let denominator =
                    f64::powi(2.0, 2 * counter) * (u32::pow(factorial(counter as u32), 2) as f64);
                let numerator = i32::pow(-1, counter as u32) as f64 * f64::powi(*r, 2 * counter);
                answer += next_derivative * numerator / denominator;
                counter += 1;
                diff = f64::abs(answer - old_answer);
            }
            answer
        }
        fn get_radial_field(&self, z: &f64, r: &f64, tol: &f64) -> f64 {
            let mut answer = 0.0;
            let mut diff = 1.0;
            let mut counter: i32 = 0;
            while diff > *tol && counter < self.get_max_depth() as i32 {
                let old_answer = answer;
                let next_derivative = self.get_nth_derivative(&(1 + 2 * counter as u32), &z);
                let denominator = f64::powi(2.0, 2 * counter + 1)
                    * ((factorial(counter as u32) * factorial((counter + 1) as u32)) as f64);
                let numerator =
                    i32::pow(-1, 1 + counter as u32) as f64 * f64::powi(*r, 2 * counter + 1);
                answer += next_derivative * numerator / denominator;
                counter += 1;
                diff = f64::abs(answer - old_answer);
            }
            answer
        }
        fn get_fields(&self, z: &f64, r: &f64, tol: &f64) -> (f64, f64) {
            (
                self.get_axial_field(z, r, tol),
                self.get_radial_field(z, r, tol),
            )
        }
    }

    fn factorial(x: u32) -> u32 {
        if x < 2 {
            1
        } else {
            x * factorial(x - 1)
        }
    }

    fn _get_loop_normalized_b(n: &u32, x: &f64, norm: &f64) -> f64 {
        let (poly_hcf, poly) = polynomials::primitive_polynomials::get_loop_poly(*n);
        let polynomial = polynomials::SolenoidPolynomial::new(poly, poly_hcf);

        let poly_res_x = polynomial.compute(&x);

        let x_denom = f64::sqrt(1.0 + x * x) * f64::powi(1.0 + x * x, (n + 1) as i32);

        let total_norm = 1.0 / f64::powi(*norm, *n as i32);
        total_norm * (poly_res_x / x_denom)
    }
    //function returns b^(n)(x) normalised by 1/(norm^n) for the annular disk.
    fn _get_annular_normalized_b(n: &u32, x: &f64, norm: &f64) -> f64 {
        let (p_hcf, p_poly, q_hcf, q_poly) =
            polynomials::primitive_polynomials::get_annular_poly(*n);

        let p_polynomial = polynomials::SolenoidPolynomial::new(p_poly, p_hcf);
        let q_polynomial = polynomials::SolenoidPolynomial::new(q_poly, q_hcf);

        let x_res_p = p_polynomial.compute(&x);
        let x_res_q = q_polynomial.compute(&x);

        fn p_denom(x: &f64, n: &u32) -> f64 {
            f64::powi(f64::sqrt(x * x + 1.0) + 1.0, *n as i32) * f64::powi(x * x + 1.0, *n as i32)
        }
        fn q_denom(x: &f64, n: &u32) -> f64 {
            f64::powi(f64::sqrt(x * x + 1.0) + 1.0, *n as i32)
                * f64::powi(x * x + 1.0, *n as i32)
                * f64::sqrt(x * x + 1.0)
        }

        let b_deriv_x = x_res_p / p_denom(&x, &n) + x_res_q / q_denom(&x, &n);

        (1.0 / f64::powi(*norm, *n as i32)) * b_deriv_x
    }
    // IDEAL LOOP PRIMITIVE.
    #[derive(Debug, PartialEq, Clone, Copy)]
    pub struct IdealWire {
        radius: f64,  //radius of the wire loop
        current: f64, //current flowing in the wire
        z0: f64,      //placement of the wire along the axis of symmetry
        max_depth: u32,
    }
    impl IdealWire {
        pub fn new(radius: f64, current: f64, z0: f64) -> IdealWire {
            IdealWire {
                radius,
                current,
                z0,
                max_depth: 7,
            }
        }
        pub fn set_radius(&mut self, radius: f64) {
            self.radius = radius;
        }
        pub fn set_current(&mut self, current: f64) {
            self.current = current;
        }
        pub fn set_z0(&mut self, z0: f64) {
            self.z0 = z0;
        }
    }
    impl Primitive for IdealWire {
        fn get_nth_derivative(&self, n: &u32, z: &f64) -> f64 {
            let x = (z - self.z0) / self.radius;
            let b0 = (self.current * MU0) / (2.0 * self.radius);

            let normed_b = _get_loop_normalized_b(&n, &x, &self.radius);
            b0 * normed_b
        }
        fn get_max_depth(&self) -> u32 {
            self.max_depth
        }
    }

    impl fmt::Display for IdealWire {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(
                f,
                "Ideal Wire: radius={},axial position={},current={} ",
                self.radius, self.z0, self.current
            )
        }
    }
    //THIN ANNULAR PRIMITIVE

    #[derive(Debug, PartialEq, Clone, Copy)]
    pub struct ThinAnnular {
        radius: f64,          //units of m
        current_density: f64, //units of A/m
        z0: f64,              //units of m
        thickness: f64,       //units of m
        max_depth: u32,       //10 hardwired terms so depth of 5
    }
    impl ThinAnnular {
        pub fn new(radius: f64, thickness: f64, current_density: f64, z0: f64) -> ThinAnnular {
            ThinAnnular {
                radius,
                thickness,
                current_density,
                z0,
                max_depth: 5,
            }
        }
        pub fn set_radius(&mut self, radius: f64) {
            self.radius = radius;
        }
        //takes a current measured in A and converts it into the current density as
        //current/thickness
        pub fn set_current(&mut self, current: f64) {
            self.current_density = current / self.thickness;
        }
        pub fn set_z0(&mut self, z0: f64) {
            self.z0 = z0;
        }
        pub fn set_thickness(&mut self, thickness: f64) {
            let current = self.thickness * self.current_density; //find the total current.
            self.thickness = thickness;
            self.set_current(current);
        }
    }
    impl Primitive for ThinAnnular {
        fn get_nth_derivative(&self, n: &u32, z: &f64) -> f64 {
            let x = (z - self.z0) / self.radius;
            let rho = (self.radius + self.thickness) / self.radius;
            let xi = x / rho;
            let b0 = MU0 * self.current_density / 2.0;
            // if zeroth derivative return the on axis field value Eq.32 in reference paper.
            // else compute the polynomials.
            if *n == 0 {
                let prefactor = MU0 * self.current_density / 2.0;
                let x_factor =
                    (1.0 / f64::sqrt(1.0 + x * x)) - f64::ln(1.0 + f64::sqrt(1.0 + x * x));
                let xi_factor =
                    (1.0 / f64::sqrt(1.0 + xi * xi)) - f64::ln(1.0 + f64::sqrt(1.0 + xi * xi));
                let rho_factor = f64::ln(rho);
                return prefactor * (x_factor - xi_factor - rho_factor);
            }
            let term1 = _get_annular_normalized_b(n, &x, &self.radius);
            let term2 = _get_annular_normalized_b(n, &xi, &(self.radius + self.thickness));
            b0 * (term1 - term2)
        }
        fn get_max_depth(&self) -> u32 {
            self.max_depth
        }
    }

    impl fmt::Display for ThinAnnular {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(
                f,
                "ThinAnnular: radius={},thickness={},axial position={},current={} ",
                self.radius, self.thickness, self.z0, self.current_density
            )
        }
    }

    // ThinSolenoid primitive.
    #[derive(Debug, PartialEq, Clone, Copy)]
    pub struct ThinSolenoid {
        radius: f64,
        current_density: f64,
        length: f64,
        z0: f64,
        max_depth: u32,
    }
    impl ThinSolenoid {
        pub fn new(radius: f64, length: f64, current_density: f64, z0: f64) -> ThinSolenoid {
            ThinSolenoid {
                radius,
                current_density,
                length,
                z0,
                max_depth: 14,
            }
        }
        pub fn set_radius(&mut self, radius: f64) {
            self.radius = radius;
        }
        //sets the current_density of the solenoid as current/length
        pub fn set_current(&mut self, current: f64) {
            self.current_density = current / self.length;
        }
        pub fn set_z0(&mut self, z0: f64) {
            self.z0 = z0;
        }
        pub fn set_length(&mut self, length: f64) {
            let current = self.current_density * self.length;
            self.length = length;
            self.set_current(current);
        }
    }
    impl Primitive for ThinSolenoid {
        fn get_nth_derivative(&self, n: &u32, z: &f64) -> f64 {
            let x = (z - self.z0) / self.radius;
            let eta = self.length / self.radius;
            let b0 = (self.current_density * MU0) / 2.0;

            if *n == 0 {
                let term1 = x / (f64::sqrt(1.0 + x * x));
                let term2 = (x - eta) / f64::sqrt(1.0 + (x - eta) * (x - eta));
                return b0 * (term1 - term2);
            }

            let term_x = (1.0 / f64::powi(self.radius, *n as i32))
                * _get_loop_normalized_b(&(n - 1), &x, &1.0);
            let term_eta = (1.0 / f64::powi(self.radius, *n as i32))
                * _get_loop_normalized_b(&(n - 1), &(x - eta), &1.0);

            b0 * (term_x - term_eta)
        }
        fn get_max_depth(&self) -> u32 {
            self.max_depth
        }
    }

    impl fmt::Display for ThinSolenoid {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(
                f,
                "ThinSolenoid: radius={},length={},axial position={},current={} ",
                self.radius, self.length, self.z0, self.current_density
            )
        }
    }

    // coil primitive.

    #[derive(Debug, PartialEq, Clone, Copy)]
    pub struct CoilSolenoid {
        radius: f64,
        current_density: f64,
        length: f64,
        thickness: f64,
        z0: f64,
        max_depth: u32,
    }
    impl CoilSolenoid {
        pub fn new(
            radius: f64,
            length: f64,
            thickness: f64,
            current_density: f64,
            z0: f64,
        ) -> CoilSolenoid {
            CoilSolenoid {
                radius,
                current_density,
                thickness,
                length,
                z0,
                max_depth: 5,
            }
        }
        pub fn set_radius(&mut self, radius: f64) {
            self.radius = radius;
        }
        // set current density as current/(thickness*length)
        pub fn set_current(&mut self, current: f64) {
            self.current_density = current / (self.thickness * self.length);
        }
        pub fn set_z0(&mut self, z0: f64) {
            self.z0 = z0;
        }
        pub fn set_length(&mut self, length: f64) {
            let current = self.current_density * self.length * self.thickness;
            self.length = length;
            self.set_current(current);
        }
        pub fn set_thickness(&mut self, thickness: f64) {
            let current = self.current_density * self.length * self.thickness;
            self.thickness = thickness;
            self.set_current(current);
        }
    }
    impl Primitive for CoilSolenoid {
        fn get_nth_derivative(&self, n: &u32, z: &f64) -> f64 {
            let x = (z - self.z0) / self.radius;
            let eta = self.length / self.radius;
            let rho = (self.radius + self.thickness) / self.radius;
            let xi = x / rho;
            let b0 = self.current_density * MU0 / 2.0;
            // the zeroth derivative is the axial equation given in Eq. 31.
            if *n == 0 {
                return _zeroth_order_coil(self.current_density, self.radius, x, rho, eta);
            }
            //The first derivative is found from the equations for the field of an Annular Eq.32 and Eq.34
            if *n == 1 {
                let term_entrance = _annular_normalized(x) - _annular_normalized(xi) - f64::ln(rho);
                let term_exit = _annular_normalized(x - eta)
                    - _annular_normalized(xi - eta / rho)
                    - f64::ln(rho);
                return b0 * (term_entrance - term_exit);
            }
            //All other orders are computed from the polynomials of an annular disk.
            let term1 =
                (1.0 / f64::powi(self.radius, *n as i32)) * _get_annular_normalized_b(n, &x, &1.0);
            let term2 = (1.0 / f64::powi(self.radius + self.thickness, *n as i32))
                * _get_annular_normalized_b(n, &(x - eta), &1.0);
            let term3 =
                (1.0 / f64::powi(self.radius, *n as i32)) * _get_annular_normalized_b(n, &xi, &1.0);
            let term4 = (1.0 / f64::powi(self.radius + self.thickness, *n as i32))
                * _get_annular_normalized_b(n, &(xi - eta / rho), &1.0);
            b0 * (term1 - term2 - (term3 - term4))
        }
        fn get_max_depth(&self) -> u32 {
            self.max_depth
        }
    }

    //returns the normalized b(x) for the annular.
    fn _annular_normalized(x: f64) -> f64 {
        1.0 / f64::sqrt(1.0 + x * x) - f64::ln(1.0 + f64::sqrt(1.0 + x * x))
    }

    // impliment the coil axis magnetic field.
    fn _zeroth_order_coil(current_density: f64, radius: f64, x: f64, rho: f64, eta: f64) -> f64 {
        let physical = MU0 * current_density * radius / 2.0;
        let term1_numerator = rho + f64::sqrt(rho * rho + x * x);
        let term1_denom = 1.0 + f64::sqrt(1.0 + x * x);
        let term1 = x * f64::ln(term1_numerator / term1_denom);

        let term2_numerator = rho + f64::sqrt(f64::powi(rho, 2) + (x - eta) * (x - eta));
        let term2_denom = 1.0 + f64::sqrt(1.0 + (x - eta) * (x - eta));
        let term2 = (x - eta) * f64::ln(term2_numerator / term2_denom);
        physical * (term1 - term2)
    }

    impl fmt::Display for CoilSolenoid {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(
                f,
                "CoilSolenoid: radius={},length={},thickness={},axial position={},current={} ",
                self.radius, self.length, self.thickness, self.z0, self.current_density
            )
        }
    }
    #[cfg(test)]
    mod test_coil {
        use super::*;
        //The axial field for a coil if given in Eq. 31.
        //
        //assuming the following physical dimensions for the coils
        //
        // radius = 0.5,thickness = 1.0, length = 2.0, z=1.0,z0=0.0,current_density=4.0
        //
        // x = 2.0
        // eta = 4.0
        // rho = 3.0
        #[test]
        fn test_zeroth_order_derivative() {
            let current = 4.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let thickness = 1.0;

            let x = (z - z0) / radius;
            let eta = length / radius;
            let rho = (radius + thickness) / radius;

            let zeroth = _zeroth_order_coil(current, radius, x, rho, eta);

            assert_eq!(
                zeroth / MU0,
                4.0 * f64::ln((3.0 + f64::sqrt(13.0)) / (1.0 + f64::sqrt(5.0)))
            );
        }
    }
    #[cfg(test)]
    mod test_thin_solenoid {
        use super::*;
        // the field on axis is analytical and can be found using Eq. 29
        //
        // B0*(x/(sqrt(1+x*x)) - (x-eta)/(sqrt(1+(x-eta)^2)))
        //
        // where x = z-z0/radius
        // eta = length/radius
        //
        // if let current = 2.0,z = 1.0, z0 = 0.0, length = 2.0, radius = 0.5
        // x = 2
        // eta = 4
        //
        // B/MU0 = 2/sqrt(5) - (2-4)/sqrt(5) = 4/sqrt(5)
        #[test]
        fn test_order_0() {
            let current = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);
            let first = solenoid.get_nth_derivative(&0, &z);
            assert_eq!(first / MU0, 4.0 / f64::sqrt(5.0));
        }

        // the first order is given by
        // (B0/r0)*(1/((1+x^2)^n+1+0.5) - 1/((1+(x-eta)^2)^n+1+0.5
        //
        // if same parameters as above reduces to
        // 1/((5)^2.5) - 1/((5)^2.5) = 0.0
        #[test]
        fn test_order_1() {
            let current = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);
            let first = solenoid.get_nth_derivative(&1, &z);
            assert_eq!(first / MU0, 0.0);
        }

        // the second order is given by
        //
        // P_1(x) = -3x
        //
        // -3x/(1+x^2)^n-0.5
        // 1/(0.25) * (-6/5^2.5 - 6/5^2.5) = -48/5^2.5
        #[test]
        fn test_order_2() {
            let current = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);
            let second = solenoid.get_nth_derivative(&2, &z);
            assert_eq!(second / MU0, -48.0 / f64::powf(5.0, 2.5));
        }
        #[test]
        fn test_order_3() {
            let current = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);
            let second = solenoid.get_nth_derivative(&3, &z);
            assert_eq!(second / MU0, 0.0);
        }
        #[test]
        fn test_order_4() {
            let current = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let length = 2.0;
            let radius = 0.5;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);
            let second = solenoid.get_nth_derivative(&4, &z);
            assert_eq!(
                second / MU0,
                16.0 * (-26.0 * 15.0 * 2.0 / f64::powf(5.0, 4.5))
            );
        }

        // test against the on-axis magnetic field equation given in Eq.29
        //
        // B_z(0,z) = {MU0*J/2.0][x/sqrt(1+x^2) - (x-eta)/sqrt(1+(x-eta)^2)
        fn analytical(current_density: f64, x: f64, eta: f64) -> f64 {
            (MU0 * current_density / 2.0)
                * (x / f64::sqrt(1.0 + x * x) - (x - eta) / f64::sqrt(1.0 + (x - eta) * (x - eta)))
        }
        #[test]
        fn test_on_axis_field() {
            let current = 2.0;
            let radius = 0.5;
            let length = 2.0;
            let z = 1.0;
            let z0 = 0.0;
            let solenoid = ThinSolenoid::new(radius, length, current, z0);

            let x = (z - solenoid.z0) / solenoid.radius;
            let eta = solenoid.length / solenoid.radius;
            let ana_answer = analytical(current, x, eta);
            assert_eq!(ana_answer, solenoid.get_axial_field(&z, &0.0, &1e-6));
        }
    }
    #[cfg(test)]
    mod test_thin_annular {
        use super::*;
        // To test the thin annular derivative function have hand computed the expressions up to
        // order 5.

        //the zeroth order is analytically given in Eq.32.
        // [MU0*J/2]*( 1/(sqrt(1+x*x)) - ln(1+sqrt(1+x*x)) + ln(1+sqrt(1+xi*xi)) -
        // 1/(sqrt(1+xi*xi)))
        // assuming z = 1.0, z0 = 0.0, radius = 1.0, thickness=1.0,J=1
        // x = (z -z0)/radius => x = 1.0
        // rho = (radius+thickness)/radius => 2.0
        // xi = (x/rho) => 0.5
        // B0/MU0 = 1/2 *((1/sqrt(2) - ln(1+sqrt(2)) - 1/sqrt(1.25) + ln(1+sqrt(1.25)) - ln(2) ) =
        // -0.5056764413
        #[test]
        fn annular_disk_order_0() {
            let radius = 1.0;
            let thickness = 1.0;
            let current_density = 1.0;
            let z0 = 0.0;
            let z = 1.0;
            let annular = ThinAnnular::new(radius, thickness, current_density, z0);
            let answer = annular.get_nth_derivative(&0, &z);
            let error = f64::abs(answer / MU0 - (-0.5056764413));
            assert!(error < 1e-9);
        }
        //the first order can be computed from the polynomial experssions given in Eq.36 using
        //polynomials given in Table 2.
        //
        //P_1(x) = -x = -1.0
        //Q_1(x) = -x^3-2x = -3.0
        //
        //P_1(xi) = -0.5
        //Q_1(xi) = -1.125
        //Using the same initial conditions stated in above test:
        //The nth deriviative is given by:
        // b_1(x) => P_1(x)/(2sqrt(2)+2) + Q_1(x)/(4+2sqrt(2))
        // b_1(xi)=> P_1(xi)/(sqrt(1.25)+1) + Q_1(xi)/(1.25(1.25+sqrt(1.25)))
        //
        // => B_0/MU0 = -1.0/(2sqrt(2)+2) + -3.0/(4+2sqrt(2))-( (1/2.0)*((-0.5)/(sqrt(1.25)+1) +
        // -1.125/(1.25*(1.25+sqrt(1.25))))) = -0.1809941811033...
        #[test]
        fn annular_disk_order_1() {
            let radius = 1.0;
            let thickness = 1.0;
            let current_density = 1.0;
            let z0 = 0.0;
            let z = 1.0;
            let annular = ThinAnnular::new(radius, thickness, current_density, z0);
            let answer = annular.get_nth_derivative(&1, &z);
            let error = f64::abs(answer / MU0 - (-0.1809941811033));
            assert!(error < 1e-9);
        }

        //test the second order derivative for the same parameters again.
        //
        // P_2(x) = x^4 + 4x^2 - 3
        // q_2(x) = 2x^4 + 2x^2 - 3
        //
        // P_2(1.0) = 2
        // Q_2(1.0) = 1
        //
        // P_2(0.5) = -1.9375
        // Q_2(0.5) = -2.375
        //
        // denominator for P is (sqrt(1+x^2)+1)^n * (x*x+1)^2
        // denominator for Q is sqrt(x*x+1) * (sqrt(1+x^2)+1)^n * (x*x+1)^2
        //
        // P(1.0) = 2/((sqrt(2)+1)^2 * 4)
        // Q(1.0) = 1/(4*sqrt(2)*(sqrt(2)+1)^2 )
        // P(0.5) = -1.9375/ (1.25^2 *(sqrt(1.25+1)^2 :)
        // Q(0.5) = -2.375/(sqrt(1.25)*(1.25^2 *(sqrt(1.25+1)^2)
        // answer =0.5*(P(1.0)+Q(1.0) - 0.25*(P(0.5)-Q(0.5))) = 0.130491663998
        #[test]
        fn annular_disk_order_2() {
            let radius = 1.0;
            let thickness = 1.0;
            let current_density = 1.0;
            let z0 = 0.0;
            let z = 1.0;
            let annular = ThinAnnular::new(radius, thickness, current_density, z0);
            let answer = annular.get_nth_derivative(&2, &z);
            let error = f64::abs(answer / MU0 - (0.130491663998));
            assert!(error < 1e-9);
        }
    }

    #[cfg(test)]
    mod test_ideal_loop {
        use super::*;

        //test against hand computed values for the derivatives using
        //equations in reference paper.
        #[test]
        fn test_field_derivative() {
            let radius = 1.0;
            let current = 3.0;
            let z0 = 0.0;
            let myloop = IdealWire::new(radius, current, z0);
            let answer = (MU0 * current) / (2.0 * radius) * (-3.0) / (f64::sqrt(2.0) * 4.0);
            let first_derivative = myloop.get_nth_derivative(&1, &1.0);
            let diff = first_derivative - answer;
            let mut boo = false;
            if diff < 1e-10 {
                boo = true;
            }
            assert!(boo);
        }
        #[test]
        fn test_field_derivative_2() {
            let radius = 1.0;
            let current = 3.0;
            let z0 = 1.0;
            let myloop = IdealWire::new(radius, current, z0);
            let answer = (MU0 * current) / (2.0 * radius) * (-3.0) / (f64::sqrt(5.0) * 25.0);
            let first_derivative = myloop.get_nth_derivative(&1, &3.0);
            let diff = first_derivative - answer;
            let mut boo = false;
            if diff < 1e-10 {
                boo = true;
            }
            assert!(boo);
        }
        #[test]
        fn test_field_derivative_3() {
            let radius = 1.0;
            let current = 3.0;
            let z0 = 1.0;
            let myloop = IdealWire::new(radius, current, z0);
            let answer = (MU0 * current) / (2.0 * radius) * (3.0 * 15.0) / (f64::sqrt(5.0) * 125.0);
            let second_derivative = myloop.get_nth_derivative(&2, &3.0);
            let diff = second_derivative - answer;
            let mut boo = false;
            if diff < 1e-10 {
                boo = true;
            }
            assert!(boo);
        }
        #[test]
        fn test_field_derivative_4() {
            let radius = 1.0;
            let current = 3.0;
            let z0 = 1.0;
            let myloop = IdealWire::new(radius, current, z0);
            let answer =
                (MU0 * current) / (2.0 * radius) * (-15.0 * 26.0) / (f64::sqrt(5.0) * 625.0);
            let second_derivative = myloop.get_nth_derivative(&3, &3.0);
            let diff = second_derivative - answer;
            let mut boo = false;
            if diff < 1e-10 {
                boo = true;
            }
            assert!(boo);
        }
        #[test]
        fn test_factorial() {
            let factorial_results = vec![1, 1, 2, 6, 24, 120, 720, 5040, 40320];
            let mut res: Vec<u32> = Vec::new();
            for i in 0..9 {
                res.push(factorial(i));
            }
            assert_eq!(factorial_results, res);
        }
        //test the on axis field against the anayltical field given in Eq. 23 of reference paper.
        #[test]
        fn test_onaxis_field_ideal_loop() {
            fn on_axis_field(z: &f64, current: &f64, radius: &f64, z0: &f64) -> f64 {
                let physical_part = 1.0 * current * MU0 / 2.0;
                let numerator = radius * radius;
                let demoninator = f64::powf(radius * radius + (z - z0) * (z - z0), 1.5);
                physical_part * numerator / demoninator
            }

            let current = 1.0;
            let radius = 1.0;
            let z0 = 0.0;

            let myloop = IdealWire::new(radius, current, z0);
            let z_positions = vec![
                -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
            ];
            let mut ana_res: Vec<f64> = Vec::new();
            let mut res: Vec<f64> = Vec::new();
            let mut diff_vec: Vec<f64> = Vec::new();
            for z in z_positions {
                let ana = on_axis_field(&z, &current, &radius, &z0);
                ana_res.push(ana);
                let comp = myloop.get_axial_field(&z, &0.0, &1e-6);
                res.push(comp);
                diff_vec.push(f64::abs(comp - ana));
            }

            diff_vec = diff_vec.into_iter().filter(|x| x < &&1e-6).collect();
            assert_eq!(diff_vec.len(), 13);
        }
        #[test]
        fn test_onaxis_field_ideal_loop2() {
            fn on_axis_field(z: &f64, current: &f64, radius: &f64, z0: &f64) -> f64 {
                let physical_part = 1.0 * current * MU0 / 2.0;
                let numerator = radius * radius;
                let demoninator = f64::powf(radius * radius + (z - z0) * (z - z0), 1.5);
                physical_part * numerator / demoninator
            }

            let current = 1.3;
            let radius = 1.9;
            let z0 = 1.7;

            let myloop = IdealWire::new(radius, current, z0);
            let z_positions = vec![
                -3.1, -2.4, -2.2, -1.8, -1.6, -0.33, 0.020, 0.5234, 1.340, 1.534, 2.034, 2.5234,
                3.0234,
            ];
            let mut ana_res: Vec<f64> = Vec::new();
            let mut res: Vec<f64> = Vec::new();
            let mut diff_vec: Vec<f64> = Vec::new();
            for z in z_positions {
                let ana = on_axis_field(&z, &current, &radius, &z0);
                ana_res.push(ana);
                let comp = myloop.get_axial_field(&z, &0.0, &1e-6);
                res.push(comp);
                diff_vec.push(f64::abs(comp - ana));
            }

            diff_vec = diff_vec.into_iter().filter(|x| x < &&1e-6).collect();
            assert_eq!(diff_vec.len(), 13);
        }
        #[test]
        fn test_onaxis_field_ideal_loop3_check_not_equal_when_anayltic_is_wrong() {
            fn on_axis_field(z: &f64, current: &f64, radius: &f64, z0: &f64) -> f64 {
                let physical_part = 3.8 * 1.0 * current * MU0 / 2.0;
                let numerator = radius * radius;
                let demoninator = f64::powf(radius * radius + (z - z0) * (z - z0), 1.5);
                physical_part * numerator / demoninator
            }

            let current = 1.3;
            let radius = 1.9;
            let z0 = 1.7;

            let myloop = IdealWire::new(radius, current, z0);
            let z_positions = vec![
                -3.1, -2.4, -2.2, -1.8, -1.6, -0.33, 0.020, 0.5234, 1.340, 1.534, 2.034, 2.5234,
                3.0234,
            ];
            let mut ana_res: Vec<f64> = Vec::new();
            let mut res: Vec<f64> = Vec::new();
            let mut diff_vec: Vec<f64> = Vec::new();
            for z in z_positions {
                let ana = on_axis_field(&z, &current, &radius, &z0);
                ana_res.push(ana);
                let comp = myloop.get_axial_field(&z, &0.0, &1e-6);
                res.push(comp);
                diff_vec.push(f64::abs(comp - ana));
            }

            diff_vec = diff_vec.into_iter().filter(|x| x < &&1e-6).collect();
            assert!(diff_vec.len() < 13);
        }

        fn on_axis_field(z: &f64, current: &f64, radius: &f64, z0: &f64) -> f64 {
            let physical_part = 1.0 * current * MU0 / 2.0;
            let numerator = radius * radius;
            let demoninator = f64::powf(radius * radius + (z - z0) * (z - z0), 1.5);
            physical_part * numerator / demoninator
        }

        fn off_axis_field(z: &f64, r: &f64, current: &f64, radius: &f64, z0: &f64) -> f64 {
            let x = (z - z0) / radius;
            let b0 = current * MU0 / (2.0 * radius);
            let term1 = 1.0 / (f64::powf(1.0 + x * x, 1.5));
            let term2 = -(12.0 * x - 3.0) * (r * r) / (4.0 * f64::powf(1.0 + x * x, 3.5));
            let term3 = (45.0 * (8.0 * x * x * x * x - 12.0 * x * x + 1.0)) * (r * r * r * r)
                / (f64::powi(2.0, 6) * f64::powf(1.0 + x * x, 5.5));
            let answer = b0 * (term1 + term2 + term3);
            answer
        }
        #[test]
        fn test_off_axis_ana_equals_on_axis() {
            let current = 1.0;
            let radius = 1.0;
            let z0 = 0.0;
            let r = 0.0;

            let z = 0.0;
            let ana = off_axis_field(&z, &r, &current, &radius, &z0);
            let ana2 = on_axis_field(&z, &current, &radius, &z0);
            assert_eq!(ana, ana2);
        }
        #[test]
        fn test_off_axis() {
            let current = 1.0;
            let radius = 1.0;
            let z0 = 0.0;
            let r = 0.001;
            let z = 0.0;
            let myloop = IdealWire::new(radius, current, z0);
            let ans = myloop.get_axial_field(&z, &r, &1e-10);
            let ana = off_axis_field(&z, &r, &current, &radius, &z0);
            let percentage_error = (ans - ans) / ((ans + ana) / 2.0);
            assert!(percentage_error < 1e-3);
        }
        #[test]
        fn test_off_axis1() {
            let current = 1.0;
            let radius = 1.0;
            let z0 = 0.0;
            let r = 0.01;
            let z = 0.0;
            let myloop = IdealWire::new(radius, current, z0);
            let ans = myloop.get_axial_field(&z, &r, &1e-6);
            let ana = off_axis_field(&z, &r, &current, &radius, &z0);
            let percentage_error = f64::abs((ans - ana) / ((ans + ana) / 2.0));
            assert!(percentage_error < 1e-3);
        }
        #[test]
        fn test_off_axis2() {
            let current = 1.0;
            let radius = 1.0;
            let z0 = 0.0;
            let r = 0.01;
            let z = 0.05;
            let myloop = IdealWire::new(radius, current, z0);
            let ans = myloop.get_axial_field(&z, &r, &1e-6);
            let ana = off_axis_field(&z, &r, &current, &radius, &z0);
            let percentage_error = f64::abs((ans - ana) / ((ans + ana) / 2.0));
            assert!(percentage_error < 1e-3);
        }
    }
}
mod polynomials {

    pub struct SolenoidPolynomial {
        coefficients: Vec<f64>,
        hcf: f64,
        //number_terms: usize,
    }

    impl SolenoidPolynomial {
        pub fn new(coefficients: Vec<f64>, hcf: f64) -> SolenoidPolynomial {
            SolenoidPolynomial { coefficients, hcf }
        }
        // compute the sum of the polynomial using Horner's Method.
        pub fn compute(&self, value: &f64) -> f64 {
            let answer = horners_method(&self.coefficients, *value);
            answer * self.hcf
        }
    }

    // given a polynomial containting a0 + x * a1 + x**2 * a2 + ....... x^n * an
    // compute the value for value
    // assumes all the coefficient terms are present even if value is 0.
    fn horners_method(coefficients: &[f64], value: f64) -> f64 {
        let length = coefficients.len();
        let mut answer: f64 = coefficients[length - 1];
        for i in 1..length {
            let j = length - 1 - i;
            answer = answer * value + coefficients[j]
        }
        answer
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn test_horners_method() {
            let mycoeffs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
            let answer = horners_method(&mycoeffs, 1.0);
            assert_eq!(answer, 15.0);
        }
        #[test]
        fn test_horners_method_2() {
            let mycoeffs = vec![9.0, 1.0, 2.0, 3.0, 4.0, 5.0];
            let answer = horners_method(&mycoeffs, 1.0);
            assert_eq!(answer, 24.0);
        }
        #[test]
        fn test_horners_method_3() {
            let mycoeffs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
            let answer = horners_method(&mycoeffs, 2.0);
            assert_eq!(answer, 258.0);
        }
        #[test]
        fn test_horners_method_6() {
            let mycoeffs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
            let answer = horners_method(&mycoeffs, -2.0);
            assert_eq!(answer, -114.0);
        }
        #[test]
        fn test_horners_method_4() {
            let mycoeffs = vec![9.0];
            let answer = horners_method(&mycoeffs, 1.0);
            assert_eq!(answer, 9.0);
        }
        #[test]
        fn test_horners_method_5() {
            let mycoeffs = vec![9.0];
            let answer = horners_method(&mycoeffs, 1.0);
            assert_eq!(answer, 9.0);
        }
        #[test]
        fn test_polynomial_1() {
            let mycoeffs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
            let mypolynomial = SolenoidPolynomial::new(mycoeffs, 2.0);
            let answer = mypolynomial.compute(&2.0);
            assert_eq!(answer, 516.0);
        }
        #[test]
        fn test_polynomial_2() {
            let mycoeffs = vec![0.0, 1.0, 0.0, 3.0, 4.0, 5.0];
            let mypolynomial = SolenoidPolynomial::new(mycoeffs, 2.0);
            let answer = mypolynomial.compute(&2.0);
            assert_eq!(answer, 500.0);
        }
        #[test]
        fn test_polynomial_3() {
            let mycoeffs = vec![0.9, 1.0, 0.0, 3.0, 4.0, 5.0];
            let mypolynomial = SolenoidPolynomial::new(mycoeffs, 2.0);
            let answer = mypolynomial.compute(&2.0);
            assert_eq!(answer, 501.8);
        }
    }

    // module containing the hardwired polynomials for the primitive shapes
    pub mod primitive_polynomials {
        // define the polynomial coefficients for an infinitesimal loop
        // values taken from Table 1 in https://ieeexplore.ieee.org/abstract/document/760416
        const INFLOOP_0: [f64; 1] = [1.0];
        const INFLOOP_0_HCF: f64 = 1.0;
        const INFLOOP_1: [f64; 2] = [0.0, 1.0];
        const INFLOOP_1_HCF: f64 = -3.0;
        const INFLOOP_2: [f64; 3] = [-1.0, 0.0, 4.0];
        const INFLOOP_2_HCF: f64 = 3.0;
        const INFLOOP_3: [f64; 4] = [0.0, -3.0, 0.0, 4.0];
        const INFLOOP_3_HCF: f64 = -15.0;
        const INFLOOP_4: [f64; 5] = [1.0, 0.0, -12.0, 0.0, 8.0];
        const INFLOOP_4_HCF: f64 = 45.0;
        const INFLOOP_5: [f64; 6] = [0.0, 5.0, 0.0, -20.0, 0.0, 8.0];
        const INFLOOP_5_HCF: f64 = -315.0;
        const INFLOOP_6: [f64; 7] = [-5.0, 0.0, 120.0, 0.0, -240.0, 0.0, 64.0];
        const INFLOOP_6_HCF: f64 = 315.0;
        const INFLOOP_7: [f64; 8] = [0.0, -35.0, 0.0, 280.0, 0.0, -336.0, 0.0, 64.0];
        const INFLOOP_7_HCF: f64 = -2835.0;
        const INFLOOP_8: [f64; 9] = [7.0, 0.0, -280.0, 0.0, 1120.0, 0.0, -896.0, 0.0, 128.0];
        const INFLOOP_8_HCF: f64 = 14175.0;
        const INFLOOP_9: [f64; 10] = [
            0.0, 63.0, 0.0, -840.0, 0.0, 2016.0, 0.0, -1152.0, 0.0, 128.0,
        ];
        const INFLOOP_9_HCF: f64 = -155925.0;
        const INFLOOP_10: [f64; 11] = [
            -21.0, 0.0, 1260.0, 0.0, -8400.0, 0.0, 13440.0, 0.0, 5760.0, 0.0, 512.0,
        ];
        const INFLOOP_10_HCF: f64 = 467775.0;
        const INFLOOP_11: [f64; 12] = [
            0.0, -231.0, 0.0, 4620.0, 0.0, -18480.0, 0.0, 21120.0, 0.0, -7040.0, 0.0, 512.0,
        ];
        const INFLOOP_11_HCF: f64 = -6081075.0;
        const INFLOOP_12: [f64; 13] = [
            33.0, 0.0, -2772.0, 0.0, 27720.0, 0.0, -73920.0, 0.0, 63360.0, 0.0, -16896.0, 0.0,
            1024.0,
        ];
        const INFLOOP_12_HCF: f64 = 42567525.0;
        const INFLOOP_13: [f64; 14] = [
            0.0, 429.0, 0.0, -12012.0, 0.0, 72072.0, 0.0, -137280.0, 0.0, 91520.0, 0.0, -19968.0,
            0.0, 1024.0,
        ];
        const INFLOOP_13_HCF: f64 = -638512875.0;
        const INFLOOP_14: [f64; 15] = [
            -429.0, 0.0, 48048.0, 0.0, -672672.0, 0.0, 2690688.0, 0.0, -3843840.0, 0.0, 2050048.0,
            0.0, -372736.0, 0.0, 16384.0,
        ];
        const INFLOOP_14_HCF: f64 = 638512875.0;
        const INFLOOP_15: [f64; 16] = [
            0.0, 6435.0, 0.0, 240240.0, 0.0, -2018016.0, 0.0, 5765760.0, 0.0, -6406400.0, 0.0,
            2795520.0, 0.0, -430080.0, 0.0, 16384.0,
        ];
        const INFLOOP_15_HCF: f64 = -10854718875.0;
        pub fn get_loop_poly(n: u32) -> (f64, Vec<f64>) {
            match n {
                0 => (INFLOOP_0_HCF, INFLOOP_0.to_vec()),
                1 => (INFLOOP_1_HCF, INFLOOP_1.to_vec()),
                2 => (INFLOOP_2_HCF, INFLOOP_2.to_vec()),
                3 => (INFLOOP_3_HCF, INFLOOP_3.to_vec()),
                4 => (INFLOOP_4_HCF, INFLOOP_4.to_vec()),
                5 => (INFLOOP_5_HCF, INFLOOP_5.to_vec()),
                6 => (INFLOOP_6_HCF, INFLOOP_6.to_vec()),
                7 => (INFLOOP_7_HCF, INFLOOP_7.to_vec()),
                8 => (INFLOOP_8_HCF, INFLOOP_8.to_vec()),
                9 => (INFLOOP_9_HCF, INFLOOP_9.to_vec()),
                10 => (INFLOOP_10_HCF, INFLOOP_10.to_vec()),
                11 => (INFLOOP_11_HCF, INFLOOP_11.to_vec()),
                12 => (INFLOOP_12_HCF, INFLOOP_12.to_vec()),
                13 => (INFLOOP_13_HCF, INFLOOP_13.to_vec()),
                14 => (INFLOOP_14_HCF, INFLOOP_14.to_vec()),
                15 => (INFLOOP_15_HCF, INFLOOP_15.to_vec()),
                _ => (0.0, vec![0.0]),
            }
        }

        // define the polynomials for the P(x) olynomial defined in table 2 of reference paper.

        const ANNULAR_P_1: [f64; 2] = [0.0, 1.0];
        const ANNULAR_P_1_HCF: f64 = -1.0;
        const ANNULAR_P_2: [f64; 5] = [-3.0, 0.0, 4.0, 0.0, 1.0];
        const ANNULAR_P_2_HCF: f64 = 1.0;
        const ANNULAR_P_3: [f64; 6] = [0.0, -15.0, 0.0, 2.0, 0.0, 2.0];
        const ANNULAR_P_3_HCF: f64 = -3.0;
        const ANNULAR_P_4: [f64; 9] = [15.0, 0.0, -100.0, 0.0, -34.0, 0.0, 12.0, 0.0, 1.0];
        const ANNULAR_P_4_HCF: f64 = 6.0;
        const ANNULAR_P_5: [f64; 10] = [0.0, 280.0, 0.0, -455.0, 0.0, -412.0, 0.0, 16.0, 0.0, 8.0];
        const ANNULAR_P_5_HCF: f64 = -15.0;
        const ANNULAR_P_6: [f64; 13] = [
            -280.0, 0.0, 4515.0, 0.0, -924.0, 0.0, -4392.0, 0.0, -660.0, 0.0, 96.0, 0.0, 4.0,
        ];
        const ANNULAR_P_6_HCF: f64 = 30.0;
        const ANNULAR_P_7: [f64; 14] = [
            0.0, -2520.0, 0.0, 11214.0, 0.0, 6867.0, 0.0, -6872.0, 0.0, -2976.0, 0.0, 48.0, 0.0,
            16.0,
        ];
        const ANNULAR_P_7_HCF: f64 = -315.0;
        const ANNULAR_P_8: [f64; 17] = [
            630.0, 0.0, -18648.0, 0.0, 26775.0, 0.0, 45288.0, 0.0, -4756.0, 0.0, -11536.0, 0.0,
            -1032.0, 0.0, 80.0, 0.0, 2.0,
        ];
        const ANNULAR_P_8_HCF: f64 = 2520.0;
        const ANNULAR_P_9: [f64; 18] = [
            0.0, 88704.0, 0.0, -763840.0, 0.0, 42768.0, 0.0, 1372437.0, 0.0, 422968.0, 0.0,
            -226528.0, 0.0, -62592.0, 0.0, 512.0, 0.0, 128.0,
        ];
        const ANNULAR_P_9_HCF: f64 = -2835.0;
        const ANNULAR_P_10: [f64; 21] = [
            -88704.0,
            0.0,
            4176480.0,
            0.0,
            -13953720.0,
            0.0,
            -13785255.0,
            0.0,
            16604060.0,
            0.0,
            13983440.0,
            0.0,
            -587520.0,
            0.0,
            -1329920.0,
            0.0,
            -80320.0,
            0.0,
            3840.0,
            0.0,
            64.0,
        ];
        const ANNULAR_P_10_HCF: f64 = 5670.0;

        //impliment the polynomials given in table 3.
        //
        const ANNULAR_Q_1: [f64; 4] = [0.0, 2.0, 0.0, 1.0];
        const ANNULAR_Q_1_HCF: f64 = -1.0;
        const ANNULAR_Q_2: [f64; 5] = [-3.0, 0.0, 2.0, 0.0, 2.0];
        const ANNULAR_Q_2_HCF: f64 = 1.0;
        const ANNULAR_Q_3: [f64; 8] = [0.0, -45.0, 0.0, -16.0, 0.0, 16.0, 0.0, 2.0];
        const ANNULAR_Q_3_HCF: f64 = -1.0;
        const ANNULAR_Q_4: [f64; 9] = [30.0, 0.0, -185.0, 0.0, 172.0, 0.0, 16.0, 0.0, 8.0];
        const ANNULAR_Q_4_HCF: f64 = 3.0;
        const ANNULAR_Q_5: [f64; 12] = [
            0.0, 1400.0, 0.0, -1575.0, 0.0, -3372.0, 0.0, -576.0, 0.0, 144.0, 0.0, 8.0,
        ];
        const ANNULAR_Q_5_HCF: f64 = -3.0;
        const ANNULAR_Q_6: [f64; 13] = [
            -560.0, 0.0, 8750.0, 0.0, 2737.0, 0.0, -10872.0, 0.0, -4896.0, 0.0, 144.0, 0.0, 48.0,
        ];
        const ANNULAR_Q_6_HCF: f64 = 15.0;
        const ANNULAR_Q_7: [f64; 16] = [
            0.0, -17640.0, 0.0, 69678.0, 0.0, 89523.0, 0.0, -34984.0, 0.0, -45296.0, 0.0, -4608.0,
            0.0, 512.0, 0.0, 16.0,
        ];
        const ANNULAR_Q_7_HCF: f64 = -45.0;
        const ANNULAR_Q_8: [f64; 17] = [
            5040.0, 0.0, -146664.0, 0.0, 138978.0, 0.0, 488367.0, 0.0, 106808.0, 0.0, -137248.0,
            0.0, -39552.0, 0.0, 512.0, 0.0, 128.0,
        ];
        const ANNULAR_Q_8_HCF: f64 = 315.0;
        const ANNULAR_Q_9: [f64; 20] = [
            0.0,
            798336.0,
            0.0,
            -6475392.0,
            0.0,
            -3152160.0,
            0.0,
            -13453605.0,
            0.0,
            9473720.0,
            0.0,
            -1364960.0,
            0.0,
            -1505920.0,
            0.0,
            -102400.0,
            0.0,
            6400.0,
            0.0,
            128.0,
        ];
        const ANNULAR_Q_9_HCF: f64 = -315.0;
        const ANNULAR_Q_10: [f64; 21] = [
            -177408.0,
            0.0,
            8264256.0,
            0.0,
            -23708784.0,
            0.0,
            -42579438.0,
            0.0,
            23440285.0,
            0.0,
            45941900.0,
            0.0,
            8256400.0,
            0.0,
            -4528000.0,
            0.0,
            -937600.0,
            0.0,
            6400.0,
            0.0,
            1280.0,
        ];
        const ANNULAR_Q_10_HCF: f64 = 2835.0;

        pub fn get_annular_poly(n: u32) -> (f64, Vec<f64>, f64, Vec<f64>) {
            match n {
                1 => (
                    ANNULAR_P_1_HCF,
                    ANNULAR_P_1.to_vec(),
                    ANNULAR_Q_1_HCF,
                    ANNULAR_Q_1.to_vec(),
                ),
                2 => (
                    ANNULAR_P_2_HCF,
                    ANNULAR_P_2.to_vec(),
                    ANNULAR_Q_2_HCF,
                    ANNULAR_Q_2.to_vec(),
                ),
                3 => (
                    ANNULAR_P_3_HCF,
                    ANNULAR_P_3.to_vec(),
                    ANNULAR_Q_3_HCF,
                    ANNULAR_Q_3.to_vec(),
                ),
                4 => (
                    ANNULAR_P_4_HCF,
                    ANNULAR_P_4.to_vec(),
                    ANNULAR_Q_4_HCF,
                    ANNULAR_Q_4.to_vec(),
                ),
                5 => (
                    ANNULAR_P_5_HCF,
                    ANNULAR_P_5.to_vec(),
                    ANNULAR_Q_5_HCF,
                    ANNULAR_Q_5.to_vec(),
                ),
                6 => (
                    ANNULAR_P_6_HCF,
                    ANNULAR_P_6.to_vec(),
                    ANNULAR_Q_6_HCF,
                    ANNULAR_Q_6.to_vec(),
                ),
                7 => (
                    ANNULAR_P_7_HCF,
                    ANNULAR_P_7.to_vec(),
                    ANNULAR_Q_7_HCF,
                    ANNULAR_Q_7.to_vec(),
                ),
                8 => (
                    ANNULAR_P_8_HCF,
                    ANNULAR_P_8.to_vec(),
                    ANNULAR_Q_8_HCF,
                    ANNULAR_Q_8.to_vec(),
                ),
                9 => (
                    ANNULAR_P_9_HCF,
                    ANNULAR_P_9.to_vec(),
                    ANNULAR_Q_9_HCF,
                    ANNULAR_Q_9.to_vec(),
                ),
                10 => (
                    ANNULAR_P_10_HCF,
                    ANNULAR_P_10.to_vec(),
                    ANNULAR_Q_10_HCF,
                    ANNULAR_Q_10.to_vec(),
                ),
                _ => (0.0, vec![0.0], 0.0, vec![0.0]),
            }
        }

        #[cfg(test)]
        mod test {
            use super::*;
            #[test]
            fn test_hcf_0() {
                let hcf = 1.0;
                assert_eq!(get_loop_poly(0).0, hcf);
            }
            #[test] //test all the hcf components
            fn test_hcf_1() {
                let hcfs = vec![
                    1.0,
                    -3.0,
                    3.0,
                    -15.0,
                    45.0,
                    -315.0,
                    315.0,
                    -2835.0,
                    14175.0,
                    -155925.0,
                    467775.0,
                    -6081075.0,
                    42567525.0,
                    -638512875.0,
                    638512875.0,
                    -10854718875.0,
                ];
                let mut constant_hcf = Vec::new();
                for i in 0..16 {
                    constant_hcf.push(get_loop_poly(i).0);
                }

                assert_eq!(constant_hcf, hcfs);
            }
        }
    }
}
