### FiniteSolenoid

This small package impliment the numerical method of R.Jackson in https://ieeexplore.ieee.org/document/760416.
It is intended for solving the off axis magnetic field due to systems with cyclindrical symmetry, i.e. coils and solenoids.
Analytic solutions exist for such systems based on elliptical functions however these are numerically slower to impliment. 
For cases where the function will be called many times within a larger function a faster solution is desired.
This package is written in Rust to take advantage of its speed and concurrency (mainly because I needed an excuse to use Rust) and provides a Python api.
The method imployed makes use of partial power series decomposition of the magneto-scalar potential.
It is designed to be 'fast enough to be used interactivley and compact enough to be embedded as a module in other codes'. 

The algorithim makes use of primitive shapes with cylindrical symmetry:
     - Circular loops
     - Annular disks
     - Thin solenoid 
     - Full Coil

Need accurate computation of the on-axis potential derivatives, and accurate summation of the series.
This is solved by taking analytical high order derivatives of the on-axis field equations for the primitive shapes. 

On-axis derivatives of the ideal loop. The axial field for a single turn, infinitely thin current carrying loop can be found in many places,

$$
B_z(0,z) = B_z(z) = \frac{\mu_0 I}{2} \frac{r_0^2}{(r_0^2 +(z-z_0)^2)^{3/2}}
$$
where r_0 and z_0 are the radius and axial position of the loop, I is the current and $\mu_0$ is the permeability of free space. This equation can be renormalized as 
$$
B_z(z) = B_0 b(x) = \frac{\mu_0 I}{2r_0} \frac{1}{(1+x^2)^{3/2}
$$
where $B_0 = \frac{\mu_0 I}{2 r_0}$, $b(x) = \frac{1}{(1+x^2)^(3/2)}$,$ x = \frac{z-z0}{r_0}$.

This splits the equation into a component containing the physical parameters and one having the functional dependance. Differention can now be performed using the "natural" scale of the field without reference to the physical parameters of the loop which can be ignored until the end of the precess. 

The normalized magnetic field equation translates very naturally into the expansion solution sums where each term normalises as $\frac{d^k}{dz^k} \rightarrow \frac{1}{r_0^k}\frac{d^k}{dx^k}$.
The expansion for any loop field reduces to generating an expansion solution for the normalised $b(x)$ and multiplying by $B_0$. 
The $n^{\mathrm{th}}$ derivative of b(x) can be put in the general form:
$$
b^(n)(x) = \frac{P_n(x)}{\sqrt{1+x^2}(1+x^2)^{n+1}}
$$
A recursion relationship exists between the polynomials $Pn$
