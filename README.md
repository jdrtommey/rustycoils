# RustyCoils
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/Test/badge.svg)](https://github.com/jdrtommey/solenoid/actions)
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/LintFormat/badge.svg)](https://github.com/jdrtommey/solenoid/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This crate is an implementation of the numerical method for finding the off-axis magnetic field due to current systems with cyclindrical symmetry, i.e. coils and solenoids, presented in https://ieeexplore.ieee.org/document/760416 by R.H. Jackson. 
Analytical solutions exist for simple systems based on elliptical integrals, however these are numerically slower to implement. Further, provided that the region of interest is sufficiently close to the symmetry axis, no more accurate. Quoting the paper, this algorithm is designed to be "fast enough to be used interactively and compact enough to be embedded as a module in other codes".

This situation arises in simulating beams of atoms in Rydberg states in coherent quantum optics experiments moving in homogenous magentic fields provided by either solenoids or systems of coils. At each time-step the magnetic field needs to be calculated for a distribution of atoms close to the center of the magnetic field. The small deviations due to a physical implimentation can lead to small forces or Zeeman energy shifts. 
Similar experiments such as Zeeman slowers purposely make use of inhomogenous magnetic fields to impart forces. I thought this could be useful for someone else as an alternative to the elliptical integral methods I could find. 

This package is written in Rust to take advantage of its speed and (eventually) parallelisation abilities, but really it's because I wanted an excuse to use Rust for a project, and I plan on providing a Python wrapper using pyo3 in a seperated repository.
 

The algorithm makes use of primitive shapes with cylindrical symmetry:

* Ideal current loops
* Annular disks
* Thin solenoids 
* Thick solenoids (Coils)

# Warning

I put this together very quickly and it is still very much a prototype. I've tested the case of an ideal wire loop which has a simple analytical solution
and this appears to be working (test script in tests/idealloop.rs) but have not tested the other primitives yet. If anyone finds, and try and use, this package in its current form please check the solutions against another method before attempting to use it for anything that matters.


# Usage 

The crate exposes a single struct 
```rust

//object containing primitives sharing a symmetry axis
mycoil = AxialSystem::default(); 
```
which defines a symmetry axis. Currently this symmetry axis can only be defined along the three cartesian axes (defaults to x). Eventually this will be arbitrary. 

```rust

mycoil.transform_x();
mycoil.transform_y();
mycoil.transform_z();

```

Individual primitive coils can be added to the AxialSystem with a unique UTF-8 identifier.

```rust
//define physical parameters
let radius = 1.0;
let thickness = 0.1;
let current = 1.0;
let length = 5.0;
let position = 2.0; //position along the symmetry axis
mycoil.add_loop("loop1",radius,position,current);
mycoil.add_annular("foo",radius,thickness,position,current);
mycoil.add_solenoid("bar",radius,length,position,length);
mycoil.add_coil("coil1",radius,length,thickness,position,current);
```
The parameters controlling these primitives can be be modified by using the functions 
```rust
//change radius of the current loop
mycoil.modify_radius("loop1",6.0);  "loop1"
//change length of the solenoid "bar"
mycoil.modify_length("bar",3.0); 
//change length of the coil "coil1"
mycoil.modify_position("coil1",3.0); 
//change thickness of the annular "foo"
mycoil.modify_thickness("foo",1.0); 
//change current of the annular "foo"
mycoil.modify_current("foo",1.0); 
```

These functions accept keywords to modify multiple primitives at once. Note these keywords can not be used as identifiers for primitives.

| Reserved word  | Meaning |
| -------------  | ------------- |
| *   | Apply to all  |
| LOOP   | Apply to current loops  |
| ANNULAR   | Apply to annulars  |
| SOLENOID   | Apply to solenoids  |
| COIL   | Apply to coils  |


```rust
//changes all the current of all primitives
mycoil.modify_current("*",6.0); 
//changes all the current of all current loop primitives
mycoil.modify_current("LOOP",6.0); 
//changes all the radius of all annular primitives
mycoil.modify_radius("ANNULAR",6.0);  
//changes all the length of all solenoid primitives
mycoil.modify_length("SOLENOID",6.0);
//changes all the length of all coil primitives
mycoil.modify_thickness("COIL",6.0);  
```

The magnetic field in each of the cartesian directions can be computed from 
```rust
(mag_x,mag_y,mag_z) = mycoil.get_field([x,y,z],1e-10);
```
where 1e-10 is the tolerance to stop including additional terms in the power expansion.

# TO-DO

Things I plan on implementing soon

- Tests of Annular and Coil primitives
- Arbitrary orientation of symmetry axes
- Parallelise the field computation using Rayon
