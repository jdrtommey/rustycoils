# RustyCoils
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/Test/badge.svg)](https://github.com/jdrtommey/rustycoils/actions)
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/LintFormat/badge.svg)](https://github.com/jdrtommey/rustycoils/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This crate implements a numerical method for determining the off-axis magnetic field of current systems with cyclindrical symmetry, i.e., solenoids and coils. The underlying algorithm is taken from [Off-axis expansion solution of Laplace's equation: Application to accurate and rapid calculation of coil magnetic fields](https://ieeexplore.ieee.org/document/760416) by R.H. Jackson. This method makes use of primitive coil/loop shapes for which the on-axis magnetic field (and all their derivatives) are analytically known. A power series of these derivatives can be used to determine the full off-axis magnetic field. Provided the position of interest is near enough the axis this provides very rapid and accurate computation of the magnetic field. Quoting the paper, this algorithm is designed to be "fast enough to be used interactively and compact enough to be embedded as a module in other codes".

Fully analytic solutions for the infinitely thin wire loop exist based on elliptical integrals and can be combined to compute the field of more complicated systems of coils. This can become slow when multiple wire loops are being used to model a solenoid or set of coils. In particular when the region of interest is close to the axis. For this scenario the approach implemented here which includes such solenoids as a basic primitive is useful. For example, this arises in atom optics experiments such as simulating beams of atoms in Rydberg states moving in (in)homogenous magentic fields provided by either solenoids or systems of coils. At each time-step the magnetic field needs to be calculated for a distribution of atoms close to the axis of the system. The small deviations due to a physical implementation can lead to small forces or Zeeman energy shifts.

This package is written in Rust to take advantage of its speed and parallelisation abilities, but really it's because I wanted an excuse to use Rust for a project, and I wrote a Python wrapper using pyo3 in a seperate repository [rustpycoils](https://github.com/jdrtommey/rustpycoils).

The algorithm makes use of primitive shapes with cylindrical symmetry:

* Ideal current loops
* Annular disks
* Thin solenoids 
* Thick solenoids (Coils)

# Warning
 
I put this together very quickly and it is still very much a prototype. I've tested the case of an ideal wire loop which has a simple analytical solution
and this appears to be working (test script in [idealloop.rs](./tests/idealloop.rs). But I have not tested the other primitives yet, beyond checking the under certain limits they give the same answer from combining many current loops. If anyone tries to use this package in its current form please check the solutions against another method before attempting to use it for anything that matters.


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
mycoil.modify_radius("loop1",6.0); 
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
[mag_x,mag_y,mag_z] = mycoil.get_field([x,y,z],1e-18);
```
where 1e-18 is the tolerance to stop including additional terms in the power expansion.

# TO-DO

Things I plan on implementing soon

- Arbitrary orientation of symmetry axes
- Combine multiple AxialSystems along different axes
