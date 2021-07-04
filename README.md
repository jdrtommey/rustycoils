# RustyCoils
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/Test/badge.svg)](https://github.com/jdrtommey/solenoid/actions)
[![Actions Status](https://github.com/jdrtommey/solenoid/workflows/LintFormat/badge.svg)](https://github.com/jdrtommey/solenoid/actions)

Implimentation of a numerical method for finding the off-axis magnetic field due to current systems with cyclindrical symmetry, i.e. coils and solenoids. The exact method used is shown in https://ieeexplore.ieee.org/document/760416. 
Analytical solutions exist for simple systems based on elliptical functions however these are numerically slower to impliment, and provided the region on interest is sufficiently close to the middle of the symmetry axis are no more accurate.
For cases where the function will be called many times within a larger function the quickest possible implimentation is desirable.
This package is written in Rust to take advantage of its speed and (eventually) concurrency, but really it's because I wanted an excuse to use Rust for a project.
There is (eventually) a python wrapper using Py03 in a seperate repository.
To quote the paper this method is taken from the algorithim is designed to be 'fast enough to be used interactivley and compact enough to be embedded as a module in other codes'. 

The algorithim makes use of primitive shapes with cylindrical symmetry:
* Ideal current loops
* Annular disks
* Thin solenoids 
* Thick solenoids (Coils)

## Warning

I put this together very quickly and it is still very much a prototype. Ive tested the case of an ideal wire loop which has a simple analytical solution
and this appears to be working (test script in tests/idealloop.rs) but have not tested the other primitives yet, if anyone finds, and try and use, this package in its current form please check the solutions against another method before attempting to use it for anything that matters.


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
(mag_x,mag_y,mag_z) = mycoil.get_field((x,y,z),1e-10);
```
where 1e-10 is the tolerance to stop including additional terms in the power expansion.

##### TO-DO

- Arbitrary symmetry direction axes
- Parallel the computation within AxialObjects and on collections of them.
