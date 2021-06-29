### FiniteSolenoid

This small package impliment the numerical method of R.Jackson in https://ieeexplore.ieee.org/document/760416.
It is intended for solving the off axis magnetic field due to systems with cyclindrical symmetry, i.e. coils and solenoids.
Analytic solutions exist for such systems based on elliptical functions however these are numerically slower to impliment. 
For cases where the function will be called many times within a larger function a faster solution is desired.
This package is written in Rust to take advantage of its speed and (eventually) concurrency, but essientially because I wanted an excuse to write a rust code with a python wrapper. 
(mainly because I needed an excuse to use Rust) and provides a Python api.
The method employed makes use of the partial power series decomposition of the magneto-scalar potential.
It is designed to be 'fast enough to be used interactivley and compact enough to be embedded as a module in other codes'. 

The algorithim makes use of primitive shapes with cylindrical symmetry:
* Circular loops
* Annular disks
* Thin solenoid 
* Full Coil

These primitives can be added to the publically exposed AxialObject struct as 

