//! # RustyCoils
//!
//! Off-axis magnetic fields for systems with cyclindrical symmetry.  
//!
//! This crate impliments a fast approximate method for finding off axis magnetic fields
//! of systems with cylindrical symmetry. Using a power series decomposition of the on-axis
//! magentic field the full magentic field for "basic" primitive shapes can be determined.
//! These can then be combined within an object with a given orientation and location in order
//! to build up a larger system.
//!
//! The method used in this crate comes directly from "Off-Axis Expansion Solution of Laplace's
//! Equation: Application to Accurate and Rapid Calculation of Coil Magentic Fields" by Robert H.
//! Jackson. The author sums up the utility of this method with the sentance "The simplicity,
//! compactness and speed of this method make it a good adjunct to other techniques and ideal as
//! a module for incorporation into more general programs". Near the axis of symmetry the method
//! can give very accurate fields at a much faster speed than other methods (Author of paper
//! determines that out to about 70% the error is less than 0.1% of the exact elliptical integral
//! solution for an ideal current loop.

#[warn(missing_docs)]
mod axialobject;
mod fieldcalc;
#[cfg(feature = "rayon")]
pub mod rayon;
pub use axialobject::{AxialError, AxialSystem};
