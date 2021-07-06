use crate::fieldcalc::primitives::{CoilSolenoid, IdealWire, Primitive, ThinAnnular, ThinSolenoid};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;

//public api struct. Groups a set of primitive types which
//share a common symmetry axis. Primitives stored in a HashMap<string,primitive>
//and are indexed and interacted with via string based identifiers. The whole object
//has a single origin and orientation and the individual primitives are placed relative
//to this origin.In V 0.1 will only provide x,y,z centered at origin
//There are currently three main function types for the AxialObject
//Add_primative("id",radius,etc...) these functions add a new primative to the object,
//Modify_physical_parameter("id",value) these functins allow for propertires of the
//primative to be changed, if "*" passed as id it will modify all primatives which
//have that field, if "COIL" passed will modify all COILS etc.
//get_field((x,y,z)) this function calls the magnetic field due to all primatives in
//the class and returns the sum.
//
#[derive(Debug, PartialEq, Clone)]
/// Structure which defines a system of primitives with a shared symmetry axis.
///
/// This allows for large number of primitive shapes to be combined to generate a more
/// complicated magentic field struture. The system is defined by an origin vector and
/// a orientation vector. Currently only orientations along the x,y,z axes and placed at
/// the origin of the global coordinate system are allowed.
///
/// individual primitives can be added to the AxialSystem and are stored in a HashMap with String
/// based keys. These keys allow individual primitives to be accessed and modifies. Functions exist
/// to modify individual physical parameters such as radius,length,thickness,position,current. The
/// magnetic field is computed currently by working out the individual magentic field of each
/// primitive individually. TO-DO: include rayon support for parallel compuation of the primitive
/// magnetic fields.
pub struct AxialSystem {
    objects: HashMap<String, Primitives>,
    origin: [f64; 3],
    orientation: [f64; 3],
}

impl fmt::Display for AxialSystem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Origin:({},{},{}),Orientation:({},{},{})",
            self.origin[0],
            self.origin[1],
            self.origin[2],
            self.orientation[0],
            self.orientation[1],
            self.orientation[2]
        )
    }
}
impl AxialSystem {
    /// Returns a new AxialSystem
    ///
    /// Currently only supports objects located at the origin
    /// with their orientation along any of the x,y,z directions
    /// and so internally calls the default method.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::AxialSystem;
    /// let axial = AxialSystem::new([0.0,0.0,0.0],[1.0,0.0,0.0]);
    /// ```
    pub fn new(origin: [f64; 3], orientation: [f64; 3]) -> AxialSystem {
        AxialSystem {
            objects: HashMap::new(),
            origin,
            orientation,
        }
    }
    /// Returns the default AxialSystem
    ///
    /// This has the shared symmetry axis located at the global origin (0,0,0)
    /// with its symmetry axis along the x axis (1,0,0).
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let axial = AxialSystem::default();
    /// ```
    pub fn default() -> AxialSystem {
        AxialSystem {
            objects: HashMap::new(),
            origin: [0.0, 0.0, 0.0],
            orientation: [1.0, 0.0, 0.0],
        }
    }
}
//Add functions
impl AxialSystem {
    /// Adds an instance of the ideal loop to the AxialSystem
    ///
    /// Provide a unique identifer for the primitive as a String
    /// method checks if the identifier is allowed due to either
    /// clashing with a reserved word or due to a primitive already
    /// sharing the name.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let result = axial.add_loop("loop1".to_string(),1.0,0.0,1.0);
    /// assert_eq!(result,Ok(()));
    /// let result2 = axial.add_loop("loop1".to_string(),2.0,1.0,1.0);
    /// assert_eq!(result2,Err(AxialError::KeyDuplicateError("loop1".to_string())));
    /// let result3 = axial.add_loop("LOOP".to_string(),2.0,1.0,1.0);
    /// assert_eq!(result3,Err(AxialError::ReservedWordError("LOOP".to_string())));
    /// ```
    pub fn add_loop(
        &mut self,
        id: String,
        radius: f64,
        origin: f64,
        current: f64,
    ) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyMissingError(_) => {}
            error => return Err(error),
        }
        let new_loop = Primitives::IdealWire(IdealWire::new(radius, current, origin));
        self.objects.insert(id, new_loop);
        Ok(())
    }
    /// Adds an instance of the annular primitive to the AxialSystem
    ///
    /// Provide a unique identifer for the primitive as a String
    /// method checks if the identifier is allowed due to either
    /// clashing with a reserved word or due to a primitive already
    /// sharing the name.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let result = axial.add_annular("annular1".to_string(),1.0,0.1,0.0,1.0);
    /// assert_eq!(result,Ok(()));
    /// let result2 = axial.add_annular("annular1".to_string(),2.0,0.1,1.0,1.0);
    /// assert_eq!(result2,Err(AxialError::KeyDuplicateError("annular1".to_string())));
    /// let result3 = axial.add_annular("ANNULAR".to_string(),2.0,0.1,1.0,1.0);
    /// assert_eq!(result3,Err(AxialError::ReservedWordError("ANNULAR".to_string())));
    /// ```
    pub fn add_annular(
        &mut self,
        id: String,
        radius: f64,
        thickness: f64,
        origin: f64,
        current: f64,
    ) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyMissingError(_) => {}
            error => return Err(error),
        }
        let new_annular = Primitives::Annular(ThinAnnular::new(radius, thickness, current, origin));
        self.objects.insert(id, new_annular);
        Ok(())
    }
    /// Adds an instance of the thin solenoid primitive to the AxialSystem
    ///
    /// Provide a unique identifer for the primitive as a String
    /// method checks if the identifier is allowed due to either
    /// clashing with a reserved word or due to a primitive already
    /// sharing the name.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let result = axial.add_thin_solenoid("solenoid1".to_string(),1.0,10.0,0.0,1.0);
    /// assert_eq!(result,Ok(()));
    /// let result2 = axial.add_thin_solenoid("solenoid1".to_string(),2.0,0.1,1.0,1.0);
    /// assert_eq!(result2,Err(AxialError::KeyDuplicateError("solenoid1".to_string())));
    /// let result3 = axial.add_thin_solenoid("SOLENOID".to_string(),2.0,0.1,1.0,1.0);
    /// assert_eq!(result3,Err(AxialError::ReservedWordError("SOLENOID".to_string())));
    /// ```
    pub fn add_thin_solenoid(
        &mut self,
        id: String,
        radius: f64,
        length: f64,
        origin: f64,
        current: f64,
    ) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyMissingError(_) => {}
            error => return Err(error),
        }
        let new_solenoid =
            Primitives::ThinSolenoid(ThinSolenoid::new(radius, length, current, origin));
        self.objects.insert(id, new_solenoid);
        Ok(())
    }
    /// Adds an instance of the coil solenoid primitive to the AxialSystem
    ///
    /// Provide a unique identifer for the primitive as a String
    /// method checks if the identifier is allowed due to either
    /// clashing with a reserved word or due to a primitive already
    /// sharing the name.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let result = axial.add_coil_solenoid("coil1".to_string(),1.0,10.0,0.1,0.0,1.0);
    /// assert_eq!(result,Ok(()));
    /// let result2 = axial.add_coil_solenoid("coil1".to_string(),2.0,10.0,0.1,1.0,1.0);
    /// assert_eq!(result2,Err(AxialError::KeyDuplicateError("coil1".to_string())));
    /// let result3 = axial.add_coil_solenoid("COIL".to_string(),2.0,10.0,0.1,1.0,1.0);
    /// assert_eq!(result3,Err(AxialError::ReservedWordError("COIL".to_string())));
    /// ```
    pub fn add_coil_solenoid(
        &mut self,
        id: String,
        radius: f64,
        length: f64,
        thickness: f64,
        origin: f64,
        current: f64,
    ) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyMissingError(_) => {}
            error => return Err(error),
        }
        let new_coil = Primitives::CoilSolenoid(CoilSolenoid::new(
            radius, length, thickness, current, origin,
        ));
        self.objects.insert(id, new_coil);
        Ok(())
    }
    /// Removes the primitive matching the provided id.
    ///
    /// # Arguments
    ///
    /// * `id` - &str containing the ID to remove
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    ///
    /// let result = axial.add_loop("loop1".to_string(),1.0,0.0,1.0);
    /// let result_wrong_id = axial.remove("loop2");
    /// assert_eq!(result_wrong_id,Err(AxialError::KeyMissingError("loop2".to_string())));
    /// ```
    pub fn remove(&mut self, id: &str) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyDuplicateError(_) => {}
            error => return Err(error),
        }
        self.objects.remove(id);
        Ok(())
    }
}

//view functions
impl AxialSystem {
    /// Returns the display string of the primitive.
    ///
    /// # Arguments
    ///
    /// * `id` - &str containing the ID to view
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    ///
    /// let result = axial.add_loop("loop1".to_string(),1.0,0.0,1.0);
    /// let result_wrong_id = axial.view("loop1");
    /// ```
    pub fn view(&self, id: &str) -> Result<String, AxialError> {
        match _is_id_valid(&self.objects, id) {
            AxialError::KeyDuplicateError(_) => {}
            error => return Err(error),
        }
        let primitive = self.objects.get(id);
        if let Some(primitive) = primitive {
            Ok((*primitive).to_string())
        } else {
            Err(AxialError::KeyMissingError(id.to_string()))
        }
    }
}
#[cfg(test)]
mod test_view {
    use super::*;
    #[test]
    fn test_view_correct_id() -> Result<(), AxialError> {
        let mut mycoil = AxialSystem::default();
        let _res = mycoil.add_loop("loop1".to_string(), 0.0, 0.0, 0.0);
        let string = mycoil.view("loop1")?;
        assert_eq!(
            Primitives::IdealWire(IdealWire::new(0.0, 0.0, 0.0)).to_string(),
            string
        );
        Ok(())
    }
    #[test]
    fn test_view_wrong_id() -> Result<(), AxialError> {
        let mut mycoil = AxialSystem::default();
        let _res = mycoil.add_loop("loop1".to_string(), 0.0, 0.0, 0.0);
        let answer = mycoil.view("loop2");
        assert_eq!(
            answer,
            Err(AxialError::KeyMissingError("loop2".to_string()))
        );
        Ok(())
    }
}
//transform functions
impl AxialSystem {
    /// Transforms the symmetry axis to point along the x axis
    ///
    /// i.e. converts the AxialSystems orientation to (1,0,0)
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// axial.transform_x();
    /// ```
    pub fn transform_x(&mut self) {
        self.orientation = [1.0, 0.0, 0.0];
    }
    /// Transforms the symmetry axis to point along the y axis
    ///
    /// i.e. converts the AxialSystems orientation to (0,1,0)
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// axial.transform_y();
    /// ```
    pub fn transform_y(&mut self) {
        self.orientation = [0.0, 1.0, 0.0];
    }
    /// Transforms the symmetry axis to point along the z axis
    ///
    /// i.e. converts the AxialSystems orientation to (0,0,1)
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// axial.transform_z();
    /// ```
    pub fn transform_z(&mut self) {
        self.orientation = [0.0, 0.0, 1.0];
    }
}

//modify functions
//these functions allow the individual physical parameters of
//underlying primitives to be modified.

impl AxialSystem {
    /// Modifies the radius of a given primitive/ set of primitives
    ///
    /// Can provide the ID of a single primitive or provide one
    /// of the possible reserved keywords to modify a set of primitives
    ///
    /// *`LOOP` change radius of all loop primitives
    /// *`ANNULAR` change radius of all annular primitives
    /// *`SOLENOID` change radius of all solenoid primitives
    /// *`COILS` change radius of all coil primitives
    ///
    /// # Arguments
    ///
    /// * `id` contains the ID of the primitive to modify
    /// * `radius` new radius of the primitive
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,0.0);
    /// let res = axial.modify_radius("loop1",2.0);
    /// assert_eq!(res,Ok(()));
    /// let res = axial.modify_radius("loop2",2.0);
    /// assert_eq!(res,Err(AxialError::KeyMissingError("loop2".to_string())));
    /// ```
    pub fn modify_radius(&mut self, id: &str, radius: f64) -> Result<(), AxialError> {
        //generate a list of the objects to be modified, and run through and modify their
        //radius via the set_radius() function

        if id == "*" {
            for (_, object) in self.objects.iter_mut() {
                match object {
                    Primitives::IdealWire(primitive) => primitive.set_radius(radius),
                    Primitives::Annular(primitive) => primitive.set_radius(radius),
                    Primitives::ThinSolenoid(primitive) => primitive.set_radius(radius),
                    Primitives::CoilSolenoid(primitive) => primitive.set_radius(radius),
                }
            }
        } else if id == "LOOP" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::IdealWire(primitive) = object {
                    primitive.set_radius(radius)
                }
            }
        } else if id == "ANNULAR" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::Annular(primitive) = object {
                    primitive.set_radius(radius)
                }
            }
        } else if id == "SOLENOID" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::ThinSolenoid(primitive) = object {
                    primitive.set_radius(radius)
                }
            }
        } else if id == "COIL" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::CoilSolenoid(primitive) = object {
                    primitive.set_radius(radius)
                }
            }
        } else {
            let primative = self.objects.get_mut(id);
            match primative {
                Some(primative) => match primative {
                    Primitives::IdealWire(primitive) => primitive.set_radius(radius),
                    Primitives::Annular(primitive) => primitive.set_radius(radius),
                    Primitives::ThinSolenoid(primitive) => primitive.set_radius(radius),
                    Primitives::CoilSolenoid(primitive) => primitive.set_radius(radius),
                },
                None => return Err(AxialError::KeyMissingError(id.to_string())),
            }
        }
        Ok(())
    }

    /// Modifies the position of a given primitive/set of primitives relative to the AxialSystem
    /// along the symmetry axis
    ///
    /// Can provide the ID of a single primitive or provide one
    /// of the possible reserved keywords to modify a set of primitives
    ///
    /// *`LOOP` change radius of all loop primitives
    /// *`ANNULAR` change radius of all annular primitives
    /// *`SOLENOID` change radius of all solenoid primitives
    /// *`COILS` change radius of all coil primitives
    ///
    /// # Arguments
    ///
    /// * `id` contains the ID of the primitive to modify
    /// * `position` new position of the primitive
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,0.0);
    /// let res = axial.modify_position("loop1",1.0);
    /// assert_eq!(res,Ok(()));
    /// ```
    pub fn modify_position(&mut self, id: &str, position: f64) -> Result<(), AxialError> {
        //generate a list of the objects to be modified, and run through and modify their
        //radius via the set_radius() function

        if id == "*" {
            for (_, object) in self.objects.iter_mut() {
                match object {
                    Primitives::IdealWire(primitive) => primitive.set_z0(position),
                    Primitives::Annular(primitive) => primitive.set_z0(position),
                    Primitives::ThinSolenoid(primitive) => primitive.set_z0(position),
                    Primitives::CoilSolenoid(primitive) => primitive.set_z0(position),
                }
            }
        } else if id == "LOOP" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::IdealWire(primitive) = object {
                    primitive.set_z0(position)
                }
            }
        } else if id == "ANNULAR" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::Annular(primitive) = object {
                    primitive.set_z0(position)
                }
            }
        } else if id == "SOLENOID" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::ThinSolenoid(primitive) = object {
                    primitive.set_z0(position)
                }
            }
        } else if id == "COIL" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::CoilSolenoid(primitive) = object {
                    primitive.set_z0(position)
                }
            }
        } else {
            let primative = self.objects.get_mut(id);
            match primative {
                Some(primative) => match primative {
                    Primitives::IdealWire(primitive) => primitive.set_z0(position),
                    Primitives::Annular(primitive) => primitive.set_z0(position),
                    Primitives::ThinSolenoid(primitive) => primitive.set_z0(position),
                    Primitives::CoilSolenoid(primitive) => primitive.set_z0(position),
                },
                None => return Err(AxialError::KeyMissingError(id.to_string())),
            }
        }
        Ok(())
    }

    /// Modifies the length of a given primitive/set of primitives
    ///
    /// Can provide the ID of a single primitive or provide one
    /// of the possible reserved keywords to modify a set of primitives
    /// If the provided ID belongs to a primitive that does not possess
    /// length as a parameter returns a AxialError::IncompatiblePrimitiveError
    ///
    /// *`LOOP` change radius of all loop primitives
    /// *`ANNULAR` change radius of all annular primitives
    /// *`SOLENOID` change radius of all solenoid primitives
    /// *`COILS` change radius of all coil primitives
    ///
    /// # Arguments
    ///
    /// * `id` contains the ID of the primitive to modify
    /// * `length` new position of the primitive
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_thin_solenoid("solenoid1".to_string(),1.0,10.0,0.0,0.0);
    /// let res = axial.modify_length("solenoid1",5.0);
    /// assert_eq!(res,Ok(()));
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,0.0);
    /// let res = axial.modify_length("loop1",5.0);
    /// assert_eq!(res,Err(AxialError::IncompatiblePrimitiveError("loop1".to_string(),"LOOP".to_string())));
    /// ```
    pub fn modify_length(&mut self, id: &str, length: f64) -> Result<(), AxialError> {
        //generate a list of the objects to be modified, and run through and modify their
        //radius via the set_radius() function

        if id == "*" {
            for (_, object) in self.objects.iter_mut() {
                match object {
                    Primitives::ThinSolenoid(primitive) => primitive.set_length(length),
                    Primitives::CoilSolenoid(primitive) => primitive.set_length(length),
                    _ => {}
                }
            }
        } else if id == "LOOP" {
            return Err(AxialError::IncompatiblePrimitiveError(
                id.to_string(),
                "LOOP".to_string(),
            ));
        } else if id == "ANNULAR" {
            return Err(AxialError::IncompatiblePrimitiveError(
                id.to_string(),
                "ANNULAR".to_string(),
            ));
        } else if id == "SOLENOID" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::ThinSolenoid(primitive) = object {
                    primitive.set_length(length)
                }
            }
        } else if id == "COIL" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::CoilSolenoid(primitive) = object {
                    primitive.set_length(length)
                }
            }
        } else {
            let primative = self.objects.get_mut(id);
            match primative {
                Some(primative) => match primative {
                    Primitives::IdealWire(_) => {
                        return Err(AxialError::IncompatiblePrimitiveError(
                            id.to_string(),
                            "LOOP".to_string(),
                        ))
                    }
                    Primitives::Annular(_) => {
                        return Err(AxialError::IncompatiblePrimitiveError(
                            id.to_string(),
                            "ANNULAR".to_string(),
                        ))
                    }
                    Primitives::ThinSolenoid(primitive) => primitive.set_length(length),
                    Primitives::CoilSolenoid(primitive) => primitive.set_length(length),
                },
                None => return Err(AxialError::KeyMissingError(id.to_string())),
            }
        }
        Ok(())
    }
    /// Modifies the thickness of a given primitive/set of primitives
    ///
    /// Can provide the ID of a single primitive or provide one
    /// of the possible reserved keywords to modify a set of primitives
    /// If the provided ID belongs to a primitive that does not possess
    /// length as a parameter returns a AxialError::IncompatiblePrimitiveError
    ///
    /// *`LOOP` change radius of all loop primitives
    /// *`ANNULAR` change radius of all annular primitives
    /// *`SOLENOID` change radius of all solenoid primitives
    /// *`COILS` change radius of all coil primitives
    ///
    /// # Arguments
    ///
    /// * `id` contains the ID of the primitive to modify
    /// * `thickness` new position of the primitive
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_annular("annular1".to_string(),1.0,1.0,0.0,0.0);
    /// let res = axial.modify_thickness("annular1",5.0);
    /// assert_eq!(res,Ok(()));
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,0.0);
    /// let res = axial.modify_thickness("loop1",5.0);
    /// assert_eq!(res,Err(AxialError::IncompatiblePrimitiveError("loop1".to_string(),"LOOP".to_string())));
    /// ```
    pub fn modify_thickness(&mut self, id: &str, thickness: f64) -> Result<(), AxialError> {
        //generate a list of the objects to be modified, and run through and modify their
        //radius via the set_radius() function

        if id == "*" {
            for (_, object) in self.objects.iter_mut() {
                match object {
                    Primitives::Annular(primitive) => primitive.set_thickness(thickness),
                    Primitives::CoilSolenoid(primitive) => primitive.set_thickness(thickness),
                    _ => {}
                }
            }
        } else if id == "LOOP" {
            return Err(AxialError::IncompatiblePrimitiveError(
                id.to_string(),
                "LOOP".to_string(),
            ));
        } else if id == "SOLENOID" {
            return Err(AxialError::IncompatiblePrimitiveError(
                id.to_string(),
                "SOLENOID".to_string(),
            ));
        } else if id == "ANNULAR" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::Annular(primitive) = object {
                    primitive.set_thickness(thickness)
                }
            }
        } else if id == "COIL" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::CoilSolenoid(primitive) = object {
                    primitive.set_thickness(thickness)
                }
            }
        } else {
            let primative = self.objects.get_mut(id);
            match primative {
                Some(primative) => match primative {
                    Primitives::IdealWire(_) => {
                        return Err(AxialError::IncompatiblePrimitiveError(
                            id.to_string(),
                            "LOOP".to_string(),
                        ))
                    }
                    Primitives::ThinSolenoid(_) => {
                        return Err(AxialError::IncompatiblePrimitiveError(
                            id.to_string(),
                            "SOLENOID".to_string(),
                        ))
                    }
                    Primitives::Annular(primitive) => primitive.set_thickness(thickness),
                    Primitives::CoilSolenoid(primitive) => primitive.set_thickness(thickness),
                },
                None => return Err(AxialError::KeyMissingError(id.to_string())),
            }
        }
        Ok(())
    }
    /// Modifies the current of a given primitive/set of primitives
    ///
    /// Can provide the ID of a single primitive or provide one
    /// of the possible reserved keywords to modify a set of primitives
    ///
    /// *`LOOP` change radius of all loop primitives
    /// *`ANNULAR` change radius of all annular primitives
    /// *`SOLENOID` change radius of all solenoid primitives
    /// *`COILS` change radius of all coil primitives
    ///
    /// # Arguments
    ///
    /// * `id` contains the ID of the primitive to modify
    /// * `thickness` new position of the primitive
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_annular("annular1".to_string(),1.0,1.0,0.0,0.0);
    /// let res = axial.modify_current("annular1",5.0);
    /// assert_eq!(res,Ok(()));
    /// ```
    pub fn modify_current(&mut self, id: &str, current: f64) -> Result<(), AxialError> {
        //generate a list of the objects to be modified, and run through and modify their
        //radius via the set_radius() function

        if id == "*" {
            for (_, object) in self.objects.iter_mut() {
                match object {
                    Primitives::IdealWire(primitive) => primitive.set_current(current),
                    Primitives::Annular(primitive) => primitive.set_current(current),
                    Primitives::ThinSolenoid(primitive) => primitive.set_current(current),
                    Primitives::CoilSolenoid(primitive) => primitive.set_current(current),
                }
            }
        } else if id == "LOOP" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::IdealWire(primitive) = object {
                    primitive.set_current(current)
                }
            }
        } else if id == "ANNULAR" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::Annular(primitive) = object {
                    primitive.set_current(current)
                }
            }
        } else if id == "SOLENOID" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::ThinSolenoid(primitive) = object {
                    primitive.set_current(current)
                }
            }
        } else if id == "COIL" {
            for (_, object) in self.objects.iter_mut() {
                if let Primitives::CoilSolenoid(primitive) = object {
                    primitive.set_current(current)
                }
            }
        } else {
            let primative = self.objects.get_mut(id);
            match primative {
                Some(primative) => match primative {
                    Primitives::IdealWire(primitive) => primitive.set_current(current),
                    Primitives::Annular(primitive) => primitive.set_current(current),
                    Primitives::ThinSolenoid(primitive) => primitive.set_current(current),
                    Primitives::CoilSolenoid(primitive) => primitive.set_current(current),
                },
                None => return Err(AxialError::KeyMissingError(id.to_string())),
            }
        }
        Ok(())
    }
}

impl AxialSystem {
    fn get_field_in_frame(&self, pos: [f64; 3], tol: &f64) -> [f64; 2] {
        let [z, r, _theta] = _convert_cartesian_to_axial(pos, self.origin, self.orientation);
        let mut field_z = 0.0;
        let mut field_r = 0.0;
        for (_, object) in self.objects.iter() {
            let [new_z, new_r] = object.get_field(&z, &r, tol);
            field_z += new_z;
            field_r += new_r;
        }
        [field_z, field_r]
    }
    /// Computes the magnetic field of the axial system
    ///
    /// Takes the position coordinates of the location for which the
    /// magentic field is desired. These coordinates are in the global
    /// space and not relative to the axial system.
    ///
    /// # Arguments
    /// * `(x,y,z)` tuple containing the x,y,z coordinates.
    /// * `tol` the tolerance at which the series expansion shuold terminate.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,1.0);
    /// let magnetic_field = axial.get_field([0.0,0.0,0.0],&1e-16);
    /// ```
    pub fn get_field(&self, pos: [f64; 3], tol: &f64) -> [f64; 3] {
        let [_z, _r, theta] = _convert_cartesian_to_axial(pos, self.origin, self.orientation);
        let [field_z, field_r] = self.get_field_in_frame(pos, tol);
        _convert_axial_to_cartesian([field_z, field_r, theta], self.origin, self.orientation)
    }
    /// Computes the magnetic field of the axial system in relative frame
    ///
    /// Takes the position coordinates of the location for which the
    /// magentic field is desired. These coordinates are in the global
    /// space and not relative to the axial system.
    ///
    /// # Arguments
    /// * `z` axial position relative to AxialSystem
    /// * `r` radial position relative to AxialSystem
    /// * `tol` the tolerance at which the series expansion shuold terminate.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use rustycoils::{AxialSystem,AxialError};
    /// let mut axial = AxialSystem::default();
    /// let res = axial.add_loop("loop1".to_string(),1.0,0.0,1.0);
    /// let magnetic_field = axial.get_field_axial(&2.0,&0.1,&1e-16);
    /// ```
    pub fn get_field_axial(&self, z: &f64, r: &f64, tol: &f64) -> [f64; 2] {
        let mut field_z = 0.0;
        let mut field_r = 0.0;
        for (_, object) in self.objects.iter() {
            let [new_z, new_r] = object.get_field(&z, &r, tol);
            field_z += new_z;
            field_r += new_r;
        }
        [field_z, field_r]
    }
}

//compute functions
//
//
#[derive(Debug, PartialEq, Clone, Copy)]
enum Primitives {
    IdealWire(IdealWire),
    ThinSolenoid(ThinSolenoid),
    Annular(ThinAnnular),
    CoilSolenoid(CoilSolenoid),
}
impl fmt::Display for Primitives {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Primitives::IdealWire(a) => write!(f, "{}", a),
            Primitives::ThinSolenoid(a) => write!(f, "{}", a),
            Primitives::Annular(a) => write!(f, "{}", a),
            Primitives::CoilSolenoid(a) => write!(f, "{}", a),
        }
    }
}
impl Primitives {
    pub fn get_field(&self, z: &f64, r: &f64, tol: &f64) -> [f64; 2] {
        match self {
            Primitives::IdealWire(wire) => wire.get_fields(z, r, tol),
            Primitives::ThinSolenoid(solenoid) => solenoid.get_fields(z, r, tol),
            Primitives::Annular(annular) => annular.get_fields(z, r, tol),
            Primitives::CoilSolenoid(coil) => coil.get_fields(z, r, tol),
        }
    }
}
#[derive(Debug, PartialEq)]
/// Errors produced by interacting with AxialSystem
pub enum AxialError {
    /// Error produced when provide an ID which is already used
    KeyDuplicateError(String),
    /// Error produced when wrong key provided to modify functions
    KeyMissingError(String),
    /// Error produced when provide an ID which matches one of the reserved words
    /// *`LOOP`
    /// *`ANNULAR`
    /// *`SOLENOID`
    /// *`COIL`
    /// *`*`
    ReservedWordError(String),
    /// Error produced when attempt to modify a primitive with wrong physical parameter.
    IncompatiblePrimitiveError(String, String),
}
impl fmt::Display for AxialError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AxialError::KeyDuplicateError(a) => {
                write!(f, "Supplied already in objects: {}.", a)
            }
            AxialError::KeyMissingError(a) => write!(f, "Could not find object ID: {}.", a),
            AxialError::ReservedWordError(a) => {
                write!(f, "Supplied ID is a reserved keyword: {}.", a)
            }
            AxialError::IncompatiblePrimitiveError(modification, primitive) => write!(
                f,
                "Can not perform modification {} on primitive {}.",
                modification, primitive
            ),
        }
    }
}
impl Error for AxialError {}

const RESERVED_WORDS: [&str; 5] = ["*", "COIL", "LOOP", "ANNULAR", "SOLENOID"];
//checks is the provided ID is allowed. Checks that the given
//string is not one the of reservered words and that the hashmap
//does not already contain the provided ID.
fn _is_id_valid(map: &HashMap<String, Primitives>, id: &str) -> AxialError {
    // check if word in the reserved words list
    for word in RESERVED_WORDS.iter() {
        if &id == word {
            return AxialError::ReservedWordError(id.to_string());
        }
    }
    // check if the hashmap already contains the given id.k

    match map.contains_key(id) {
        true => AxialError::KeyDuplicateError(id.to_string()),
        false => AxialError::KeyMissingError(id.to_string()),
    }
}

// functions takes a vector containing the real coordinates, and the vectors defining the axial
// space and converts the coordinates into the frame of the AxialObject.
// currently assumes that the origin is (0.0,0.0,0.0) and orientation is
// (1.0,0.0,0.0)/(0.0,1.0,0.0)/(0.0,0.0,1.0)
// hacky solution until impliment a general transformation.
#[allow(unused_variables)]
fn _convert_cartesian_to_axial(
    coordinates: [f64; 3],
    origin: [f64; 3],
    orientation: [f64; 3],
) -> [f64; 3] {
    let mut z = 0.0;
    let mut perp1 = 0.0;
    let mut perp2 = 0.0;
    if orientation[0] > 0.1 {
        z = coordinates[0];
        perp1 = coordinates[1];
        perp2 = coordinates[2];
    }
    if orientation[1] > 0.1 {
        z = coordinates[1];
        perp1 = coordinates[2];
        perp2 = coordinates[0];
    }
    if orientation[2] > 0.1 {
        z = coordinates[2];
        perp1 = coordinates[0];
        perp2 = coordinates[1];
    }

    let r = (perp1.powi(2) + perp2.powi(2)).sqrt();
    let theta = perp1.atan2(perp2);
    [z, r, theta]
}
//takes z,r,theta coordinates and converts into x,y,z
//again a hacky solution when only x,y,z at origin are allowed.
#[allow(unused_variables)]
fn _convert_axial_to_cartesian(
    coordinates: [f64; 3],
    origin: [f64; 3],
    orientation: [f64; 3],
) -> [f64; 3] {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    if orientation[0] > 0.1 {
        x = coordinates[0];
        y = coordinates[1] * (f64::sin(coordinates[2]));
        z = coordinates[1] * (f64::cos(coordinates[2]));
    }
    if orientation[1] > 0.1 {
        y = coordinates[0];
        z = coordinates[1] * (f64::sin(coordinates[2]));
        x = coordinates[1] * (f64::cos(coordinates[2]));
    }
    if orientation[2] > 0.1 {
        z = coordinates[0];
        x = coordinates[1] * (f64::sin(coordinates[2]));
        y = coordinates[1] * (f64::cos(coordinates[2]));
    }
    [x, y, z]
}

#[cfg(test)]
mod test_transformations {
    use super::*;
    #[test]
    fn check_temporary_converter() {
        let x = 1.0;
        let y = 2.0;
        let z = 4.0;

        let axial_coordinate =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert_eq!(
            axial_coordinate,
            [4.0, f64::sqrt(x * x + y * y), x.atan2(y)]
        );
        let axial_coordinate =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(
            axial_coordinate,
            [1.0, f64::sqrt(z * z + y * y), y.atan2(z)]
        );
        let axial_coordinate =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(
            axial_coordinate,
            [2.0, f64::sqrt(x * x + z * z), z.atan2(x)]
        );
    }
    #[test]
    //assume have a coil along the y axis and have a magnetic field (B_z,B_r)
    //want to check that once back in correct frame to converts.
    //assume that located along x axis.
    fn check_magentic_field_converts() {
        let b_z = 1.0;
        let b_r = 10.0;
        let theta = 0.0;

        let [bx, by, bz] =
            _convert_axial_to_cartesian([b_z, b_r, theta], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((by - b_z).abs() < 1e-8);
        assert!((bz - 0.0).abs() < 1e-8);
        assert!((bx - 10.0).abs() < 1e-8);
    }
    #[test]
    fn check_get_coordinates_back() {
        let x = 1.0;
        let y = 2.0;
        let z = 3.0;

        let axial_coordinates =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let new_carts =
            _convert_axial_to_cartesian(axial_coordinates, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert!((new_carts[0] - x).abs() < 1e-6);
        assert!((new_carts[1] - y).abs() < 1e-6);
        assert!((new_carts[2] - z).abs() < 1e-6);
        let axial_coordinates =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let new_carts =
            _convert_axial_to_cartesian(axial_coordinates, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((new_carts[0] - x).abs() < 1e-6);
        assert!((new_carts[1] - y).abs() < 1e-6);
        assert!((new_carts[2] - z).abs() < 1e-6);
        let axial_coordinates =
            _convert_cartesian_to_axial([x, y, z], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let new_carts =
            _convert_axial_to_cartesian(axial_coordinates, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((new_carts[0] - x).abs() < 1e-6);
        assert!((new_carts[1] - y).abs() < 1e-6);
        assert!((new_carts[2] - z).abs() < 1e-6);
    }
}

#[cfg(test)]
mod test_add_functions {
    use super::*;
    #[test]
    fn _test_is_id_valid_reserved_words() {
        let mymap: HashMap<String, Primitives> = HashMap::new();
        assert_eq!(
            AxialError::ReservedWordError("COIL".to_string()),
            _is_id_valid(&mymap, &"COIL".to_string())
        );
        assert_eq!(
            AxialError::ReservedWordError("LOOP".to_string()),
            _is_id_valid(&mymap, &"LOOP".to_string())
        );
        assert_eq!(
            AxialError::ReservedWordError("ANNULAR".to_string()),
            _is_id_valid(&mymap, &"ANNULAR".to_string())
        );
        assert_eq!(
            AxialError::ReservedWordError("SOLENOID".to_string()),
            _is_id_valid(&mymap, &"SOLENOID".to_string())
        );
        assert_eq!(
            AxialError::ReservedWordError("*".to_string()),
            _is_id_valid(&mymap, &"*".to_string())
        );
    }
    #[test]
    fn _test_is_id_valid_hash_map_cotains() {
        let mut mymap: HashMap<String, Primitives> = HashMap::new();
        let myobj = Primitives::IdealWire(IdealWire::new(0.0, 0.0, 0.0));
        mymap.insert("foo".to_string(), myobj);
        assert_eq!(
            AxialError::KeyDuplicateError("foo".to_string()),
            _is_id_valid(&mymap, &"foo".to_string())
        );
        assert_eq!(
            AxialError::KeyMissingError("bar".to_string()),
            _is_id_valid(&mymap, &"bar".to_string())
        );
    }

    #[test]
    fn test_add_coil() {
        let mut myaxial = AxialSystem::default();
        let radius = 1.0;
        let x0 = 1.0;
        let current = 1.0;
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Err(AxialError::KeyDuplicateError("loop1".to_string())));
        let res = myaxial.add_loop("loop2".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
    }
    #[test]
    fn test_add_and_remove() {
        let mut myaxial = AxialSystem::default();
        let radius = 1.0;
        let x0 = 1.0;
        let current = 1.0;
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
        let res = myaxial.remove("loop1");
        assert_eq!(res, Ok(()));
    }
    #[test]
    fn test_add_and_remove_reserved() {
        let mut myaxial = AxialSystem::default();
        let radius = 1.0;
        let x0 = 1.0;
        let current = 1.0;
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
        let res = myaxial.remove("loop2");
        assert_eq!(res, Err(AxialError::KeyMissingError("loop2".to_string())));
        let res = myaxial.remove("*");
        assert_eq!(res, Err(AxialError::ReservedWordError("*".to_string())));
    }
}

#[cfg(test)]
mod test_modify_functions {
    use super::*;
    #[test]
    fn test_modify_length() {
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("LOOP1".to_string(), 1.0, 1.0, 1.0);

        let res = myaxial.modify_radius("LOOP1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("LOOP1", 2.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "LOOP1".to_string(),
                "LOOP".to_string()
            ))
        );
        let res = myaxial.modify_length("LOOP2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("LOOP2".to_string())));
    }
    #[test]
    fn test_modify_coil() {
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_coil_solenoid("coil1".to_string(), 1.0, 1.0, 1.0, 1.0, 1.0);

        let res = myaxial.modify_thickness("coil1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("coil1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_radius("coil1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("coil1", 2.0);
        assert_eq!(res, Ok(()));

        let res = myaxial.modify_thickness("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_length("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_radius("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_position("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));

        let res = myaxial.modify_thickness("COIL", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("COIL", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("COIL", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_radius("COIL", 3.0);
        assert_eq!(res, Ok(()));
    }
    #[test]
    fn test_modify_loop() {
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_loop("loop1".to_string(), 1.0, 1.0, 1.0);

        let res = myaxial.modify_thickness("loop1", 2.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "loop1".to_string(),
                "LOOP".to_string()
            ))
        );
        let res = myaxial.modify_length("loop1", 2.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "loop1".to_string(),
                "LOOP".to_string()
            ))
        );
        let res = myaxial.modify_radius("loop1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("loop1", 2.0);
        assert_eq!(res, Ok(()));

        let res = myaxial.modify_thickness("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_length("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_radius("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_position("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));

        let res = myaxial.modify_thickness("LOOP", 3.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "LOOP".to_string(),
                "LOOP".to_string()
            ))
        );
        let res = myaxial.modify_position("LOOP", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("LOOP", 3.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "LOOP".to_string(),
                "LOOP".to_string()
            ))
        );
        let res = myaxial.modify_radius("LOOP", 3.0);
        assert_eq!(res, Ok(()));
    }
    #[test]
    fn test_modify_solenoid() {
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_thin_solenoid("solenoid1".to_string(), 1.0, 1.0, 1.0, 1.0);

        let res = myaxial.modify_length("solenoid1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_thickness("solenoid1", 2.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "solenoid1".to_string(),
                "SOLENOID".to_string()
            ))
        );
        let res = myaxial.modify_radius("solenoid1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("solenoid1", 2.0);
        assert_eq!(res, Ok(()));

        let res = myaxial.modify_thickness("solenoid2", 2.0);
        assert_eq!(
            res,
            Err(AxialError::KeyMissingError("solenoid2".to_string()))
        );
        let res = myaxial.modify_length("solenoid2", 2.0);
        assert_eq!(
            res,
            Err(AxialError::KeyMissingError("solenoid2".to_string()))
        );
        let res = myaxial.modify_radius("solenoid2", 2.0);
        assert_eq!(
            res,
            Err(AxialError::KeyMissingError("solenoid2".to_string()))
        );
        let res = myaxial.modify_position("solenoid2", 2.0);
        assert_eq!(
            res,
            Err(AxialError::KeyMissingError("solenoid2".to_string()))
        );

        let res = myaxial.modify_thickness("SOLENOID", 3.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "SOLENOID".to_string(),
                "SOLENOID".to_string()
            ))
        );
        let res = myaxial.modify_position("SOLENOID", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("SOLENOID", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_radius("SOLENOID", 3.0);
        assert_eq!(res, Ok(()));
    }

    #[test]
    fn test_modify_annular() {
        let mut myaxial = AxialSystem::default();
        let _res = myaxial.add_annular("annular1".to_string(), 1.0, 1.0, 1.0, 1.0);

        let res = myaxial.modify_thickness("annular1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("annular1", 2.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "annular1".to_string(),
                "ANNULAR".to_string()
            ))
        );
        let res = myaxial.modify_radius("annular1", 2.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("annular1", 2.0);
        assert_eq!(res, Ok(()));

        let res = myaxial.modify_thickness("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_length("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_radius("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));
        let res = myaxial.modify_position("coil2", 2.0);
        assert_eq!(res, Err(AxialError::KeyMissingError("coil2".to_string())));

        let res = myaxial.modify_thickness("ANNULAR", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_position("ANNULAR", 3.0);
        assert_eq!(res, Ok(()));
        let res = myaxial.modify_length("ANNULAR", 3.0);
        assert_eq!(
            res,
            Err(AxialError::IncompatiblePrimitiveError(
                "ANNULAR".to_string(),
                "ANNULAR".to_string()
            ))
        );
        let res = myaxial.modify_radius("ANNULAR", 3.0);
        assert_eq!(res, Ok(()));
    }
}
