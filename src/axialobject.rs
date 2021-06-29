use crate::fieldcalc::primitives::{CoilSolenoid, IdealWire, Primitive, ThinAnnular, ThinSolenoid};
use std::collections::HashMap;
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
pub struct AxialObject {
    objects: HashMap<String, Primitives>,
    origin: (f64, f64, f64),
    orientation: (f64, f64, f64),
}

impl AxialObject {
    pub fn new(origin: (f64, f64, f64), orientation: (f64, f64, f64)) -> AxialObject {
        AxialObject {
            objects: HashMap::new(),
            origin,
            orientation,
        }
    }
    pub fn default() -> AxialObject {
        AxialObject {
            objects: HashMap::new(),
            origin: (0.0, 0.0, 0.0),
            orientation: (1.0, 0.0, 0.0),
        }
    }
}
//Add functions
impl AxialObject {
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
    //removes an object matching the provided id.
    pub fn remove(&mut self, id: String) -> Result<(), AxialError> {
        match _is_id_valid(&self.objects, &id) {
            AxialError::KeyDuplicateError(_) => {}
            error => return Err(error),
        }
        self.objects.remove(&id);
        Ok(())
    }
}

//transform functions
impl AxialObject {
    pub fn transform_x(&mut self) {
        self.orientation = (1.0, 0.0, 0.0);
    }
    pub fn transform_y(&mut self) {
        self.orientation = (0.0, 1.0, 0.0);
    }
    pub fn transform_z(&mut self) {
        self.orientation = (0.0, 0.0, 1.0);
    }
}

//modify functions
//these functions allow the individual physical parameters of
//underlying primitives to be modified.

impl AxialObject {
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

impl AxialObject {
    pub fn get_field(&self, (x, y, z): &(f64, f64, f64), tol: &f64) -> (f64, f64) {
        let (z, r) = _convert_cartesian_to_axial((*x, *y, *z), self.origin, self.orientation);
        let mut field_z = 0.0;
        let mut field_r = 0.0;
        for (_, object) in self.objects.iter() {
            let (new_z, new_r) = object.get_field(&z, &r, tol);
            field_z += new_z;
            field_r += new_r;
        }
        (field_z, field_r)
    }
}

//compute functions
//
//
#[derive(Debug, PartialEq)]
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
    pub fn get_field(&self, z: &f64, r: &f64, tol: &f64) -> (f64, f64) {
        match self {
            Primitives::IdealWire(wire) => wire.get_fields(z, r, tol),
            Primitives::ThinSolenoid(solenoid) => solenoid.get_fields(z, r, tol),
            Primitives::Annular(annular) => annular.get_fields(z, r, tol),
            Primitives::CoilSolenoid(coil) => coil.get_fields(z, r, tol),
        }
    }
}
#[derive(Debug, PartialEq)]
pub enum AxialError {
    KeyDuplicateError(String),
    KeyMissingError(String),
    ReservedWordError(String),
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
// in future will make general in origin and orientation.
#[allow(unused_variables)]
fn _convert_cartesian_to_axial(
    coordinates: (f64, f64, f64),
    origin: (f64, f64, f64),
    orientation: (f64, f64, f64),
) -> (f64, f64) {
    let z_direction = coordinates.0 * orientation.0
        + coordinates.1 * orientation.1
        + coordinates.2 * orientation.2;
    let r_direction = f64::sqrt(
        f64::powi((1.0 - orientation.0) * coordinates.0, 2)
            + f64::powi((1.0 - orientation.1) * coordinates.1, 2)
            + f64::powi((1.0 - orientation.2) * coordinates.2, 2),
    );
    (z_direction, r_direction)
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
            _convert_cartesian_to_axial((x, y, z), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0));
        assert_eq!(axial_coordinate, (4.0, f64::sqrt(x * x + y * y)));
        let axial_coordinate =
            _convert_cartesian_to_axial((x, y, z), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0));
        assert_eq!(axial_coordinate, (1.0, f64::sqrt(z * z + y * y)));
        let axial_coordinate =
            _convert_cartesian_to_axial((x, y, z), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0));
        assert_eq!(axial_coordinate, (2.0, f64::sqrt(x * x + z * z)));
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
        let mut myaxial = AxialObject::default();
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
        let mut myaxial = AxialObject::default();
        let radius = 1.0;
        let x0 = 1.0;
        let current = 1.0;
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
        let res = myaxial.remove("loop1".to_string());
        assert_eq!(res, Ok(()));
    }
    #[test]
    fn test_add_and_remove_reserved() {
        let mut myaxial = AxialObject::default();
        let radius = 1.0;
        let x0 = 1.0;
        let current = 1.0;
        let res = myaxial.add_loop("loop1".to_string(), radius, x0, current);
        assert_eq!(res, Ok(()));
        let res = myaxial.remove("loop2".to_string());
        assert_eq!(res, Err(AxialError::KeyMissingError("loop2".to_string())));
        let res = myaxial.remove("*".to_string());
        assert_eq!(res, Err(AxialError::ReservedWordError("*".to_string())));
    }
}

#[cfg(test)]
mod test_modify_functions {
    use super::*;
    #[test]
    fn test_modify_length() {
        let mut myaxial = AxialObject::default();
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
        let mut myaxial = AxialObject::default();
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
        let mut myaxial = AxialObject::default();
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
        let mut myaxial = AxialObject::default();
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
        let mut myaxial = AxialObject::default();
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
