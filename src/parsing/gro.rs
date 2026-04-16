use crate::view_rs::To3dViewMolecule;
use lin_alg::f32::Vec3;
use moleucle_3dview_rs::{
    molecule::{Atom, Bond},
    Molecule,
};
use std::collections::HashMap;
use std::fs::File;
use std::fmt::Write;
use std::io::{self, BufRead, BufReader};

// Gromacs GRO coordinates and box vectors are represented in nanometers.
#[allow(dead_code)]
pub const GROMACS_LENGTH_UNIT: &str = "nm";
pub const ANGSTROM_TO_NM: f32 = 0.1;

// --- GRO Structures ---

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GroFixed5([u8; 5]);

impl GroFixed5 {
    fn from_field(field: &str) -> Self {
        let mut out = [b' '; 5];
        let bytes = field.as_bytes();
        let n = bytes.len().min(5);
        out[..n].copy_from_slice(&bytes[..n]);
        Self(out)
    }

    fn from_left_aligned(value: &str) -> Self {
        let mut out = [b' '; 5];
        let trimmed = value.trim();
        let bytes = trimmed.as_bytes();
        let n = bytes.len().min(5);
        out[..n].copy_from_slice(&bytes[..n]);
        Self(out)
    }

    fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or("     ")
    }

    pub fn trimmed(&self) -> &str {
        self.as_str().trim()
    }
}

#[derive(Debug, Clone)]
pub struct GroAtomRecord {
    pub res_num: i32,
    pub res_name: GroFixed5,
    pub atom_name: GroFixed5,
    pub atom_id: usize,
    // GRO coordinate unit is nanometer (nm).
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl GroAtomRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        let res_num = line.get(0..5)?.trim().parse().ok()?;
        let res_name = GroFixed5::from_field(line.get(5..10)?);
        let atom_name = GroFixed5::from_field(line.get(10..15)?);
        let atom_id = line.get(15..20)?.trim().parse().ok()?;
        let x = line.get(20..28)?.trim().parse().ok()?;
        let y = line.get(28..36)?.trim().parse().ok()?;
        let z = line.get(36..44)?.trim().parse().ok()?;

        Some(Self {
            res_num,
            res_name,
            atom_name,
            atom_id,
            x,
            y,
            z,
        })
    }

    pub fn set_res_name(&mut self, value: &str) {
        self.res_name = GroFixed5::from_left_aligned(value);
    }

    pub fn to_line(&self) -> String {
        format!(
            "{:>5}{}{}{:>5}{:>8.3}{:>8.3}{:>8.3}",
            self.res_num,
            self.res_name.as_str(),
            self.atom_name.as_str(),
            self.atom_id,
            self.x,
            self.y,
            self.z
        )
    }
}

#[derive(Debug, Clone, Default)]
pub struct GroFile {
    pub title: String,
    pub atom_count_line: String,
    pub atoms: Vec<GroAtomRecord>,
    pub box_line: String,
}

impl GroFile {
    pub fn load(content: &str) -> Self {
        let mut iter = content.lines();
        let title = iter.next().unwrap_or("").to_string();
        let atom_count_line = iter.next().unwrap_or("").to_string();
        let declared_atom_count = atom_count_line.trim().parse::<usize>().unwrap_or(0);

        let mut atoms = Vec::with_capacity(declared_atom_count);
        let mut box_line = String::new();

        for line in iter {
            if box_line.is_empty() && !line.trim().is_empty() {
                if let Some(atom) = GroAtomRecord::from_line(line) {
                    atoms.push(atom);
                } else {
                    box_line = line.to_string();
                }
            } else if !box_line.is_empty() {
                box_line = line.to_string();
            }
        }

        Self {
            title,
            atom_count_line,
            atoms,
            box_line,
        }
    }

    pub fn load_from_reader<R: BufRead>(reader: R) -> io::Result<Self> {
        let mut iter = reader.lines();
        let title = iter.next().transpose()?.unwrap_or_default();
        let atom_count_line = iter.next().transpose()?.unwrap_or_default();
        let declared_atom_count = atom_count_line.trim().parse::<usize>().unwrap_or(0);

        let mut atoms = Vec::with_capacity(declared_atom_count);
        let mut box_line = String::new();

        for line in iter {
            let line = line?;
            if box_line.is_empty() && !line.trim().is_empty() {
                if let Some(atom) = GroAtomRecord::from_line(&line) {
                    atoms.push(atom);
                } else {
                    box_line = line;
                }
            } else if !box_line.is_empty() {
                box_line = line;
            }
        }

        Ok(Self {
            title,
            atom_count_line,
            atoms,
            box_line,
        })
    }

    pub fn load_from_path(path: &std::path::Path) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::load_from_reader(reader)
    }

    pub fn dump(&self) -> String {
        let mut out = String::new();
        writeln!(out, "{}", self.title).unwrap();
        if self.atom_count_line.trim().parse::<usize>().ok() == Some(self.atoms.len()) {
            writeln!(out, "{}", self.atom_count_line).unwrap();
        } else {
            writeln!(out, "{:>5}", self.atoms.len()).unwrap();
        }
        for atom in &self.atoms {
            writeln!(out, "{}", atom.to_line()).unwrap();
        }
        writeln!(out, "{}", self.box_line).unwrap();
        out
    }

    pub fn atoms(&self) -> impl Iterator<Item = &GroAtomRecord> {
        self.atoms.iter()
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut GroAtomRecord> {
        self.atoms.iter_mut()
    }

    fn infer_element(atom_name: &str) -> String {
        let mut chars = atom_name.trim().chars().filter(|c| c.is_ascii_alphabetic());
        let first = chars.next();
        let second = chars.next();

        match (first, second) {
            (Some(a), Some(b)) => format!("{}{}", a, b).to_uppercase(),
            (Some(a), None) => a.to_string().to_uppercase(),
            _ => "X".to_string(),
        }
    }

    fn covalent_radius_nm(element: &str) -> f32 {
        let angstrom = match element.trim().to_uppercase().as_str() {
            "H" => 0.31,
            "C" => 0.76,
            "N" => 0.71,
            "O" => 0.66,
            "F" => 0.57,
            "P" => 1.07,
            "S" => 1.05,
            "CL" => 1.02,
            "BR" => 1.20,
            "I" => 1.39,
            _ => 0.77,
        };
        angstrom * ANGSTROM_TO_NM
    }

    fn infer_single_bonds_from_distance(atoms: &[Atom]) -> Vec<Bond> {
        let mut bonds = Vec::new();

        // Bond inference is computed in nm (same as GRO and viewer coordinates).
        const EXTRA_TOLERANCE: f32 = 0.45 * ANGSTROM_TO_NM;
        const MIN_DISTANCE_2: f32 = (0.4 * ANGSTROM_TO_NM) * (0.4 * ANGSTROM_TO_NM);
        const MAX_COVALENT_RADIUS: f32 = 1.39 * ANGSTROM_TO_NM;
        const CELL_SIZE: f32 = MAX_COVALENT_RADIUS * 2.0 + EXTRA_TOLERANCE;
        const CELL_SIZE_INV: f32 = 1.0 / CELL_SIZE;

        let mut spatial_index: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

        for (i, atom) in atoms.iter().enumerate() {
            let cell = (
                (atom.position.x * CELL_SIZE_INV).floor() as i32,
                (atom.position.y * CELL_SIZE_INV).floor() as i32,
                (atom.position.z * CELL_SIZE_INV).floor() as i32,
            );

            let ri = Self::covalent_radius_nm(&atom.element);

            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(candidate_indices) =
                            spatial_index.get(&(cell.0 + dx, cell.1 + dy, cell.2 + dz))
                        {
                            for &j in candidate_indices {
                                let delta = atoms[j].position - atom.position;
                                let distance =
                                    delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                                if distance < MIN_DISTANCE_2 {
                                    continue;
                                }

                                let rj = Self::covalent_radius_nm(&atoms[j].element);
                                let max_bond_distance = ri + rj + EXTRA_TOLERANCE;
                                if distance <= max_bond_distance * max_bond_distance {
                                    bonds.push(Bond {
                                        atom_a: j,
                                        atom_b: i,
                                        order: 1,
                                    });
                                }
                            }
                        }
                    }
                }
            }

            spatial_index.entry(cell).or_default().push(i);
        }

        bonds
    }
}

impl GroFile {
    pub fn to_molecule_with_metadata(&self, include_metadata: bool) -> Molecule {
        // Viewer molecule coordinates are interpreted in nm.
        // GRO is also nm, so pass through coordinates without scaling.

        let atoms = self
            .atoms()
            .enumerate()
            .map(|(idx, atom)| Atom {
                position: Vec3::new(
                    atom.x,
                    atom.y,
                    atom.z,
                ),
                element: Self::infer_element(atom.atom_name.trimmed()),
                id: idx,
                name: if include_metadata {
                    Some(atom.atom_name.trimmed().to_string())
                } else {
                    None
                },
                res_name: if include_metadata {
                    Some(atom.res_name.trimmed().to_string())
                } else {
                    None
                },
                chain_id: Some('A'),
                res_seq: Some(atom.res_num),
                occupancy: None,
                temp_factor: None,
                charge: None,
            })
            .collect::<Vec<_>>();

        let bonds = Self::infer_single_bonds_from_distance(&atoms);

        Molecule { atoms, bonds }
    }
}

impl To3dViewMolecule for GroFile {
    fn to_molecule(&self) -> Molecule {
        self.to_molecule_with_metadata(true)
    }
}