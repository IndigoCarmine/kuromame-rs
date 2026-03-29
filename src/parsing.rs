use crate::view_rs::To3dViewMolecule;
use moleucle_3dview_rs::{
    Molecule,
    molecule::{Atom, Bond},
};
use lin_alg::f32::Vec3;
use std::collections::HashMap;
use std::fmt::Write;

// --- PDB Structures ---

#[derive(Debug, Clone, PartialEq)]
pub struct AtomRecord {
    pub serial: usize,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: i32,
    pub i_code: char,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub occupancy: f32,
    pub temp_factor: f32,
    pub element: String,
    pub charge: String,
}

impl AtomRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        if line.len() < 54 {
            return None;
        }

        let parse_int =
            |range: std::ops::Range<usize>| -> Option<i32> { line.get(range)?.trim().parse().ok() };
        let parse_float =
            |range: std::ops::Range<usize>| -> Option<f32> { line.get(range)?.trim().parse().ok() };
        let parse_str = |range: std::ops::Range<usize>| -> String {
            line.get(range)
                .map(|s| s.trim().to_string())
                .unwrap_or_default()
        };
        let parse_char = |index: usize| -> char {
            line.get(index..index + 1)
                .and_then(|s| s.chars().next())
                .unwrap_or(' ')
        };

        Some(AtomRecord {
            serial: parse_int(6..11)? as usize,
            name: parse_str(12..16),
            alt_loc: parse_char(16),
            res_name: parse_str(17..20),
            chain_id: parse_char(21),
            res_seq: parse_int(22..26)?,
            i_code: parse_char(26),
            x: parse_float(30..38)?,
            y: parse_float(38..46)?,
            z: parse_float(46..54)?,
            occupancy: parse_float(54..60).unwrap_or(0.0),
            temp_factor: parse_float(60..66).unwrap_or(0.0),
            element: parse_str(76..78),
            charge: parse_str(78..80),
        })
    }

    pub fn to_line(&self) -> String {
        format!(
            "ATOM  {:5} {:<4}{}{:>3} {}{:4}{:<4}{:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:>2}",
            self.serial,
            self.name,
            self.alt_loc,
            self.res_name,
            self.chain_id,
            self.res_seq,
            self.i_code,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.temp_factor,
            self.element,
            self.charge
        )
    }
}

#[derive(Debug, Clone)]
pub struct ConectRecord {
    pub serial: usize,
    pub bonded: Vec<usize>,
}

impl ConectRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        let serial = line.get(6..11)?.trim().parse().ok()?;
        let mut bonded = Vec::new();

        for range in [(11, 16), (16, 21), (21, 26), (26, 31)] {
            if let Some(s) = line.get(range.0..range.1) {
                if let Ok(v) = s.trim().parse::<usize>() {
                    bonded.push(v);
                }
            }
        }

        Some(ConectRecord { serial, bonded })
    }

    pub fn to_line(&self) -> String {
        let mut s = format!("CONECT{:5}", self.serial);
        for b in &self.bonded {
            write!(s, "{:5}", b).unwrap();
        }
        s
    }
}

#[derive(Debug, Clone)]
pub enum PdbLine {
    Atom(AtomRecord),
    Conect(ConectRecord),
    Other(String),
}

#[derive(Debug, Clone, Default)]
pub struct PdbFile {
    pub lines: Vec<PdbLine>,
}

impl PdbFile {
    pub fn load(content: &str) -> Self {
        let mut lines = Vec::new();
        for line in content.lines() {
            if line.starts_with("ATOM") {
                if let Some(atom) = AtomRecord::from_line(line) {
                    lines.push(PdbLine::Atom(atom));
                    continue;
                }
            } else if line.starts_with("CONECT") {
                if let Some(conect) = ConectRecord::from_line(line) {
                    lines.push(PdbLine::Conect(conect));
                    continue;
                }
            }
            lines.push(PdbLine::Other(line.to_string()));
        }
        Self { lines }
    }

    pub fn dump(&mut self) -> String {
        self.update_resseq();
        let mut out = String::new();
        for line in &self.lines {
            match line {
                PdbLine::Atom(a) => writeln!(out, "{}", a.to_line()).unwrap(),
                PdbLine::Conect(c) => writeln!(out, "{}", c.to_line()).unwrap(),
                PdbLine::Other(s) => writeln!(out, "{}", s).unwrap(),
            }
        }
        out
    }

    pub fn atoms(&self) -> impl Iterator<Item = &AtomRecord> {
        self.lines.iter().filter_map(|l| match l {
            PdbLine::Atom(a) => Some(a),
            _ => None,
        })
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomRecord> {
        self.lines.iter_mut().filter_map(|l| match l {
            PdbLine::Atom(a) => Some(a),
            _ => None,
        })
    }

    pub fn update_resseq(&mut self) {
        // Collect unique residue names from atoms
        let mut resnames: Vec<String> = self
            .atoms()
            .map(|a| a.res_name.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        resnames.sort();

        // Create a mapping from resname to resseq
        let resseq_map: HashMap<String, i32> = resnames
            .into_iter()
            .enumerate()
            .map(|(i, name)| (name, (i + 1) as i32))
            .collect();

        // Update all atoms with new resseq
        for atom in self.atoms_mut() {
            if let Some(&new_seq) = resseq_map.get(&atom.res_name) {
                atom.res_seq = new_seq;
            }
        }
    }

    pub fn find_connected_hydrogen(&self, atom: &AtomRecord) -> Vec<AtomRecord> {
        let mut connected_serials = std::collections::HashSet::new();
        let target_serial = atom.serial;

        // Scan all CONECT records
        for line in &self.lines {
            if let PdbLine::Conect(conect) = line {
                // Case 1: The record is for the target atom
                if conect.serial == target_serial {
                    for serial in &conect.bonded {
                        connected_serials.insert(*serial);
                    }
                }
                // Case 2: The target atom is in the bonded list of another atom
                else if conect.bonded.contains(&target_serial) {
                    connected_serials.insert(conect.serial);
                }
            }
        }

        // Find hydrogen atoms in connected serials
        let mut hydrogens = Vec::new();
        for atom_rec in self.atoms() {
            if connected_serials.contains(&atom_rec.serial) && atom_rec.element == "H" {
                hydrogens.push(atom_rec.clone());
            }
        }
        hydrogens
    }
}

impl To3dViewMolecule for PdbFile {
    fn to_molecule(&self) -> Molecule {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();

        let mut serial_to_index = HashMap::new();

        // 1. Collect atoms
        for (i, record) in self.atoms().enumerate() {
            serial_to_index.insert(record.serial, i);
            atoms.push(Atom {
                position: Vec3::new(record.x, record.y, record.z),
                element: record.element.clone(),
                id: i,
                name: Some(record.name.clone()),
                res_name: Some(record.res_name.clone()),
                chain_id: Some(record.chain_id),
                res_seq: Some(record.res_seq),
                occupancy: Some(record.occupancy),
                temp_factor: Some(record.temp_factor),
                charge: Some(record.charge.clone()),
            });
        }

        // 2. Collect bonds from CONECT records
        for line in &self.lines {
            if let PdbLine::Conect(c) = line {
                if let Some(&idx_a) = serial_to_index.get(&c.serial) {
                    for &bonded_serial in &c.bonded {
                        if let Some(&idx_b) = serial_to_index.get(&bonded_serial) {
                            if idx_a < idx_b {
                                bonds.push(Bond {
                                    atom_a: idx_a,
                                    atom_b: idx_b,
                                    order: 1,
                                });
                            }
                        }
                    }
                }
            }
        }

        Molecule { atoms, bonds }
    }
}

// --- Mol2 Structures ---

#[derive(Debug, Clone)]
pub struct Mol2AtomRecord {
    pub atom_id: usize,
    pub atom_name: String,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_type: String,
    pub subst_id: i32,
    pub subst_name: String,
    pub charge: f32,
    pub status_bit: String,
}

impl Mol2AtomRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 9 {
            return None;
        }

        Some(Mol2AtomRecord {
            atom_id: parts[0].parse().ok()?,
            atom_name: parts[1].to_string(),
            x: parts[2].parse().ok()?,
            y: parts[3].parse().ok()?,
            z: parts[4].parse().ok()?,
            atom_type: parts[5].to_string(),
            subst_id: parts[6].parse().ok()?,
            subst_name: parts[7].to_string(),
            charge: parts[8].parse().ok()?,
            status_bit: parts.get(9).unwrap_or(&"").to_string(),
        })
    }

    pub fn to_line(&self) -> String {
        format!(
            "{:>7} {:<8}{:>10.4}{:>10.4}{:>10.4} {:<5}{:>6}     {:<4}{:>10.4}",
            self.atom_id,
            self.atom_name,
            self.x,
            self.y,
            self.z,
            self.atom_type,
            self.subst_id,
            self.subst_name,
            self.charge
        )
    }
}

#[derive(Debug, Clone)]
pub struct Mol2BondRecord {
    pub bond_id: usize,
    pub origin_atom_id: usize,
    pub target_atom_id: usize,
    pub bond_type: String,
    pub status_bit: String,
}

impl Mol2BondRecord {
    pub fn from_line(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return None;
        }

        Some(Mol2BondRecord {
            bond_id: parts[0].parse().ok()?,
            origin_atom_id: parts[1].parse().ok()?,
            target_atom_id: parts[2].parse().ok()?,
            bond_type: parts[3].to_string(),
            status_bit: parts.get(4).unwrap_or(&"").to_string(),
        })
    }

    pub fn to_line(&self) -> String {
        format!(
            "{:>6}{:>6}{:>6}{:>6}",
            self.bond_id, self.origin_atom_id, self.target_atom_id, self.bond_type
        )
    }
}

#[derive(Debug, Clone)]
pub enum Mol2Line {
    Atom(Mol2AtomRecord),
    Bond(Mol2BondRecord),
    SectionHeader(String),
    MoleculeLine(String),
    Other(String),
    Empty,
}

#[derive(Debug, Clone, Default)]
pub struct Mol2File {
    pub lines: Vec<Mol2Line>,
}

impl Mol2File {
    pub fn load(content: &str) -> Self {
        let mut lines = Vec::new();
        let mut current_section = "";

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                lines.push(Mol2Line::Empty);
                continue;
            }

            if line.starts_with("@<TRIPOS>") {
                current_section = &line[9..];
                lines.push(Mol2Line::SectionHeader(line.to_string()));
                continue;
            }

            match current_section {
                "MOLECULE" => lines.push(Mol2Line::MoleculeLine(line.to_string())),
                "ATOM" => {
                    if let Some(atom) = Mol2AtomRecord::from_line(line) {
                        lines.push(Mol2Line::Atom(atom));
                    } else {
                        lines.push(Mol2Line::Other(line.to_string()));
                    }
                }
                "BOND" => {
                    if let Some(bond) = Mol2BondRecord::from_line(line) {
                        lines.push(Mol2Line::Bond(bond));
                    } else {
                        lines.push(Mol2Line::Other(line.to_string()));
                    }
                }
                _ => lines.push(Mol2Line::Other(line.to_string())),
            }
        }
        Self { lines }
    }

    pub fn dump(&self) -> String {
        let mut out = String::new();
        for line in &self.lines {
            match line {
                Mol2Line::Atom(a) => writeln!(out, "{}", a.to_line()).unwrap(),
                Mol2Line::Bond(b) => writeln!(out, "{}", b.to_line()).unwrap(),
                Mol2Line::SectionHeader(s) | Mol2Line::MoleculeLine(s) | Mol2Line::Other(s) => {
                    writeln!(out, "{}", s).unwrap()
                }
                Mol2Line::Empty => writeln!(out).unwrap(),
            }
        }
        out
    }

    pub fn atoms(&self) -> impl Iterator<Item = &Mol2AtomRecord> {
        self.lines.iter().filter_map(|l| match l {
            Mol2Line::Atom(a) => Some(a),
            _ => None,
        })
    }

    pub fn bonds(&self) -> impl Iterator<Item = &Mol2BondRecord> {
        self.lines.iter().filter_map(|l| match l {
            Mol2Line::Bond(b) => Some(b),
            _ => None,
        })
    }
}

impl To3dViewMolecule for Mol2File {
    fn to_molecule(&self) -> Molecule {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();

        // Mol2 IDs are usually 1-based, but we map them to 0-based index
        let mut id_to_index = HashMap::new();

        for (i, record) in self.atoms().enumerate() {
            id_to_index.insert(record.atom_id, i);
            atoms.push(Atom {
                position: Vec3::new(record.x, record.y, record.z),
                element: record
                    .atom_type
                    .split('.')
                    .next()
                    .unwrap_or("?")
                    .to_uppercase(),
                id: i,
                name: None,
                res_name: None,
                chain_id: None,
                res_seq: None,
                occupancy: None,
                temp_factor: None,
                charge: None,
            });
        }

        for record in self.bonds() {
            if let (Some(&idx_a), Some(&idx_b)) = (
                id_to_index.get(&record.origin_atom_id),
                id_to_index.get(&record.target_atom_id),
            ) {
                let order = match record.bond_type.as_str() {
                    "2" => 2,
                    "3" => 3,
                    _ => 1,
                };
                bonds.push(Bond {
                    atom_a: idx_a,
                    atom_b: idx_b,
                    order,
                });
            }
        }

        Molecule { atoms, bonds }
    }
}
