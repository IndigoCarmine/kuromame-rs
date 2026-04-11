use crate::view_rs::To3dViewMolecule;
use moleucle_3dview_rs::{
    Molecule,
    molecule::{Atom, Bond},
};
pub use moleucle_3dview_rs::AtomRecord;
use lin_alg::f32::Vec3;
use std::collections::HashMap;
use std::fmt::Write;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

// --- PDB Structures ---

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
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
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

    pub fn from_molecule(molecule: &Molecule) -> Self {
        let mut lines: Vec<PdbLine> = Vec::new();
        lines.reserve(molecule.atoms.len() + molecule.bonds.len());

        for (idx, atom) in molecule.atoms.iter().enumerate() {
            let serial = idx + 1;
            let name = atom
                .name
                .as_deref()
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .unwrap_or_else(|| atom.element.clone());

            let res_name = atom
                .res_name
                .as_deref()
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .unwrap_or_else(|| "LIG".to_string());

            let element = atom
                .element
                .trim()
                .chars()
                .take(2)
                .collect::<String>()
                .to_uppercase();

            lines.push(PdbLine::Atom(AtomRecord {
                serial,
                name,
                alt_loc: ' ',
                res_name,
                chain_id: atom.chain_id.unwrap_or('A'),
                res_seq: atom.res_seq.unwrap_or(1),
                i_code: ' ',
                x: atom.position.x,
                y: atom.position.y,
                z: atom.position.z,
                occupancy: atom.occupancy.unwrap_or(1.0),
                temp_factor: atom.temp_factor.unwrap_or(0.0),
                element,
                charge: atom.charge.clone().unwrap_or_default(),
            }));
        }

        let mut bonded_by_atom: HashMap<usize, Vec<usize>> = HashMap::new();
        for bond in &molecule.bonds {
            let a = bond.atom_a + 1;
            let b = bond.atom_b + 1;
            bonded_by_atom.entry(a).or_default().push(b);
            bonded_by_atom.entry(b).or_default().push(a);
        }

        let mut serials: Vec<usize> = bonded_by_atom.keys().copied().collect();
        serials.sort_unstable();
        for serial in serials {
            if let Some(mut bonded) = bonded_by_atom.remove(&serial) {
                bonded.sort_unstable();
                bonded.dedup();
                lines.push(PdbLine::Conect(ConectRecord { serial, bonded }));
            }
        }

        Self { lines }
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

    fn covalent_radius_angstrom(element: &str) -> f32 {
        match element.trim().to_uppercase().as_str() {
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
        }
    }

    fn infer_single_bonds_from_distance(atoms: &[Atom]) -> Vec<Bond> {
        let mut bonds = Vec::new();

        // Generous tolerance for coordinate noise and coarse input structures.
        const EXTRA_TOLERANCE: f32 = 0.45;
        const MIN_DISTANCE_2: f32 = 0.4 * 0.4; // Minimum distance squared to avoid false positives
        const MAX_COVALENT_RADIUS: f32 = 1.39;
        const CELL_SIZE: f32 = MAX_COVALENT_RADIUS * 2.0 + EXTRA_TOLERANCE;
        const CELL_SIZE_INV: f32 = 1.0 / CELL_SIZE;

        let mut spatial_index: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

        for (i, atom) in atoms.iter().enumerate() {
            let cell = (
                (atom.position.x * CELL_SIZE_INV).floor() as i32,
                (atom.position.y * CELL_SIZE_INV).floor() as i32,
                (atom.position.z * CELL_SIZE_INV).floor() as i32,
            );

            let ri = Self::covalent_radius_angstrom(&atom.element);

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

                                let rj = Self::covalent_radius_angstrom(&atoms[j].element);
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
        const GRO_TO_VIEW_SCALE: f32 = 10.0;

        let atoms = self
            .atoms()
            .enumerate()
            .map(|(idx, atom)| Atom {
                position: Vec3::new(
                    atom.x * GRO_TO_VIEW_SCALE,
                    atom.y * GRO_TO_VIEW_SCALE,
                    atom.z * GRO_TO_VIEW_SCALE,
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

        Molecule { atoms, bonds: Vec::new() }
    }
}

impl To3dViewMolecule for GroFile {
    fn to_molecule(&self) -> Molecule {
        self.to_molecule_with_metadata(true)
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
