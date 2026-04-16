use crate::parsing::GroFile;
use std::collections::BTreeSet;

#[derive(Debug, Clone)]
pub struct TopAtomRecord {
    pub nr: usize,
    pub atom_type: String,
    pub resi: i32,
    pub res: String,
    pub atom: String,
    pub cgnr: i32,
    pub charge: f32,
    pub mass: f32,
    pub comment: Option<String>,
}

impl TopAtomRecord {
    fn split_comment(line: &str) -> (&str, Option<String>) {
        if let Some((head, tail)) = line.split_once(';') {
            (head.trim_end(), Some(tail.trim().to_string()))
        } else {
            (line.trim_end(), None)
        }
    }

    pub fn from_line(line: &str) -> Option<Self> {
        let (data, comment) = Self::split_comment(line);
        let parts: Vec<&str> = data.split_whitespace().collect();
        if parts.len() < 8 {
            return None;
        }

        Some(Self {
            nr: parts[0].parse().ok()?,
            atom_type: parts[1].to_string(),
            resi: parts[2].parse().ok()?,
            res: parts[3].to_string(),
            atom: parts[4].to_string(),
            cgnr: parts[5].parse().ok()?,
            charge: parts[6].parse().ok()?,
            mass: parts[7].parse().ok()?,
            comment,
        })
    }

    pub fn to_line(&self) -> String {
        let mut line = format!(
            "{:>6} {:<8}{:>6} {:<5} {:<5}{:>6} {:>12.6} {:>11.5}",
            self.nr, self.atom_type, self.resi, self.res, self.atom, self.cgnr, self.charge, self.mass
        );
        if let Some(comment) = &self.comment {
            if !comment.is_empty() {
                line.push_str(" ; ");
                line.push_str(comment);
            }
        }
        line
    }

    pub fn set_res_name(&mut self, value: &str) {
        self.res = value.trim().to_string();
    }
}

#[derive(Debug, Clone)]
pub struct TopBondRecord {
    pub ai: usize,
    pub aj: usize,
    pub funct: i32,
    pub r: Option<f32>,
    pub k: Option<f32>,
    pub comment: Option<String>,
}

impl TopBondRecord {
    fn split_comment(line: &str) -> (&str, Option<String>) {
        if let Some((head, tail)) = line.split_once(';') {
            (head.trim_end(), Some(tail.trim().to_string()))
        } else {
            (line.trim_end(), None)
        }
    }

    pub fn from_line(line: &str) -> Option<Self> {
        let (data, comment) = Self::split_comment(line);
        let parts: Vec<&str> = data.split_whitespace().collect();
        if parts.len() < 3 {
            return None;
        }

        Some(Self {
            ai: parts[0].parse().ok()?,
            aj: parts[1].parse().ok()?,
            funct: parts[2].parse().ok()?,
            r: parts.get(3).and_then(|value| value.parse().ok()),
            k: parts.get(4).and_then(|value| value.parse().ok()),
            comment,
        })
    }

    pub fn to_line(&self) -> String {
        let mut line = format!("{:>6}{:>7}{:>6}", self.ai, self.aj, self.funct);
        if let Some(r) = self.r {
            line.push_str(&format!("{:>13.4e}", r));
        }
        if let Some(k) = self.k {
            line.push_str(&format!("{:>13.4e}", k));
        }
        if let Some(comment) = &self.comment {
            if !comment.is_empty() {
                line.push_str(" ; ");
                line.push_str(comment);
            }
        }
        line
    }
}

#[derive(Debug, Clone)]
pub enum TopLine {
    SectionHeader(String),
    Atom(TopAtomRecord),
    Bond(TopBondRecord),
    Other(String),
    Empty,
}

#[derive(Debug, Clone, Default)]
pub struct TopFile {
    pub lines: Vec<TopLine>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TopGroComparison {
    pub atom_count_match: bool,
    pub atom_order_match: bool,
    pub bond_count_match: bool,
    pub bond_connectivity_match: bool,
}

impl TopGroComparison {
    pub fn matches(&self) -> bool {
        self.atom_count_match
            && self.atom_order_match
            && self.bond_count_match
            && self.bond_connectivity_match
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TopGroRelaxedComparison {
    pub atom_count_match: bool,
    pub atom_name_order_match: bool,
}

impl TopGroRelaxedComparison {
    pub fn matches(&self) -> bool {
        self.atom_count_match && self.atom_name_order_match
    }
}

impl TopFile {
    pub fn load(content: &str) -> Self {
        let mut lines = Vec::new();
        let mut current_section = String::new();

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                lines.push(TopLine::Empty);
                continue;
            }

            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                current_section = trimmed
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .trim()
                    .to_ascii_lowercase();
                lines.push(TopLine::SectionHeader(line.to_string()));
                continue;
            }

            match current_section.as_str() {
                "atoms" => {
                    if let Some(atom) = TopAtomRecord::from_line(line) {
                        lines.push(TopLine::Atom(atom));
                    } else {
                        lines.push(TopLine::Other(line.to_string()));
                    }
                }
                "bonds" => {
                    if let Some(bond) = TopBondRecord::from_line(line) {
                        lines.push(TopLine::Bond(bond));
                    } else {
                        lines.push(TopLine::Other(line.to_string()));
                    }
                }
                _ => lines.push(TopLine::Other(line.to_string())),
            }
        }

        Self { lines }
    }

    pub fn dump(&self) -> String {
        let mut out = String::new();
        for line in &self.lines {
            match line {
                TopLine::SectionHeader(text) | TopLine::Other(text) => {
                    out.push_str(text);
                    out.push('\n');
                }
                TopLine::Atom(atom) => {
                    out.push_str(&atom.to_line());
                    out.push('\n');
                }
                TopLine::Bond(bond) => {
                    out.push_str(&bond.to_line());
                    out.push('\n');
                }
                TopLine::Empty => out.push('\n'),
            }
        }
        out
    }

    pub fn atoms(&self) -> impl Iterator<Item = &TopAtomRecord> {
        self.lines.iter().filter_map(|line| match line {
            TopLine::Atom(atom) => Some(atom),
            _ => None,
        })
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut TopAtomRecord> {
        self.lines.iter_mut().filter_map(|line| match line {
            TopLine::Atom(atom) => Some(atom),
            _ => None,
        })
    }

    pub fn bonds(&self) -> impl Iterator<Item = &TopBondRecord> {
        self.lines.iter().filter_map(|line| match line {
            TopLine::Bond(bond) => Some(bond),
            _ => None,
        })
    }

    pub fn compare_with_gro(&self, gro: &GroFile) -> TopGroComparison {
        let top_atoms: Vec<&TopAtomRecord> = self.atoms().collect();
        let gro_atoms: Vec<_> = gro.atoms().collect();

        let atom_count_match = top_atoms.len() == gro_atoms.len();
        let atom_order_match = atom_count_match
            && top_atoms.iter().zip(&gro_atoms).all(|(top_atom, gro_atom)| {
                top_atom.nr == gro_atom.atom_id
                    && top_atom.atom.trim() == gro_atom.atom_name.trimmed()
                    && top_atom.res.trim() == gro_atom.res_name.trimmed()
                    && top_atom.resi == gro_atom.res_num
            });

        let top_bonds: BTreeSet<(usize, usize)> = self
            .bonds()
            .map(|bond| normalized_pair(bond.ai, bond.aj))
            .collect();
        let gro_bonds: BTreeSet<(usize, usize)> = gro
            .to_molecule_with_metadata(true)
            .bonds
            .into_iter()
            .map(|bond| normalized_pair(bond.atom_a + 1, bond.atom_b + 1))
            .collect();

        TopGroComparison {
            atom_count_match,
            atom_order_match,
            bond_count_match: top_bonds.len() == gro_bonds.len(),
            bond_connectivity_match: top_bonds == gro_bonds,
        }
    }

    pub fn matches_gro(&self, gro: &GroFile) -> bool {
        self.compare_with_gro(gro).matches()
    }

    pub fn compare_with_gro_relaxed(&self, gro: &GroFile) -> TopGroRelaxedComparison {
        let top_atoms: Vec<&TopAtomRecord> = self.atoms().collect();
        let gro_atoms: Vec<_> = gro.atoms().collect();

        let atom_count_match = top_atoms.len() == gro_atoms.len();
        let atom_name_order_match = atom_count_match
            && top_atoms
                .iter()
                .zip(&gro_atoms)
                .all(|(top_atom, gro_atom)| top_atom.atom.trim() == gro_atom.atom_name.trimmed());

        TopGroRelaxedComparison {
            atom_count_match,
            atom_name_order_match,
        }
    }

    pub fn matches_gro_relaxed(&self, gro: &GroFile) -> bool {
        self.compare_with_gro_relaxed(gro).matches()
    }

    pub fn can_sync_resnames_from_gro(&self, gro: &GroFile) -> Result<(), String> {
        // Check 1: atom count and atom-name order must match.
        let relaxed = self.compare_with_gro_relaxed(gro);
        if !relaxed.matches() {
            return Err(format!(
                "Check 1 failed (atom count/order): {:?}",
                relaxed
            ));
        }

        // Check 2: inferred topology connectivity must match.
        let strict = self.compare_with_gro(gro);
        if !strict.bond_count_match || !strict.bond_connectivity_match {
            return Err(format!(
                "Check 2 failed (bond topology): bond_count_match={}, bond_connectivity_match={}",
                strict.bond_count_match, strict.bond_connectivity_match
            ));
        }

        Ok(())
    }

    pub fn sync_resnames_from_gro(&mut self, gro: &GroFile) -> Result<usize, String> {
        self.can_sync_resnames_from_gro(gro)?;

        let gro_res_names: Vec<String> = gro
            .atoms()
            .map(|atom| atom.res_name.trimmed().to_string())
            .collect();

        let mut updated = 0usize;
        for (top_atom, gro_res_name) in self.atoms_mut().zip(gro_res_names.iter()) {
            let trimmed = top_atom.res.trim();
            if trimmed != gro_res_name {
                top_atom.set_res_name(gro_res_name);
                updated += 1;
            }
        }

        Ok(updated)
    }
}

fn normalized_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}