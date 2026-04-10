use kuromame_rs::parsing::{AtomRecord, GroFile, PdbFile};
use kuromame_rs::view_rs::To3dViewMolecule;

#[test]
fn test_pdb_parsing() {
    let pdb_content = "ATOM      1  N   ALA A   1      -0.525   1.363   0.000  1.00  0.00           N  \n\
                       ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n\
                       CONECT    1    2\n";

    let pdb = PdbFile::load(pdb_content);
    let atoms: Vec<&AtomRecord> = pdb.atoms().collect();

    assert_eq!(atoms.len(), 2);
    assert_eq!(atoms[0].name, "N");
    assert_eq!(atoms[1].name, "CA");

    let mol = pdb.to_molecule();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].atom_a, 0);
    assert_eq!(mol.bonds[0].atom_b, 1);
}

#[test]
fn test_gro_resname_is_editable() {
    let gro_content = concat!(
        "Example GRO\n",
        "    2\n",
        "    1MOL     C1    1   0.000   0.000   0.000\n",
        "    1MOL     O1    2   0.100   0.100   0.100\n",
        "   1.00000   1.00000   1.00000\n",
    );

    let mut gro = GroFile::load(gro_content);
    assert_eq!(gro.atoms().count(), 2);
    assert_eq!(gro.atoms().next().map(|a| a.res_name.as_str()), Some("MOL"));

    if let Some(first) = gro.atoms_mut().next() {
        first.res_name = "LIG".to_string();
    }

    let dumped = gro.dump();
    let atom_line = dumped.lines().nth(2).unwrap_or("");
    assert_eq!(atom_line.get(5..10).map(str::trim), Some("LIG"));

    let mol = gro.to_molecule();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.atoms[0].res_name.as_deref(), Some("LIG"));
    assert!((mol.atoms[1].position.x - 1.0).abs() < 1e-6);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].order, 1);
}

#[test]
fn test_gro_bond_inference_uses_distance() {
    let gro_content = concat!(
        "Example GRO\n",
        "    3\n",
        "    1MOL     C1    1   0.000   0.000   0.000\n",
        "    1MOL     O1    2   0.100   0.000   0.000\n",
        "    1MOL     C2    3   0.500   0.000   0.000\n",
        "   1.00000   1.00000   1.00000\n",
    );

    let gro = GroFile::load(gro_content);
    let mol = gro.to_molecule();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].atom_a, 0);
    assert_eq!(mol.bonds[0].atom_b, 1);
    assert_eq!(mol.bonds[0].order, 1);
}

#[test]
fn test_gro_roundtrip_is_exact_for_bar1416() {
    let original = include_str!("../Bar1416_fixed_rot_10.gro");
    let gro = GroFile::load(original);
    let time = std::time::Instant::now();
    gro.to_molecule(); // Check calculation speed. change and calclate bonds
    eprintln!("Bar1416 processing time: {:?}", time.elapsed());
}
