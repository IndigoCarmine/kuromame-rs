use kuromame_rs::parsing::{AtomRecord, GroFile, PdbFile, TopFile};
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
    assert_eq!(gro.atoms().next().map(|a| a.res_name.trimmed()), Some("MOL"));

    if let Some(first) = gro.atoms_mut().next() {
        first.set_res_name("LIG");
    }

    let dumped = gro.dump();
    let atom_line = dumped.lines().nth(2).unwrap_or("");
    assert_eq!(atom_line.get(5..10).map(str::trim), Some("LIG"));

    let mol = gro.to_molecule();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.atoms[0].res_name.as_deref(), Some("LIG"));
    assert!((mol.atoms[1].position.x - 0.1).abs() < 1e-6);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].atom_a, 0);
    assert_eq!(mol.bonds[0].atom_b, 1);
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
}

#[test]
fn test_top_matches_gro_for_bar148() {
    let top_content = include_str!("../Bar148.top");
    let gro_content = include_str!("../Bar148.gro");

    let top = TopFile::load(top_content);
    let gro = GroFile::load(gro_content);

    assert_eq!(top.atoms().count(), 119);
    assert_eq!(gro.atoms().count(), 119);

    let comparison = top.compare_with_gro(&gro);
    assert!(comparison.matches(), "{:?}", comparison);

    let relaxed = top.compare_with_gro_relaxed(&gro);
    assert!(relaxed.matches(), "{:?}", relaxed);

    let dumped = top.dump();
    let reparsed = TopFile::load(&dumped);
    assert_eq!(reparsed.atoms().count(), top.atoms().count());
    assert_eq!(reparsed.bonds().count(), top.bonds().count());
}

#[test]
fn test_top_resname_can_be_rewritten_to_aaa() {
    let top_content = include_str!("../Bar148.top");

    let mut top = TopFile::load(top_content);
    let original_bond_count = top.bonds().count();
    assert!(top.atoms().next().is_some());

    for atom in top.atoms_mut() {
        atom.set_res_name("AAA");
    }

    let dumped = top.dump();

    //save file for debugging
    std::fs::write("debug.top", &dumped).unwrap();
    

    let reparsed = TopFile::load(&dumped);

    assert!(reparsed.atoms().all(|atom| atom.res == "AAA"));
    assert_eq!(reparsed.atoms().count(), 119);
    assert_eq!(reparsed.bonds().count(), original_bond_count);

}

#[test]
fn test_top_roundtrip_is_exact_when_unmodified() {
    let top_content = include_str!("../Bar148.top");
    let top = TopFile::load(top_content);
    let dumped = top.dump();

    assert_eq!(dumped, top_content);
}

#[test]
fn test_top_resname_sync_from_gro_passes_two_checks_and_updates() {
    let top_content = include_str!("../Bar148.top");
    let mut top = TopFile::load(top_content);
    let mut gro = GroFile::load(include_str!("../Bar148.gro"));

    for atom in gro.atoms_mut() {
        atom.set_res_name("AAA");
    }

    let updated = top.sync_resnames_from_gro(&gro).expect("sync should pass");

    assert_eq!(updated, top.atoms().count());
    assert!(top.atoms().all(|atom| atom.res == "AAA"));
}

#[test]
fn test_top_resname_sync_from_gro_fails_when_atom_order_mismatch() {
    let top_content = include_str!("../Bar148.top");
    let mut top = TopFile::load(top_content);
    let gro_text = include_str!("../Bar148.gro");
    let bad_gro_text = gro_text.replacen("C1", "ZZ", 1);
    let gro = GroFile::load(&bad_gro_text);

    let result = top.sync_resnames_from_gro(&gro);
    assert!(result.is_err());
}
