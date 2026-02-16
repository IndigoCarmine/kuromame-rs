use kuromame_rs::parsing::{AtomRecord, PdbFile};
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
