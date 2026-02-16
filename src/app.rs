use crate::parsing::{AtomRecord, PdbFile};
use crate::view_rs::To3dViewMolecule;
use egui;
use graphics::{EngineUpdates, Scene, winit::event::WindowEvent};
use moleucle_3dview_rs::OrbitalCamera;
use moleucle_3dview_rs::{
    CameraController, Molecule, MoleculeViewer, SelectedAtomRender, viewer::ViewerEvent,
};
use petgraph::algo::astar;
use petgraph::graph::UnGraph;
use rfd::FileDialog;
use std::path::PathBuf;

pub struct KuromameApp {
    // Viewer state
    pub viewer: MoleculeViewer<SelectedAtomRender>,
    pub controller: CameraController<OrbitalCamera>,

    // Data state
    pub pdb_file: Option<PdbFile>,
    pub current_file_path: Option<PathBuf>,

    // Selection state
    pub with_hbond_chk: bool,

    // UI state
    pub status_msg: String,
    pub show_edit_dialog: bool,
    pub new_res_name: String,
}

impl KuromameApp {
    pub fn new() -> Self {
        let mut viewer = MoleculeViewer::new();
        viewer.additional_render = Some(Box::new(SelectedAtomRender::new()));

        Self {
            viewer,
            controller: CameraController::new(),
            pdb_file: None,
            current_file_path: None,
            with_hbond_chk: false,
            status_msg: "Ready".to_string(),
            show_edit_dialog: false,
            new_res_name: "".to_string(),
        }
    }

    pub fn open_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("PDB Files", &["pdb", "ent", "cif"])
            .add_filter("MOL2 Files", &["mol2"])
            .pick_file()
        {
            self.load_file(path);
        }
    }

    pub fn load_file(&mut self, path: PathBuf) {
        if let Ok(content) = std::fs::read_to_string(&path) {
            if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
                match ext.to_lowercase().as_str() {
                    "pdb" | "ent" => {
                        let pdb = PdbFile::load(&content);
                        let mol = pdb.to_molecule();
                        self.viewer.set_molecule(mol);
                        self.pdb_file = Some(pdb);
                        self.current_file_path = Some(path);
                        self.status_msg = "Loaded PDB".to_string();
                    }
                    "mol2" => {
                        // Placeholder for Mol2
                        // let mol2 = Mol2File::load(&content);
                        self.status_msg = "MOL2 loading not fully wired yet".to_string();
                    }
                    _ => {
                        self.status_msg = "Unsupported file type".to_string();
                    }
                }
                if let Some(renderer) = &mut self.viewer.additional_render {
                    *renderer = Box::new(SelectedAtomRender::new());
                }
            }
        } else {
            self.status_msg = "Failed to read file".to_string();
        }
    }

    pub fn update_scene(&mut self, scene: &mut Scene) {
        self.viewer.update_scene(scene);
        self.controller.update_scene_camera(scene);
    }

    pub fn handle_event(
        &mut self,
        event: &WindowEvent,
        scene: &Scene,
    ) -> (Option<ViewerEvent>, EngineUpdates) {
        self.controller.handle_event(event, scene, &self.viewer)
    }

    pub fn render_ui(&mut self, ctx: &egui::Context) {
        let mut open_edit_dialog = self.show_edit_dialog;

        if open_edit_dialog {
            egui::Window::new("Edit Residue Name")
                .open(&mut open_edit_dialog)
                .show(ctx, |ui| {
                    ui.label("Enter new residue name (3 letters):");
                    ui.text_edit_singleline(&mut self.new_res_name);
                    if ui.button("Apply").clicked() {
                        self.apply_res_name_change();
                        self.show_edit_dialog = false;
                    }
                });
            self.show_edit_dialog = open_edit_dialog;
        }

        egui::SidePanel::right("control_panel")
            .resizable(true)
            .default_width(300.0)
            .show(ctx, |ui| {
                ui.heading("Controls");
                if ui.button("Open Molecule File (PDB/MOL2)").clicked() {
                    self.open_file();
                }

                if let Some(path) = &self.current_file_path {
                    ui.label(format!("Loaded: {:?}", path.file_name().unwrap()));
                } else {
                    ui.label("No file loaded");
                }

                ui.checkbox(&mut self.with_hbond_chk, "Select with connected hydrogens");

                ui.separator();
                ui.label("Selected atoms:");

                let mut selected_indices = Vec::new();
                if let Some(renderer) = &self.viewer.additional_render {
                    selected_indices = renderer.selected_atoms.clone();
                }

                ui.label(format!("Count: {}", selected_indices.len()));

                egui::ScrollArea::vertical()
                    .max_height(400.0)
                    .show(ui, |ui| {
                        if let Some(pdb) = &self.pdb_file {
                            let atoms_vec: Vec<&AtomRecord> = pdb.atoms().collect();

                            for &idx in &selected_indices {
                                if let Some(atom) = atoms_vec.get(idx) {
                                    ui.label(format!(
                                        "{} {} {} {} chain {}",
                                        atom.serial,
                                        atom.name,
                                        atom.res_name,
                                        atom.res_seq,
                                        atom.chain_id
                                    ));
                                } else {
                                    ui.label(format!("Index {} out of bounds", idx));
                                }
                            }
                        } else {
                            ui.label("No PDB loaded to map indices");
                        }
                    });

                if ui.button("Clear Selection").clicked() {
                    if let Some(renderer) = &mut self.viewer.additional_render {
                        *renderer = Box::new(SelectedAtomRender::new());
                        self.viewer.dirty = true;
                    }
                }

                let can_select_path = selected_indices.len() == 2 && !self.with_hbond_chk;
                if ui
                    .add_enabled(
                        can_select_path,
                        egui::Button::new("Select atoms between (2 atoms)"),
                    )
                    .clicked()
                {
                    self.select_shortest_path(selected_indices[0], selected_indices[1]);
                }

                if ui
                    .add_enabled(
                        !selected_indices.is_empty(),
                        egui::Button::new("Change selected residues' resname..."),
                    )
                    .clicked()
                {
                    self.show_edit_dialog = true;
                    self.new_res_name = "ALA".to_string(); // default
                }

                if ui.button("Export PDB (Save)").clicked() {
                    self.export_pdb();
                }

                ui.label(&self.status_msg);
            });
    }

    fn select_shortest_path(&mut self, start: usize, end: usize) {
        let mut path_found = false;

        if let Some(mol) = &self.viewer.molecule {
            let mut graph = UnGraph::<usize, ()>::new_undirected();
            let nodes: Vec<_> = (0..mol.atoms.len()).map(|i| graph.add_node(i)).collect();

            for bond in &mol.bonds {
                graph.add_edge(nodes[bond.atom_a], nodes[bond.atom_b], ());
            }

            if let Some((_, path)) = astar(
                &graph,
                nodes[start],
                |finish| finish == nodes[end],
                |_| 1,
                |_| 0,
            ) {
                // Collect indices to add
                let mut indices_to_add = Vec::new();
                for node_idx in path {
                    let idx = graph[node_idx];
                    indices_to_add.push(idx);
                }

                // If H-bond check is on, collect connected hydrogens
                if self.with_hbond_chk {
                    let mut hydrogens = Vec::new();
                    for &idx in &indices_to_add {
                        Self::collect_connected_hydrogens(idx, mol, &mut hydrogens);
                    }
                    // Avoid duplicates
                    for h in hydrogens {
                        if !indices_to_add.contains(&h) {
                            indices_to_add.push(h);
                        }
                    }
                }

                // Update renderer
                if let Some(renderer) = &mut self.viewer.additional_render {
                    for idx in indices_to_add {
                        if !renderer.selected_atoms.contains(&idx) {
                            renderer.selected_atoms.push(idx);
                        }
                    }
                }
                path_found = true;
            }
        }

        if path_found {
            self.viewer.dirty = true;
        }
    }

    fn collect_connected_hydrogens(atom_idx: usize, mol: &Molecule, out: &mut Vec<usize>) {
        for bond in &mol.bonds {
            let neighbor = if bond.atom_a == atom_idx {
                Some(bond.atom_b)
            } else if bond.atom_b == atom_idx {
                Some(bond.atom_a)
            } else {
                None
            };

            if let Some(n_idx) = neighbor {
                if let Some(atom) = mol.atoms.get(n_idx) {
                    if atom.element == "H" && !out.contains(&n_idx) {
                        out.push(n_idx);
                    }
                }
            }
        }
    }

    fn apply_res_name_change(&mut self) {
        let new_name = self.new_res_name.trim().to_uppercase();
        if new_name.len() > 3 {
            self.status_msg = "Residue name too long".to_string();
            return;
        }
        let new_name = format!("{:>3}", new_name); // Pad to 3 chars

        if let Some(pdb) = &mut self.pdb_file {
            let mut indices_to_update = Vec::new();
            if let Some(renderer) = &self.viewer.additional_render {
                indices_to_update = renderer.selected_atoms.clone();
            }

            // Update PDB atoms
            let mut atoms_vec: Vec<&mut AtomRecord> = pdb.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.res_name = new_name.clone();
                }
            }

            // Update viewer
            // We need to reload the molecule from updated PDB
            let mol = pdb.to_molecule();
            self.viewer.set_molecule(mol);

            // Clear selection?
            if let Some(renderer) = &mut self.viewer.additional_render {
                *renderer = Box::new(SelectedAtomRender::new());
            }
            self.status_msg = "Residue names updated".to_string();
        }
    }

    fn export_pdb(&self) {
        if let Some(pdb) = &self.pdb_file {
            if let Some(path) = FileDialog::new()
                .set_file_name("edited.pdb")
                .add_filter("PDB", &["pdb"])
                .save_file()
            {
                let content = pdb.dump();
                if let Ok(_) = std::fs::write(path, content) {
                    // self.status_msg = "Saved".to_string(); // Need &mut self
                }
            }
        }
    }
}
