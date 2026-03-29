use eframe::egui::{self, PointerButton, Sense};
use lin_alg::f32::Vec2;
use moleucle_3dview_rs::{
    camera, Camera, CameraController, Molecule, MoleculeViewer, OffscreenRenderer,
    RenderStyle, SelectedAtomRender,
};
use crate::parsing::{AtomRecord, PdbFile};
use crate::view_rs::To3dViewMolecule;
use rfd::FileDialog;
use std::path::PathBuf;

pub struct KuromameApp {
    // Viewer state
    pub viewer: MoleculeViewer<SelectedAtomRender>,
    pub controller: CameraController<camera::OrbitalCamera>,
    pub offscreen: OffscreenRenderer,
    pub render_state: Option<egui_wgpu::RenderState>,

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
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut viewer = MoleculeViewer::new();
        viewer.additional_render = Some(Box::new(SelectedAtomRender::new()));

        Self {
            viewer,
            controller: CameraController::new(),
            offscreen: OffscreenRenderer::new(),
            render_state: cc.wgpu_render_state.clone(),
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
                        if let Ok(mol) = Molecule::from_mol2(&path) {
                            self.viewer.set_molecule(mol);
                            self.current_file_path = Some(path);
                            self.status_msg = "Loaded MOL2".to_string();
                        } else {
                            self.status_msg = "Failed to load MOL2 file".to_string();
                        }
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

    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped_paths: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .filter_map(|file| file.path.clone())
                .collect()
        });

        if let Some(path) = dropped_paths.into_iter().next() {
            self.load_file(path);
        }
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

                let selected_indices = self.viewer.selected_atoms();

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
        if let Some(mol) = &self.viewer.molecule {
            // Find all atoms on all simple paths between start and end
            let atoms_on_path = Self::find_atoms_between_dfs(mol, start, end);
            
            let mut selected_indices = atoms_on_path.clone();

            // If H-bond check is on, collect connected hydrogens for each atom on path
            if self.with_hbond_chk {
                let mut hydrogens = Vec::new();
                for &idx in &atoms_on_path {
                    Self::collect_connected_hydrogens(idx, mol, &mut hydrogens);
                }
                // Avoid duplicates
                for h in hydrogens {
                    if !selected_indices.contains(&h) {
                        selected_indices.push(h);
                    }
                }
            }

            // Update renderer with all selected atoms
            if let Some(renderer) = &mut self.viewer.additional_render {
                for idx in selected_indices {
                    renderer.toggle_atom(idx);
                }
            }
            self.viewer.dirty = true;
        }
    }

    fn find_atoms_between_dfs(mol: &Molecule, start: usize, end: usize) -> Vec<usize> {
        if start == end {
            return vec![start];
        }

        // 1. Build adjacency map from bonds
        let mut adj: std::collections::HashMap<usize, std::collections::HashSet<usize>> =
            std::collections::HashMap::new();

        for bond in &mol.bonds {
            adj.entry(bond.atom_a)
                .or_insert_with(std::collections::HashSet::new)
                .insert(bond.atom_b);
            adj.entry(bond.atom_b)
                .or_insert_with(std::collections::HashSet::new)
                .insert(bond.atom_a);
        }

        // 2. DFS to find all simple paths and collect all atoms on any path
        let mut all_path_atoms = std::collections::HashSet::new();

        fn dfs(
            current: usize,
            target: usize,
            adj: &std::collections::HashMap<usize, std::collections::HashSet<usize>>,
            visited: &mut std::collections::HashSet<usize>,
            path: &mut Vec<usize>,
            all_path_atoms: &mut std::collections::HashSet<usize>,
        ) {
            if current == target {
                // Found a path - add all atoms in this path
                all_path_atoms.extend(path.iter());
                return;
            }

            if let Some(neighbors) = adj.get(&current) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor);
                        path.push(neighbor);
                        dfs(neighbor, target, adj, visited, path, all_path_atoms);
                        path.pop();
                        visited.remove(&neighbor);
                    }
                }
            }
        }

        let mut visited = std::collections::HashSet::new();
        visited.insert(start);
        let mut path = vec![start];
        dfs(start, end, &adj, &mut visited, &mut path, &mut all_path_atoms);

        // 3. Return as sorted vector
        let mut result: Vec<usize> = all_path_atoms.into_iter().collect();
        result.sort();
        result
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
                    // Check if element starts with "H" (matching Python's "H" in name check)
                    if atom.element.starts_with("H") && !out.contains(&n_idx) {
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
            let indices_to_update = self.viewer.selected_atoms();

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

            // Clear selection
            if let Some(renderer) = &mut self.viewer.additional_render {
                *renderer = Box::new(SelectedAtomRender::new());
            }
            self.status_msg = "Residue names updated".to_string();
        }
    }

    fn export_pdb(&mut self) {
        if let Some(pdb) = &mut self.pdb_file {
            if let Some(path) = FileDialog::new()
                .set_file_name("edited.pdb")
                .add_filter("PDB", &["pdb"])
                .save_file()
            {
                let content = pdb.dump();
                if std::fs::write(path, content).is_ok() {
                    // Success
                }
            }
        }
    }
}

impl eframe::App for KuromameApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_dropped_files(ctx);

        egui::TopBottomPanel::top("help")
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("LMB: pick atom  RMB drag: orbit  MMB/Shift+RMB drag: pan  Wheel: dolly");
                ui.label("Drop a PDB/MOL2 file anywhere in the window to load it");
                ui.label(format!("Selected atoms: {:?}", self.viewer.selected_atoms()));
                
                ui.horizontal(|ui| {
                    ui.label("Render Style:");
                    let mut style = self.offscreen.render_style();
                    ui.selectable_value(&mut style, RenderStyle::BallStick, "BallStick");
                    ui.selectable_value(&mut style, RenderStyle::Wireframe, "Wireframe");
                    self.offscreen.set_render_style(style);
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let Some(render_state) = &self.render_state else {
                ui.heading("WGPU backend is unavailable");
                ui.label("Start with the wgpu backend enabled in eframe.");
                return;
            };

            let available = ui.available_size_before_wrap();
            let width = available.x.max(1.0) as u32;
            let height = available.y.max(1.0) as u32;

            // Update camera aspect ratio
            self.controller
                .camera
                .set_aspect(width as f32 / height as f32);

            // Ensure rendering resources are available
            if let Err(err) = self.offscreen.ensure_resources(render_state, width, height) {
                ui.colored_label(egui::Color32::RED, format!("Offscreen init failed: {err}"));
                return;
            }

            // Render the molecule
            let selected = self.viewer.selected_atoms();
            let view_proj = self.controller.camera.view_projection().data;
            if let Err(err) = self.offscreen.render_frame(
                render_state,
                self.viewer.molecule.as_ref(),
                &selected,
                view_proj,
            ) {
                ui.colored_label(egui::Color32::RED, format!("Offscreen render failed: {err}"));
                return;
            }

            let Some(texture_id) = self.offscreen.texture_id() else {
                ui.colored_label(egui::Color32::RED, "No texture id registered");
                return;
            };

            // Display the rendered texture
            let response = ui.add(
                egui::Image::from_texture(egui::load::SizedTexture::new(
                    texture_id,
                    egui::vec2(width as f32, height as f32),
                ))
                .sense(Sense::click_and_drag()),
            );

            // Handle mouse wheel for dolly
            if response.hovered() {
                let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
                if scroll.abs() > f32::EPSILON {
                    self.controller.camera.dolly(scroll * 0.02);
                }
            }

            // Handle mouse drag for orbit/pan
            if response.hovered() {
                let (delta, sec_down, mid_down, shift_down) = ctx.input(|i| {
                    (
                        i.pointer.delta(),
                        i.pointer.button_down(PointerButton::Secondary),
                        i.pointer.button_down(PointerButton::Middle),
                        i.modifiers.shift,
                    )
                });

                if sec_down || mid_down {
                    if mid_down || shift_down {
                        // Pan
                        self.controller
                            .camera
                            .pan(Vec2::new(delta.x * 0.01, delta.y * 0.01));
                    } else {
                        // Orbit
                        self.controller.camera.orbit(delta.x * 0.005, delta.y * 0.005);
                    }
                }
            }

            // Handle mouse click for atom picking
            if response.clicked_by(PointerButton::Primary) {
                if let Some(pointer) = response.interact_pointer_pos() {
                    let local = pointer - response.rect.min;
                    let (ray_origin, ray_dir) = self.controller.camera.ray_from_screen(
                        local.x,
                        local.y,
                        response.rect.width().max(1.0),
                        response.rect.height().max(1.0),
                    );

                    if let Some(event) = self.viewer.pick(ray_origin, ray_dir) {
                        if let moleucle_3dview_rs::viewer::ViewerEvent::AtomClicked(i) = event {
                            if let Some(renderer) = &mut self.viewer.additional_render {
                                renderer.toggle_atom(i);
                                self.viewer.dirty = true;
                            }
                        }
                    }
                }
            }
        });

        // Render the UI panel
        self.render_ui(ctx);

        // Request continuous repaint
        ctx.request_repaint();
    }

    fn on_exit(&mut self, _ctx: Option<&eframe::glow::Context>) {
        if let Some(render_state) = &self.render_state {
            self.offscreen.free_egui_texture(render_state);
        }
    }
}
