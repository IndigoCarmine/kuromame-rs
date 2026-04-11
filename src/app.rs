use eframe::egui::{self, PointerButton, Sense};
use lin_alg::f32::{Vec2, Vec3};
use moleucle_3dview_rs::{
    camera, Atom, Camera, CameraController, Molecule, MoleculeViewer, OffscreenRenderer,
    RenderStyle, SelectedAtomRender,
};
use crate::parsing::{AtomRecord, GroFile, PdbFile};
use rfd::FileDialog;
use std::path::PathBuf;

fn hsl_to_rgb(h: f32, s: f32, l: f32) -> (f32, f32, f32) {
    if s <= 0.0 {
        return (l, l, l);
    }

    let q = if l < 0.5 {
        l * (1.0 + s)
    } else {
        l + s - (l * s)
    };
    let p = 2.0 * l - q;

    fn hue_to_rgb(p: f32, q: f32, mut t: f32) -> f32 {
        if t < 0.0 {
            t += 1.0;
        }
        if t > 1.0 {
            t -= 1.0;
        }
        if t < 1.0 / 6.0 {
            return p + (q - p) * 6.0 * t;
        }
        if t < 1.0 / 2.0 {
            return q;
        }
        if t < 2.0 / 3.0 {
            return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
        }
        p
    }

    (
        hue_to_rgb(p, q, h + 1.0 / 3.0),
        hue_to_rgb(p, q, h),
        hue_to_rgb(p, q, h - 1.0 / 3.0),
    )
}

fn color_by_res_name(atom: &Atom, is_selected: bool) -> (f32, f32, f32) {
    if is_selected {
        return (1.0, 0.0, 0.0);
    }

    let key = atom
        .res_name
        .as_deref()
        .filter(|s| !s.trim().is_empty())
        .unwrap_or(atom.element.as_str());

    // Deterministic hash so the same residue name gets the same color each run.
    let hash = key
        .bytes()
        .fold(2166136261u32, |acc, b| (acc ^ (b as u32)).wrapping_mul(16777619));
    let hue = (hash % 360) as f32 / 360.0;
    hsl_to_rgb(hue, 0.65, 0.52)
}

pub struct KuromameApp {
    // Viewer state
    pub viewer: MoleculeViewer<SelectedAtomRender>,
    pub controller: CameraController<camera::OrbitalCamera>,
    pub offscreen: OffscreenRenderer,
    pub render_state: Option<egui_wgpu::RenderState>,

    // Data state
    pub pdb_file: Option<PdbFile>,
    pub gro_file: Option<GroFile>,
    pub current_file_path: Option<PathBuf>,

    // Selection state
    pub with_hbond_chk: bool,
    pub selected_atom_indices: Vec<usize>,

    // UI state
    pub status_msg: String,
    pub show_edit_dialog: bool,
    pub new_res_name: String,
}

impl KuromameApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut viewer = MoleculeViewer::new();
        viewer.set_color_fn(color_by_res_name);
        viewer.additional_render = Some(Box::new(SelectedAtomRender::new()));

        Self {
            viewer,
            controller: CameraController::new(),
            offscreen: OffscreenRenderer::new(),
            render_state: cc.wgpu_render_state.clone(),
            pdb_file: None,
            gro_file: None,
            current_file_path: None,
            with_hbond_chk: false,
            selected_atom_indices: Vec::new(),
            status_msg: "Ready".to_string(),
            show_edit_dialog: false,
            new_res_name: String::new(),
        }
    }

    fn sync_selection_to_renderer(&mut self) {
        self.viewer
            .set_selected_atoms(self.selected_atom_indices.clone());
        self.viewer.dirty = true;
    }

    fn toggle_selected_atom(&mut self, atom_index: usize) -> bool {
        let was_selected = self.selected_atom_indices.contains(&atom_index);
        if was_selected {
            self.selected_atom_indices.retain(|&i| i != atom_index);
        } else {
            self.selected_atom_indices.push(atom_index);
        }

        !was_selected
    }

    fn add_connected_hydrogens(&mut self, atom_index: usize) {
        let Some(mol) = &self.viewer.molecule else {
            return;
        };

        let mut targets = Self::collect_connected_hydrogens(atom_index, mol);
        targets.sort_unstable();
        targets.dedup();

        for idx in targets {
            if !self.selected_atom_indices.contains(&idx) {
                self.selected_atom_indices.push(idx);
            }
        }
    }

    fn remove_connected_hydrogens(&mut self, atom_index: usize) {
        let Some(mol) = &self.viewer.molecule else {
            return;
        };

        let mut targets = Self::collect_connected_hydrogens(atom_index, mol);
        targets.sort_unstable();
        targets.dedup();

        self.selected_atom_indices.retain(|idx| !targets.contains(idx));
    }

    pub fn open_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("PDB Files", &["pdb", "ent", "cif"])
            .add_filter("MOL2 Files", &["mol2"])
            .add_filter("GRO Files", &["gro"])
            .pick_file()
        {
            self.load_file(path);
        }
    }

    pub fn load_file(&mut self, path: PathBuf) {
        let Some(ext) = path.extension().and_then(|s| s.to_str()) else {
            self.status_msg = "Unsupported file type".to_string();
            return;
        };

        match ext.to_lowercase().as_str() {
            "pdb" | "ent" => {
                match Molecule::from_pdb(&path) {
                    Ok(mol) => {
                        self.set_molecule_and_frame(mol);
                        self.pdb_file = std::fs::read_to_string(&path).ok().map(|content| PdbFile::load(&content));
                        self.gro_file = None;
                        self.current_file_path = Some(path);
                        self.status_msg = "Loaded PDB".to_string();
                    }
                    Err(_) => {
                        self.status_msg = "Failed to load PDB file".to_string();
                    }
                }
            }
            "mol2" => {
                if let Ok(mol) = Molecule::from_mol2(&path) {
                    let pdb_from_mol2 = PdbFile::from_molecule(&mol);
                    self.set_molecule_and_frame(mol);
                    self.pdb_file = Some(pdb_from_mol2);
                    self.gro_file = None;
                    self.current_file_path = Some(path);
                    self.status_msg = "Loaded MOL2".to_string();
                } else {
                    self.status_msg = "Failed to load MOL2 file".to_string();
                }
            }
            "gro" => {
                const LARGE_GRO_THRESHOLD_BYTES: u64 = 20 * 1024 * 1024;
                let large_gro = std::fs::metadata(&path)
                    .map(|m| m.len() >= LARGE_GRO_THRESHOLD_BYTES)
                    .unwrap_or(false);

                match GroFile::load_from_path(&path) {
                    Ok(gro) => {
                        let mol = gro.to_molecule_with_metadata(!large_gro);
                        self.set_molecule_and_frame(mol);
                        self.gro_file = if large_gro { None } else { Some(gro) };
                        self.pdb_file = None;
                        self.current_file_path = Some(path);
                        self.status_msg = if large_gro {
                            "Loaded GRO (compact mode: reduced memory, GRO edit/export disabled)"
                                .to_string()
                        } else {
                            "Loaded GRO".to_string()
                        };
                    }
                    Err(_) => {
                        self.status_msg = "Failed to load GRO file".to_string();
                    }
                }
            }
            _ => {
                self.status_msg = "Unsupported file type".to_string();
            }
        }

        if let Some(renderer) = &mut self.viewer.additional_render {
            *renderer = Box::new(SelectedAtomRender::new());
        }
        self.selected_atom_indices.clear();
    }

    fn set_molecule_and_frame(&mut self, molecule: Molecule) {
        self.viewer.set_molecule(molecule);

        let Some(mol) = self.viewer.molecule.as_ref() else {
            return;
        };

        if mol.atoms.is_empty() {
            return;
        }

        let mut min = mol.atoms[0].position;
        let mut max = mol.atoms[0].position;

        for atom in &mol.atoms[1..] {
            min.x = min.x.min(atom.position.x);
            min.y = min.y.min(atom.position.y);
            min.z = min.z.min(atom.position.z);
            max.x = max.x.max(atom.position.x);
            max.y = max.y.max(atom.position.y);
            max.z = max.z.max(atom.position.z);
        }

        let center = (min + max) * 0.5;
        let radius = (max - min).magnitude() * 0.5;
        let camera_radius = if radius > 0.0 {
            let fit_distance = radius / (self.controller.camera.fov_y() * 0.5).tan();
            fit_distance.max(10.0) * 1.2
        } else {
            10.0
        };

        self.controller.camera.look_at(
            center + Vec3::new(0.0, 0.0, camera_radius),
            center,
            Vec3::new(0.0, 1.0, 0.0),
        );
        self.controller.camera.far = camera_radius + radius * 2.0 + 10.0;
    }

    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped_path = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .find_map(|file| file.path.clone())
        });

        if let Some(path) = dropped_path {
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
                if ui.button("Open Molecule File (PDB/MOL2/GRO)").clicked() {
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

                let selected_indices = self.selected_atom_indices.clone();

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
                        } else if let Some(gro) = &self.gro_file {
                            let atoms_vec: Vec<_> = gro.atoms().collect();

                            for &idx in &selected_indices {
                                if let Some(atom) = atoms_vec.get(idx) {
                                    ui.label(format!(
                                        "{} {} {} {}",
                                        atom.atom_id,
                                        atom.atom_name.trimmed(),
                                        atom.res_name.trimmed(),
                                        atom.res_num,
                                    ));
                                } else {
                                    ui.label(format!("Index {} out of bounds", idx));
                                }
                            }
                        } else {
                            ui.label("No structure loaded to map indices");
                        }
                    });

                if ui.button("Clear Selection").clicked() {
                    self.selected_atom_indices.clear();
                    self.sync_selection_to_renderer();
                }

                let can_select_path = selected_indices.len() == 2;
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

                if ui.button("Export (Save)").clicked() {
                    self.export_structure();
                }

                ui.label(&self.status_msg);
            });
    }

    fn select_shortest_path(&mut self, start: usize, end: usize) {
        let Some(mol) = &self.viewer.molecule else {
            return;
        };

        // Find all atoms on all simple paths between start and end
        let atoms_on_path = Self::find_atoms_between_dfs(mol, start, end);

        // Toggle only the atoms on the path first.
        for idx in atoms_on_path.iter().copied() {
            self.toggle_selected_atom(idx);
        }

        if self.with_hbond_chk {
            for idx in atoms_on_path {
                if self.selected_atom_indices.contains(&idx) {
                    self.add_connected_hydrogens(idx);
                }
            }
        }

        self.sync_selection_to_renderer();
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
                .or_default()
                .insert(bond.atom_b);
            adj.entry(bond.atom_b)
                .or_default()
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

    fn collect_connected_hydrogens(atom_idx: usize, mol: &Molecule) -> Vec<usize> {
        let mut hydrogens = Vec::new();
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
                    if atom.element.starts_with("H") && !hydrogens.contains(&n_idx) {
                        hydrogens.push(n_idx);
                    }
                }
            }
        }
        hydrogens
    }

    fn apply_res_name_change(&mut self) {
        let new_name = self.new_res_name.trim().to_uppercase();
        if new_name.len() > 3 {
            self.status_msg = "Residue name too long".to_string();
            return;
        }
        let new_name = format!("{:>3}", new_name); // Pad to 3 chars
        let indices_to_update = self.selected_atom_indices.clone();

        if indices_to_update.is_empty() {
            self.status_msg = "No atoms selected".to_string();
            return;
        }

        if let Some(pdb) = &mut self.pdb_file {
            // Update PDB atoms
            let mut atoms_vec: Vec<&mut AtomRecord> = pdb.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.res_name = new_name.clone();
                }
            }
        }

        if let Some(gro) = &mut self.gro_file {
            // In GRO, the 2nd field is residue name (resname).
            let mut atoms_vec: Vec<_> = gro.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.set_res_name(&new_name);
                }
            }
        }

        // Update viewer atoms in-place so we keep the original molecule graph
        // (bond topology from loader/inference) and avoid disappearing geometry.
        if let Some(mol) = &mut self.viewer.molecule {
            for &idx in &indices_to_update {
                if let Some(atom) = mol.atoms.get_mut(idx) {
                    atom.res_name = Some(new_name.clone());
                }
            }
        }
        self.viewer.dirty = true;

        // Clear selection
        self.selected_atom_indices.clear();
        self.sync_selection_to_renderer();
        self.status_msg = "Residue names updated".to_string();
    }

    fn export_structure(&mut self) {
        if self.gro_file.is_none() && self.pdb_file.is_none() {
            if let Some(mol) = &self.viewer.molecule {
                self.pdb_file = Some(PdbFile::from_molecule(mol));
            }
        }

        if let Some(gro) = &self.gro_file {
            if let Some(path) = FileDialog::new()
                .set_file_name("edited.gro")
                .add_filter("GRO", &["gro"])
                .save_file()
            {
                let content = gro.dump();
                if std::fs::write(path, content).is_ok() {
                    self.status_msg = "Exported GRO".to_string();
                } else {
                    self.status_msg = "Failed to export GRO".to_string();
                }
            }
        } else if let Some(pdb) = &mut self.pdb_file {
            if let Some(path) = FileDialog::new()
                .set_file_name("edited.pdb")
                .add_filter("PDB", &["pdb"])
                .save_file()
            {
                let content = pdb.dump();
                if std::fs::write(path, content).is_ok() {
                    self.status_msg = "Exported PDB".to_string();
                } else {
                    self.status_msg = "Failed to export PDB".to_string();
                }
            }
        } else {
            self.status_msg = "No molecule loaded to export".to_string();
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
                ui.label("Drop a PDB/MOL2/GRO file anywhere in the window to load it");
                ui.label(format!("Selected atoms: {:?}", self.selected_atom_indices));
                
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
            let view_proj = self.controller.camera.view_projection().data;
            if let Err(err) = self.offscreen.render_frame(
                render_state,
                self.viewer.molecule.as_ref(),
                &self.selected_atom_indices,
                view_proj,
                self.viewer.color_fn,
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

            // Track hovered atom residue name for tooltip rendering.
            let mut hovered_resname: Option<String> = None;

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

            // Handle hover picking for residue name tooltip.
            if response.hovered() {
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
                            if let Some(mol) = &self.viewer.molecule {
                                if let Some(atom) = mol.atoms.get(i) {
                                    let resname = atom
                                        .res_name
                                        .as_deref()
                                        .map(str::trim)
                                        .filter(|s| !s.is_empty())
                                        .unwrap_or("-");
                                    hovered_resname = Some(format!("resname: {}", resname));
                                }
                            }
                        }
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
                            let selected_now = self.toggle_selected_atom(i);
                            if self.with_hbond_chk {
                                if selected_now {
                                    self.add_connected_hydrogens(i);
                                } else {
                                    self.remove_connected_hydrogens(i);
                                }
                            }
                            self.sync_selection_to_renderer();
                        }
                    }
                }
            }

            if let Some(label) = hovered_resname {
                response.clone().on_hover_text(label);
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
