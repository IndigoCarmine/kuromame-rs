use eframe::egui::{self, PointerButton, Sense};
use lin_alg::f32::{Vec2, Vec3};
use moleucle_3dview_rs::{
    camera, Atom, Camera, CameraController, Molecule, MoleculeViewer, OffscreenRenderer,
    SelectedAtomRender,
};
use crate::parsing::{AtomRecord, GroFile, Mol2File, PdbFile, TopFile};
use crate::view_rs::To3dViewMolecule;
use rfd::FileDialog;
use std::path::PathBuf;

#[path = "app_ui.rs"]
mod app_ui;

struct LoadedDataState {
    pdb_file: Option<PdbFile>,
    gro_file: Option<GroFile>,
    top_file: Option<TopFile>,
    current_file_path: Option<PathBuf>,
    loaded_summary: String,
    is_modified: bool,
}

impl LoadedDataState {
    fn clear_structures(&mut self) {
        self.pdb_file = None;
        self.gro_file = None;
        self.top_file = None;
    }
}

struct SelectionState {
    with_hbond_chk: bool,
    selected_atom_indices: Vec<usize>,
}

struct UiState {
    status_msg: String,
    show_edit_dialog: bool,
    new_res_name: String,
    hovered_atom_info: String,
}

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

    data: LoadedDataState,
    selection: SelectionState,
    ui: UiState,
}

impl KuromameApp {
    fn apply_visual_theme(ctx: &egui::Context) {
        let primary = egui::Color32::from_rgb(19, 161, 152);
        let secondary = egui::Color32::from_rgb(241, 98, 69);

        let mut style = (*ctx.style()).clone();
        style.spacing.item_spacing = egui::vec2(10.0, 10.0);
        style.spacing.button_padding = egui::vec2(12.0, 8.0);

        style.text_styles.insert(
            egui::TextStyle::Heading,
            egui::FontId::new(24.0, egui::FontFamily::Proportional),
        );
        style.text_styles.insert(
            egui::TextStyle::Body,
            egui::FontId::new(18.0, egui::FontFamily::Proportional),
        );
        style.text_styles.insert(
            egui::TextStyle::Button,
            egui::FontId::new(18.0, egui::FontFamily::Proportional),
        );
        style.text_styles.insert(
            egui::TextStyle::Monospace,
            egui::FontId::new(16.0, egui::FontFamily::Monospace),
        );
        style.text_styles.insert(
            egui::TextStyle::Small,
            egui::FontId::new(15.0, egui::FontFamily::Proportional),
        );

        style.visuals.widgets.active.bg_fill = primary;
        style.visuals.widgets.hovered.bg_fill = secondary;
        style.visuals.widgets.active.fg_stroke.color = egui::Color32::WHITE;
        style.visuals.widgets.hovered.fg_stroke.color = egui::Color32::WHITE;
        style.visuals.selection.bg_fill = primary;
        style.visuals.hyperlink_color = secondary;

        ctx.set_style(style);
    }

    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut fonts = egui::FontDefinitions::default();
        fonts.font_data.insert(
            "material_icons".to_string(),
            egui::FontData::from_static(material_icons::FONT).into(),
        );
        if let Some(family) = fonts.families.get_mut(&egui::FontFamily::Proportional) {
            family.push("material_icons".to_string());
        }
        if let Some(family) = fonts.families.get_mut(&egui::FontFamily::Monospace) {
            family.push("material_icons".to_string());
        }
        cc.egui_ctx.set_fonts(fonts);
        Self::apply_visual_theme(&cc.egui_ctx);

        let mut viewer = MoleculeViewer::new();
        viewer.set_color_fn(color_by_res_name);
        viewer.additional_render = Some(Box::new(SelectedAtomRender::new()));

        Self {
            viewer,
            controller: CameraController::new(),
            offscreen: OffscreenRenderer::new(),
            render_state: cc.wgpu_render_state.clone(),
            data: LoadedDataState {
                pdb_file: None,
                gro_file: None,
                top_file: None,
                current_file_path: None,
                loaded_summary: "No file loaded".to_string(),
                is_modified: false,
            },
            selection: SelectionState {
                with_hbond_chk: false,
                selected_atom_indices: Vec::new(),
            },
            ui: UiState {
                status_msg: "Ready".to_string(),
                show_edit_dialog: false,
                new_res_name: String::new(),
                hovered_atom_info: "Hover an atom for details".to_string(),
            },
        }
    }

    fn sync_selection_to_renderer(&mut self) {
        self.viewer
            .set_selected_atoms(self.selection.selected_atom_indices.clone());
        self.viewer.dirty = true;
    }

    fn post_load_cleanup(&mut self) {
        if let Some(renderer) = &mut self.viewer.additional_render {
            *renderer = Box::new(SelectedAtomRender::new());
        }
        self.selection.selected_atom_indices.clear();
    }

    fn set_status(&mut self, msg: impl Into<String>) {
        self.ui.status_msg = msg.into();
    }

    fn set_loaded_summary(&mut self, summary: impl Into<String>) {
        self.data.loaded_summary = summary.into();
    }

    fn mark_modified(&mut self) {
        self.data.is_modified = true;
    }

    fn mark_clean(&mut self) {
        self.data.is_modified = false;
    }

    fn open_resname_dialog(&mut self) {
        self.ui.show_edit_dialog = true;
        self.ui.new_res_name = "ALA".to_string();
    }

    fn clear_selection(&mut self) {
        self.selection.selected_atom_indices.clear();
        self.sync_selection_to_renderer();
    }

    fn toggle_hbond_selection(&mut self) {
        self.selection.with_hbond_chk = !self.selection.with_hbond_chk;
    }

    fn handle_keyboard_shortcuts(&mut self, ctx: &egui::Context) {
        let shortcuts = ctx.input(|i| {
            let ctrl = i.modifiers.ctrl;
            let shift = i.modifiers.shift;
            (
                ctrl && !shift && i.key_pressed(egui::Key::O),
                ctrl && shift && i.key_pressed(egui::Key::O),
                ctrl && !shift && i.key_pressed(egui::Key::R),
                ctrl && !shift && i.key_pressed(egui::Key::S),
                ctrl && !shift && i.key_pressed(egui::Key::H),
                ctrl && !shift && i.key_pressed(egui::Key::B),
                ctrl && shift && i.key_pressed(egui::Key::A),
            )
        });

        if shortcuts.0 {
            self.open_file();
        }
        if shortcuts.1 {
            self.open_top_and_gro_for_resname_sync();
        }
        if shortcuts.2 {
            self.open_resname_dialog();
        }
        if shortcuts.3 {
            self.export_structure();
        }
        if shortcuts.4 {
            self.toggle_hbond_selection();
        }
        if shortcuts.5 && self.selection.selected_atom_indices.len() == 2 {
            self.select_shortest_path(
                self.selection.selected_atom_indices[0],
                self.selection.selected_atom_indices[1],
            );
        }
        if shortcuts.6 {
            self.clear_selection();
        }
    }

    fn hovered_atom_info(&self, atom_index: usize) -> Option<String> {
        let mut atom_name: Option<String> = None;
        let mut res_name: Option<String> = None;

        if let Some(mol) = &self.viewer.molecule {
            if let Some(atom) = mol.atoms.get(atom_index) {
                atom_name = atom
                    .name
                    .as_ref()
                    .map(|name| name.trim().to_string())
                    .filter(|name| !name.is_empty());

                if atom_name.is_none() && !atom.element.trim().is_empty() {
                    atom_name = Some(atom.element.trim().to_string());
                }

                res_name = atom
                    .res_name
                    .as_ref()
                    .map(|name| name.trim().to_string())
                    .filter(|name| !name.is_empty());
            }
        }

        if let Some(pdb) = &self.data.pdb_file {
            if let Some(atom) = pdb.atoms().nth(atom_index) {
                if atom_name.is_none() && !atom.name.trim().is_empty() {
                    atom_name = Some(atom.name.trim().to_string());
                }
                if res_name.is_none() && !atom.res_name.trim().is_empty() {
                    res_name = Some(atom.res_name.trim().to_string());
                }
            }
        }

        if let Some(gro) = &self.data.gro_file {
            if let Some(atom) = gro.atoms().nth(atom_index) {
                if atom_name.is_none() {
                    let name = atom.atom_name.trimmed();
                    if !name.is_empty() {
                        atom_name = Some(name.to_string());
                    }
                }
                if res_name.is_none() {
                    let name = atom.res_name.trimmed();
                    if !name.is_empty() {
                        res_name = Some(name.to_string());
                    }
                }
            }
        }

        if let Some(top) = &self.data.top_file {
            if let Some(atom) = top.atoms().nth(atom_index) {
                if atom_name.is_none() && !atom.atom.trim().is_empty() {
                    atom_name = Some(atom.atom.trim().to_string());
                }
                if res_name.is_none() && !atom.res.trim().is_empty() {
                    res_name = Some(atom.res.trim().to_string());
                }
            }
        }

        if atom_name.is_none() && res_name.is_none() {
            return None;
        }

        Some(format!(
            "Index={} AtomName={} Resname={}",
            atom_index,
            atom_name.unwrap_or_else(|| "-".to_string()),
            res_name.unwrap_or_else(|| "-".to_string())
        ))
    }

    fn sync_viewer_resnames_from_loaded_files(&mut self) {
        let viewer_atom_count = self
            .viewer
            .molecule
            .as_ref()
            .map(|mol| mol.atoms.len())
            .unwrap_or(0);

        if viewer_atom_count == 0 {
            return;
        }

        let top_resnames = self.data.top_file.as_ref().map(|top| {
            top.atoms()
                .map(|atom| atom.res.trim().to_string())
                .collect::<Vec<_>>()
        });

        let gro_resnames = self.data.gro_file.as_ref().map(|gro| {
            gro.atoms()
                .map(|atom| atom.res_name.trimmed().to_string())
                .collect::<Vec<_>>()
        });

        let pdb_resnames = self.data.pdb_file.as_ref().map(|pdb| {
            pdb.atoms()
                .map(|atom| atom.res_name.trim().to_string())
                .collect::<Vec<_>>()
        });

        let resnames: Vec<String> = if let Some(names) = top_resnames {
            if names.len() == viewer_atom_count {
                names
            } else {
                gro_resnames
                    .filter(|names| names.len() == viewer_atom_count)
                    .or_else(|| pdb_resnames.filter(|names| names.len() == viewer_atom_count))
                    .unwrap_or_default()
            }
        } else {
            gro_resnames
                .filter(|names| names.len() == viewer_atom_count)
                .or_else(|| pdb_resnames.filter(|names| names.len() == viewer_atom_count))
                .unwrap_or_default()
        };

        if resnames.is_empty() {
            return;
        }

        if let Some(mol) = &mut self.viewer.molecule {
            for (atom, name) in mol.atoms.iter_mut().zip(resnames.into_iter()) {
                atom.res_name = Some(name);
            }
            self.viewer.dirty = true;
        }
    }

    fn toggle_selected_atom(&mut self, atom_index: usize) -> bool {
        let was_selected = self.selection.selected_atom_indices.contains(&atom_index);
        if was_selected {
            self.selection
                .selected_atom_indices
                .retain(|&i| i != atom_index);
        } else {
            self.selection.selected_atom_indices.push(atom_index);
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
            if !self.selection.selected_atom_indices.contains(&idx) {
                self.selection.selected_atom_indices.push(idx);
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

        self.selection
            .selected_atom_indices
            .retain(|idx| !targets.contains(idx));
    }

    pub fn open_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("PDB Files", &["pdb", "ent", "cif"])
            .add_filter("MOL2 Files", &["mol2"])
            .add_filter("GRO Files", &["gro"])
            .add_filter("TOP Files", &["top"])
            .pick_file()
        {
            self.load_file(path);
        }
    }

    pub fn open_top_and_gro_for_resname_sync(&mut self) {
        let top_path = FileDialog::new()
            .add_filter("TOP Files", &["top"])
            .set_title("Select TOP file")
            .pick_file();
        let gro_path = FileDialog::new()
            .add_filter("GRO Files", &["gro"])
            .set_title("Select GRO file")
            .pick_file();

        match (top_path, gro_path) {
            (Some(top), Some(gro)) => self.load_top_and_gro_for_resname_sync(top, gro),
            _ => {
                self.set_status("TOP/GRO pair selection cancelled");
            }
        }
    }

    fn load_top_and_gro_for_resname_sync(&mut self, top_path: PathBuf, gro_path: PathBuf) {
        let top_name = top_path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        let gro_name = gro_path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        let top_content = match std::fs::read_to_string(&top_path) {
            Ok(content) => content,
            Err(_) => {
                self.set_status("Failed to read TOP file");
                return;
            }
        };

        let mut top = TopFile::load(&top_content);
        let gro = match GroFile::load_from_path(&gro_path) {
            Ok(gro) => gro,
            Err(_) => {
                self.set_status("Failed to read GRO file");
                return;
            }
        };

        match top.sync_resnames_from_gro(&gro) {
            Ok(updated) => {
                let mol = gro.to_molecule_with_metadata(true);
                self.set_molecule_and_frame(mol);
                self.data.clear_structures();
                self.data.top_file = Some(top);
                self.data.gro_file = Some(gro);
                self.data.current_file_path = Some(top_path);
                self.set_loaded_summary(format!("TOP+GRO: {} + {}", top_name, gro_name));
                self.mark_clean();
                self.set_status(format!(
                    "Loaded TOP+GRO. Check1/Check2 passed; synced {} TOP resnames from GRO",
                    updated
                ));
            }
            Err(err) => {
                self.set_status(format!("TOP/GRO consistency check failed: {}", err));
            }
        }

        self.post_load_cleanup();
    }

    pub fn load_file(&mut self, path: PathBuf) {
        let Some(ext) = path.extension().and_then(|s| s.to_str()) else {
            self.set_status("Unsupported file type");
            return;
        };

        match ext.to_lowercase().as_str() {
            "pdb" | "ent" => self.load_pdb_file(path),
            "mol2" => self.load_mol2_file(path),
            "gro" => self.load_gro_file(path),
            "top" => self.load_top_file(path),
            _ => {
                self.set_status("Unsupported file type");
            }
        }
        self.post_load_cleanup();
    }

    fn load_pdb_file(&mut self, path: PathBuf) {
        let file_name = path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        match std::fs::read_to_string(&path) {
            Ok(content) => {
                let pdb = PdbFile::load(&content);
                let mol = pdb.to_molecule();
                self.set_molecule_and_frame(mol);
                self.data.clear_structures();
                self.data.pdb_file = Some(pdb);
                self.data.current_file_path = Some(path);
                self.set_loaded_summary(format!("PDB: {}", file_name));
                self.mark_clean();
                self.set_status("Loaded PDB");
            }
            Err(_) => self.set_status("Failed to load PDB file"),
        }
    }

    fn load_mol2_file(&mut self, path: PathBuf) {
        let file_name = path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        match std::fs::read_to_string(&path) {
            Ok(content) => {
                let mol2 = Mol2File::load(&content);
                let mol = mol2.to_molecule();
                let pdb_from_mol2 = PdbFile::from_molecule(&mol);
                self.set_molecule_and_frame(mol);
                self.data.clear_structures();
                self.data.pdb_file = Some(pdb_from_mol2);
                self.data.current_file_path = Some(path);
                self.set_loaded_summary(format!("MOL2: {}", file_name));
                self.mark_clean();
                self.set_status("Loaded MOL2");
            }
            Err(_) => self.set_status("Failed to load MOL2 file"),
        }
    }

    fn load_gro_file(&mut self, path: PathBuf) {
        let file_name = path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        const LARGE_GRO_THRESHOLD_BYTES: u64 = 20 * 1024 * 1024;
        let large_gro = std::fs::metadata(&path)
            .map(|m| m.len() >= LARGE_GRO_THRESHOLD_BYTES)
            .unwrap_or(false);

        match GroFile::load_from_path(&path) {
            Ok(gro) => {
                let mol = gro.to_molecule_with_metadata(!large_gro);
                self.set_molecule_and_frame(mol);
                self.data.clear_structures();
                self.data.gro_file = if large_gro { None } else { Some(gro) };
                self.data.current_file_path = Some(path);
                self.set_loaded_summary(format!("GRO: {}", file_name));
                self.mark_clean();
                if large_gro {
                    self.set_status(
                        "Loaded GRO (compact mode: reduced memory, GRO edit/export disabled)",
                    );
                } else {
                    self.set_status("Loaded GRO");
                }
            }
            Err(_) => self.set_status("Failed to load GRO file"),
        }
    }

    fn load_top_file(&mut self, path: PathBuf) {
        let file_name = path.file_name().and_then(|name| name.to_str()).unwrap_or("unknown").to_string();
        match std::fs::read_to_string(&path).ok().map(|content| TopFile::load(&content)) {
            Some(top) => {
                let can_sync_with_loaded_gro = self
                    .data
                    .gro_file
                    .as_ref()
                    .map(|gro| {
                        let comp = top.compare_with_gro(gro);
                        comp.atom_count_match && comp.atom_order_match
                    })
                    .unwrap_or(false);

                self.data.pdb_file = None;
                self.data.top_file = Some(top);
                self.data.current_file_path = Some(path);
                self.set_loaded_summary(format!("TOP: {}", file_name));
                self.mark_clean();

                if can_sync_with_loaded_gro {
                    let top_resnames: Vec<String> = self
                        .data
                        .top_file
                        .as_ref()
                        .map(|top_file| {
                            top_file
                                .atoms()
                                .map(|atom| atom.res.trim().to_string())
                                .collect()
                        })
                        .unwrap_or_default();

                    if let Some(gro) = &mut self.data.gro_file {
                        for (atom, resname) in gro.atoms_mut().zip(top_resnames.iter()) {
                            atom.set_res_name(resname);
                        }
                    }
                    self.sync_viewer_resnames_from_loaded_files();
                    self.set_loaded_summary(format!("TOP+GRO: {}", file_name));
                    self.set_status("Loaded TOP and synchronized residue names to current GRO/view");
                } else if self.data.gro_file.is_some() {
                    self.set_status(
                        "Loaded TOP, but atom order/count does not match loaded GRO (view not updated)",
                    );
                } else {
                    self.set_status("Loaded TOP (view unchanged: load GRO to visualize residue changes)");
                }
            }
            None => self.set_status("Failed to load TOP file"),
        }
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
        let dropped_paths: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .filter_map(|file| file.path.clone())
                .collect()
        });

        if dropped_paths.is_empty() {
            return;
        }

        let mut top_path: Option<PathBuf> = None;
        let mut gro_path: Option<PathBuf> = None;

        for path in &dropped_paths {
            if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
                match ext.to_ascii_lowercase().as_str() {
                    "top" => top_path = Some(path.clone()),
                    "gro" => gro_path = Some(path.clone()),
                    _ => {}
                }
            }
        }

        if let (Some(top), Some(gro)) = (top_path, gro_path) {
            self.load_top_and_gro_for_resname_sync(top, gro);
            return;
        }

        if let Some(path) = dropped_paths.into_iter().next() {
            self.load_file(path);
        }
    }

    pub fn render_ui(&mut self, ctx: &egui::Context) {
        app_ui::render_edit_dialog(self, ctx);
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

        if self.selection.with_hbond_chk {
            for idx in atoms_on_path {
                if self.selection.selected_atom_indices.contains(&idx) {
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
        let new_name = self.ui.new_res_name.trim().to_uppercase();
        if new_name.len() > 3 {
            self.set_status("Residue name too long");
            return;
        }
        let new_name = format!("{:>3}", new_name); // Pad to 3 chars
        let indices_to_update = self.selection.selected_atom_indices.clone();

        if indices_to_update.is_empty() {
            self.set_status("No atoms selected");
            return;
        }

        if let Some(pdb) = &mut self.data.pdb_file {
            // Update PDB atoms
            let mut atoms_vec: Vec<&mut AtomRecord> = pdb.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.res_name = new_name.clone();
                }
            }
        }

        if let Some(gro) = &mut self.data.gro_file {
            // In GRO, the 2nd field is residue name (resname).
            let mut atoms_vec: Vec<_> = gro.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.set_res_name(&new_name);
                }
            }
        }

        if let Some(top) = &mut self.data.top_file {
            // Keep TOP in sync with current selection indices when TOP+GRO are loaded together.
            let mut atoms_vec: Vec<_> = top.atoms_mut().collect();
            for &idx in &indices_to_update {
                if let Some(atom) = atoms_vec.get_mut(idx) {
                    atom.set_res_name(&new_name);
                }
            }
        }

        // Keep renderer metadata aligned with currently loaded structural data.
        self.sync_viewer_resnames_from_loaded_files();
        self.mark_modified();

        // Clear selection
        self.selection.selected_atom_indices.clear();
        self.sync_selection_to_renderer();
        self.set_status("Residue names updated");
    }

    fn export_structure(&mut self) {
        if self.data.top_file.is_none() && self.data.gro_file.is_none() && self.data.pdb_file.is_none() {
            if let Some(mol) = &self.viewer.molecule {
                self.data.pdb_file = Some(PdbFile::from_molecule(mol));
            }
        }

        if let Some(path) = FileDialog::new().save_file() {
            let saved = if let Some(top) = &self.data.top_file {
                let content = top.dump();
                std::fs::write(&path, content).is_ok()
            } else if let Some(gro) = &self.data.gro_file {
                let content = gro.dump();
                std::fs::write(&path, content).is_ok()
            } else if let Some(pdb) = &mut self.data.pdb_file {
                let content = pdb.dump();
                std::fs::write(&path, content).is_ok()
            } else {
                false
            };

            if saved {
                self.mark_clean();
                self.set_status("Exported structure");
            } else {
                self.set_status("Failed to export structure");
            }
        }
    }

    fn pick_atom_index(&self, ctx: &egui::Context, response: &egui::Response) -> Option<usize> {
        let pointer = ctx.input(|i| i.pointer.hover_pos())?;
        if !response.rect.contains(pointer) {
            return None;
        }
        let local = pointer - response.rect.min;
        let (ray_origin, ray_dir) = self.controller.camera.ray_from_screen(
            local.x,
            local.y,
            response.rect.width().max(1.0),
            response.rect.height().max(1.0),
        );

        let event = self.viewer.pick(ray_origin, ray_dir)?;
        if let moleucle_3dview_rs::viewer::ViewerEvent::AtomClicked(i) = event {
            Some(i)
        } else {
            None
        }
    }
}

impl eframe::App for KuromameApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_keyboard_shortcuts(ctx);
        self.handle_dropped_files(ctx);

        app_ui::render_top_info_panel(self, ctx);
        app_ui::render_edit_dialog(self, ctx);

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
                &self.selection.selected_atom_indices,
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
            let mut hovered_atom_info: Option<String> = None;

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
                if let Some(i) = self.pick_atom_index(ctx, &response) {
                    hovered_atom_info = self.hovered_atom_info(i);
                }
            }

            // Handle mouse click for atom picking
            if response.clicked_by(PointerButton::Primary) {
                if let Some(i) = self.pick_atom_index(ctx, &response) {
                    let selected_now = self.toggle_selected_atom(i);
                    if self.selection.with_hbond_chk {
                        if selected_now {
                            self.add_connected_hydrogens(i);
                        } else {
                            self.remove_connected_hydrogens(i);
                        }
                    }
                    self.sync_selection_to_renderer();
                }
            }

            if let Some(label) = hovered_atom_info {
                response.clone().on_hover_text(label.clone());
                self.ui.hovered_atom_info = label;
            } else {
                self.ui.hovered_atom_info = "Hover an atom for details".to_string();
            }
        });

        app_ui::render_bottom_status_bar(self, ctx);

        // Request continuous repaint
        ctx.request_repaint();
    }

    fn on_exit(&mut self, _ctx: Option<&eframe::glow::Context>) {
        if let Some(render_state) = &self.render_state {
            self.offscreen.free_egui_texture(render_state);
        }
    }
}
