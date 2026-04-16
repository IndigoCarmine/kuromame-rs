use eframe::egui;
use moleucle_3dview_rs::RenderStyle;

use super::{AtomRecord, KuromameApp};

pub fn render_edit_dialog(app: &mut KuromameApp, ctx: &egui::Context) {
    let mut open_edit_dialog = app.ui.show_edit_dialog;

    if open_edit_dialog {
        egui::Window::new("Edit Residue Name")
            .open(&mut open_edit_dialog)
            .show(ctx, |ui| {
                ui.label("Enter new residue name (3 letters):");
                ui.text_edit_singleline(&mut app.ui.new_res_name);
                if ui.button("Apply").clicked() {
                    app.apply_res_name_change();
                    app.ui.show_edit_dialog = false;
                }
            });
        app.ui.show_edit_dialog = open_edit_dialog;
    }
}

pub fn render_top_help_panel(app: &mut KuromameApp, ctx: &egui::Context) {
    egui::TopBottomPanel::top("help")
        .resizable(false)
        .show(ctx, |ui| {
            ui.label("LMB: pick atom  RMB drag: orbit  MMB/Shift+RMB drag: pan  Wheel: dolly");
            ui.label("Drop a PDB/MOL2/GRO file anywhere in the window to load it");
            ui.label(format!(
                "Selected atoms: {:?}",
                app.selection.selected_atom_indices
            ));

            ui.horizontal(|ui| {
                ui.label("Render Style:");
                let mut style = app.offscreen.render_style();
                ui.selectable_value(&mut style, RenderStyle::BallStick, "BallStick");
                ui.selectable_value(&mut style, RenderStyle::Wireframe, "Wireframe");
                app.offscreen.set_render_style(style);
            });
        });
}

pub fn render_side_panel(app: &mut KuromameApp, ctx: &egui::Context) {
    egui::SidePanel::right("control_panel")
        .resizable(true)
        .default_width(300.0)
        .show(ctx, |ui| {
            ui.heading("Controls");
            if ui.button("Open Molecule File (PDB/MOL2/GRO)").clicked() {
                app.open_file();
            }
            if ui.button("Open TOP + GRO (sync TOP resname from GRO)").clicked() {
                app.open_top_and_gro_for_resname_sync();
            }

            if let Some(path) = &app.data.current_file_path {
                ui.label(format!("Loaded: {:?}", path.file_name().unwrap()));
            } else {
                ui.label("No file loaded");
            }

            ui.checkbox(
                &mut app.selection.with_hbond_chk,
                "Select with connected hydrogens",
            );

            ui.separator();
            ui.label("Selected atoms:");

            let selected_indices = app.selection.selected_atom_indices.clone();

            ui.label(format!("Count: {}", selected_indices.len()));

            egui::ScrollArea::vertical()
                .max_height(400.0)
                .show(ui, |ui| {
                    if let Some(pdb) = &app.data.pdb_file {
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
                    } else if let Some(gro) = &app.data.gro_file {
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
                app.selection.selected_atom_indices.clear();
                app.sync_selection_to_renderer();
            }

            let can_select_path = selected_indices.len() == 2;
            if ui
                .add_enabled(
                    can_select_path,
                    egui::Button::new("Select atoms between (2 atoms)"),
                )
                .clicked()
            {
                app.select_shortest_path(selected_indices[0], selected_indices[1]);
            }

            if ui
                .add_enabled(
                    !selected_indices.is_empty(),
                    egui::Button::new("Change selected residues' resname..."),
                )
                .clicked()
            {
                app.ui.show_edit_dialog = true;
                app.ui.new_res_name = "ALA".to_string();
            }

            if ui.button("Export (Save)").clicked() {
                app.export_structure();
            }

            ui.label(&app.ui.status_msg);
        });
}
