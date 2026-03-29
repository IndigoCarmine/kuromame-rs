use kuromame_rs::app::KuromameApp;

fn main() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1280.0, 820.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Kuromame - Molecule Editor",
        options,
        Box::new(|cc| Ok(Box::new(KuromameApp::new(cc)))),
    )
}
