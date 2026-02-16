use graphics::{EngineUpdates, EntityUpdate, GraphicsSettings, Scene, UiSettings, run};
use kuromame_rs::app::KuromameApp;
use moleucle_3dview_rs::viewer::ViewerEvent;

fn main() {
    let mut app = KuromameApp::new();

    // Sync initial camera state
    app.controller.camera.position = nalgebra::Point3::new(0.0, 0.0, -20.0);

    let mut scene = Scene::default();
    app.update_scene(&mut scene);

    run(
        app,
        scene,
        UiSettings::default(),
        GraphicsSettings::default(),
        // Render Handler
        |app, scene, _dt| {
            let mut updates = EngineUpdates::default();
            if app.viewer.dirty {
                app.viewer.update_scene(scene);
                updates.meshes = true;
                updates.entities = EntityUpdate::All;
            }
            app.controller.update_scene_camera(scene);
            updates.camera = true;
            updates
        },
        // Device Event Handler
        |_state, _event, _scene, _is_synthetic, _dt| EngineUpdates::default(),
        // Window Event Handler
        |app, event, scene, _dt| {
            let (picked, updates) = app.handle_event(&event, scene);

            if let Some(event) = picked {
                match event {
                    ViewerEvent::AtomClicked(i) => {
                        println!("Atom {} Clicked", i);
                        if let Some(renderer) = &mut app.viewer.additional_render {
                            renderer.add_atom(i);
                            app.viewer.dirty = true;
                        }
                    }
                    _ => {}
                }
            }
            updates
        },
        // GUI Handler
        |app, ctx, _scene| {
            app.render_ui(ctx);
            EngineUpdates::default()
        },
    );
}
