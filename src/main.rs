#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use eframe::egui;
use egui_wgpu::{WgpuSetup, WgpuSetupCreateNew, wgpu};
use kuromame_rs::app::KuromameApp;

fn main() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default().with_inner_size([1280.0, 820.0]),
        wgpu_options: egui_wgpu::WgpuConfiguration {
            wgpu_setup: WgpuSetup::CreateNew(WgpuSetupCreateNew {
                // Prefer Vulkan and fall back to DX12 when Vulkan runtime/driver is unavailable.
                instance_descriptor: wgpu::InstanceDescriptor {
                    backends: wgpu::Backends::VULKAN | wgpu::Backends::DX12,
                    ..Default::default()
                },
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            }),
            ..Default::default()
        },
        ..Default::default()
    };

    eframe::run_native(
        "Kuromame - Molecule Editor",
        options,
        Box::new(|cc| Ok(Box::new(KuromameApp::new(cc)))),
    )
}
