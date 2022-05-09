use cgmath::Vector3;
use criterion::{criterion_group, criterion_main, Criterion};
use mini_raytracer::components::{Camera, Light, Rgba, Sphere};
use mini_raytracer::World;
use std::time::Duration;

const WIDTH: i32 = 500;
const HEIGHT: i32 = 500;
const REFLECTION_PASSES: u32 = 3;

pub fn world_draw(c: &mut Criterion) {
    let mut group = c.benchmark_group("world.draw");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    group.bench_function("world.draw (sequential)", |b| {
        let world = generate_world();
        let mut frame = generate_frame();

        b.iter(|| world.draw(frame.as_mut_slice()))
    });
    group.bench_function("world.draw (parallel)", |b| {
        let world = generate_world();
        let mut frame = generate_frame();

        b.iter(|| world.draw_parallel(frame.as_mut_slice()))
    });
    group.finish();
}

criterion_group!(benches, world_draw);
criterion_main!(benches);

pub fn generate_frame() -> Vec<u8> {
    let n = (WIDTH * HEIGHT * 4) as usize;
    vec![0; n]
}

pub fn generate_world() -> World {
    let spheres = vec![
        Sphere::new(
            Vector3::new(0.0, 1.0, 4.0),
            Rgba::new(255.0, 255.0, 255.0, 255.0),
            0.7,
            // very shiny
            1000.0,
            // 70% reflective
            0.7,
        ),
        Sphere::new(
            Vector3::new(0.0, -1.0, 3.5),
            Rgba::new(255.0, 0.0, 0.0, 255.0),
            0.95,
            // shiny
            500.0,
            // 20% reflective
            0.2,
        ),
        Sphere::new(
            Vector3::new(2.0, 0.0, 4.5),
            Rgba::new(0.0, 0.0, 255.0, 255.0),
            0.95,
            // shiny
            500.0,
            // 30% reflective
            0.3,
        ),
        Sphere::new(
            Vector3::new(-2.0, 0.0, 4.5),
            Rgba::new(0.0, 255.0, 0.0, 255.0),
            0.95,
            // somewhat shiny
            10.0,
            //40% reflective
            0.4,
        ),
        Sphere::new(
            Vector3::new(0.0, -5001.2, 0.0),
            Rgba::new(100.0, 100.0, 100.0, 255.0),
            5000.0,
            // very shiny
            1000.0,
            // 50% reflective
            0.5,
        ),
    ];

    let lights = vec![
        Light::Ambient(0.2),
        Light::Point(0.6, Vector3::new(2.0, 1.0, 0.0)),
        Light::Directional(0.2, Vector3::new(1.0, 4.0, 4.0)),
    ];

    let background_color = Rgba::new(0.0, 0.0, 0.0, 255.0);

    let camera = Camera::new(Vector3::new(0.0, 0.0, 0.0));

    World::new(
        WIDTH.into(),
        HEIGHT.into(),
        REFLECTION_PASSES,
        camera,
        spheres,
        lights,
        background_color,
    )
}
