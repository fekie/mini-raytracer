use cgmath::{Vector2, Vector3};

pub struct Sphere {
    center: Vector3<f64>,
    radius: f64,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, radius: f64) -> Self {
        Self { center, radius }
    }
}

pub struct Viewport {
    width: u32,
    height: u32,
}

impl Viewport {
    pub fn new(width: u32, height: u32) -> Self {
        Self { width, height }
    }
}

pub struct Canvas {
    width: u32,
    height: u32,
}

impl Canvas {
    pub fn new(width: u32, height: u32) -> Self {
        Self { width, height }
    }
}
