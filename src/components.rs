use cgmath::{InnerSpace, Vector2, Vector3};

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vector3<f64>,
    pub color: Rgba,
    pub radius: f64,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, color: Rgba, radius: f64) -> Self {
        Self {
            center,
            color,
            radius,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Rgba {
    pub r: f64,
    pub g: f64,
    pub b: f64,
    pub a: f64,
}

impl Rgba {
    pub fn new(r: f64, g: f64, b: f64, a: f64) -> Self {
        Self { r, g, b, a }
    }

    pub fn to_u8_array(&self) -> [u8; 4] {
        [self.r as u8, self.g as u8, self.b as u8, self.a as u8]
    }
}

pub struct Viewport {
    pub width: f64,
    pub height: f64,
    pub depth: f64,
}

impl Viewport {
    pub fn new(width: f64, height: f64, depth: f64) -> Self {
        Self {
            width,
            height,
            depth,
        }
    }
}

pub struct Canvas {
    pub width: f64,
    pub height: f64,
}

impl Canvas {
    pub fn new(width: f64, height: f64) -> Self {
        Self { width, height }
    }
}
