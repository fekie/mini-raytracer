use cgmath::{InnerSpace, Vector2, Vector3};

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vector3<f64>,
    pub color: Rgba,
    pub radius: f64,
    pub specular: f64,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, color: Rgba, radius: f64, specular: f64) -> Self {
        Self {
            center,
            color,
            radius,
            specular,
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

    /// does not mutate rgba
    pub fn multiply(&self, f: f64) -> Self {
        Rgba::new(self.r * f, self.g * f, self.b * f, self.a)
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

/// Contains all the possible lights:
///
/// Direction (contains direction value)
///
/// Point (contains position and default intensity value)
///
/// Ambient (contains intensity value)
pub enum Light {
    Directional(f64, Vector3<f64>),
    Point(f64, Vector3<f64>),
    Ambient(f64),
}
