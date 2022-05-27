use cgmath::Vector3;

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vector3<f64>,
    pub color: Rgba,
    pub radius: f64,
    pub specular: f64,
    pub reflective: f64,
}

impl Sphere {
    pub fn new(
        center: Vector3<f64>,
        color: Rgba,
        radius: f64,
        specular: f64,
        reflective: f64,
    ) -> Self {
        Self {
            center,
            color,
            radius,
            specular,
            reflective,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Camera {
    pub position: Vector3<f64>,
}

impl Camera {
    pub fn new(position: Vector3<f64>) -> Self {
        Self { position }
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

    /// Finds the difference of two [`Rgba`]s.
    /// The result is normalized and will be a value between 0.0 and 1.0.
    pub fn difference(&self, other: Self) -> f64 {
        let denom = (255 * 3) as f64;
        let mut raw_difference = 0.0;
        raw_difference += (self.r - other.r).abs();
        raw_difference += (self.g - other.g).abs();
        raw_difference += (self.b - other.b).abs();
        raw_difference += (self.a - other.a).abs();
        raw_difference / denom
    }
}

impl std::ops::Add for Rgba {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            r: self.r + other.r,
            g: self.g + other.g,
            b: self.b + other.b,
            a: self.a + other.a,
        }
    }
}

#[derive(Clone, Debug, Copy)]
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

#[derive(Clone, Debug, Copy)]
pub struct Canvas {
    pub width: f64,
    pub height: f64,
}

impl Canvas {
    pub fn new(width: f64, height: f64) -> Self {
        Self { width, height }
    }
}

#[derive(Clone, Debug, Copy)]
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
