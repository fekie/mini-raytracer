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

#[derive(Debug, Clone, Copy)]
pub struct CanvasCoords {
    pub coords: Vector2<f64>,
}

impl CanvasCoords {
    pub fn new(coords: Vector2<f64>) -> Self {
        Self { coords }
    }

    pub fn to_viewport_coords(&self, canvas: &Canvas, viewport: &Viewport) -> Vector3<f64> {
        let vx = self.coords.x * (viewport.width / canvas.width);
        let vy = self.coords.y * (viewport.height / canvas.height);
        Vector3::new(vx, vy, viewport.depth)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct FrameCoords {
    pub coords: Vector2<i32>,
}

impl FrameCoords {
    pub fn new(coords: Vector2<i32>) -> Self {
        Self { coords }
    }

    pub fn to_canvas_coords(&self, canvas: &Canvas) -> CanvasCoords {
        let cx = self.coords.x as f64 - (canvas.width as f64 / 2.0);
        let cy = (canvas.height as f64 / 2.0) - self.coords.y as f64;
        CanvasCoords::new(Vector2::new(cx, cy))
    }

    pub fn to_viewport_coords(&self, canvas: &Canvas, viewport: &Viewport) -> Vector3<f64> {
        let canvas_coords = self.to_canvas_coords(canvas);
        canvas_coords.to_viewport_coords(canvas, viewport)
    }
}

#[cfg(test)]
mod tests {
    use crate::components::{Canvas, CanvasCoords, FrameCoords, Viewport};
    use cgmath::{Vector2, Vector3};

    #[test]
    fn frame_to_canvas() {
        let frame_coords = FrameCoords::new(Vector2::new(50, 50));
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        assert_eq!(canvas_coords.coords, Vector2::new(-100.0, 100.0));

        let frame_coords = FrameCoords::new(Vector2::new(150, 150));
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        assert_eq!(canvas_coords.coords, Vector2::new(0.0, 0.0));

        // 299 is the highest x and y value in frame coords with a width and height of 300 because the index starts at 0
        let frame_coords = FrameCoords::new(Vector2::new(299, 299));
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        assert_eq!(canvas_coords.coords, Vector2::new(149.0, -149.0));

        let frame_coords = FrameCoords::new(Vector2::new(0, 299));
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        assert_eq!(canvas_coords.coords, Vector2::new(-150.0, -149.0));

        let frame_coords = FrameCoords::new(Vector2::new(299, 0));
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        assert_eq!(canvas_coords.coords, Vector2::new(149.0, 150.0));
    }

    #[test]
    fn canvas_to_viewport() {
        let canvas = Canvas::new(300.0, 300.0);
        let viewport = Viewport::new(300.0, 300.0, 1.0);
        let canvas_coords = CanvasCoords::new(Vector2::new(100.0, 100.0));
        let viewport_coords = canvas_coords.to_viewport_coords(&canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(100.0, 100.0, 1.0));

        // "squeezes" the viewport by a factor of 2
        let canvas = Canvas::new(600.0, 600.0);
        let viewport = Viewport::new(300.0, 300.0, 1.0);
        let canvas_coords = CanvasCoords::new(Vector2::new(100.0, 100.0));
        let viewport_coords = canvas_coords.to_viewport_coords(&canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(50.0, 50.0, 1.0));

        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let canvas_coords = CanvasCoords::new(Vector2::new(40.0, 40.0));
        let viewport_coords = canvas_coords.to_viewport_coords(&canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(0.4, 0.4, 1.0));

        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let canvas_coords = CanvasCoords::new(Vector2::new(-40.0, -40.0));
        let viewport_coords = canvas_coords.to_viewport_coords(&canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(-0.4, -0.4, 1.0));
    }

    #[test]
    fn frame_to_viewport() {
        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let frame_coords = FrameCoords::new(Vector2::new(90, 10));
        let canvas_coords = frame_coords.to_canvas_coords(&canvas);
        let viewport_coords = canvas_coords.to_viewport_coords(&canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(0.4, 0.4, 1.0));
    }
}
