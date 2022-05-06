use cgmath::{InnerSpace, Vector2, Vector3};
use components::{Canvas, Rgba, Sphere, Viewport};

pub mod components;

/// Returns a color if the ray hits a sphere.
/// Returns None if the ray does not hit
pub fn trace_ray(
    origin: Vector3<f64>,
    direction: Vector3<f64>,
    spheres: &Vec<Sphere>,
    t_min: f64,
    t_max: f64,
) -> Option<Rgba> {
    let mut closest_t = f64::MAX;
    let mut closest_sphere = None;

    for sphere in spheres {
        let (t1, t2) = intersect_ray_sphere(origin, direction, *sphere);
        if (t1 >= t_min) && (t1 <= t_max) && (t1 < closest_t) {
            closest_t = t1;
            closest_sphere = Some(sphere);
        }
        if (t2 >= t_min) && (t2 <= t_max) && (t2 < closest_t) {
            closest_t = t2;
            closest_sphere = Some(sphere);
        }
    }

    match closest_sphere {
        Some(sphere) => Some(sphere.color),
        None => None,
    }
}

/// Loops through all the spheres and finds any intersections between the ray and points on the spheres.
/// It will return two values, t1 and t2, which are the scalars needed to show the intersection.
///
/// `P = O + t1(D)`
///
/// ` P = O + t2(D)`
pub fn intersect_ray_sphere(
    origin: Vector3<f64>,
    direction: Vector3<f64>,
    sphere: Sphere,
) -> (f64, f64) {
    let r = sphere.radius;
    let co = origin - sphere.center;
    let a = direction.dot(direction);
    let b = 2.0 * co.dot(direction);
    let c = co.dot(co) - r * r;

    // quadratic formula!
    let discriminant = (b * b) - (4.0 * a * c);
    if discriminant < 0.0 {
        return (f64::MAX, f64::MAX);
    }

    let t1 = (-b + discriminant.sqrt()) / (2.0 * a);
    let t2 = (-b - discriminant.sqrt()) / (2.0 * a);
    (t1, t2)
}

pub fn compute_lighting_intensity(
    surface_point: Vector3<f64>,
    surface_normal: Vector3<f64>,
    ambient_light_intensity: f64,
    directional_lights: Vec<DirectionalLight>,
    point_lights: Vec<PointLight>,
) {
}

pub fn frame_to_canvas_coords(frame_coords: Vector2<f64>, canvas: &Canvas) -> Vector2<f64> {
    let cx = frame_coords.x - (canvas.width as f64 / 2.0);
    let cy = (canvas.height as f64 / 2.0) - frame_coords.y;
    Vector2::new(cx, cy)
}

pub fn canvas_to_viewport_coords(
    canvas_coords: Vector2<f64>,
    canvas: &Canvas,
    viewport: &Viewport,
) -> Vector3<f64> {
    let vx = canvas_coords.x * (viewport.width / canvas.width);
    let vy = canvas_coords.y * (viewport.height / canvas.height);
    Vector3::new(vx, vy, viewport.depth)
}

pub fn frame_to_viewport_coords(
    frame_coords: Vector2<f64>,
    canvas: &Canvas,
    viewport: &Viewport,
) -> Vector3<f64> {
    let canvas_coords = frame_to_canvas_coords(frame_coords, canvas);
    canvas_to_viewport_coords(canvas_coords, canvas, viewport)
}

#[cfg(test)]
mod tests {
    use crate::components::{Canvas, Viewport};
    use crate::{canvas_to_viewport_coords, frame_to_canvas_coords, frame_to_viewport_coords};
    use cgmath::{Vector2, Vector3};

    #[test]
    fn frame_to_canvas() {
        let frame_coords = Vector2::new(50.0, 50.0);
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_to_canvas_coords(frame_coords, &canvas);
        assert_eq!(canvas_coords, Vector2::new(-100.0, 100.0));

        let frame_coords = Vector2::new(150.0, 150.0);
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_to_canvas_coords(frame_coords, &canvas);
        assert_eq!(canvas_coords, Vector2::new(0.0, 0.0));

        // 299 is the highest x and y value in frame coords with a width and height of 300 because the index starts at 0
        let frame_coords = Vector2::new(299.0, 299.0);
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_to_canvas_coords(frame_coords, &canvas);
        assert_eq!(canvas_coords, Vector2::new(149.0, -149.0));

        let frame_coords = Vector2::new(0.0, 299.0);
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_to_canvas_coords(frame_coords, &canvas);
        assert_eq!(canvas_coords, Vector2::new(-150.0, -149.0));

        let frame_coords = Vector2::new(299.0, 0.0);
        let canvas = Canvas::new(300.0, 300.0);
        let canvas_coords = frame_to_canvas_coords(frame_coords, &canvas);
        assert_eq!(canvas_coords, Vector2::new(149.0, 150.0));
    }

    #[test]
    fn canvas_to_viewport() {
        let canvas = Canvas::new(300.0, 300.0);
        let viewport = Viewport::new(300.0, 300.0, 1.0);
        let canvas_coords = Vector2::new(100.0, 100.0);
        let viewport_coords = canvas_to_viewport_coords(canvas_coords, &canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(100.0, 100.0, 1.0));

        // "squeezes" the viewport by a factor of 2
        let canvas = Canvas::new(600.0, 600.0);
        let viewport = Viewport::new(300.0, 300.0, 1.0);
        let canvas_coords = Vector2::new(100.0, 100.0);
        let viewport_coords = canvas_to_viewport_coords(canvas_coords, &canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(50.0, 50.0, 1.0));

        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let canvas_coords = Vector2::new(40.0, 40.0);
        let viewport_coords = canvas_to_viewport_coords(canvas_coords, &canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(0.4, 0.4, 1.0));

        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let canvas_coords = Vector2::new(-40.0, -40.0);
        let viewport_coords = canvas_to_viewport_coords(canvas_coords, &canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(-0.4, -0.4, 1.0));
    }

    #[test]
    fn frame_to_viewport() {
        let canvas = Canvas::new(100.0, 100.0);
        let viewport = Viewport::new(1.0, 1.0, 1.0);
        let frame_coords = Vector2::new(90.0, 10.0);
        let viewport_coords = frame_to_viewport_coords(frame_coords, &canvas, &viewport);
        assert_eq!(viewport_coords, Vector3::new(0.4, 0.4, 1.0));
    }
}
