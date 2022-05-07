use cgmath::{num_traits::Pow, InnerSpace, Vector2, Vector3};
use components::{Canvas, Light, Rgba, Sphere, Viewport};

pub mod components;

/// Returns a color if the ray hits a sphere.
/// Returns None if the ray does not hit
pub fn trace_ray(
    camera: Vector3<f64>,
    direction: Vector3<f64>,
    spheres: &Vec<Sphere>,
    lights: &Vec<Light>,
    t_min: f64,
    t_max: f64,
) -> Option<Rgba> {
    let (closest_sphere, closest_t) =
        closest_sphere_intersection(camera, direction, spheres, t_min, t_max);

    let surface_point = camera + (closest_t * direction);

    closest_sphere.map(|sphere| {
        sphere.color.multiply(compute_lighting_intensity(
            surface_point,
            sphere,
            camera,
            lights,
        ))
    })
}

/// Returns the sphere, and the value t, if an intersection was found.
/// t will == f64::MAX if no intersection is found.
pub fn closest_sphere_intersection(
    origin: Vector3<f64>,
    direction: Vector3<f64>,
    spheres: &Vec<Sphere>,
    t_min: f64,
    t_max: f64,
) -> (Option<Sphere>, f64) {
    let mut closest_t = f64::MAX;
    let mut closest_sphere = None;

    for sphere in spheres {
        let (t1, t2) = intersect_ray_sphere(origin, direction, *sphere);
        if (t1 >= t_min) && (t1 <= t_max) && (t1 < closest_t) {
            closest_t = t1;
            closest_sphere = Some(*sphere);
        }
        if (t2 >= t_min) && (t2 <= t_max) && (t2 < closest_t) {
            closest_t = t2;
            closest_sphere = Some(*sphere);
        }
    }

    (closest_sphere, closest_t)
}

/// Loops through all the spheres and finds any intersections between the ray and points on the spheres.
/// It will return two values, t1 and t2, which are the scalars needed to show the intersection.
///
/// `P = O + t1(D)`
///
/// ` P = O + t2(D)`
///
/// Origin will usually be the camera, but can also be a different value when calulating for any blocking spheres between two
/// arbitrary points. It is not referring to the 3d space origin (0,0,0).
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

// Computes the lighting intensity using lighting diffusion and specular illumination
pub fn compute_lighting_intensity(
    surface_point: Vector3<f64>,
    sphere: Sphere,
    camera: Vector3<f64>,
    lights: &Vec<Light>,
) -> f64 {
    let mut i = 0.0;

    let surface_normal = (surface_point - sphere.center).normalize();

    for light in lights {
        if let Light::Ambient(_i) = light {
            i += _i;
            continue;
        };

        // direction is the vector starting from the surface point, and pointing to the light
        let (intensity, direction) = match *light {
            Light::Point(intensity, position) => (intensity, position - surface_point),
            Light::Directional(intensity, direction) => (intensity, direction),
            _ => unreachable!(),
        };

        let normal_dot_direction = surface_normal.dot(direction);

        // if the direction is below 0 then that means we're lighting the back of the surface!

        // we're done calculating diffuse illumination
        if normal_dot_direction > 0.0 {
            i += intensity
                * (normal_dot_direction / (surface_normal.magnitude() * direction.magnitude()));
        };

        //TODO: simplify this
        // r is going to be the direction of the light (starting from the surface point) reflected over the normal of the surface
        // we're going to break the direction into two components, dn and dp, where dn is parallel to the normal
        // direction = dn + dp
        // this means that r = dn - dp
        // dn is the vector projection of the direction onto the surface normal.
        // normally we would divide the dot product by |n|, but |n| = 1 so no division is neccessary

        // broken down code example for calculating r
        // let dn = surface_normal * (surface_normal.dot(direction));
        // let dp = direction - dn;
        // let r = dn - dp;

        // compact version
        let r = (2.0 * surface_normal * normal_dot_direction) - direction;

        // v is the "view vector" that points from the surface point to the camera
        let v = camera - surface_point;
        let specular_multiplier = (r.dot(v) / (r.magnitude() * v.magnitude())).pow(sphere.specular);

        i += intensity * (specular_multiplier);
    }

    i
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
