use cgmath::{InnerSpace, Vector3};
use components::{Rgba, Sphere};

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
