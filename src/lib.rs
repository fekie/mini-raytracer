use cgmath::{num_traits::Pow, InnerSpace, Vector2, Vector3};
use components::{Camera, Canvas, Light, Rgba, Sphere, Viewport};
use rayon::prelude::*;
use std::collections::HashMap;

pub mod components;

/// [`BigChunk`] carries the four corners of an 8x8 chunk of the screen.
/// If the four corners are different by a threshold of h, subdivide the chunk into 4x4 chunks.
/// When subdivided, this chunk becomes a [`MediumChunk`].
#[derive(Clone, Debug, Copy)]
pub struct BigChunk {
    top_left: Vector2<usize>,
    top_left_color: Option<Rgba>,
    top_right: Vector2<usize>,
    top_right_color: Option<Rgba>,
    bottom_left: Vector2<usize>,
    bottom_left_color: Option<Rgba>,
    bottom_right: Vector2<usize>,
    bottom_right_color: Option<Rgba>,
}

impl BigChunk {
    /// Only needs the x and y values of the top left corner of the chunk.
    pub fn new(x: usize, y: usize) -> Self {
        Self {
            top_left: Vector2::new(x, y),
            top_left_color: None,
            top_right: Vector2::new(x + 7, y),
            top_right_color: None,
            bottom_left: Vector2::new(x, y + 7),
            bottom_left_color: None,
            bottom_right: Vector2::new(x + 7, y + 7),
            bottom_right_color: None,
        }
    }

    /// Chunks will go from left to right, top to bottom.
    pub fn from_frame_dimensions(width: usize, height: usize) -> Vec<BigChunk> {
        assert_eq!(width, height);
        assert_eq!(width % 8, 0);

        let chunks_per_row = width / 8;
        let big_chunk_amount = chunks_per_row * chunks_per_row;
        let mut big_chunks = Vec::with_capacity(big_chunk_amount);

        for y in (0..width).step_by(8) {
            for x in (0..width).step_by(8) {
                let big_chunk = Self::new(x, y);
                big_chunks.push(big_chunk);
            }
        }

        big_chunks
    }

    /// Index 0 is top left, index 1 is top right
    /// index 2 is bottom left, index 3 is bottom right.
    pub fn subdivide(&self) -> Vec<MediumChunk> {
        todo!()
    }

    // TODO: make this only enabled by a feature or something.
    pub fn debug(big_chunks: &Vec<BigChunk>, canvas_width: usize, frame: &mut [u8]) {
        for big_chunk in big_chunks {
            let index = ((big_chunk.top_left.y * canvas_width) + big_chunk.top_left.x) * 4;
            frame[index] = 255;
            frame[index + 1] = 255;
            frame[index + 2] = 255;
        }
    }
}

/// [`MediumChunk`] carries the four corners of a 4x4 chunk of the screen.
/// Subdivides to [`SmallChunk`] (2x2)
#[derive(Clone, Debug, Copy)]
pub struct MediumChunk {
    top_left: Vector2<u32>,
    top_left_color: Option<Rgba>,
    top_right: Vector2<u32>,
    top_right_color: Option<Rgba>,
    bottom_left: Vector2<u32>,
    bottom_left_color: Option<Rgba>,
    bottom_right: Vector2<u32>,
    bottom_right_color: Option<Rgba>,
}

impl MediumChunk {
    /// Index 0 is top left, index 1 is top right
    /// index 2 is bottom left, index 3 is bottom right.
    pub fn subdivide(&self) -> Vec<SmallChunk> {
        todo!()
    }
}

/// [`SmallChunk`] carries the four corners of a 2x2 chunk of the screen.
/// Cannot be subdivided further.
#[derive(Clone, Debug, Copy)]
pub struct SmallChunk {
    top_left: Vector2<u32>,
    top_left_color: Option<Rgba>,
    top_right: Vector2<u32>,
    top_right_color: Option<Rgba>,
    bottom_left: Vector2<u32>,
    bottom_left_color: Option<Rgba>,
    bottom_right: Vector2<u32>,
    bottom_right_color: Option<Rgba>,
}

#[derive(Clone, Debug, Copy)]
/// [`ThickRow`] starts at start_y and ends at end_y. It includes all of the x values. End_y is non-inclusive.
pub struct ThickRow {
    start_y: u32,
    end_y: u32,
}

impl ThickRow {
    pub fn new(start_y: u32, end_y: u32) -> Self {
        Self { start_y, end_y }
    }

    pub fn split_frame(row_size: u32, height: u32) -> Vec<Self> {
        assert_eq!(height % row_size, 0);
        let mut frame_chunks = Vec::new();

        for start_y in (0..height).step_by(row_size as usize) {
            let end_y = start_y + row_size;
            frame_chunks.push(Self::new(start_y, end_y));
        }

        frame_chunks
    }
}

pub struct World {
    width: f64,
    #[allow(dead_code)]
    height: f64,
    reflection_passes: u32,
    camera: Camera,
    canvas: Canvas,
    viewport: Viewport,
    pub spheres: Vec<Sphere>,
    pub lights: Vec<Light>,
    background_color: Rgba,
}

impl World {
    pub fn new(
        width: f64,
        height: f64,
        reflection_passes: u32,
        camera: Camera,
        spheres: Vec<Sphere>,
        lights: Vec<Light>,
        background_color: Rgba,
    ) -> Self {
        Self {
            width,
            height,
            reflection_passes,
            camera,
            canvas: Canvas::new(width, height),
            viewport: Viewport::new(1.0, 1.0, 1.0),
            spheres,
            lights,
            background_color,
        }
    }

    #[inline]
    pub fn draw(&self, frame: &mut [u8]) {
        // iterates through each pixel, we will need to do raytracing to find the color of each pixel
        for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
            let x = (i % self.width as usize) as f64;
            let y = (i / self.width as usize) as f64;

            let frame_coords = Vector2::new(x, y);

            // convert the canvas coords to the viewport coords
            let direction =
                Self::frame_to_viewport_coords(frame_coords, &self.canvas, &self.viewport);

            let rgba = Self::trace_ray(
                self.camera.position,
                self.camera,
                direction,
                self.background_color,
                &self.spheres,
                &self.lights,
                1.0,
                f64::MAX,
                self.reflection_passes,
            );

            pixel.copy_from_slice(&rgba.to_u8_array());
        }
    }

    #[inline]
    pub fn draw_parallel(&self, frame: &mut [u8]) {
        // TODO: make column_size and row_size changable
        // splits frame computation in 4 threads
        let thick_rows = ThickRow::split_frame(125, self.height as u32);

        let width = self.width as u32;
        let canvas = self.canvas;
        let viewport = self.viewport;
        let camera = self.camera;
        let background_color = self.background_color;
        let all_spheres = self.spheres.clone();
        let all_lights = self.lights.clone();
        let reflection_passes = self.reflection_passes;

        let mut new_frame = Vec::with_capacity(frame.len());
        let mut meow = thick_rows
            .par_iter()
            .enumerate()
            .map(|(i, thick_row)| {
                let mut pixels = Vec::new();

                for y in (thick_row.start_y)..(thick_row.end_y) {
                    for x in 0..(width) {
                        let frame_coords = Vector2::new(x as f64, y as f64);

                        // convert the canvas coords to the viewport coords
                        let direction =
                            Self::frame_to_viewport_coords(frame_coords, &canvas, &viewport);

                        let rgba = Self::trace_ray(
                            camera.position,
                            camera,
                            direction,
                            background_color,
                            &all_spheres,
                            &all_lights,
                            1.0,
                            f64::MAX,
                            reflection_passes,
                        );

                        pixels.push(rgba.r as u8);
                        pixels.push(rgba.g as u8);
                        pixels.push(rgba.b as u8);
                        pixels.push(rgba.a as u8);
                    }
                }
                (i, pixels)
            })
            .collect::<HashMap<usize, Vec<u8>>>();

        for i in 0..meow.len() {
            new_frame.append(meow.get_mut(&i).unwrap())
        }

        frame.copy_from_slice(&new_frame);
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

    /// Reflects a ray over a normal. The second value MUST be a normal as this uses an optimization of the vector projection
    /// algorithm that only works when the length of the value being projected on has a magnitude of 1.
    pub fn reflect_ray_over_normal(ray: Vector3<f64>, normal: Vector3<f64>) -> Vector3<f64> {
        (2.0 * normal * normal.dot(ray)) - ray
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
        let canvas_coords = Self::frame_to_canvas_coords(frame_coords, canvas);
        Self::canvas_to_viewport_coords(canvas_coords, canvas, viewport)
    }

    #[allow(clippy::too_many_arguments)]
    /// Returns a color if the ray hits a sphere.
    /// Returns background_color if the ray does not hit
    pub fn trace_ray(
        origin: Vector3<f64>,
        camera: Camera,
        direction: Vector3<f64>,
        background_color: Rgba,
        all_spheres: &Vec<Sphere>,
        all_lights: &Vec<Light>,
        t_min: f64,
        t_max: f64,
        recursion_depth: u32,
    ) -> Rgba {
        let (closest_sphere, closest_t) =
            Self::closest_sphere_intersection(origin, direction, all_spheres, t_min, t_max);

        let surface_point = origin + (closest_t * direction);

        if closest_sphere.is_none() {
            return background_color;
        }

        // this will always not panic because we already checked to see if it was None
        let unwrapped_sphere = closest_sphere.unwrap();

        let local_color = unwrapped_sphere
            .color
            .multiply(Self::compute_lighting_intensity(
                surface_point,
                camera,
                unwrapped_sphere,
                all_spheres,
                all_lights,
            ));

        // if we hit the recursion limit, or the object is not reflective, we can go ahead and return the color
        let reflective = unwrapped_sphere.reflective;
        if (recursion_depth == 0) || (reflective <= 0.0) {
            local_color
        } else {
            let surface_normal = (surface_point - unwrapped_sphere.center).normalize();
            // this reflection will point in the direction of where we need to look for intersections
            let reflection = Self::reflect_ray_over_normal(-direction, surface_normal);
            let reflected_color = Self::trace_ray(
                surface_point,
                camera,
                reflection,
                background_color,
                all_spheres,
                all_lights,
                0.00001,
                f64::MAX,
                recursion_depth - 1,
            );

            // we "mix" the reflected color into our local color
            // a reflectiveness of 0.2 means that the new color includes 20% of the relected color
            local_color.multiply(1.0 - reflective) + reflected_color.multiply(reflective)
        }
    }

    // Computes the lighting intensity using lighting diffusion and specular illumination
    pub fn compute_lighting_intensity(
        surface_point: Vector3<f64>,
        camera: Camera,
        sphere: Sphere,
        all_spheres: &Vec<Sphere>,
        all_lights: &Vec<Light>,
    ) -> f64 {
        let mut i = 0.0;

        let surface_normal = (surface_point - sphere.center).normalize();

        for light in all_lights {
            if let Light::Ambient(_i) = light {
                i += _i;
                continue;
            };

            // direction is the vector starting from the surface point, and pointing to the light
            // t_max is the scalar multiplier we will use when checking for shadows
            // point lights have a scalar value of 1, because we dont want to look for shadows if a sphere is behind the light
            // directional lights have a scalar value of f64::MAX
            let (default_intensity, direction, t_max) = match *light {
                Light::Point(intensity, position) => (intensity, position - surface_point, 1.0),
                Light::Directional(intensity, direction) => (intensity, direction, f64::MAX),
                _ => unreachable!(),
            };

            // shadow check
            // we dont want t_min == 0 so that so that the sphere doesnt intersect itself
            let (shadow_sphere, _shadow_t) = Self::closest_sphere_intersection(
                surface_point,
                direction,
                all_spheres,
                0.0000001,
                t_max,
            );

            if shadow_sphere.is_some() {
                continue;
            }

            let normal_dot_direction = surface_normal.dot(direction);

            // if the direction is below 0 then that means we're lighting the back of the surface!
            if normal_dot_direction > 0.0 {
                // calculate and add diffuse illumination
                i += default_intensity
                    * (normal_dot_direction / (surface_normal.magnitude() * direction.magnitude()));
            };

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
            // reflect_ray_over_normal would also work here, but we can reuse normal_dot_direction if we do it this way
            let r = (2.0 * surface_normal * normal_dot_direction) - direction;

            // v is the "view vector" that points from the surface point to the camera
            let v = camera.position - surface_point;
            let specular_multiplier =
                (r.dot(v) / (r.magnitude() * v.magnitude())).pow(sphere.specular);

            i += default_intensity * (specular_multiplier);
        }

        i
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
            let (t1, t2) = Self::intersect_ray_sphere(origin, direction, *sphere);
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
}

#[cfg(test)]
mod tests {
    use crate::BigChunk;

    #[test]
    fn big_chunks_from_frame_dimensions() {
        let width = 600;
        let height = 600;
        let big_chunks = BigChunk::from_frame_dimensions(width, height);
        assert_eq!(big_chunks.len(), 5625)
    }
}
