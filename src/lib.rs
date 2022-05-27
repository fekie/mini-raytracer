use cgmath::{num_traits::Pow, InnerSpace, Vector2, Vector3};
use components::{Camera, Canvas, Light, Rgba, Sphere, Viewport};
use rayon::prelude::*;
use std::{collections::HashMap, os};

pub mod components;

/// Higher value means that the chunks will subdivide less often,
/// leading to a faster but lower quality rendering.
const ALLOWED_COLOR_DIFFERENCE: f64 = 0.01;

const LARGE_CHUNK_WIDTH: usize = 8;
const MEDIUM_CHUNK_WIDTH: usize = 4;
const SMALL_CHUNK_WIDTH: usize = 2;

/// Carries a mutable slice of the frame, as well as the offset for the top left corner.
#[derive(Debug)]
pub enum Chunk {
    // 8x8
    Large([u8; 256], Vector2<usize>),
    // 4x4
    Medium([u8; 256], Vector2<usize>),
    // 2x2
    // A bit hacky, will need to refactor somehow.
    // Vec<usize> contains the indices that will be subdivided further.
    // It is this way because I forgot that I would not be subdividing every chunk if one needs to be subdivided
    // because that would happen to every chunk that lies partially off of the sphere.
    Small([u8; 256], Vector2<usize>, Vec<usize>),
}

impl Chunk {
    /// Chunks will go from left to right, top to bottom.
    pub fn from_frame_dimensions(width: usize, height: usize) -> Vec<Self> {
        assert_eq!(width, height);
        assert_eq!(width % 8, 0);

        let chunks_per_row = width / LARGE_CHUNK_WIDTH;
        let big_chunk_amount = chunks_per_row * chunks_per_row;
        let mut big_chunks = Vec::with_capacity(big_chunk_amount);

        for y in (0..width).step_by(LARGE_CHUNK_WIDTH) {
            for x in (0..width).step_by(LARGE_CHUNK_WIDTH) {
                let big_chunk = Self::Large([0; 256], Vector2::new(x, y));
                big_chunks.push(big_chunk);
            }
        }

        big_chunks
    }

    // TODO: refactor so I don't need section_indices.
    /// Finds the color of the four corner pixels.
    /// If the color difference between them is too great, it will subdivide into smaller chunks.
    /// On subdivision, it finds the color of the four corner pixels, if needed it will subdivide again.
    /// This function will not subdivide below small chunks.
    pub fn subdivide(&mut self, section_indices: Option<Vec<usize>>) {
        match self {
            Self::Large(chunk_frame, offset) => {
                *self = Self::Medium(*chunk_frame, *offset);
            }
            Self::Medium(chunk_frame, offset) => {
                *self = Self::Small(*chunk_frame, *offset, section_indices.unwrap());
            }
            Self::Small(_, _, _) => {
                panic!("Can't divide a small chunk into smaller chunks!");
            }
        }
    }

    /// Initially takes a [`Chunk::Large`], then recursively calls render() until all pixels have been rendered.
    #[allow(clippy::erasing_op)]
    #[allow(clippy::identity_op)]
    pub fn render(&mut self, world: &World) {
        match self {
            // Assumes that no pixels have been rendered.
            Self::Large(chunk_frame, offset) => {
                // We will find the colors of the four corners.
                // If the color difference between them is too big, we subdivide it and call render again.
                let top_left_starting_index = ((0 * LARGE_CHUNK_WIDTH) + 0) * 4;
                let top_right_starting_index = ((0 * LARGE_CHUNK_WIDTH) + 7) * 4;
                let bottom_left_starting_index = ((7 * LARGE_CHUNK_WIDTH) + 0) * 4;
                let bottom_right_starting_index = ((7 * LARGE_CHUNK_WIDTH) + 7) * 4;

                let floating_offset = Vector2::new(offset.x as f64, offset.y as f64);
                let top_left_frame_position = Vector2::new(0.0, 0.0) + floating_offset;
                let top_right_frame_position = Vector2::new(7.0, 0.0) + floating_offset;
                let bottom_left_frame_position = Vector2::new(0.0, 7.0) + floating_offset;
                let bottom_right_frame_position = Vector2::new(7.0, 7.0) + floating_offset;

                let top_left_color = world.trace_ray_parameter_shortcut(top_left_frame_position);
                let top_right_color = world.trace_ray_parameter_shortcut(top_right_frame_position);
                let bottom_left_color =
                    world.trace_ray_parameter_shortcut(bottom_left_frame_position);
                let bottom_right_color =
                    world.trace_ray_parameter_shortcut(bottom_right_frame_position);

                Self::set_chunk_frame_pixel_to_color_by_starting_index(
                    chunk_frame,
                    top_left_starting_index,
                    top_left_color,
                );
                Self::set_chunk_frame_pixel_to_color_by_starting_index(
                    chunk_frame,
                    top_right_starting_index,
                    top_right_color,
                );
                Self::set_chunk_frame_pixel_to_color_by_starting_index(
                    chunk_frame,
                    bottom_left_starting_index,
                    bottom_left_color,
                );
                Self::set_chunk_frame_pixel_to_color_by_starting_index(
                    chunk_frame,
                    bottom_right_starting_index,
                    bottom_right_color,
                );

                let rgbas = [
                    top_left_color,
                    top_right_color,
                    bottom_left_color,
                    bottom_right_color,
                ];

                if Self::needs_subdivision(rgbas) {
                    // subdivide
                    //dbg!("I need to subdivide to a medium chunk!");
                    self.subdivide(None);
                    self.render(world);
                } else {
                    // TODO: linear interpolation here
                }
            }
            // Assumes that one pixel in each of the corners of the 8x8 square is rendered.
            // This means that one pixel in each 4x4 chunk is rendered.
            Self::Medium(chunk_frame, offset) => {
                // Corresponds to top_left, top_right, bottom_left, and bottom_right in that order.
                let sub_chunk_offsets = [
                    Vector2::new(0, 0),
                    Vector2::new(4, 0),
                    Vector2::new(0, 4),
                    Vector2::new(4, 4),
                ];

                // Corresponds to the offsets of the corners of the sub chunks.
                let sub_chunk_corner_offsets: [Vector2<usize>; 4] = [
                    Vector2::new(0, 0),
                    Vector2::new(3, 0),
                    Vector2::new(0, 3),
                    Vector2::new(3, 3),
                ];

                // Find the colors in these pixels so we arent re-doing work we don't have to.
                let already_rendered_rgbas = {
                    // This refers to the 8x8 chunk.
                    let top_left_starting_index = ((0 * LARGE_CHUNK_WIDTH) + 0) * 4;
                    let top_right_starting_index = ((0 * LARGE_CHUNK_WIDTH) + 7) * 4;
                    let bottom_left_starting_index = ((7 * LARGE_CHUNK_WIDTH) + 0) * 4;
                    let bottom_right_starting_index = ((7 * LARGE_CHUNK_WIDTH) + 7) * 4;
                    [
                        Self::color_from_chunk_frame_by_starting_index(
                            chunk_frame,
                            top_left_starting_index,
                        ),
                        Self::color_from_chunk_frame_by_starting_index(
                            chunk_frame,
                            top_right_starting_index,
                        ),
                        Self::color_from_chunk_frame_by_starting_index(
                            chunk_frame,
                            bottom_left_starting_index,
                        ),
                        Self::color_from_chunk_frame_by_starting_index(
                            chunk_frame,
                            bottom_right_starting_index,
                        ),
                    ]
                };

                // TODO: refactor this out
                let mut small_subdivision_indices = Vec::new();

                // Iterates through each of the four sub_chunks
                for i in 0..4 {
                    let mut sub_chunk_corner_colors = [Rgba::new(0.0, 0.0, 0.0, 0.0); 4];

                    // total offset for the top left corner of the 4x4 chunk
                    for k in 0..4 {
                        // We can skip finding the color of the pixel if we already know it.
                        if k == i {
                            //dbg!("meow");
                            continue;
                        }

                        let total_offset = sub_chunk_offsets[i] + sub_chunk_corner_offsets[k];

                        let color = world.trace_ray_parameter_shortcut(Vector2::new(
                            total_offset.x as f64,
                            total_offset.y as f64,
                        ));

                        let starting_index =
                            ((total_offset.y * MEDIUM_CHUNK_WIDTH) + total_offset.x) * 4;

                        Self::set_chunk_frame_pixel_to_color_by_starting_index(
                            chunk_frame,
                            starting_index,
                            color,
                        );

                        sub_chunk_corner_colors[i] = color;
                    }

                    if Self::needs_subdivision(sub_chunk_corner_colors) {
                        //dbg!("I need to subdivide to a small chunk!");
                        small_subdivision_indices.push(i);
                    } else {
                        // linear interpolation here
                    }
                }

                if !small_subdivision_indices.is_empty() {
                    *self = Self::Small(*chunk_frame, *offset, small_subdivision_indices.clone());
                    for i in 0..4 {
                        if small_subdivision_indices.contains(&i) {
                            continue;
                        }
                        // linear interpolation
                    }
                } else {
                    // linear interpolation here
                }

                self.render(world);
            }
            // only needs to render the subdivided parts
            Self::Small(chunk_frame, offset, section_indices) => {
                let medium_chunk_offsets = [
                    Vector2::new(0, 0),
                    Vector2::new(4, 0),
                    Vector2::new(0, 4),
                    Vector2::new(4, 4),
                ];

                let small_chunk_offsets = [
                    Vector2::new(0, 0),
                    Vector2::new(2, 0),
                    Vector2::new(0, 2),
                    Vector2::new(2, 2),
                ];

                let small_chunk_corner_offsets = [
                    Vector2::new(0, 0),
                    Vector2::new(1, 0),
                    Vector2::new(0, 1),
                    Vector2::new(1, 1),
                ];

                for section_index in section_indices {
                    let chunk_offset: Vector2<usize> =
                        medium_chunk_offsets[*section_index] + small_chunk_offsets[*section_index];

                    for corner_offset in small_chunk_corner_offsets {
                        let frame_coords = chunk_offset + corner_offset;

                        let color = world.trace_ray_parameter_shortcut(Vector2::new(
                            frame_coords.x as f64,
                            frame_coords.y as f64,
                        ));

                        let starting_index =
                            ((frame_coords.y * SMALL_CHUNK_WIDTH) + frame_coords.x) * 4;

                        Self::set_chunk_frame_pixel_to_color_by_starting_index(
                            chunk_frame,
                            starting_index,
                            color,
                        );
                    }
                }
            }
        }
    }

    /// Translates pixels in chunks to pixels in the frame buffer,
    #[allow(clippy::erasing_op)]
    #[allow(clippy::identity_op)]
    pub fn apply_to_frame(&self, frame: &mut [u8], frame_width: usize) {
        match self {
            Self::Large(chunk_frame, offset) => {
                Self::translate_chunk_frame_to_frame(
                    frame,
                    frame_width,
                    chunk_frame,
                    LARGE_CHUNK_WIDTH,
                    *offset,
                );
            }
            Self::Medium(chunk_frame, offset) => {
                Self::translate_chunk_frame_to_frame(
                    frame,
                    frame_width,
                    chunk_frame,
                    LARGE_CHUNK_WIDTH,
                    *offset,
                );
            }
            Self::Small(chunk_frame, offset, _) => {
                Self::translate_chunk_frame_to_frame(
                    frame,
                    frame_width,
                    chunk_frame,
                    LARGE_CHUNK_WIDTH,
                    *offset,
                );
            }
        }
    }

    fn needs_subdivision(rgbas: [Rgba; 4]) -> bool {
        for i in 1..4 {
            let diff = rgbas[0].difference(rgbas[i]);
            if diff > ALLOWED_COLOR_DIFFERENCE {
                return true;
            }
        }

        for i in 2..4 {
            let diff = rgbas[1].difference(rgbas[i]);
            if diff > ALLOWED_COLOR_DIFFERENCE {
                return true;
            }
        }

        let diff = rgbas[2].difference(rgbas[3]);
        if diff > ALLOWED_COLOR_DIFFERENCE {
            return true;
        }

        false
    }

    fn translate_chunk_frame_to_frame(
        frame: &mut [u8],
        frame_width: usize,
        chunk_frame: &[u8; 256],
        chunk_width: usize,
        offset: Vector2<usize>,
    ) {
        for y in 0..8 {
            for x in 0..8 {
                let chunk_frame_starting_index = ((y * chunk_width) + x) * 4;
                let rgba_slice =
                    &chunk_frame[chunk_frame_starting_index..chunk_frame_starting_index + 4];

                let frame_starting_index = (((y + offset.y) * frame_width) + (x + offset.x)) * 4;
                let pixel = &mut frame[frame_starting_index..frame_starting_index + 4];
                pixel.copy_from_slice(rgba_slice);
            }
        }
    }

    fn color_from_chunk_frame_by_starting_index(
        chunk_frame: &mut [u8; 256],
        starting_index: usize,
    ) -> Rgba {
        let r = chunk_frame[starting_index];
        let g = chunk_frame[starting_index + 1];
        let b = chunk_frame[starting_index + 2];
        let a = chunk_frame[starting_index + 3];

        Rgba::new(r as f64, g as f64, b as f64, a as f64)
    }

    fn set_chunk_frame_pixel_to_color_by_starting_index(
        chunk_frame: &mut [u8; 256],
        starting_index: usize,
        rgba: Rgba,
    ) {
        chunk_frame[starting_index] = rgba.r as u8;
        chunk_frame[starting_index + 1] = rgba.g as u8;
        chunk_frame[starting_index + 2] = rgba.b as u8;
        chunk_frame[starting_index + 3] = rgba.a as u8;
    }

    // TODO: make this only enabled by a feature or something.
    /// Draws all big chunks blue.
    pub fn debug(chunks: &Vec<Self>, canvas_width: usize, frame: &mut [u8]) {
        for chunk in chunks {
            match chunk {
                // outlines the large chunks in blue
                Self::Large(_chunk_frame, offset) => {
                    for y in offset.y..(offset.y + 8) {
                        for x in offset.x..(offset.x + 8) {
                            // filter out insides
                            if (x > offset.x)
                                && (x < (offset.x + 7))
                                && (y > offset.y)
                                && (y < offset.y + 7)
                            {
                                continue;
                            }

                            let frame_index = ((y * canvas_width) + x) * 4;

                            // multiplying by 8 instead of the canvas width because the width of the chunk is 8x8
                            // this is just here for future reference
                            let _chunk_frame_index = (((y - offset.y) * 8) + (x - offset.x)) * 4;

                            frame[frame_index] = 0;
                            frame[frame_index + 1] = 0;
                            frame[frame_index + 2] = 255;
                            frame[frame_index + 3] = 255;
                        }
                    }
                }
                Self::Medium(_chunk_frame, offset) => {
                    for y in offset.y..(offset.y + 8) {
                        for x in offset.x..(offset.x + 8) {
                            // filter out inside of top left subchunk
                            if (x > offset.x)
                                && (x < (offset.x + 4))
                                && (y > offset.y)
                                && (y < offset.y + 4)
                            {
                                continue;
                            }

                            if (x > offset.x + 4)
                                && (x < offset.x + 8)
                                && (y > offset.y)
                                && (y < offset.y + 4)
                            {
                                continue;
                            }

                            if (x > offset.x)
                                && (x < (offset.x + 4))
                                && (y > offset.y + 4)
                                && (y < offset.y + 8)
                            {
                                continue;
                            }

                            if (x > offset.x + 4)
                                && (x < (offset.x + 8))
                                && (y > offset.y + 4)
                                && (y < offset.y + 8)
                            {
                                continue;
                            }

                            let frame_index = ((y * canvas_width) + x) * 4;

                            // multiplying by 8 instead of the canvas width because the width of the chunk is 8x8
                            // this is just here for future reference
                            let _chunk_frame_index = (((y - offset.y) * 4) + (x - offset.x)) * 4;

                            frame[frame_index] = 0;
                            frame[frame_index + 1] = 255;
                            frame[frame_index + 2] = 0;
                            frame[frame_index + 3] = 255;
                        }
                    }
                }
                Self::Small(_chunk_frame, offset, section_indices) => {
                    let medium_chunk_offsets = [
                        Vector2::new(0, 0),
                        Vector2::new(4, 0),
                        Vector2::new(0, 4),
                        Vector2::new(4, 4),
                    ];

                    let section_offsets = [
                        Vector2::new(0, 0),
                        Vector2::new(2, 0),
                        Vector2::new(0, 2),
                        Vector2::new(2, 2),
                    ];

                    for k in 0..4 {
                        if !section_indices.contains(&k) {
                            let med_chunk_offset = medium_chunk_offsets[k];
                            for y in med_chunk_offset.y..(med_chunk_offset.y + 4) {
                                for x in med_chunk_offset.x..(med_chunk_offset.x + 4) {
                                    if (x > offset.x)
                                        && (x < (offset.x + 4))
                                        && (y > offset.y)
                                        && (y < offset.y + 4)
                                    {
                                        continue;
                                    }

                                    let frame_index = ((y * canvas_width) + x) * 4;

                                    frame[frame_index] = 0;
                                    frame[frame_index + 1] = 255;
                                    frame[frame_index + 2] = 0;
                                    frame[frame_index + 3] = 255;
                                }
                            }
                        }
                    }

                    for i in section_indices {
                        let med_chunk_offset = medium_chunk_offsets[*i];
                        for section_offset in section_offsets {
                            let frame_coords = *offset + med_chunk_offset + section_offset;

                            let frame_index =
                                ((frame_coords.y * canvas_width) + frame_coords.x) * 4;

                            frame[frame_index] = 255;
                            frame[frame_index + 1] = 0;
                            frame[frame_index + 2] = 0;
                            frame[frame_index + 3] = 255;
                        }
                    }
                }
            }
        }
    }
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

    // TODO: choose a better name for this.
    /// The same as trace_ray(), except it takes a world parameter.
    /// The other function exists for when you do not have access to a shared mutable reference to a [`World`] you can pass in.
    /// The original trace_ray() may be able to be factored out.
    pub fn trace_ray_parameter_shortcut(&self, frame_position: Vector2<f64>) -> Rgba {
        World::trace_ray(
            self.camera.position,
            self.camera,
            World::frame_to_viewport_coords(frame_position, &self.canvas, &self.viewport),
            self.background_color,
            &self.spheres,
            &self.lights,
            1.0,
            f64::MAX,
            self.reflection_passes,
        )
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
