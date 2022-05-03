use cgmath::Vector2;
use mini_raytracer::{Canvas, Sphere, Viewport};
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::WindowBuilder;
use winit_input_helper::WinitInputHelper;

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;

struct World {
    canvas: Canvas,
    viewport: Viewport,
    spheres: Vec<Sphere>,
}

impl World {
    pub fn new(width: u32, height: u32, spheres: Vec<Sphere>) -> Self {
        Self {
            canvas: Canvas::new(width, height),
            viewport: Viewport::new(width, height),
            spheres,
        }
    }

    fn draw(&self, frame: &mut [u8]) {
        // iterates through each pixel, we will need to do raytracing to find the color of each pixel
        for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
            let x = (i % WIDTH as usize) as i16;
            let y = (i / WIDTH as usize) as i16;

            // convert the frame coords to canvas coords
            let canvas_coords =
                self.frame_coords_to_canvas_coords(Vector2::new(x as f64, y as f64));

            // convert the canvas coords to the viewport coords
            let viewport_coords = self.viewport_coords_to_canvas_coords(canvas_coords);

            // do the math here

            let rgba = [0x48, 0xb2, 0xe8, 0xff];

            pixel.copy_from_slice(&rgba);
        }
    }

    pub fn canvas_coords_to_frame_coords(&self, coords: Vector2<f64>) -> Vector2<f64> {
        unimplemented!()
    }

    /// Sets frame coords (coords with the origin in the top left), to canvas coords (coords with the origin in the center)
    pub fn frame_coords_to_canvas_coords(&self, coords: Vector2<f64>) -> Vector2<f64> {
        unimplemented!()
    }

    pub fn canvas_coords_to_viewport_coords(&self, coords: Vector2<f64>) -> Vector2<f64> {
        unimplemented!()
    }

    pub fn viewport_coords_to_canvas_coords(&self, coords: Vector2<f64>) -> Vector2<f64> {
        unimplemented!()
    }
}

fn main() -> Result<(), Error> {
    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window = {
        let size = LogicalSize::new(WIDTH as f64, HEIGHT as f64);
        WindowBuilder::new()
            .with_title("Basic Ray Collisions")
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .unwrap()
    };

    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(WIDTH, HEIGHT, surface_texture)?
    };

    let spheres = vec![];
    let mut world = World::new(WIDTH, HEIGHT, spheres);

    event_loop.run(move |event, _, control_flow| {
        if let Event::RedrawRequested(_) = event {
            world.draw(pixels.get_frame());
            if pixels
                .render()
                .map_err(|e| panic!("pixels.render() failed: {}", e))
                .is_err()
            {
                *control_flow = ControlFlow::Exit;
                return;
            }
        }

        // Handle input events
        if input.update(&event) {
            // Close events
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
                return;
            }

            // Resize the window
            if let Some(size) = input.window_resized() {
                pixels.resize_surface(size.width, size.height);
            }

            window.request_redraw();
        }
    });
}
