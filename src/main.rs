#![feature(trait_alias)]

extern crate minifb;

mod vectors;
mod fluid_sim;
use crate::fluid_sim::*;

use std::{thread, time};
use std::time::{Instant, Duration};

use minifb::{
    Key,
    WindowOptions,
    Window,
};

const W: usize = 64;
const H: usize = 64;


fn main() {
    let mut sim: fluid_c = fluid_c::new(W);

    let renderField = true;

    if renderField {
        let cell = 19;

        let mut opts = WindowOptions::default();
        let mut win = Window::new("Fluid Simulation", (W+2)*cell, (H+2)*cell, opts).unwrap();

        while win.is_open() {
            // let buffer = sim.render();
            // win.update_with_buffer(&buffer, W + 2, H + 2).unwrap();
            let mut now = Instant::now();
            let buffer = sim.renderVelocity(cell, 2); // render velocity buffer
            let mut after = Instant::now();
            println!("rendering took {:?}", after.duration_since(now));

            now = Instant::now();
            win.update_with_buffer(&buffer, (W+2)*cell, (W+2)*cell).unwrap();
            after = Instant::now();
            println!("updating buffer took {:?}", after.duration_since(now));

            now = Instant::now();
            sim.advance(0.01);
            after = Instant::now();
            println!("dynamics took {:?}\n\n", after.duration_since(now));
        }
    }
    else {
        let mut opts = WindowOptions::default();
        opts.scale = minifb::Scale::X8;
        let mut win = Window::new("Fluid Simulation", W+2, H+2, opts).unwrap();

        while win.is_open() {
            let mut now = Instant::now();

            let buffer = sim.render();

            let mut after = Instant::now();
            println!("rendering took {:?}", after.duration_since(now));

            now = Instant::now();

            win.update_with_buffer(&buffer, W + 2, H + 2).unwrap();

            after = Instant::now();
            println!("updating buffer took {:?}", after.duration_since(now));

            now = Instant::now();

            sim.advance(0.01);

            after = Instant::now();
            println!("dynamics took {:?}\n\n", after.duration_since(now));
        }
    }
}