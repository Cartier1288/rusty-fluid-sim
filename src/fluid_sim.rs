use crate::vectors::vec::*;
use crate::vectors::vec_field::*;
use std::mem::swap;

const SPARSE_SOLVE_ACCURACY: usize = 15;

pub struct fluid_c {
    pub N: usize,           // side-length of (square) container
    W: usize,
    H: usize,
    v_field: vec_field<f32, 2>,     //  
    density: vec_field<f32, 1>       // 
}

// code taken directly from the minifb rust page:
// https://docs.rs/minifb/latest/minifb/struct.Window.html#method.update_with_buffer 
fn from_u8_rgb(r: u8, g: u8, b: u8) -> u32 {
    let (r, g, b) = (r as u32, g as u32, b as u32);
    return (r << 16) | (g << 8) | b;
}

fn from_one_u8_rgb(n: u8) -> u32 {
    let n32 = n as u32;
    return (n32 << 16) | (n32 << 8) | n32;
}

#[inline(always)]
fn draw_line(p1: vec<f32,2>, p2: vec<f32,2>, buffer: &mut Vec<u32>, W: usize, H: usize, color: u32) {
    // there's probably a better way to avoid skipping pixels than just picking whichever stride 
    // is the largest, but this is what came to my head fastest, and it works for now :D
    //print!("({},{}) -> ({},{})", p1[0], p1[1], p2[0], p2[1]);

    let diffX = (p1[0] - p2[0]).abs();
    let diffY = (p1[1] - p2[1]).abs();

    if diffX > diffY { // let y be the dependent variable
        let end = (p2[0] - p1[0]).round().abs() as usize;
        let c = (p2[0] - p1[0]).abs() / (p2[0] - p1[0]);

        let slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);

        for i in 0..=end {
            let mut x: f32 = c*(i as f32);
            let y: f32 = (p1[1] + slope*x).round();
            x += p1[0];

            let idx: f32 = (y * (W as f32)) + x;

            buffer[idx as usize] = color; 
        }
    }
    else if diffX != 0.0 || diffY != 0.0 {
        let end = (p2[1] - p1[1]).abs().round() as usize;
        let c = (p2[1] - p1[1]).abs() / (p2[1] - p1[1]);

        let slope = (p2[0] - p1[0]) / (p2[1] - p1[1]);

        for i in 0..=end {
            let mut y: f32 = c*(i as f32);
            let x: f32 = (p1[0] + slope*y).round();
            y += p1[1];

            let idx: f32 = (y * (W as f32)) + x;

            buffer[idx as usize] = color; 
        }
    }
}

impl fluid_c {
    pub fn new(N:usize) -> fluid_c {
        let W = N + 2;
        let H = N + 2;

        let mut field: vec_field<f32, 2> = vec_field::new(W, H, vec::new([1.0, 4.0]));
        let mut density: vec_field<f32, 1> = vec_field::new(W, H, vec::new([0.0]));

        for c in 1..=N {
            for r in 1..=N {
                if c > (N/4) && c < 3*(N/4) && r > (N/4) && r < 3*(N/4) {
                    // field.set(r,c,
                    //     vec::new([0.5, 2.0]));
                    density.set(r,c,
                        vec::new([100.0]));
                }
                else {
                    field.set(r,c,
                        -field.at(r,c));
                }
            }
        }

        return fluid_c {
            N,
            W,
            H,
            v_field: field, // for now we are simply using 2 dimensional velocity-field
            density: density
        };
    }

    pub fn resetFields(&mut self) {
        self.v_field = vec_field::new(self.W, self.H, vec::new([0.0, 0.0]));
        self.density = vec_field::new(self.W, self.H, vec::new([0.0]));
    }

    // dt: time step
    pub fn advance(&mut self, dt: f32) {
        self.advanceVelocity(dt);
        self.advanceDensity(dt);
    }

    pub fn render(&self) -> Vec<u32> {
        let len = self.density.len();
        let mut v: Vec<u32> = Vec::with_capacity(len);

        for i in 0..len {
            //let norm: u8 = ((((self.density.atIndex(i)[0] / f32::MAX) + 1.0) / 2.0) * (u8::MAX as f32)) as u8;
            let norm: u8 = (self.density.atIndex(i)[0] / 100.0 * u8::MAX as f32) as u8;
            v.push(from_u8_rgb(norm, norm, norm));
        }

        return v;
    }

    pub fn renderVelocity(&self, cell:usize, border:usize) -> Vec<u32> {
        let len = self.v_field.len();
        let mut v: Vec<u32> = vec![0;len*cell*cell];

        let blue = from_u8_rgb(0,0,255);
        let red = from_u8_rgb(255,0,0);
        let white = from_u8_rgb(255,255,255);

        let center = cell / 2;
        let left = border;
        let right = cell - border - 1;

        for x in 0..self.W {
            for y in 0..self.H {
                let top = (x*cell) + (((y*cell) + border) * (cell*self.W));

                let mid = top + center*(cell*self.W);
                let bot = top + (cell-1)*(cell*self.W);

                let mut vel: vec<f32,2> = self.v_field.at(y,x).clone();

                let norm: u32 = from_one_u8_rgb((self.density.at(y,x)[0] / 120.0 * u8::MAX as f32) as u8);

                for dy in border..(cell-border) {
                    let depth = top + (dy*(cell*self.W));
                    for dx in border..(cell-border) {
                        v[depth + dx] = norm;
                    }
                }


                // if vel[0] > 0.0 { // rightward force
                //     v[mid + 2] = red;
                //     if vel[1] > 0.0 { // downward force
                //         v[bot+right] = red;
                //     }
                //     else if vel[1] < 0.0 { // upward force
                //         v[top+right] = red;
                //     }
                // }
                // else if vel[0] < 0.0 { // leftward force
                //     v[mid] = red;
                //     if vel[1] > 0.0 { // downward force
                //         v[bot+left] = red;
                //     }
                //     else if vel[1] < 0.0 { // upward force
                //         v[top+left] = red;
                //     }
                // }

                // if vel[1] > 0.0 { // downward force
                //     v[bot+center] = red;
                // }
                // else if vel[1] < 0.0 { // upward force
                //     v[top+center] = red;
                // }

                // if vel[0] == 0.0 && vel[1] == 0.0 { // no force
                //     v[mid+center] = red;
                // }

                // let toptop = (x*cell) + ((y*cell) * (cell*self.W));
                // for i in 0..cell {
                //     v[toptop + i] = blue;
                //     v[toptop + i + 4*(cell*self.W)] = blue;
                //     v[toptop + i*(cell*self.W)] = blue;
                //     v[toptop + i*(cell*self.W) + 4*(self.W) - 1] = blue;
                // }
            }
        }

        for x in 0..self.W {
            for y in 0..self.H {
                let mut vel: vec<f32,2> = self.v_field.at(y,x).clone();

                let p1:vec<f32,2> = vec::new([(x*cell + center) as f32, (y*cell + center + 2) as f32]);
                let mut p2:vec<f32,2> = (&vel * center as f32 / 2.0) + &p1;
                
                // ensure we aren't going out of bounds
                // obviously this will distort the lines if we aren't careful
                p2[0] = p2[0].min((self.W*cell) as f32);
                p2[0] = p2[0].max(0.0);
                p2[1] = p2[1].min((self.H*cell) as f32);
                p2[1] = p2[1].max(0.0);

                draw_line(p1, p2, &mut v, self.W * cell, self.H * cell, red);
            }
        }

        return v;
    }

    fn enforceBound<const M: usize>(&self, negate: bool, next: &mut vec_field<f32, M>) {
        for i in 1..=self.N {
            let mut left = next.at(i,1).clone();
            let mut right = next.at(i,self.N).clone();
            let mut top = next.at(1,i).clone();
            let mut bottom = next.at(self.N,i).clone();

            if negate {
                left[0] = -left[0];         // if left, go right
                right[0] = -right[0];       // if right, go left
                top[1] = -top[1];           // if up, go down
                bottom[1] = -bottom[1];     // if down, go up
            }

            next.set(       i,0,    left);
            next.set(i,self.N+1,    right);
            next.set(       0,i,    top);
            next.set(self.N+1,i,    bottom);
        }

        // corners are the average of their two adjacent cells
        next.set(0,0,               (next.at(0,1) + next.at(1,0)) * 0.5);       // top left
        next.set(self.N+1,self.N+1, (next.at(self.N,self.N+1) + next.at(self.N+1,self.N)) * 0.5);   // bottom right
        next.set(0,self.N+1,        (next.at(0,self.N) + next.at(1,self.N+1)) * 0.5);     // top right
        next.set(self.N+1,0,        (next.at(self.N,0) + next.at(self.N+1,1)) * 0.5);     // bottom left
    }

    fn sparse_solve<const M: usize>(&self, boundrev: bool, a: f32, d: f32, current: &vec_field<f32, M>, next: &mut vec_field<f32, M>, steps: usize) {
        for _ in 0..steps {
            for r in 1..=self.N {
                let depth: usize = r * self.W;
                for c in 1..=self.N {
                    let start: usize = depth + c;
                
                    next.setAtIndex(start,
                        (current.atIndex(start)
                         + ((next.at(r-1, c) + next.at(r+1, c)
                             + next.at(r, c-1) + next.at(r, c+1)) * a))
                        / d
                    );
                }
            }

            self.enforceBound(boundrev, next);
        }
    }

    fn addSource<const M: usize>(from: &vec_field<f32, M>, to: &mut vec_field<f32, M>, dt: f32) {
        assert!(from.W == to.W && from.H == to.H);

        let size = to.W * to.H;

        for i in 0..size {
            to.setAtIndex(i,
                to.atIndex(i) + (from.atIndex(i) * dt)
            )
        }
    }

    fn diffuse<const M: usize>(&self, boundrev: bool, current: &vec_field<f32, M>, next: &mut vec_field<f32, M>, diff: f32, dt: f32) {
        let a: f32 = dt * diff * self.N as f32 * self.N as f32;
        let d: f32 = 1.0 + 4.0*a;

        self.sparse_solve(false, a, d, current, next, SPARSE_SOLVE_ACCURACY);
    }

    /*
     * The idea with advection is to look at each cell and consider the cells that would reach its center.
     * 
     * That is, at cell (r,c), we look backward dt0 * v[0,0] (approximate our distance function), and linearly 
     * interpolate our density based off of the surrounding cells.
     */
    fn advect<const M: usize>(&self, boundrev: bool, prev: &vec_field<f32, M>, next: &mut vec_field<f32, M>, dt: f32) {
        let dt0 = dt * self.N as f32;

        let mut p = vec::new([1.0, 1.0]);

        for c in 1..=self.N {
            p[0] = c as f32;
            for r in 1..=self.N {
                p[1] = r as f32;

                // see how far back we go with this timestep
                let mut d = &p - (self.v_field.at(r, c) * dt0);

                d.apply(|x: &mut f32| {
                    if *x < 0.5 {
                        *x = 0.5;
                    }
                    else if *x > (self.N as f32) + 0.5 {
                        *x = (self.N as f32) + 0.5;
                    }
                });

                let xFlr = d[0].floor() as usize; let xNext = xFlr + 1;
                let yFlr = d[1].floor() as usize; let yNext = yFlr + 1;

                let diffX = d[0] - (xFlr as f32); let diffXInv = 1.0 - diffX;
                let diffY = d[1] - (yFlr as f32); let diffYInv = 1.0 - diffY;

                // println!("({}, {}) acted on by -({}, {}) goes to ({}, {}): \n xFlr {}, xNext {}, yFlr {}, yNext {}, \n diffX {}, diffXInv {}, diffY {}, diffYInv {}",
                //     c, r, self.v_field.at(r,c)[0], self.v_field.at(r,c)[1], d[0], d[1],
                //     xFlr, xNext, yFlr, yNext, diffX, diffXInv, diffY, diffYInv);

                // println!("
                //     (({} * {} + {} * {}) * {}) +
                //     (({} * {} + {} * {}) * {})
                // ",
                //     prev.at(yFlr, xFlr)[0], diffYInv, prev.at(yNext, xFlr)[0], diffY, diffXInv,
                //     prev.at(yFlr, xNext)[0], diffYInv, prev.at(yNext, xNext)[0], diffY, diffX
                // );

                next.set(r,c,
                    // interpolate density, based on points around vec-d
                    (((prev.at(yFlr, xFlr) * diffYInv) + (prev.at(yNext, xFlr) * diffY)) * diffXInv) +
                    (((prev.at(yFlr, xNext) * diffYInv) + (prev.at(yNext, xNext) * diffY)) * diffX)
                )
            }
        }

        self.enforceBound(boundrev, next);
    }

    fn project(&self, current: &vec_field<f32, 2>, next: &mut vec_field<f32, 2>) {
        let h = 1.0 / self.N as f32;

        let mut p = vec_field::new(current.W, current.H, vec::new([0.0]));
        let mut div = vec_field::new(current.W, current.H, vec::new([0.0]));

        for c in 1..=self.N {
            for r in 1..=self.N {
                div.set(r,c, vec::new([
                    -(h/2.0) * (next.at(r,c+1)[0] - next.at(r,c-1)[0] +
                                   next.at(r+1,c)[1] - next.at(r-1,c)[1])
                ]));
            }
        }
        self.enforceBound(false, &mut div);

        self.sparse_solve(false, 1.0, 4.0, &div, &mut p, SPARSE_SOLVE_ACCURACY);

        for r in 1..=self.N {
            let depth: usize = r * self.W;
            for c in 1..=self.N {
                let start: usize = depth + c;

                let v: vec<f32, 2> = vec::new([
                    (0.5/h) * (p.at(r,c+1)[0] - p.at(r, c-1)[0]),
                    (0.5/h) * (p.at(r+1,c)[0] - p.at(r-1, c)[0])
                ]);

                next.set(r,c,
                    next.at(r,c) - v
                );
            }
        }
        self.enforceBound(true, next);
    }

    fn advanceDensity(&mut self, dt: f32) {
        // let mut from = vec_field::new(self.W, self.H, vec::new([0.0]));
        let mut next = self.density.clone();

        // fluid_c::addSource(&from, &mut next, dt);
        // swap(&mut self.density, &mut next);
        
        self.diffuse(false, &self.density, &mut next, 0.003, dt);
        swap(&mut self.density, &mut next);

        self.advect(false, &self.density, &mut next, dt);
        swap(&mut self.density, &mut next);
    }

    fn advanceVelocity(&mut self, dt: f32) {
        // let mut from = vec_field::new(self.W, self.H, vec::new([0.0, 0.0]));
        let mut next = self.v_field.clone();
        self.enforceBound(true, &mut next);

        // fluid_c::addSource(&from, &mut next, dt);
        // swap(&mut self.v_field, &mut next);

        self.diffuse(true, &self.v_field, &mut next, 0.003, dt);

        self.project(&self.v_field, &mut next);
        swap(&mut self.v_field, &mut next);

        self.advect(true, &self.v_field, &mut next, dt);

        self.project(&self.v_field, &mut next);
        swap(&mut self.v_field, &mut next);
    }
}