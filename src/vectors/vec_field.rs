use crate::vectors::vec::*;

pub struct vec_field<T: VecTrait<T>, const N: usize> {
    pub W: usize,               // width of field
    pub H: usize,               // height of field
    field: Vec<vec<T, N>>       // field is of length N * W * H
                                // order: ( v1_1, v1_2, ..., v1_N, ..., v(W*H)_1, ..., v(W*H)_N )
}

impl<T: VecTrait<T>, const N: usize> vec_field<T, N> {
    pub fn new(W:usize, H:usize, default: vec<T, N>) -> vec_field<T, N> {
        let field: Vec<vec<T, N>> = vec![default; (W+2) * (H+2)];

        return vec_field {
            W,
            H,
            field
        };
    }

    pub fn len(&self) -> usize {
        return self.field.len();
    }

    pub fn clone(&self) -> vec_field<T, N> {
        return vec_field {
            W: self.W,
            H: self.H,
            field: self.field.clone()
        };
    }

    // returns an immutable reference to the vector field
    pub fn getField(&self) -> &Vec<vec<T, N>> {
        return &self.field;
    }

    // returns the vector at (r,c)
    pub fn at(&self, r:usize, c:usize) -> &vec<T, N> {
        assert!(r < self.H && c < self.W);
        return &self.field[(self.W * r) + c];
    }

    // returns a mutable reference at (r,c)
    pub fn atMut(&mut self, r:usize, c:usize) -> &mut vec<T, N> {
        assert!(r < self.H && c < self.W);
        return &mut self.field[(self.W * r) + c];
    }

    pub fn atIndex(&self, i:usize) -> &vec<T, N> {
        assert!(i < self.field.len());
        return &self.field[i];
    }

    // sets a single value of a vector of the field, assumes understanding
    // of the layout of field
    pub fn setAtIndex(&mut self, i: usize, v: vec<T, N>) {
        assert!(i < self.field.len());
        self.field[i] = v;
    }

    // sets a vector at (r,c)
    pub fn set(&mut self, r:usize, c:usize, v:vec<T, N>) {
        assert!(r < self.H && c < self.W);
        self.field[(self.W * r) + c] = v;
    }

    pub fn setByRef(&mut self, r:usize, c:usize, v:&vec<T, N>) {
        assert!(r < self.H && c < self.W);
        self.field[(self.W * r) + c] = v.clone();
    }

    pub fn add(&mut self, V: &vec_field<T, N>, C: Option<f32>) {
        // make sure the fields have the same dimensions
        assert!(V.W == self.W && V.H == self.H);
        
        // some constant to multiply the field by, with a default value of 1
        let c = C.unwrap_or(1.0);

        for r in 1..self.H {
            let depth = r * self.W;

            for c in 1..self.W {
                self.field[depth + c] = &self.field[depth + c] + V.at(r, c);
            }
        }
    }
}