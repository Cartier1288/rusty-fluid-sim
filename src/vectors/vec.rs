use std::ops::Add;
use std::ops::Sub;
use std::ops::Mul;
use std::ops::Div;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Neg;
use std::ops::MulAssign;
use std::ops::DivAssign;
use std::ops::AddAssign;

pub trait VecTrait<T>: Clone + Copy + Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + MulAssign + DivAssign + AddAssign {}
impl<T> VecTrait<T> for T where T: Clone + Copy + Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + MulAssign + DivAssign + AddAssign {}

pub struct vec<T: VecTrait<T>, const N: usize> {
    val: [T; N]
}

impl<T: VecTrait<T>, const N: usize> Clone for vec<T, N> {
    fn clone(&self) -> vec<T, N> {
        return vec {
            val: self.val.clone()
        }
    }
}

impl<T: VecTrait<T>, const N: usize> Index<usize> for vec<T, N> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        return &self.val[i];
    }
}

impl<T: VecTrait<T>, const N: usize> IndexMut<usize> for vec<T, N> {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut T {
        return &mut self.val[i];
    }
}

// == <Mul operations> ===========================================================================//
impl<T: VecTrait<T>, const N: usize>
Mul<T> for vec<T, N> {
    type Output = vec<T, N>;

    fn mul(self, _rhs: T) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] *= _rhs;
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize>
Mul<T> for &vec<T, N> {
    type Output = vec<T, N>;

    fn mul(self, _rhs: T) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] *= _rhs;
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize>
MulAssign<T> for vec<T, N> {
    fn mul_assign(&mut self, _rhs: T) {
        for i in 0..N {
            self[i] *= _rhs;
        }
    }
}
// == </Mul operations> ==========================================================================//

// == <Div operations> ===========================================================================//
impl<T: VecTrait<T>, const N: usize>
Div<T> for vec<T, N> {
    type Output = vec<T, N>;

    fn div(self, _rhs: T) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] /= _rhs;
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize>
Div<T> for &vec<T, N> {
    type Output = vec<T, N>;

    fn div(self, _rhs: T) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] /= _rhs;
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize>
DivAssign<T> for vec<T, N> {
    fn div_assign(&mut self, _rhs: T) {
        for i in 0..N {
            self.val[i] /= _rhs;
        }
    }
}
// == </Div operations> ==========================================================================//

// == <Add operations> ===========================================================================//
impl<T: VecTrait<T>, const N: usize> Add<vec<T, N>> for vec<T, N> {
    type Output = vec<T, N>;

    fn add(self, _rhs: vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] + _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Add<vec<T, N>> for &vec<T, N> {
    type Output = vec<T, N>;

    fn add(self, _rhs: vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] + _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Add<&vec<T, N>> for vec<T, N> {
    type Output = vec<T, N>;

    fn add(self, _rhs: &vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] += _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Add<&vec<T, N>> for &vec<T, N> {
    type Output = vec<T, N>;

    fn add(self, _rhs: &vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] += _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize>
AddAssign<&vec<T,N>> for vec<T,N> {
    fn add_assign(&mut self, _rhs: &vec<T,N>) {
        for i in 0..N {
            self.val[i] += _rhs[i];
        }
    }
}

impl<T: VecTrait<T>, const N: usize>
AddAssign<T> for vec<T, N> {
    fn add_assign(&mut self, _rhs: T) {
        for i in 0..N {
            self.val[i] += _rhs;
        }
    }
}
// == </Add operations> ==========================================================================//

// == <Sub operations> ===========================================================================//
impl<T: VecTrait<T>, const N: usize> Sub<vec<T, N>> for vec<T, N> {
    type Output = vec<T, N>;

    fn sub(self, _rhs: vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] - _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Sub<vec<T, N>> for &vec<T, N> {
    type Output = vec<T, N>;

    fn sub(self, _rhs: vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] - _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Sub<&vec<T, N>> for vec<T, N> {
    type Output = vec<T, N>;

    fn sub(self, _rhs: &vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] - _rhs[i];
        }

        return new;
    }
}

impl<T: VecTrait<T>, const N: usize> Sub<&vec<T, N>> for &vec<T, N> {
    type Output = vec<T, N>;

    fn sub(self, _rhs: &vec<T, N>) -> Self::Output {
        let mut new: vec<T, N> = self.clone();

        for i in 0..N {
            new.val[i] = self[i] - _rhs[i];
        }

        return new;
    }
}
// == </Sub operations> ==========================================================================//

// == <Neg operations> ===========================================================================//
impl<const N: usize>
Neg for vec<f32, N> {
    type Output = vec<f32, N>;

    fn neg(self) -> Self::Output {
        let mut new: vec<f32, N> = self.clone();

        for i in 0..N {
            new.val[i] = -1.0 * self.val[i];
        }

        return new;
    }
}

impl<const N: usize>
Neg for &vec<f32, N> {
    type Output = vec<f32, N>;

    fn neg(self) -> Self::Output {
        let mut new: vec<f32, N> = self.clone();

        for i in 0..N {
            new.val[i] = -1.0 * self.val[i];
        }

        return new;
    }
}

// == </Neg operations> ==========================================================================//

impl<T: VecTrait<T>, const N: usize> vec<T, N> {
    pub fn new(v: [T;N]) -> vec<T, N> {
        return vec {
            val: v
        }
    }

    pub fn apply<F>(&mut self, f: F)
        where F: Fn(&mut T) {
        
        for x in self.val.iter_mut() {
            f(x);
        }
    }
}
