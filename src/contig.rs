use std::alloc::{alloc, dealloc, Layout};
use std::ptr::{self, NonNull};
use std::slice;

#[allow(dead_code)]
struct Vector<T> {
    _phantom: std::marker::PhantomData<T>,
}

#[allow(dead_code)]
struct ConcreteType {
    field0: f64,
    vec_field: Vector<f64>,
    field1: f64,
}

#[allow(dead_code)]
pub struct ReferenceType<'a> {
    pub field0: &'a mut f64,
    pub vec_vield: &'a mut [f64],
    pub field1: &'a mut f64,
    raw_buffer: NonNull<u8>,
}

impl ReferenceType<'_> {
    pub fn new(length: usize) -> Self {
        let layout = Layout::new::<f64>();
        let layout = layout
            .extend(Layout::array::<f64>(length).unwrap())
            .unwrap()
            .0;
        let layout = layout.extend(Layout::new::<f64>()).unwrap().0;

        let raw_buffer;
        unsafe {
            raw_buffer = NonNull::new(alloc(layout)).unwrap();
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_struct() {
        let _vec = Vector::<f64> {
            _phantom: std::marker::PhantomData,
        };
        let _test = ConcreteType {
            field0: 1.0,
            vec_field: _vec,
            field1: 2.0,
        };
        assert_eq!(_test.field0, 1.0);
    }
}
