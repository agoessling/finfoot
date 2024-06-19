use nalgebra::DVector;

pub mod dopri5;

pub type DerivativeFunc = dyn Fn(f64, &DVector<f64>) -> DVector<f64>;

#[derive(Debug)]
pub enum InputError {
    TimeSpan,
    StepSize,
}

#[derive(Debug)]
pub enum Error {
    Input(InputError),
    Convergence,
}

impl From<InputError> for Error {
    fn from(err: InputError) -> Error {
        Error::Input(err)
    }
}
