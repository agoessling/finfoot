pub mod dopri5;

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
