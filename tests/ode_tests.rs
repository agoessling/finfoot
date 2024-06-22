use nalgebra::DVector;
use paste::paste;
use speculoos::prelude::*;

use finfoot::ode::dopri5;
use test_util::{all_problems, OdeProblem};

fn assert_dvector_close(a: &DVector<f64>, b: &DVector<f64>, tolerance: &DVector<f64>, name: &str) {
    assert!(a.len() == b.len() && b.len() == tolerance.len());

    for i in 0..a.len() {
        assert_that!(a[i])
            .named(&format!("{name} element {i}"))
            .is_close_to(b[i], tolerance[i]);
    }
}

fn test_problem(problem: &OdeProblem) {
    const CONFIG: dopri5::Config = dopri5::Config {
        rel_tol: 1e-4,
        abs_tol: 1e-6,
    };

    let input = dopri5::Input {
        t_span: problem.t_span,
        y0: &problem.y0,
        h0: (problem.t_span[1] - problem.t_span[0]) / 100.0,
        f: &problem.f,
    };
    let output = dopri5::integrate(&input, &CONFIG);
    assert_that!(output).named(&problem.name).is_ok();
    assert_dvector_close(
        &output.unwrap().y,
        &problem.yf,
        &problem.tolerance,
        &problem.name,
    );
}

macro_rules! generate_tests {
    ($($name:ident,)*) => {
        $(
            paste! {
                #[test]
                fn [<test_ $name>]() {
                    test_problem(&all_problems()[stringify!($name)]);
                }
            }
        )*
    };
}

generate_tests! {
    exponential,
    harmonic_oscillator,
    van_der_pol_oscillator,
    lorentz_attractor,
    robertson_equations,
    coupled_oscillators,
}
