use nalgebra::{dvector, DVector};
use speculoos::prelude::*;

use finfoot::ode::{dopri5, DerivativeFunc};

struct OdeProblem {
    name: String,
    t_span: [f64; 2],
    y0: DVector<f64>,
    f: Box<DerivativeFunc>,
    yf: DVector<f64>,
    tolerance: DVector<f64>,
}

fn calc_tolerance(y: &DVector<f64>, rel_tol: f64, abs_tol: f64) -> DVector<f64> {
    (rel_tol * y).abs().map(|x| x.max(abs_tol.abs()))
}

fn all_problems() -> Vec<OdeProblem> {
    let mut problems: Vec<OdeProblem> = Vec::new();

    // Exponential decay.
    let name = String::from("exponential");
    let t_span = [0.0, 1.0];
    let tau = 1.0;
    let y0 = dvector![1.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| -tau * y);
    let yf = dvector![f64::exp(-tau * t_span[1])];
    let tolerance = calc_tolerance(&yf, 1e-4, 1e-6);

    problems.push(OdeProblem {
        name,
        t_span,
        y0,
        f,
        yf,
        tolerance,
    });

    // Harmonic oscillator.
    let name = String::from("harmonic oscillator");
    let t_span = [0.0, 1.0];
    let omega = 2.0 * std::f64::consts::PI;
    let y0 = dvector![1.0, 0.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| {
        let x = y[0];
        let dx = y[1];
        dvector![dx, -omega.powi(2) * x]
    });
    let yf = dvector![f64::cos(omega * t_span[1]), -f64::sin(omega * t_span[1])];
    let tolerance = calc_tolerance(&yf, 1e-4, 2e-4);

    problems.push(OdeProblem {
        name,
        t_span,
        y0,
        f,
        yf,
        tolerance,
    });

    // Van der Pol oscillator.
    let name = String::from("van der pol oscillator");
    let t_span = [0.0, 15.0];
    let mu = 5.0;
    let y0 = dvector![1.0, 0.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| {
        let x = y[0];
        let dx = y[1];
        dvector![dx, mu * (1.0 - x.powi(2)) * dx - x]
    });
    let yf = dvector![-1.7479415, 0.16720312]; // Calculated in ode_solutions.py
    let tolerance = calc_tolerance(&yf, 1e-4, 1e-6);

    problems.push(OdeProblem {
        name,
        t_span,
        y0,
        f,
        yf,
        tolerance,
    });

    // Lorentz attractor.
    let name = String::from("lorentz attractor");
    let t_span = [0.0, 5.0];
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let y0 = dvector![1.0, 1.0, 1.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| {
        let (x, y, z) = (y[0], y[1], y[2]);
        dvector![sigma * (y - x), x * (rho - z) - y, x * y - beta * z]
    });
    let yf = dvector![-6.51226597, -6.97427953, 23.92417887]; // Calculated in ode_solutions.py
    let tolerance = calc_tolerance(&yf, 1e-4, 7e-3);

    problems.push(OdeProblem {
        name,
        t_span,
        y0,
        f,
        yf,
        tolerance,
    });

    // Robertson equations.
    let name = String::from("robertson equations");
    let t_span = [0.0, 30.0];
    let y0 = dvector![1.0, 0.0, 0.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| {
        let (y1, y2, y3) = (y[0], y[1], y[2]);
        dvector![
            -0.04 * y1 + 1e4 * y2 * y3,
            0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2.powi(2),
            3e7 * y2.powi(2),
        ]
    });
    let yf = dvector![7.44269918e-01, 1.03732456e-05, 2.55719708e-01]; // Calculated in ode_solutions.py
    let tolerance = calc_tolerance(&yf, 1e-4, 4e-4);

    problems.push(OdeProblem {
        name,
        t_span,
        y0,
        f,
        yf,
        tolerance,
    });

    problems
}

fn assert_dvector_close(
    a: &DVector<f64>,
    b: &DVector<f64>,
    tolerance: &DVector<f64>,
    name: &str,
) -> () {
    assert!(a.len() == b.len() && b.len() == tolerance.len());

    for i in 0..a.len() {
        assert_that!(a[i])
            .named(&format!("{} element {}", name, i))
            .is_close_to(b[i], tolerance[i]);
    }
}

#[test]
fn test_dopri5() -> () {
    const CONFIG: dopri5::Config = dopri5::Config {
        rel_tol: 1e-4,
        abs_tol: 1e-6,
    };

    for problem in all_problems() {
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
}
