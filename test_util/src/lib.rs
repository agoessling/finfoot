use std::collections::HashMap;

use nalgebra::{dvector, DVector};

use finfoot::ode::DerivativeFunc;

pub struct OdeProblem {
    pub name: String,
    pub t_span: [f64; 2],
    pub y0: DVector<f64>,
    pub f: Box<DerivativeFunc>,
    pub yf: DVector<f64>,
    pub tolerance: DVector<f64>,
}

fn calc_tolerance(y: &DVector<f64>, rel_tol: f64, abs_tol: f64) -> DVector<f64> {
    (rel_tol * y).abs().map(|x| x.max(abs_tol.abs()))
}

pub fn all_problems() -> HashMap<String, OdeProblem> {
    let mut problems = HashMap::new();

    // Exponential decay.
    let name = String::from("exponential");
    let t_span = [0.0, 1.0];
    let tau = 1.0;
    let y0 = dvector![1.0];
    let f: Box<DerivativeFunc> = Box::new(move |_, y| -tau * y);
    let yf = dvector![f64::exp(-tau * t_span[1])];
    let tolerance = calc_tolerance(&yf, 1e-4, 1e-6);

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    // Harmonic oscillator.
    let name = String::from("harmonic_oscillator");
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

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    // Van der Pol oscillator.
    let name = String::from("van_der_pol_oscillator");
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

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    // Lorentz attractor.
    let name = String::from("lorentz_attractor");
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

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    // Robertson equations.
    let name = String::from("robertson_equations");
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

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    // Coupled oscillators with full coupling
    let name = String::from("coupled_oscillators");
    let t_span = [0.0, 10.0];
    let n = 100; // Number of oscillators
    let k = 1.0; // Coupling constant

    let mut y0 = DVector::zeros(2 * n); // Initial positions and velocities
    for i in 0..n {
        y0[i] = (i as f64).sin(); // Initial positions
        y0[n + i] = (i as f64).cos(); // Initial velocities
    }

    let f: Box<DerivativeFunc> = Box::new(move |_, y| {
        let mut dydt = DVector::zeros(2 * n);
        for i in 0..n {
            dydt[i] = y[n + i]; // Velocity
            let mut sum = 0.0;
            for j in 0..n {
                sum += y[i] - y[j];
            }
            dydt[n + i] = -k * sum; // Acceleration
        }
        dydt
    });

    #[rustfmt::skip]
    let yf = dvector![
        -8.97722333e-02,  6.59104964e-01,  7.66020429e-01,  1.32676575e-01,
        -6.58633178e-01, -8.80382293e-01, -3.28695657e-01,  4.89208582e-01,
         8.21353038e-01,  3.62365631e-01, -4.65762735e-01, -9.01654658e-01,
        -5.44553115e-01,  2.77224382e-01,  8.08139393e-01,  5.60071104e-01,
        -2.38907643e-01, -8.54219474e-01, -7.20149528e-01,  4.00389043e-02,
         7.27432084e-01,  7.10043892e-01,  3.86095190e-03, -7.41855398e-01,
        -8.41496985e-01, -2.03453793e-01,  5.85660210e-01,  8.00337248e-01,
         2.43204243e-01, -5.73513290e-01, -8.98929017e-01, -4.33857220e-01,
         3.94117236e-01,  8.23758454e-01,  4.60056280e-01, -3.62603184e-01,
        -8.87870622e-01, -6.32817573e-01,  1.68061366e-01,  7.78441792e-01,
         6.37142756e-01, -1.25926060e-01, -8.09202705e-01, -7.84485783e-01,
        -7.44999186e-02,  6.67997159e-01,  7.60357061e-01,  1.17664519e-01,
        -6.69191908e-01, -8.76780049e-01, -3.14244325e-01,  5.01222514e-01,
         8.19884017e-01,  3.48764268e-01, -4.78991409e-01, -9.02348262e-01,
        -5.32073952e-01,  2.91403027e-01,  8.10981738e-01,  5.48963911e-01,
        -2.53752472e-01, -8.59153672e-01, -7.10636616e-01,  5.52527991e-02,
         7.34359377e-01,  7.02315662e-01, -1.14175021e-02, -7.50637136e-01,
        -8.35708117e-01, -1.88416578e-01,  5.96120626e-01,  7.96603607e-01,
         2.28709237e-01, -5.85443019e-01, -8.97325331e-01, -4.20194541e-01,
         4.07277504e-01,  8.24316822e-01,  4.47499387e-01, -3.76730589e-01,
        -8.90579867e-01, -6.21617791e-01,  1.82873147e-01,  7.83247689e-01,
         6.27524249e-01, -1.41125760e-01, -8.16009064e-01, -7.76641066e-01,
        -5.92165222e-02,  6.76667751e-01,  7.54443146e-01,  1.02603324e-01,
        -6.79553190e-01, -8.72915303e-01, -2.99706781e-01,  5.13067105e-01,
         8.18145792e-01,  3.35041343e-01, -4.92082241e-01, -9.02771350e-01,
         8.42557161e-01,  4.70684936e+00,  4.22553393e+00, -1.58869945e-01,
        -4.41536156e+00, -4.63054215e+00, -6.06575682e-01,  3.95692164e+00,
         4.86429142e+00,  1.28130206e+00, -3.49786253e+00, -5.07926048e+00,
        -2.00896181e+00,  2.89021506e+00,  5.11398949e+00,  2.61783353e+00,
        -2.30329853e+00, -5.12494059e+00, -3.25288793e+00,  1.59170285e+00,
         4.95473734e+00,  3.74425713e+00, -9.26827846e-01, -4.76394361e+00,
        -4.23926363e+00,  1.64823752e-01,  4.39922090e+00,  4.57084261e+00,
         5.21900664e-01, -4.02502638e+00, -4.88951476e+00, -1.27675786e+00,
         3.49169229e+00,  5.03174462e+00,  1.92748211e+00, -2.96705060e+00,
        -5.15184271e+00, -2.61820643e+00,  2.30444473e+00,  5.09024800e+00,
         3.17794869e+00, -1.67429402e+00, -5.00535057e+00, -3.75266292e+00,
         9.32053671e-01,  4.74169238e+00,  4.17368895e+00, -2.49736890e-01,
        -4.46170782e+00, -4.58975719e+00, -5.16157003e-01,  4.01384352e+00,
         4.83538278e+00,  1.19314138e+00, -3.56422074e+00, -5.06280679e+00,
        -1.92482366e+00,  2.96468143e+00,  5.11032005e+00,  2.53940195e+00,
        -2.38438263e+00, -5.13412885e+00, -3.18173271e+00,  1.67778177e+00,
         4.97659939e+00,  3.68180245e+00, -1.01617872e+00, -4.79804190e+00,
        -4.18675952e+00,  2.55658218e-01,  4.44487293e+00,  4.52933994e+00,
         4.31400656e-01, -4.08131843e+00, -4.85984421e+00, -1.18840367e+00,
         3.55749769e+00,  5.01450005e+00,  1.84304215e+00, -3.04105224e+00,
        -5.14736925e+00, -2.53937075e+00,  2.38516147e+00,  5.09863521e+00,
         3.10629521e+00, -1.76011031e+00, -5.02643056e+00, -3.68962577e+00,
         1.02125190e+00,  4.77504325e+00,  4.12052982e+00, -3.40531760e-01,
        -4.50666205e+00, -4.54754007e+00, -4.25582756e-01,  4.06950134e+00,
         4.80495264e+00,  1.10460060e+00, -3.62946817e+00, -5.04477268e+00,
    ];
    let tolerance = calc_tolerance(&yf, 2e-4, 1e-5);

    problems.insert(
        name.clone(),
        OdeProblem {
            name,
            t_span,
            y0,
            f,
            yf,
            tolerance,
        },
    );

    problems
}
