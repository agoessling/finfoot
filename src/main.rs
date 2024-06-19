use nalgebra::DVector;

use finfoot::ode::dopri5;

fn main() {
    let van_der_pol = |_: f64, y: &DVector<f64>| {
        let x = y[0];
        let dx = y[1];
        DVector::from_vec(vec![dx, 5.0 * (1.0 - x * x) * dx - x])
    };

    const TIMESTEP: f64 = 30.0 - 0.1;
    let mut t: f64 = 0.0;
    let mut y = DVector::from_vec(vec![1.0, -0.5]);
    let mut h = TIMESTEP;
    let mut num_calls = 0;

    while t < 30.0 {
        println!("{:.6e}, {:.6e}, {:.6e}, {}", t, y[0], y[1], num_calls);

        let input = dopri5::Input {
            t_span: [t, t + TIMESTEP],
            y0: &y,
            h0: h,
            f: &van_der_pol,
        };

        const CONFIG: dopri5::Config = dopri5::Config {
            rel_tol: 1e-6,
            abs_tol: 1e-6,
        };

        let output = dopri5::integrate(&input, &CONFIG).unwrap();

        t += TIMESTEP;
        y = output.y;
        h = output.h;
        num_calls = output.num_calls;
    }
}
