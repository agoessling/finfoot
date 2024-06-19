use nalgebra::DVector;

use super::{Error, InputError};

#[derive(Debug)]
pub struct Input<'a, F>
where
    F: Fn(f64, &DVector<f64>) -> DVector<f64>,
{
    pub t_span: [f64; 2],
    pub y0: &'a DVector<f64>,
    pub h0: f64,
    pub f: &'a F,
}

#[derive(Debug)]
pub struct Config {
    pub rel_tol: f64,
    pub abs_tol: f64,
}

#[derive(Debug)]
pub struct Output {
    pub y: DVector<f64>,
    pub h: f64,
    pub num_calls: usize,
}

pub fn integrate<F>(input: &Input<F>, config: &Config) -> Result<Output, Error>
where
    F: Fn(f64, &DVector<f64>) -> DVector<f64>,
{
    validate_input(input)?;

    let mut t = input.t_span[0];
    let mut y = input.y0.clone();
    let mut h = input.h0;
    let mut k1: Option<DVector<f64>> = None;
    let mut num_calls = 0;
    let mut num_failures = 0;

    loop {
        h = h.min(input.t_span[1] - t);
        let t_next = t + h;

        let step_output = dopri5_step(t, &y, input.f, h, &k1);
        num_calls += step_output.num_calls;

        // h step size control.
        let error = step_output.error.abs();
        let allowed_error = (config.rel_tol * step_output.y.abs()).map(|x| x.max(config.abs_tol));

        const MIN_ERROR_RATIO: f64 = 1e-5; // (1/10)^5, 10x decrease in h.
        const MAX_ERROR_RATIO: f64 = 1e5; // 10^5, 10x increase in h.

        let error_ratio = allowed_error
            .zip_map(&error, |a, b| {
                (a / b).clamp(MIN_ERROR_RATIO, MAX_ERROR_RATIO)
            })
            .amin();

        h = 0.9 * h * error_ratio.powf(1.0 / 5.0);

        // Discard step if error is too high.
        if error_ratio < 1.0 {
            num_failures += 1;
            if num_failures > 10 {
                return Err(Error::Convergence);
            }

            continue;
        }
        num_failures = 0;

        // Propagate state.
        t = t_next;
        y = step_output.y;
        k1 = Some(step_output.k7); // First same as last property. [FSAL]

        // Terminate integration.
        if t >= input.t_span[1] {
            break;
        }
    }

    Ok(Output { y, h, num_calls })
}

fn validate_input<F>(input: &Input<F>) -> Result<(), InputError>
where
    F: Fn(f64, &DVector<f64>) -> DVector<f64>,
{
    if input.t_span[0] > input.t_span[1] {
        return Err(InputError::TimeSpan);
    }
    if input.h0 <= 0.0 {
        return Err(InputError::StepSize);
    }
    Ok(())
}

struct StepOutput {
    y: DVector<f64>,
    error: DVector<f64>,
    k7: DVector<f64>,
    num_calls: usize,
}

fn dopri5_step<F>(t: f64, y: &DVector<f64>, f: &F, h: f64, k1: &Option<DVector<f64>>) -> StepOutput
where
    F: Fn(f64, &DVector<f64>) -> DVector<f64>,
{
    const C_COEFF: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    const A_COEFF: [[f64; 6]; 6] = [
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            0.0,
            0.0,
        ],
        [
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
            0.0,
        ],
        [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const B_COEFF: [[f64; 7]; 2] = [
        [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
            0.0,
        ],
        [
            5179.0 / 57600.0,
            0.0,
            7571.0 / 16695.0,
            393.0 / 640.0,
            -92097.0 / 339200.0,
            187.0 / 2100.0,
            1.0 / 40.0,
        ],
    ];

    let num_calls = match k1 {
        Some(_) => 6,
        None => 7,
    };

    // Lazily initialize k1 if it is not provided.
    let default_k1;
    let k1 = match k1 {
        Some(k1) => k1,
        None => {
            default_k1 = f(t, y);
            &default_k1
        }
    };

    let k2 = f(t + C_COEFF[1] * h, &(y + (h * A_COEFF[0][0]) * k1));

    let k3 = f(
        t + C_COEFF[2] * h,
        &(y + (h * A_COEFF[1][0]) * k1 + (h * A_COEFF[1][1]) * &k2),
    );

    let k4 = f(
        t + C_COEFF[3] * h,
        &(y + (h * A_COEFF[2][0]) * k1 + (h * A_COEFF[2][1]) * &k2 + (h * A_COEFF[2][2]) * &k3),
    );

    let k5 = f(
        t + C_COEFF[4] * h,
        &(y + (h * A_COEFF[3][0]) * k1
            + (h * A_COEFF[3][1]) * &k2
            + (h * A_COEFF[3][2]) * &k3
            + (h * A_COEFF[3][3]) * &k4),
    );

    let k6 = f(
        t + C_COEFF[5] * h,
        &(y + (h * A_COEFF[4][0]) * k1
            + (h * A_COEFF[4][1]) * &k2
            + (h * A_COEFF[4][2]) * &k3
            + (h * A_COEFF[4][3]) * &k4
            + (h * A_COEFF[4][4]) * &k5),
    );

    // With Dormand Prince, a_7 == b_1
    // Purposefully skip a_72 as it is zero.
    let fifth_order = y
        + (h * A_COEFF[5][0]) * k1
        + (h * A_COEFF[5][2]) * &k3
        + (h * A_COEFF[5][3]) * &k4
        + (h * A_COEFF[5][4]) * &k5
        + (h * A_COEFF[5][5]) * &k6;

    let k7 = f(t + C_COEFF[6] * h, &fifth_order);

    let fourth_order = y
        + (h * B_COEFF[1][0]) * k1
        + (h * B_COEFF[1][2]) * &k3
        + (h * B_COEFF[1][3]) * &k4
        + (h * B_COEFF[1][4]) * &k5
        + (h * B_COEFF[1][5]) * &k6
        + (h * B_COEFF[1][6]) * &k7;

    let error = &fifth_order - fourth_order;

    StepOutput {
        y: fifth_order,
        error,
        k7,
        num_calls,
    }
}
