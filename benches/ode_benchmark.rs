use criterion::{black_box, criterion_group, criterion_main, Criterion};

use finfoot::ode::dopri5;
use test_util::all_problems;

macro_rules! generate_ode_benchmarks {
    ( $( $name:ident ),* $(,)?) => {
        $(
            fn $name(c: &mut Criterion) {
                const CONFIG: dopri5::Config = dopri5::Config {
                    rel_tol: 1e-4,
                    abs_tol: 1e-6,
                };

                let problem = &all_problems()[stringify!($name)];
                let input = dopri5::Input {
                    t_span: problem.t_span,
                    y0: &problem.y0,
                    h0: (problem.t_span[1] - problem.t_span[0]) / 100.0,
                    f: &problem.f,
                };
                c.bench_function(stringify!($name), |b| {
                    b.iter(|| dopri5::integrate(black_box(&input), &CONFIG))
                });
            }
        )*

        criterion_group!(
            ode_benches,
            $( $name ),*
        );
    };
}

generate_ode_benchmarks! {
    van_der_pol_oscillator,
    robertson_equations,
    coupled_oscillators,
}

criterion_main!(ode_benches);
