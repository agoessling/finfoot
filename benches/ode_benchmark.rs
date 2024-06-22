use criterion::{black_box, criterion_group, criterion_main, Criterion};

use finfoot::ode::dopri5;
use test_util::all_problems;

pub fn van_der_pol_oscillator(c: &mut Criterion) {
    const CONFIG: dopri5::Config = dopri5::Config {
        rel_tol: 1e-4,
        abs_tol: 1e-6,
    };

    let problem = &all_problems()["van_der_pol_oscillator"];
    let input = dopri5::Input {
        t_span: problem.t_span,
        y0: &problem.y0,
        h0: (problem.t_span[1] - problem.t_span[0]) / 100.0,
        f: &problem.f,
    };
    c.bench_function("van der pol oscillator", |b| {
        b.iter(|| dopri5::integrate(black_box(&input), &CONFIG))
    });
}

pub fn robertson_equations(c: &mut Criterion) {
    const CONFIG: dopri5::Config = dopri5::Config {
        rel_tol: 1e-4,
        abs_tol: 1e-6,
    };

    let problem = &all_problems()["robertson_equations"];
    let input = dopri5::Input {
        t_span: problem.t_span,
        y0: &problem.y0,
        h0: (problem.t_span[1] - problem.t_span[0]) / 100.0,
        f: &problem.f,
    };
    c.bench_function("robertson equations", |b| {
        b.iter(|| dopri5::integrate(black_box(&input), &CONFIG))
    });
}

pub fn coupled_oscillators(c: &mut Criterion) {
    const CONFIG: dopri5::Config = dopri5::Config {
        rel_tol: 1e-4,
        abs_tol: 1e-6,
    };

    let problem = &all_problems()["coupled_oscillators"];
    let input = dopri5::Input {
        t_span: problem.t_span,
        y0: &problem.y0,
        h0: (problem.t_span[1] - problem.t_span[0]) / 100.0,
        f: &problem.f,
    };
    c.bench_function("coupled oscillators", |b| {
        b.iter(|| dopri5::integrate(black_box(&input), &CONFIG))
    });
}

criterion_group!(
    ode_benches,
    van_der_pol_oscillator,
    robertson_equations,
    coupled_oscillators
);
criterion_main!(ode_benches);
