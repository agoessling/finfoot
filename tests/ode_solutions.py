import dataclasses
from typing import Callable

import numpy as np
import scipy as sp


@dataclasses.dataclass(frozen=True)
class OdeProblem:
    name: str
    t_span: tuple[float, float]
    y0: list[float]
    f: Callable[[float, list[float]], list[float]]


def all_problems() -> list[OdeProblem]:
    problems: list[OdeProblem] = []

    # Exponential decay.
    tau = 1.0
    f = lambda _, y: [-tau * y[0]]
    problems.append(
        OdeProblem(name='exponential', t_span=(0.0, 1.0), y0=[1.0], f=f))

    # Harmonic oscillator.
    omega = 2.0 * np.pi
    f = lambda _, y: [y[1], -omega**2 * y[0]]
    problems.append(
        OdeProblem(name='harmonic oscillator',
                   t_span=(0.0, 1.0),
                   y0=[1.0, 0.0],
                   f=f))

    # Van der Pol oscillator.
    mu = 5.0
    f = lambda _, y: [y[1], mu * (1.0 - y[0]**2) * y[1] - y[0]]
    problems.append(
        OdeProblem(name='van der pol oscillator',
                   t_span=(0.0, 15.0),
                   y0=[1.0, 0.0],
                   f=f))

    # Lorentz attractor.
    sigma = 10.0
    rho = 28.0
    beta = 8.0 / 3.0
    f = lambda _, y: [sigma * (y[1] - y[0]), y[0] * (rho - y[2]) - y[1], y[0] * y[1] - beta * y[2]]
    problems.append(
        OdeProblem(name='lorentz attractor',
                   t_span=(0.0, 5.0),
                   y0=[1.0, 1.0, 1.0],
                   f=f))

    # Robertson equations.
    f = lambda _, y: [
        -0.04 * y[0] + 1e4 * y[1] * y[2],
        0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1]**2,
        3e7 * y[1]**2,
    ];
    problems.append(
        OdeProblem(name='robertson equations',
                   t_span=(0.0, 30.0),
                   y0=[1.0, 0.0, 0.0],
                   f=f))

    return problems


def main() -> None:
    for problem in all_problems():
        solution = sp.integrate.solve_ivp(problem.f,
                                          problem.t_span,
                                          problem.y0,
                                          rtol=1e-6,
                                          atol=1e-6,
                                          method='RK45')

        assert solution.success
        assert solution.t[-1] == problem.t_span[1]

        print(problem)
        print(f'yf: {solution.y[:, -1]}')
        print('')


if __name__ == '__main__':
    main()
