from .utils import *
import numpy as np
import numpy.linalg as np_lin
import sys

# For the following functions, one must define a mapping class to pass to the solver.
# Seeing mappings.py for details.


def solve_scalar_kovacs(N, mapper, eps1, eps2, eps4, tolerance):
    if mapper.a != 2.0:
        sys.exit("Inner boundary must be at horizon")
    P, S, grid, first_deriv, second_deriv = ChebyshevHelper(N, mapper)
    sol = np.ones(N)
    previous_guess = 0.8 / grid
    while np.sum(np.abs(sol)) > tolerance:
        LHS = (
            np.diag(grid / (2.0 + grid)) @ second_deriv
            + np.diag((2 * (3.0 + grid) / ((2.0 + grid) ** 2.0))) @ first_deriv
            + np.diag((12.0 * eps2) / (grid**6.0))
            + np.diag((36.0 * eps4 * previous_guess * previous_guess) / (grid**6.0))
        )
        RHS = (
            -48.0 * eps1 / (grid**6.0)
            - np.diag(grid / (2.0 + grid)) @ second_deriv @ previous_guess
            - np.diag((2 * (3.0 + grid) / ((2.0 + grid) ** 2.0)))
            @ first_deriv
            @ previous_guess
            - (12.0 * eps2 * previous_guess) / (grid**6.0)
            - (12.0 * eps4 * previous_guess * previous_guess * previous_guess)
            / (grid**6.0)
        )

        LHS[0, :] = first_deriv[0, :]
        RHS[0] = -first_deriv[0, :] @ previous_guess

        LHS[-1, :] = mapper.b * first_deriv[-1, :]
        LHS[-1, -1] += 1.0
        RHS[-1] = -mapper.b * first_deriv[-1, :] @ previous_guess - previous_guess[-1]

        sol = np.matmul(np_lin.inv(LHS), RHS)
        previous_guess += sol
        print(f"Magnitude of correction: {np.sum(np.abs(sol))}")
    return grid, previous_guess, P, S, first_deriv, second_deriv


def solve_scalar_sbvp(N, mapper, eps1, eps2, eps4, tolerance):
    if mapper.a > 2.0:
        sys.exit("Inner boundary must be within horizon.")
    P, S, grid, first_deriv, second_deriv = ChebyshevHelper(N, mapper)
    sol = np.ones(N)
    previous_guess = 0.8 / grid
    while np.sum(np.abs(sol)) > tolerance:
        LHS = (
            np.diag((1.0 - 2.0 / grid)) @ second_deriv
            - np.diag((2 * (1.0 - grid) / (grid**2.0))) @ first_deriv
            + np.diag((12.0 * eps2) / (grid**6.0))
            + np.diag((36.0 * eps4 * previous_guess * previous_guess) / (grid**6.0))
        )
        RHS = (
            -48.0 * eps1 / (grid**6.0)
            - np.diag((1.0 - 2.0 / grid)) @ second_deriv @ previous_guess
            + np.diag((2 * (1.0 - grid) / (grid**2.0))) @ first_deriv @ previous_guess
            - (12.0 * eps2 * previous_guess) / (grid**6.0)
            - (12.0 * eps4 * previous_guess * previous_guess * previous_guess)
            / (grid**6.0)
        )

        LHS[-1, :] = mapper.b * first_deriv[-1, :]
        LHS[-1, -1] += 1.0
        RHS[-1] = -mapper.b * first_deriv[-1, :] @ previous_guess - previous_guess[-1]

        sol = np.matmul(np_lin.inv(LHS), RHS)
        previous_guess += sol
        print(f"Magnitude of correction: {np.sum(np.abs(sol))}")
    return grid, previous_guess, P, S, first_deriv, second_deriv
