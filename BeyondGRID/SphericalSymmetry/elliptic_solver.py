from .utils import *
import numpy as np
import numpy.linalg as np_lin

def solve_scalar_kovacs(N, b, C, eps1, eps2, eps4, tolerance):
    P, S, grid, first_deriv, second_deriv = ChebyshevHelper(N, 2, b, C)
    sol = np.ones(N)
    previous_guess = 0.8/grid
    while np.sum(np.abs(sol)) > tolerance:
        LHS = np.diag(grid/(2.+grid))@second_deriv+np.diag((2*(3.+grid)/((2.+grid)**2.)))@first_deriv+np.diag((12.*eps2)/(grid**6.))+np.diag((36.*eps4*previous_guess*previous_guess)/(grid**6.))
        RHS = -48.*eps1/(grid**6.)-np.diag(grid/(2.+grid))@second_deriv@previous_guess-np.diag((2*(3.+grid)/((2.+grid)**2.)))@first_deriv@previous_guess-(12.*eps2*previous_guess)/(grid**6.)-(12.*eps4*previous_guess*previous_guess*previous_guess)/(grid**6.)

        LHS[0, :] = first_deriv[0,:]
        RHS[0] = -first_deriv[0,:]@previous_guess


        LHS[-1, :] = b*first_deriv[-1,:]
        LHS[-1, -1] += 1.
        RHS[-1] = -b*first_deriv[-1,:]@previous_guess - previous_guess[-1]

        sol = np.matmul(np_lin.inv(LHS), RHS)
        previous_guess += sol
        print(f"Magnitude of correction: {np.sum(np.abs(sol))}")
    return grid, previous_guess, P, S, first_deriv, second_deriv

def solve_scalar_sbvp(N, a, b, C, eps1, eps2, eps4, tolerance):
    P, S, grid, first_deriv, second_deriv = ChebyshevHelper(N, a, b, C)
    sol = np.ones(N)
    previous_guess = 0.8/grid
    while np.sum(np.abs(sol)) > tolerance:
        LHS = np.diag((1.-2./grid))@second_deriv-np.diag((2*(1.-grid)/(grid**2.)))@first_deriv+np.diag((12.*eps2)/(grid**6.))+np.diag((36.*eps4*previous_guess*previous_guess)/(grid**6.))
        RHS = -48.*eps1/(grid**6.)-np.diag((1.-2./grid))@second_deriv@previous_guess+np.diag((2*(1.-grid)/(grid**2.)))@first_deriv@previous_guess-(12.*eps2*previous_guess)/(grid**6.)-(12.*eps4*previous_guess*previous_guess*previous_guess)/(grid**6.)

        LHS[-1, :] = b*first_deriv[-1,:]
        LHS[-1, -1] += 1.
        RHS[-1] = -b*first_deriv[-1,:]@previous_guess - previous_guess[-1]

        sol = np.matmul(np_lin.inv(LHS), RHS)
        previous_guess += sol
        print(f"Magnitude of correction: {np.sum(np.abs(sol))}")
    return grid, previous_guess, P, S, first_deriv, second_deriv