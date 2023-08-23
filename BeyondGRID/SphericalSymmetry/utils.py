import numpy as np


def collocation_points(N):
    result = np.zeros(N)
    for i in range(N):
        result[i] = np.cos(np.pi * i / (N - 1))
    return -result


def c_n(n, N):
    if n == 0 or n == N - 1:
        return 2
    return 1


def T(x, k):
    return np.cos(k * np.arccos(x))


def physical_matrix(N, x):
    result = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            result[i, j] = T(x[i], j)
    return result


def spectral_matrix(N, x):
    result = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            result[i, j] = 2 * T(x[j], i) / ((N - 1) * c_n(i, N) * c_n(j, N))
    return result


def D_tilde(N):
    result = np.zeros((N, N))
    result[N - 2, N - 1] = 2 * (N - 1) / c_n(N - 2, N)
    for i in np.flip(np.arange(N - 2)):
        result[i, i + 1] = 2 * (i + 1.0) / c_n(i, N)
        result[i, :] += result[i + 2, :] / c_n(i, N)
    return result


def D(P, S, D_tilde, mapper, x):
    return np.diag(mapper.jacobian(x)) @ P @ D_tilde @ S


def D_2(D):
    return np.matmul(D, D)


def interpolate(u, S, x_out):
    u_tilde = S @ u
    result = np.zeros(x_out.size)
    for i in range(u.size):
        result += u_tilde[i] * T(x_out, i)
    return result


def ChebyshevHelper(N, mapper):
    x_unscaled = collocation_points(N)
    P = physical_matrix(N, x_unscaled)
    S = spectral_matrix(N, x_unscaled)
    grid = mapper.rescale_grid(x_unscaled)
    D_til = D_tilde(N)
    first_deriv = D(P, S, D_til, mapper, x_unscaled)
    second_deriv = D_2(first_deriv)
    return P, S, grid, first_deriv, second_deriv
