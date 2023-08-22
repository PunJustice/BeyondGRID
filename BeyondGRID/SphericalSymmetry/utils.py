import numpy as np

def collocation_points(N):
    result = np.zeros(N)
    for i in range(N):
        result[i] = np.cos(np.pi*i/(N-1))
    return -result

def rescale_grid(x, a, b, C):
    return (B(a,b,C)/(x-A(a,b,C)))-C

def unscale_grid(x, a, b, C):
    return A(a, b, C) + B(a, b, C)/(x+C)

def c_n(n, N):
    if n==0 or n == N-1:
        return 2
    return 1

def T(x, k):
    return np.cos(k*np.arccos(x))

def physical_matrix(N, x):
    result = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            result[i, j] = T(x[i], j)
    return result

def spectral_matrix(N, x):
    result = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            result[i, j] = 2*T(x[j], i)/((N-1)*c_n(i, N)*c_n(j,N))
    return result

def D_tilde(N):
    result = np.zeros((N,N))
    result[N-2, N-1] = 2*(N-1)/c_n(N-2,N)
    for i in np.flip(np.arange(N-2)):
        result[i, i+1] = 2*(i+1.)/c_n(i,N)
        result[i,:] += result[i+2,:]/c_n(i,N)
    return result

def D(P, S, D_tilde, A, B, C, x):
    return np.diag((-B/((x+C)**2.)))@P@D_tilde@S

def D_2(D):
    return np.matmul(D,D)

def A(a, b, C):
    return (a+b+2.*C)/(b-a)

def B(a, b, C):
    return b-b*A(a, b, C)+C-C*A(a, b, C)

def interpolate(u, S, x_out):
    u_tilde = S@u
    result = np.zeros(x_out.size)
    for i in np.range(u.size):
        result += u_tilde[i]*T(x_out, i)
    return result

def ChebyshevHelper(N, a, b, C):
    x_unscaled = collocation_points(N)
    P = physical_matrix(N, x_unscaled)
    S = spectral_matrix(N, x_unscaled)
    grid = rescale_grid(x_unscaled, a, b, C)
    D_til = D_tilde(N)
    first_deriv = D(P, S, D_til, A(a,b,C), B(a,b,C), C, grid)
    second_deriv = D_2(first_deriv)
    return P, S, grid, first_deriv, second_deriv