"""Trim aircraft longitudinally."""
from numpy import array, cos, deg2rad, sin
from scipy.optimize import minimize
from src.projects.f104 import model


def trim_alpha_de_nonlinear(speed, altitude, gamma, n=1, tol=1e-1):
    """trim nonlinear aircraft with angle of attack and elevator."""
    w = array([0])

    def obj(x):
        out = x[2]
        return out

    def alpha_stab(x):
        u = array([0, x[1], 0, x[2]])
        x = array([speed * cos(x[0]), 0, speed * sin(x[0]), 0, x[0] + deg2rad(gamma), 0, 0, 0, 0, 0, 0, altitude])
        x_dot, y = model(x, u, w)
        const = array([x_dot[0], x_dot[7], y[1] - n])
        return const

    lim = ([-5/57.3, 20/57.3], [-30/57.3, 30/57.3], [0, 1])
    x0 = array([10/57.3, -5/57.3, 0.5])
    u_out = minimize(obj, x0, bounds=lim, tol=tol,
                     constraints=({'type': 'eq', 'fun': alpha_stab}),
                     options=({'maxiter': 400}))
    return u_out['x']
