"""Trim aircraft longitudinally."""
from numpy import abs, append, array, cos, deg2rad, isnan, nan, sin, sum, zeros
from scipy.optimize import minimize
from src.projects.F104.f104 import model, outputs


# def trim(model, outputs, x_dot_0, x_0, u_0, y_0, x_var_i, u_var_i, x_lim, u_lim, tol=1e-1, maxiter=500):
#     """nonlinear trim."""
#     x_dot_fixed = isnan(x_dot_0) != nan
#     nanned_x = isnan(x_0)
#     not_nanned_x = [not elem for elem in nanned_x]
#     x_fixed = not_nanned_x and x_var_i != 1
#     nanned_u = isnan(u_0)
#     not_nanned_u = [not elem for elem in nanned_u]
#     u_fixed = not_nanned_u and u_var_i != 1
#     y_fixed = y_0 != nan
#
#     x_iter = x_0
#     u_iter = u_0
#
#     def obj(xi):
#         return abs(sum(xi))
#
#     def constraint(xi):
#         x_iter[x_var_i] = xi[sum(u_var_i)-1:-1]
#         u_iter[u_var_i] = xi[0:sum(u_var_i)]
#         x_iter[isnan(x_iter)] = 0
#         u_iter[isnan(u_iter)] = 0
#         x_dot_iter = model(x_iter, u_iter)
#         y_iter = outputs(x_dot_iter, x_iter, u_iter)
#         delta_x_dot = x_dot_0[x_dot_fixed] - x_dot_iter[x_dot_fixed]
#         delta_x = x_0[x_fixed] - x_iter[x_fixed]
#         delta_u = u_0[u_fixed] - u_iter[u_fixed]
#         delta_y = y_0[y_fixed] - y_iter[y_fixed]
#         c = append(append(delta_x_dot, delta_x), append(delta_u, delta_y))
#         return c
#
#     lim = u_lim + x_lim
#     x0 = append(u_0[u_var_i], x_0[x_var_i])
#     u_out = minimize(obj, x0, bounds=lim, tol=tol,
#                      constraints=({'type': 'eq', 'fun': constraint}),
#                      options=({'maxiter': maxiter}))
#
#     u_out = u_out['x'][0:len(u_var_i)]
#     x_out = u_out['x'][len(u_var_i):-1]
#
#     x_0[x_var_i] = u_out['x'][len(u_var_i):-1]
#     u_0[u_var_i] = u_out['x'][0:len(u_var_i)]
#
#     x_dot_out = model(x_0, u_0)
#     y_out = outputs(x_dot_out, x_0, u_0)
#     return x_dot_out, x_out, u_out, y_out


def trim_alpha_de_nonlinear(speed, altitude, gamma, n=1, tol=1e-1):
    """trim nonlinear aircraft with angle of attack and elevator."""

    def obj(x):
        out = x[2]
        return out

    def alpha_stab(x):
        u = array([0, x[1], 0, x[2]])
        x = array([speed * cos(x[0]), 0, speed * sin(x[0]), 0, x[0] + deg2rad(gamma), 0, 0, 0, 0, 0, 0, altitude])
        x_dot = model(x, u)
        y = outputs(x_dot, x, u)
        const = array([x_dot[0], x_dot[7], y[1] - n])
        return const

    lim = ([-5/57.3, 20/57.3], [-30/57.3, 30/57.3], [0, 1])
    x0 = array([10/57.3, -5/57.3, 0.5])
    u_out = minimize(obj, x0, bounds=lim, tol=tol,
                     constraints=({'type': 'eq', 'fun': alpha_stab}),
                     options=({'maxiter': 400}))
    return u_out['x']
