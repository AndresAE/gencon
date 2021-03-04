from control import append, connect, feedback, forced_response, ss, tf, parallel, tf2ss, series, matlab
from numpy import array, interp, linspace, linalg, real, imag, roots
from matplotlib import pyplot as plt
from src.analysis.criteria import root_locus_alt

a = array([[-1.9311e-2, 8.8157, -3.217e+1, -5.7499e-1],
           [-2.5389e-4, -1.0189, 0, 9.0505e-1],
           [0, 0, 0, -1],
           [2.9465e-12, 8.2225e-1, 0, -1.0774]])
b = array([[1.737e-1], [-2.1499e-3], [0], [-1.7555e-1]])
c = array([[0, 57.2957, 0, 0], [0, 0, 0, 57.2957]])
d = array([[0], [0]])

aircraft = ss(a, b, c, d)

actuator = tf2ss(tf([20.2], [1, 20.2]))
k_alpha = -0.5
k_q = -0.25
k_i = 1.5
k_p = 0.5
lp_filter = tf2ss(tf([[[k_alpha*10], [0]]], [[[1, 10], [1]]]))
ol_airplane = series(actuator, aircraft)
sys = feedback(ol_airplane, lp_filter)
q_sys = tf2ss(tf([[[0], [k_q]]], [[[1], [1]]]))
sas = feedback(sys, q_sys)
root_locus_alt(sas)

prop = tf2ss(tf([k_p/k_i], [1]))
integral = tf2ss(tf([1], [1, 0]))
ff = series(parallel(prop, integral), tf2ss(tf([k_i], [1])))

ol = series(ff, sas)
cas = feedback(ol, tf([[[1], [0]]], [[[1], [1]]]))

t_final = 6
dt = 0.01
t = array([0, 0.99, 1, 20])
u = array([0, 0, 1, 1])
t_out = linspace(0, t_final, int(t_final/dt+1))
u_t = interp(t_out, t, u)

t, y, x = forced_response(sas, t_out, u_t)

plt.plot(t, y[0, :])
plt.show()
