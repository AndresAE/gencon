from control import forced_response, ss, ss2tf
from numpy import array, interp, linspace
from matplotlib import pyplot as plt

a = array([[-1.9311e-2, 8.8157, -3.217e+1, 5.4799e-1],
           [-2.5389e-4, -1.0189, 0, 9.0505e-1],
           [0, 0, 0, -1],
           [2.9465e-12, 8.2225e-1, 0, -1.0774]])
b = array([[1.737], [-2.1499e-3], [0], [-1.7555e-1]])
c = array([[0, 57.2957, 0, 0], [0, 0, 0, 57.2957]])
d = array([[0], [0]])

sys = ss(a, b, c, d)
tfs = ss2tf(a, b, c, d)

t_final = 10
dt = 0.01
t = array([0, 0.99, 1, 10])
u = array([0, 0, 1, 1])
t_out = linspace(0, t_final, int(t_final/dt+1))
u_t = interp(t_out, t, u)

t, y, x = forced_response(sys, t_out, u_t)

plt.plot(t, y[1, :])
plt.show()
