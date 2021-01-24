from control import feedback, tf, ss, step_response
from matplotlib import pyplot as plt

num_u = [-319.59, 26420.82, 49830.26]
num_alpha = [-18.3, -1822.4, -69.76, -103.02]
num_theta = [-1819.74, -1641.68, -143.55]
num_q = [-1819.74, -1641.68, -143.55, 0]
den = [190.635, 750.53, 1434.27, 81.27, 58.87]

long = tf([[num_u], [num_alpha], [num_theta], [num_q]], [[den], [den], [den], [den]])
sys = ss(long)
claw = ss([], [], [], [0, 0, 0, -1])
sys3 = feedback(sys, claw)

t, yout = step_response(sys3, 10)

plt.plot(t, yout[1, :])
plt.show()