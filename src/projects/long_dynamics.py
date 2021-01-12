from control import bode, feedback, nichols, step_response, tf, ss, forced_response, ss2tf, series
from matplotlib import pyplot as plt
from numpy import polymul, log10
from src.common.Pilot import Pilot
gain = 0.2
delay = 0.3

# plane = tf(4.85*polymul([1/0.0098, 1], [1/1.371, 1]), polymul([1/(0.063**2), 2*0.0714/0.063, 1], [1/(4.27**2), 2*0.493/4.27, 1]))
# pilot = Pilot(gain, delay).delay_model()
# actuator = tf(15, [1, 15])
# sys = feedback(series(actuator, plane), pilot, -1)
# mag, phase, omega = bode(sys, Plot=False)
#
# plt.subplot(2, 1, 1)
# plt.plot(omega, 20*log10(mag))
# plt.xscale("log")
# plt.grid()
# plt.subplot(2, 1, 2)
# plt.plot(omega, phase*57.3)
# plt.xscale("log")
# plt.grid()
# plt.show()

sys = ss([[1, 2], [3, 4]], [[2], [1]], [[1, 0], [0, 1]], [[0], [0]])
sys2 = ss2tf(sys)
t_out, y_out = step_response(sys)
plt.plot(t_out, y_out[0, :])
plt.show()
