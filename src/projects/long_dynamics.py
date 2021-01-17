from control import bode, feedback, nichols, tf, series
from matplotlib import pyplot as plt
from numpy import polymul
from common.Pilot import Pilot
from src.analysis.criteria import nichols
gain = 0.2
delay = 0.1

plane = tf(4.85*polymul([1/0.0098, 1], [1/1.371, 1]), polymul([1/(0.063**2), 2*0.0714/0.063, 1], [1/(4.27**2), 2*0.493/4.27, 1]))
pilot = Pilot(gain, delay).delay_model()
actuator = tf(15, [1, 15])
sys = feedback(series(actuator, plane), pilot, -0.1)
mag, phase, omega = bode(sys, plot=False)
nichols(sys)

# plt.subplot(2, 1, 1)
# plt.plot(omega, 20*log10(mag))
# plt.xscale("log")
# plt.grid()
# plt.subplot(2, 1, 2)
# plt.plot(omega, phase*57.3)
# plt.xscale("log")
# plt.grid()
# plt.show()
