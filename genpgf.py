import math
from matplotlib import pyplot as plt
import numpy as np
from nmsol import *
#general settings
y_capping = 30
t_left=-1.2
t_right=1.2
stepsize = 0.01
# Euler
euler_1 = NumericalSols(IVP1, stepsize, t_left, t_right, EULER_METHOD,y_capping)
euler_1.draw("red","IVP1")
plt.legend()
# plt.savefig('euler1.pgf')
# plt.cla()
euler_2 = NumericalSols(IVP2, stepsize, t_left, t_right, EULER_METHOD,100)
euler_2.draw("blue","IVP2")
plt.legend()
plt.savefig('euler.pgf')
plt.cla()
# Improved Euler
impeuler_1 = NumericalSols(IVP1, stepsize, t_left, t_right, IMP_EULER_METHOD,y_capping)
impeuler_1.draw("red","IVP1")
plt.legend()
# plt.savefig('impeuler1.pgf')
# plt.cla()
impeuler_2 = NumericalSols(IVP2, stepsize, t_left, t_right, IMP_EULER_METHOD,y_capping)
impeuler_2.draw("blue","IVP2")
plt.legend()
plt.savefig('impeuler.pgf')
plt.cla()
# Runge-Kutta
runge1 = NumericalSols(IVP1, stepsize, t_left, t_right, RUNGE_4TH_METHOD,100)
runge1.draw("red","IVP1")
plt.legend()
# plt.savefig('runge1.pgf')
# plt.cla()
runge2 = NumericalSols(IVP2, stepsize, t_left, t_right, RUNGE_4TH_METHOD,y_capping)
runge2.draw("blue","IVP2")
plt.legend()
plt.savefig('runge.pgf')
plt.cla()
# ABM method
abm1 = NumericalSols(IVP1, stepsize, t_left, t_right, ABM_METHOD,y_capping)
abm1.draw("red","IVP1")
plt.legend()
# plt.savefig('abm1.pgf')
# plt.cla()
abm2 = NumericalSols(IVP2, stepsize, t_left, t_right, ABM_METHOD,100)
abm2.draw("blue","IVP2")
plt.legend()
plt.savefig('abm.pgf')
plt.cla()
# Power Series
power1 = NumericalSols(IVP1, stepsize, t_left, t_right, POWER_SERIES_METHOD,100)
power1.draw("red","IVP1")
plt.legend()
# plt.savefig('power1.pgf')
# plt.cla()
power2 = NumericalSols(IVP2, stepsize, t_left, t_right, POWER_SERIES_METHOD,100)
power2.draw("blue","IVP2")
plt.legend()
plt.savefig('power.pgf')
plt.show()
plt.cla()

# overall IVP1
plotDiffMethods(IVP1, stepsize, t_left, t_right, y_capping)
plt.title("Overall IVP1")
plt.savefig('overallIVP1.pgf')
plt.show()
plt.cla()
# overall IVP2
plotDiffMethods(IVP2, stepsize, t_left, t_right, y_capping)
plt.title("Overall IVP2")
plt.savefig('overallIVP2.pgf')
plt.show()
plt.cla()