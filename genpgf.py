import math
from matplotlib import pyplot as plt
import numpy as np
from nmsol import *
#general settings
y_capping = 2000
t_left=-5
t_right=1
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