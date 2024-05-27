import math
from matplotlib import pyplot as plt
import numpy as np
from nmsol import *
#general settings
y_capping = 2000
t_left=-1
t_right=10
stepsize = 0.01
# Euler
# euler_1 = NumericalSols(IVP1, stepsize, t_left, t_right, EULER_METHOD,y_capping)
# euler_1.draw()
# plt.show() 

# euler_2 = NumericalSols(IVP2, stepsize, t_left, t_right, EULER_METHOD,y_capping)
# euler_2.draw()
# plt.show()

# Improved Euler
# impeuler_1 = NumericalSols(IVP1, stepsize, t_left, t_right, IMP_EULER_METHOD,y_capping)
# impeuler_1.draw()
# plt.show() 

# impeuler_2 = NumericalSols(IVP2, stepsize, t_left, t_right, IMP_EULER_METHOD,y_capping)
# impeuler_2.draw()
# plt.show()
# Runge-Kutta
runge1 = NumericalSols(IVP1, stepsize, t_left, t_right, RUNGE_4TH_METHOD,y_capping)
runge1.draw()
plt.show()

runge2 = NumericalSols(IVP2, stepsize, t_left, t_right, RUNGE_4TH_METHOD,y_capping)
runge2.draw()
plt.show()
# ABM method
abm1 = NumericalSols(IVP1, stepsize, t_left, t_right, ABM_METHOD,y_capping)
abm1.draw()
plt.show()

abm2 = NumericalSols(IVP2, stepsize, t_left, t_right, ABM_METHOD,y_capping)
abm2.draw()
plt.show()
