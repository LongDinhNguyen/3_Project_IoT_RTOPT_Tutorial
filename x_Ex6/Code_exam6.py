
from cvxpy import *
import scipy as sp
import scipy.io as spio
import ecos
import matplotlib.pyplot as plt
import time

# sp.random.seed(2035)
coverage_r = 1000 #m

num_uav = 1
num_d2d = 10

hUAV = 100
eta = 0.5  # EH efficiency
P_UAV = 5000
Pc_UAV = 4000
rmin = 0.4*log(2)
tau_ini = Parameter(value=0.5)


# mat = spio.loadmat('XModel_Ex6_4D2D.mat', squeeze_me=True)
# mat = spio.loadmat('XModel_Ex6_6D2D.mat', squeeze_me=True)
# mat = spio.loadmat('XModel_Ex6_8D2D.mat', squeeze_me=True)
mat = spio.loadmat('XModel_Ex6_10D2D.mat', squeeze_me=True)
h_d2d = mat['h_d2d']
h_uav_d2d = mat['h_uav_d2d']
# ## the size of h_d2d is (num_d2d, num_d2d)
# ## the size of h_uav_d2d is (1,num_d2d)

max_d2d_gains_diff = sp.copy(h_d2d[:, :])
sp.fill_diagonal(max_d2d_gains_diff, 0)
d2d_to_d2d_gains_diff = max_d2d_gains_diff[:num_d2d, :num_d2d]
d2d_to_d2d_gains_diag = sp.subtract(h_d2d, d2d_to_d2d_gains_diff)


# # ############################################################
# # This code is used to solve the opt problem
# # ############################################################
t0 = time.time()

pow_co = Variable(num_d2d)
objective = Minimize( sum(pow_co)*sp.subtract(1, tau_ini.value) )

constraints = []
constraints.append(d2d_to_d2d_gains_diag * pow_co >= (exp(rmin / (1-tau_ini.value)) - 1) * (d2d_to_d2d_gains_diff*pow_co + 1))
constraints.append(pow_co*(1-tau_ini.value) <= tau_ini.value * eta * P_UAV * h_uav_d2d)

prob = Problem(objective, constraints)
prob.solve(solver=ECOS, verbose=True)
# print prob.status
# prob.solve(solver=scs)

t1 = time.time()
time_sol = (t1 - t0)

pow_sol = pow_co.value

print "opt_power of UAV: ", sp.sum(pow_sol)
print "Time executive: ", time_sol

# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, time_sol_vec_Mon)
# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, EE_sol_vec_Mon)
# plt.show()
