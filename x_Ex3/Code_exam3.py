
from cvxpy import *
import scipy as sp
import ecos
import matplotlib.pyplot as plt
import time

sp.random.seed(17)
# ########################################################
# Create a Rectangle patch
plt.figure(figsize=(8, 6))
# rect = plt.patches.Rectangle((0, 0), 5000, 200, linewidth=1, edgecolor='r', facecolor='none')
rect = plt.Rectangle(([0, 0]), 1000, 4000, linewidth=3, color='r', fill=False)
ax = plt.gca()
ax.add_patch(rect)


num_uav = 4
num_row = 5
d_k = 300    # setup time (sec)
L_k = 1200    # battery duration (sec)
V_nk = 10   # 10m/s
C_n = 2000   # distance of row (m)

#
# # ############################################################
# # This code is used to solve the opt problem
# # ############################################################
t0 = time.time()

u_ue = Variable((num_row, num_uav), boolean=True)

obj1, obj2 = 0, 0
for k in xrange(num_uav):
    obj1 += sum(u_ue[:, k])*(C_n/V_nk)
    obj2 += d_k*sum(u_ue[:, k])

obj_opt = Minimize( obj1 + obj2 )
constraints = []
constraints.append(sum(u_ue) == num_row)
for k in xrange(num_uav):
    constraints.append(sum(u_ue[:, k])*(C_n/V_nk) <= L_k)
    # constraints.append(sum(u_ue[:, k]) <= 1)

for n in xrange(num_row):
    constraints.append(sum(u_ue[n, :]) <= 1)

prob = Problem(obj_opt, constraints)
prob.solve(solver=ECOS_BB)
# prob.solve(solver=ECOS_BB, verbose=True)
t1 = time.time()
time_sol = (t1 - t0)

u_ue_sol = u_ue.value
print "Optimal of flying time: ", obj_opt.value
print "Number of served uav: ", u_ue_sol
print "Time executive: ", time_sol

# plt.show()
