
from cvxpy import *
import scipy as sp
import ecos
import scs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time

sp.random.seed(25)

num_uav = 1
num_ue = 100
hon_size = 5000
ver_size = 200
x_max = hon_size
x_min = 0
y_max = ver_size
y_min = 0


# ########################################################
# ########################################################
# Create a Rectangle patch
plt.figure(figsize=(8, 6))
# rect = plt.patches.Rectangle((0, 0), 5000, 200, linewidth=1, edgecolor='r', facecolor='none')
rect = plt.Rectangle(([0, 0]), 5000, 200, linewidth=3, color='r', fill=False)
ax = plt.gca()
ax.add_patch(rect)

class UE_Loc(object):
    def __init__(self, id, hon_size, ver_size, dis_ref):
        x_ref = sp.multiply(id, dis_ref) + sp.multiply(dis_ref, sp.random.rand(1))
        y_ref = sp.multiply(ver_size, sp.random.rand(1))
        self.tx_x = x_ref
        self.tx_y = y_ref
        plt.scatter(self.tx_x, self.tx_y, s=20, c='red')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))
    def coordinate_ue(self):
        x_coor = self.tx_x
        y_coor = self.tx_y
        return x_coor, y_coor

dis_ref = sp.divide(hon_size, num_ue)
coor_ue = []
x_ue, y_ue = [], []
for p in xrange(num_ue):
    coor_ue.append(UE_Loc(p, hon_size, ver_size, dis_ref))
    x_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[0])
    y_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[1])


H_u = 50
D_k = 400
C_k = 200
V_u = 10

# # ############################################################
# # This code is used to solve the opt problem
# # ############################################################
t0 = time.time()

x_u = Variable(num_ue+1)
y_u = Variable(num_ue+1)
sig_u = Variable(num_ue)

obj_opt = Minimize( sum(sig_u) + sp.multiply(num_ue, sp.divide(D_k, C_k)) )
constraints = [x_min <= x_u, x_u <= x_max, y_min <= y_u, y_u <= y_max, x_u[0] == 0, y_u[0] == 0]

constraints.append(square(vstack([x_u[1:]]) - x_ue[:]) + square(vstack([y_u[1:]]) - y_ue[:]) <= sp.power(H_u, 2))
# x_u[0] == 0, y_u[0] == 0
for i in xrange(num_ue+1):
    if i >= 1:
        # constraints.append(square(x_u[i] - x_ue[i - 1]) + square(y_u[i] - y_ue[i - 1]) <= sp.power(H_u, 2))
        constraints.append(sig_u[i-1] >= norm(hstack([x_u[i] - x_u[i-1], y_u[i] - y_u[i-1]])) / V_u)

prob = Problem(obj_opt, constraints)
# prob.solve(solver=SCS)
prob.solve(solver=ECOS)
# prob.solve(solver=ECOS, verbose=True)
t1 = time.time()
time_sol = (t1 - t0)

x_uav_sol = x_u.value
y_uav_sol = y_u.value

print "optimal value: ", obj_opt.value
print "x_coordinate of UAV: ", x_uav_sol
print "y_coordinate of UAV: ", y_uav_sol
print "Time executive: ", time_sol


# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, time_sol_vec_Mon)
# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, EE_sol_vec_Mon)
plt.show()
