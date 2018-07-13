
from cvxpy import *
import scipy as sp
import ecos
import matplotlib.pyplot as plt
import time

sp.random.seed(35)

coverage_r = 500 #m
plt.figure(figsize=(8, 6))
coverage = plt.Circle((0, 0), coverage_r, color='green', fill=False)
ax = plt.gca()
ax.add_patch(coverage)

class MBS(object):
    def __init__(self, HBS, XBS, YBS):
        self.h = HBS
        plt.scatter(XBS, YBS, s=100, c='black', marker='D')

class block(object):
    def __init__(self, h):
        self.h = h
        plt.scatter(0, 0, s=10000, c='none', marker='^', edgecolors='black')

class UE_Loc(object):
    def __init__(self, id, coverage_r, origin_x=300, origin_y=0, low_tx=0.2):
        tx_d = sp.multiply(200, sp.random.uniform(low=low_tx))
        tx_angle = sp.multiply(2, sp.multiply(sp.pi, sp.random.rand(1,1)))
        self.tx_x = sp.add(origin_x, sp.multiply(tx_d, sp.cos(tx_angle)))
        self.tx_y = sp.add(origin_y, sp.multiply(tx_d, sp.sin(tx_angle)))
        plt.scatter(self.tx_x, self.tx_y, s=20, c='red')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))
    def coordinate_ue(self):
        x_coor = self.tx_x
        y_coor = self.tx_y
        return x_coor, y_coor


num_uav = 1
num_ue = 30
HBS = 30
XBS = -300
YBS = 0
MBS(HBS, XBS, YBS)
block(30)

coor_ue = []
x_ue, y_ue = [], []
for p in xrange(num_ue):
    coor_ue.append(UE_Loc(p, coverage_r, low_tx=0.2))
    x_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[0])
    y_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[1])

# print x_ue, y_ue

h_min = 47
h_max = 122
x_max = coverage_r
x_min = -coverage_r
y_max = coverage_r
y_min = -coverage_r

# # ############################################################
# # This code is used to solve the opt problem
# # ############################################################
t0 = time.time()

x_u = Variable(1)
y_u = Variable(1)
h_u = Variable(1)

R_k2 = []
for i in xrange(num_ue):
    R_k2.append(square(x_u - x_ue[i]) + square(y_u - y_ue[i]) + square(h_u))

# R_k2 = sum(square(norm(x_u - x_ue))) + sum(square(norm(y_u - y_ue))) + square(h_u)
R_BU = square(XBS - x_u) + square(YBS - y_u) + square(HBS - h_u)

obj_opt = Minimize( sum(R_k2) + R_BU )
constraints = [h_min <= h_u, h_u <= h_max, x_min <= x_u, x_u <= x_max, y_min <= y_u, y_u <= y_max]

prob = Problem(obj_opt, constraints)
prob.solve(solver=ECOS)
# prob.solve(solver=ECOS_BB, verbose=True)
t1 = time.time()
time_sol = (t1 - t0)

x_uav_sol = x_u.value
y_uav_sol = y_u.value
h_uav_sol = h_u.value

print "x_coordinate of UAV: ", x_uav_sol
print "y_coordinate of UAV: ", y_uav_sol
print "height of UAV: ", h_uav_sol
print "Time executive: ", time_sol

for p in xrange(num_uav):
    uav_x = x_u.value
    uav_y = y_u.value
    plt.scatter(uav_x, uav_y, s=50, c='black')
    plt.annotate(1, (uav_x + 20, uav_y + 20))

    conver = plt.Circle((uav_x, uav_y), h_max, color='green', fill=False)
    ax = plt.gca()
    ax.add_patch(conver)

# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, time_sol_vec_Mon)
# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, EE_sol_vec_Mon)
plt.show()
