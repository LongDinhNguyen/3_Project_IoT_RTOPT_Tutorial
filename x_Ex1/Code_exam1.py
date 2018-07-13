
from cvxpy import *
import scipy as sp
import ecos
import matplotlib.pyplot as plt
import time

# sp.random.seed(17)
# sp.random.seed(35)
sp.random.seed(55)
coverage_r = 300 #m
plt.figure(figsize=(8, 6))
coverage = plt.Circle((0, 0), coverage_r, color='green', fill=False)
ax = plt.gca()
ax.add_patch(coverage)

class UAV(object):
    def __init__(self, h):
        self.h = h
        plt.scatter(0, 0, s=200, c='black', marker='D')

class UE_Loc(object):
    def __init__(self, id, coverage_r, origin_x=0, origin_y=0, low_tx=0.1):
        tx_d = sp.multiply(coverage_r, sp.random.uniform(low=low_tx))
        tx_angle = sp.multiply(2, sp.multiply(sp.pi, sp.subtract(sp.random.rand(), 1)))
        self.tx_x = sp.add(origin_x, sp.multiply(tx_d, sp.sin(tx_angle)))
        self.tx_y = sp.add(origin_y, sp.multiply(tx_d, sp.cos(tx_angle)))
        plt.scatter(self.tx_x, self.tx_y, s=20, c='red')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))
    def coordinate_ue(self):
        x_coor = self.tx_x
        y_coor = self.tx_y
        return x_coor, y_coor

class UE_Loc_spec(object):
    def __init__(self, id, coverage_r, origin_x=0, origin_y=0, low_tx=0.1):
        tx_d = sp.multiply(coverage_r, sp.random.uniform(low=low_tx))
        tx_angle = sp.multiply(0.5, sp.multiply(sp.pi, (sp.random.rand())))
        self.tx_x = sp.add(origin_x, sp.multiply(tx_d, sp.sin(tx_angle)))
        self.tx_y = sp.add(origin_y, sp.multiply(tx_d, sp.cos(tx_angle)))
        plt.scatter(self.tx_x, self.tx_y, s=20, c='red')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))
    def coordinate_ue(self):
        x_coor = self.tx_x
        y_coor = self.tx_y
        return x_coor, y_coor


num_uav = 1
num_ue = 30

coor_ue = []
x_ue, y_ue = [], []
for p in range(0, num_ue/2+1):
    coor_ue.append(UE_Loc_spec(p, coverage_r, low_tx=0.1))
    x_ue.append(UE_Loc_spec.coordinate_ue(coor_ue[p])[0])
    y_ue.append(UE_Loc_spec.coordinate_ue(coor_ue[p])[1])

for p in range(num_ue/2+1, num_ue):
    coor_ue.append(UE_Loc(p, coverage_r, low_tx=0.1))
    x_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[0])
    y_ue.append(UE_Loc.coordinate_ue(coor_ue[p])[1])

# print x_ue, y_ue

h_min = 47
h_max = 122
x_max = coverage_r
x_min = -coverage_r
y_max = coverage_r
y_min = -coverage_r
#
# # ############################################################
# # This code is used to solve the opt problem
# # ############################################################
t0 = time.time()

x_u = Variable(1)
y_u = Variable(1)
h_u = Variable(1)
u_ue = Variable(num_ue, boolean=True)

obj_opt = Maximize( sum(u_ue) )
constraints = [h_min <= h_u, h_u <= h_max, x_min <= x_u, x_u <= x_max, y_min <= y_u, y_u <= y_max]

# for i in xrange(num_ue):
    # constraints.append(square(x_u - x_ue[i]) + square(y_u - y_ue[i]) + square(h_u) \
    # <= 1e5 + sp.multiply(sqrt(sp.power(1000,2) + sp.power(1000,2) + sp.power(122,2)), sp.subtract(1, u_ue[i])))
    # constraints.append(square(x_u - x_ue[i]) + square(y_u - y_ue[i]) <= sp.power(h_max,2) + sp.multiply(2*sp.power(1000,2), sp.subtract(1, u_ue[i])))

constraints.append(square(x_u - x_ue) + square(y_u - y_ue) <= sp.power(h_max, 2) + sp.multiply(2*sp.power(coverage_r, 2), sp.subtract(1, u_ue)))


prob = Problem(obj_opt, constraints)
prob.solve(solver=ECOS_BB)
# prob.solve(solver=ECOS_BB, verbose=True)
t1 = time.time()
time_sol = (t1 - t0)

u_ue_sol = u_ue.value
x_uav_sol = x_u.value
y_uav_sol = y_u.value
h_uav_sol = h_u.value

print "x_coordinate of UAV: ", x_uav_sol
print "y_coordinate of UAV: ", y_uav_sol
print "height of UAV: ", h_uav_sol

print "Number of served UEs: ", u_ue_sol

print "Time executive: ", time_sol


for p in xrange(num_uav):
    uav_x = x_u[p].value
    uav_y = y_u[p].value
    UE_Loc(p, uav_x, uav_y)
    conver = plt.Circle((uav_x, uav_y), h_max, color='green', fill=False)
    ax = plt.gca()
    ax.add_patch(conver)
    plt.scatter(uav_x, uav_y, s=50, c='black', marker='D')

plt.show()
