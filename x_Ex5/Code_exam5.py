
from cvxpy import *
import scipy as sp
import ecos
import matplotlib.pyplot as plt
import time

# sp.random.seed(2035)
coverage_r = 1000 #m
plt.figure(figsize=(8, 6))
coverage = plt.Circle((0, 0), coverage_r, color='green', fill=False)
ax = plt.gca()
ax.add_patch(coverage)

# ##############################################
# ####### Set up location of UAVs ##############
class UAV_loc(object):
    def __init__(self, id, coverage_r, origin_x=0, origin_y=0, low_tx=0.2):
        # uav_d = sp.add(sp.multiply(sp.subtract(coverage_r, 100), sp.random.rand(1, 1)), 50)
        uav_d = sp.add(sp.multiply(sp.subtract(coverage_r, 100), sp.random.uniform(low=low_tx)), 50)
        uav_angle = sp.multiply(2, sp.multiply(sp.pi, sp.random.rand(1, 1)))
        self.uav_x = sp.add(origin_x, sp.multiply(uav_d, sp.cos(uav_angle)))
        self.uav_y = sp.add(origin_y, sp.multiply(uav_d, sp.sin(uav_angle)))
        plt.scatter(self.uav_x, self.uav_y, s=40, c='red', marker='D')
        plt.annotate(id, (self.uav_x + 10, self.uav_y + 10))

    def coordinate_uav(self):
        xUAV = self.uav_x
        yUAV = self.uav_y
        return xUAV, yUAV

# #############################################
# ####### Set up location of UEs ##############
class UE_loc(object):
    def __init__(self, id, coverage_r, origin_x=0, origin_y=0, low_tx=0.5):
        tx_d = sp.add(sp.multiply(sp.subtract(coverage_r, 200), sp.random.uniform(low=low_tx)), 50)
        tx_angle = sp.multiply(2, sp.multiply(sp.pi, sp.random.rand(1, 1)))
        self.tx_x = sp.add(origin_x, sp.multiply(tx_d, sp.cos(tx_angle)))
        self.tx_y = sp.add(origin_y, sp.multiply(tx_d, sp.sin(tx_angle)))
        plt.scatter(self.tx_x, self.tx_y, s=20, c='blue', marker='o')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))

    def coordinate_ue(self):
        x_coor = self.tx_x
        y_coor = self.tx_y
        return x_coor, y_coor

# #############################################
# ####### Set up pathloss system ##############
class Path_loss(object):
    def __init__(self, fc, c_vel, alp_g, mu_los, mu_nlos, a, b, noise_var, hUAV, xUAV, yUAV, xUE, yUE):
        dist = sp.sqrt( sp.add(sp.square(sp.subtract(yUAV, yUE)), sp.square(sp.subtract(xUAV, xUE))) )
        R_dist = sp.sqrt( sp.add(sp.square(dist), sp.square(hUAV)) )
        temp1 = sp.multiply(10, sp.log10(sp.power(fc*4*sp.pi*R_dist/c_vel, alp_g)))
        temp2 = sp.multiply(sp.subtract(mu_los, mu_nlos), sp.divide(1, (1+a*sp.exp(-b*sp.arctan(hUAV/dist)-a))))
        temp3 = sp.add(sp.add(temp1, temp2), mu_nlos)
        self.pathloss = sp.divide(sp.real(sp.power(10, -sp.divide(temp3, 10))), noise_var)

    def pl_uav_ue(self):
        p_loss = self.pathloss
        return p_loss


num_uav = 4
num_ue = 2

hUAV = 100
fc = sp.multiply(2.5, sp.power(10, 9))
c_vel = sp.multiply(3, sp.power(10, 8))
alp_g = 2.6
mu_los = 1.6
mu_nlos = 23
a = 12.08
b = 0.11
Wband = sp.multiply(20, sp.power(10, 6))
noise_var = sp.multiply(sp.multiply(sp.multiply(290, sp.multiply(1.38, sp.power(10, -23))), Wband), sp.multiply(10, 0.9))
tau = 200
tau_u = num_ue
frac_dl = sp.subtract(1, sp.divide(tau_u, tau))

coor_uav = []
x_uav, y_uav = [], []
for p in xrange(num_uav):
    coor_uav.append(UAV_loc(p, coverage_r))
    x_uav.append(UAV_loc.coordinate_uav(coor_uav[p])[0])
    y_uav.append(UAV_loc.coordinate_uav(coor_uav[p])[1])

coor_ue = []
x_ue, y_ue = [], []
for p in xrange(num_ue):
    coor_ue.append(UE_loc(p, coverage_r, low_tx=0.2))
    x_ue.append(UE_loc.coordinate_ue(coor_ue[p])[0])
    y_ue.append(UE_loc.coordinate_ue(coor_ue[p])[1])

val_pl = []
chan_array = []
for j in xrange(num_uav):
    pl_value = []
    chan_val = []
    for i in xrange(num_ue):
        xUAV, yUAV, xUE, yUE = x_uav[j], y_uav[j], x_ue[i], y_ue[i]
        val_pl_tem = Path_loss(fc, c_vel, alp_g, mu_los, mu_nlos, a, b, noise_var, hUAV, xUAV, yUAV, xUE, yUE)
        pl_value.append(Path_loss.pl_uav_ue(val_pl_tem))

        chan_temp = sp.divide(sp.add(sp.randn(1), 1j*sp.randn(1)), sp.sqrt(2))
        chan_val.append(sp.multiply(Path_loss.pl_uav_ue(val_pl_tem), chan_temp))
    val_pl.append(pl_value)
    chan_array.append(chan_val)
# ## the size of val_pl is (num_uav, num_ue)
# ## the size of chan_array is (num_uav, num_ue)




Pf = 200
Pr = 100
P_cm = 200
P_0m = 200
P_cir = 9000
Pow_cir = P_cir + num_uav*(P_cm + P_0m)

Chan = chan_array
F_tot = sp.conj(Chan)*sp.linalg.inv(sp.transpose(Chan)*sp.conj(Chan))



#
# # # ############################################################
# # # This code is used to solve the opt problem
# # # ############################################################
# t0 = time.time()
#
# x_u = Variable(1)
# y_u = Variable(1)
# h_u = Variable(1)
#
# R_k2 = []
# for i in xrange(num_ue):
#     R_k2.append(square(x_u - x_ue[i]) + square(y_u - y_ue[i]) + square(h_u))
# # R_k2 = sum(square(x_u - x_ue)) + sum(square(y_u - y_ue)) + square(h_u)
# R_U2 = square(XBS - x_u) + square(YBS - y_u) + square(HBS - h_u)
#
# obj_opt = Minimize( sum(R_k2) + R_U2 )
# constraints = [h_min <= h_u, h_u <= h_max, x_min <= x_u, x_u <= x_max, y_min <= y_u, y_u <= y_max]
#
# prob = Problem(obj_opt, constraints)
# prob.solve(solver=ECOS_BB)
# # prob.solve(solver=ECOS_BB, verbose=True)
# t1 = time.time()
# time_sol = (t1 - t0)
#
# x_uav_sol = x_u.value
# y_uav_sol = y_u.value
# h_uav_sol = h_u.value
#
# print "x_coordinate of UAV: ", x_uav_sol
# print "y_coordinate of UAV: ", y_uav_sol
# print "height of UAV: ", h_uav_sol
# print "Time executive: ", time_sol
#
#
# for p in xrange(num_uav):
#     uav_x = x_u.value
#     uav_y = y_u.value
#     plt.scatter(uav_x, uav_y, s=50, c='black')
#     plt.annotate(1, (uav_x + 20, uav_y + 20))
#
#     conver = plt.Circle((uav_x, uav_y), h_max, color='green', fill=False)
#     ax = plt.gca()
#     ax.add_patch(conver)

# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, time_sol_vec_Mon)
# plt.figure(figsize=(8, 6))
# plt.clf()
# plt.plot(range_num_d2d_pairs, EE_sol_vec_Mon)
plt.show()
