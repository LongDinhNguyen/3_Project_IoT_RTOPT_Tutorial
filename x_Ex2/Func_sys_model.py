from cvxpy import *
import scipy as sp
import ecos
import matplotlib.pyplot as plt
import time

# sp.random.seed(2030)

# ############################################################
# This code is used to set up the UAV model
# ############################################################

coverage_r = 500 #m
plt.figure(figsize=(10, 10))
coverage = plt.Circle((0, 0), coverage_r, color='green', fill=False)
ax = plt.gca()
ax.add_patch(coverage)

class UAV(object):
    def __init__(self, h):
        self.h = h
        plt.scatter(0, 0, s=200, c='black', marker='D')

class UE_Loc(object):
    def __init__(self, id, coverage_r, max_dist, origin_x=0, origin_y=0, low_tx=0.2, low_rx=0.5):
    # def __init__(self, id, coverage_r, max_dist, origin_x=0, origin_y=0, low_tx=0.3, low_rx=0.2):
        tx_d = sp.multiply(coverage_r, sp.random.uniform(low=low_tx))
        tx_angle = sp.multiply(2, sp.multiply(sp.pi, sp.subtract(sp.random.rand(), 1)))
        self.tx_x = sp.add(origin_x, sp.multiply(tx_d, sp.sin(tx_angle)))
        self.tx_y = sp.add(origin_y, sp.multiply(tx_d, sp.cos(tx_angle)))
        plt.scatter(self.tx_x, self.tx_y, s=20, c='red')
        plt.annotate(id, (self.tx_x + 10, self.tx_y + 10))

# plt.show()