from cvxpy import *
import scipy as sp
import ecos
import scs


# a = sp.array([1, 2, 3, 4, 5, 6, 7, 8])
# ind = sp.array([0, 1])
# print a[ind]
#
# b = a[0:3]
# print b

# a = sp.array([2, 3])
# # b = a.reshape(1, 2)
# print norm(a, p=2).value

# print sp.divide(sp.randn(1) + 1j*sp.randn(1), sp.sqrt(2))

a = sp.rand(3, 3)
print a

print sp.sum(a[:, 0])

# print vstack([a-b])