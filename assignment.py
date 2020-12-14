from math import *
import numpy as np
from numpy.linalg import inv



# CONSTANTS
height = 450
length = 500
E1, E2 = 141000, 200000
A1, A2 = 300, 400
alpha1, alpha2 = 0.0000016, 0.000012
deltaT = 75


def _c(x) -> float:
	return np.cos(x)

def _s(x) -> float:
	return np.sin(x)

def Tt_k(c, s):
	# Transform matrix
	value = np.array([[c, -c],[s, -s],[-c, c],[-s, s]])
	return value

def T(c, s):
	value = np.array([[c, s, 0, 0],[0, 0, c, s]])
	# value.astype('float32')
	return value
	
beta = (np.arctan(height/length))
length2 = height/s2


# Element 1 values 
c1 = _c(0)
s1 = _s(0)

# Element 2 values 
c2 = _c(beta)
s2 = _s(beta)
Tt_k1 = Tt_k(c1, s1)
T_t1 = T(c1, s1)

Tt_k2 = Tt_k(c2, s2)
T_t2 = T(c2, s2)

# Global stiffness matrix of the first road
K_1 = np.matmul(Tt_k1, T_t1)
# Stifness for the second element
K_2 = np.matmul(Tt_k2, T_t2)


# Global stiffness matrix for claculating 
# displacement vectors

_global = [
	[K_1[0, 0], K_1[0, 1], 0, 0, K_1[0, 2], K_1[0, 3]],

	[K_1[1, 0], K_1[1, 1], 0, 0, K_1[1, 2], K_1[1, 3]],

	[0, 0, K_2[0, 0], K_2[0, 1], K_2[0, 2], K_2[0, 3]],

	[0, 0, K_2[1, 0], K_2[1, 1], K_2[1, 2], K_2[1, 3]],

	[K_1[2, 0], K_1[2, 1], K_2[2, 1], K_2[2, 2], K_1[2, 2]+K_2[2,2], K_2[2, 3]+K_1[2, 3]],

	[K_1[3, 0], K_1[3, 1], K_2[3, 1], K_2[3, 2], K_1[3, 2]+K_2[3,2], K_2[3, 3]+K_1[3, 3]],
	]

# component of global matrix 
# calculates for u3 and v3
_x = [[K_1[2, 2]+K_2[2,2], K_2[2, 3]+K_1[2, 3]],
	  [K_1[3, 2]+K_2[3,2], K_2[3, 3]+K_1[3, 3]],]


Ft1 = alpha1*deltaT*length
Ft2 = alpha2*deltaT*length2

ft2x = Ft2*(c1)
ft2y = Ft2*(s1)

Forces = [[Ft1 + ft2x], [-ft2y]]

# Forces = [[-Ft1],[0],[0],[0],[Ft1 + ft2x],[- ft2y]]

m = 10^-9
inv_global = np.linalg.inv(_global + np.eye(np.asmatrix(_global).shape[1])*m) 


inv_x = np.linalg.inv(_x) 


# Prints -----------------------------
# print(inv_global.dot(np.asmatrix(Forces)))
print(inv_x.dot(np.asmatrix(Forces)))
# print(np.asmatrix(Forces).dot(inv_global))
# print(np.asmatrix(K_1))
print()
# print(np.asmatrix(K_2))
print(np.asmatrix(_global))




	