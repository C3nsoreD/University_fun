from math import *
import numpy as np
from numpy.linalg import inv

# CONSTANTS
height = 450
length = 500
E1, E2 = 1.41e5, 2.0e5
A1, A2 = 300, 400
alpha1, alpha2 = 1.6e-6, 1.2e-5
deltaT = 75


def _c(x) -> float:
	return np.cos(x * np.pi/180)

def _s(x) -> float:
	return np.sin(x * np.pi/180)

def Tt_k(c, s):
	# Transform matrix
	value = np.array([[c, -c],[s, -s],[-c, c],[-s, s]])
	return value

def T(c, s):
	value = np.array([[c, s, 0, 0],[0, 0, c, s]])
	# value.astype('float32')
	return value
	
beta = -np.degrees((np.arctan(height/length)))
beta = 180 + beta
# print(beta)

# Element 1 values 
c1 = _c(0)
s1 = _s(0)

# Element 2 values 
c2 = _c(beta)
s2 = _s(beta)

length2 = height/s2
Tt_k1 = Tt_k(c1, s1)
T_t1 = T(c1, s1)

Tt_k2 = Tt_k(c2, s2)
T_t2 = T(c2, s2)

# Global stiffness matrix of the first road
K_1 = np.matmul(Tt_k1, T_t1) * ((E1 * A1) / length)
# Stifness for the second element
K_2 = np.matmul(Tt_k2, T_t2) * (E2 * A2) / length2

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


Ft1 = alpha1*deltaT*E1*A1
Ft2 = alpha2*deltaT*E2*A2


ft2x = Ft2*(_c(beta))
ft2y = Ft2*(_s(beta))

# Forces = [[Ft1 + ft2x], [-ft2y]]

Forces_1 = [[Ft1], [0], [-Ft1], [0]]
Forces_2 = [[ft2x], [ft2y], [-ft2x], [-ft2y]]

total_f = [Forces_1[0], Forces_1[1], Forces_2[0], Forces_2[1], [Forces_1[2][0]+Forces_2[2][0]], [Forces_1[3][0]+Forces_2[3][0]]]
# m = 0.0000000000000000001
# inv_global = np.linalg.inv(_global + np.eye(np.asmatrix(_global).shape[1])*m) 

inv_x = np.linalg.inv(_x) 
# print(total_f)
node_3 = inv_x.dot(np.asmatrix(total_f[2:][2:]))

q_total = [[0], [0], [0], [0], [node_3[0]], [node_3[1]]]


zz = np.matmul(_global, q_total)

r1 = zz[0] - total_f[0]
r1 = zz[0] - total_f[0]
r2 = zz[1] - total_f[1]
r3 = zz[2] - total_f[2]
r4 = zz[3] - total_f[3]
r6 = zz[5] - total_f[5]

# print(zz)

# Problem 2

_x_2 = K_1[2, 2] + K_2[2, 2]

q5 = total_f[2]/_x_2

_q_total = [[0], [0], [0], [0], q5, [0]]
_zz = np.matmul(_global, _q_total)

# print(np.asmatrix(_global).round(4))
print(f"_zz {np.asmatrix(_zz).round(2)}")

_r1 = _zz[0] - total_f[0]
_r2 = _zz[1] - total_f[1]
_r3 = _zz[2] - total_f[2]
_r4 = _zz[3] - total_f[3]
_r6 = _zz[5] - total_f[5]

# print(_r1)
print(q5)
print(_r1)
print(_r2)
print(_r3)
print(_r4)
print(_r6)

print()
print(r1)
print(r2)
print(r3)
print(r4)
print(r6)
# Prints -----------------------------
# print(np.asmatrix(_global).round(3))
# print(s1, s2, c1, c2)
# print(inv_global.dot(np.asmatrix(Forces)))
print()
print(f"Results of displacement {node_3} \n")
# print((q_total))
# print(np.asmatrix(total_f))
# print("individual nodes")
# print(node_3_1)
# print(np.asmatrix(Forces).dot(inv_global))
# print(np.asmatrix(K_1))
# print(np.degrees(beta))
# print(length2)
# print(f'Ft1: {Ft1} ft2y: {ft2y} ft2x: {ft2x} Ft2 {Ft2}')
# print(np.asmatrix(K_2))

# print(total_f[2:][2:])