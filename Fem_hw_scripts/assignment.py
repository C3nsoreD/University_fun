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
# beta = 180 + beta
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
__x = [[K_1[2, 2]+K_2[2,2], K_2[2, 3]+K_1[2, 3]],
	  [K_1[3, 2]+K_2[3,2], K_2[3, 3]+K_1[3, 3]],]
_x = _global[2:][2:]

# Force constants for elements
Ft1 = alpha1*deltaT*E1*A1
Ft2 = alpha2*deltaT*E2*A2

# the components of element 2
ft2x = Ft2*(_c(beta))
ft2y = Ft2*(_s(beta))

# Individual force vectors for each beam element 
Forces_1 = [[Ft1], [0], [-Ft1], [0]]
Forces_2 = [[-ft2x], [-ft2y], [ft2x], [ft2y]]

# Total force vectors
total_f = [Forces_1[0], Forces_1[1], Forces_2[0], Forces_2[1], [Forces_1[2][0]+Forces_2[2][0]], Forces_2[3]]
# m = 0.0000000000000000001
# inv_global = np.linalg.inv(_global + np.eye(np.asmatrix(_global).shape[1])*m) 

inv_x = np.linalg.inv(__x)

node_3 = inv_x.dot(np.asmatrix(total_f[2:][2:]))

q_total = [[0], [0], [0], [0], [node_3[0]], [node_3[1]]]


zz = np.matmul(_global, q_total)

### Stress calculations
"""
Calculation of stress components
"""
# element 1
_e1 = node_3[0]/length - alpha1*deltaT
sigma1 = _e1 * E1
N1 = sigma1 * A1
# element 2
_e2 = (node_3[0]*c2 + node_3[1]*s2)/length2 - alpha2*deltaT
sigma2 = _e2 * E2
N2 = sigma2 * A2

# reactions
reactions = [zz[i] - total_f[i] for i in range(len(total_f))]

_x_2 = K_1[2, 2] + K_2[2, 2]

q5 = total_f[4]/_x_2

_q_total = [[0], [0], [0], [0], [q5], [0]]
_zz = np.matmul(_global, _q_total)

# print(np.asmatrix(_global).round(4))
# print(f"_zz {np.asmatrix(_zz).round(2)}")

reactions_2 = [(_zz[i] - total_f[i]) for i in range(len(total_f))]

## print(_r1)

### Stress calculations
_e1 = q5/length - alpha1*deltaT
sigma1 = _e1 * E1
N1 = sigma1 * A1
# element 2
_e2 = (-q5*c2)/length2 - alpha2*deltaT
sigma2 = _e2 * E2
N2 = sigma2 * A2

print()
print("####### Problem 2")
print(f"displacement for p2: {q5}")

print("REACTIONS for 2")
print([round(float(i), 3) for i in np.asmatrix(reactions_2)])

print(
	"The resulting stress for part a: \n"
	f"e1:{_e1} sigma 1{sigma1} N1:{N1}"
)
print(
	"The resulting stress for part a: \n"
	f"e2:{_e2} sigma 2 {sigma2} N2:{N2}"
)
print()
print()



print("#### Problem 1")
print("Global stiffness matrix")
print(np.asmatrix(__x).round(2))
print()
print("REACTIONS")
print([round(float(i), 3) for i in np.asmatrix(reactions)])
print()
print("Total Force matrix") 
# print(# np.asmatrix(inv_x))
print(np.asmatrix(total_f))

print()
print(f"Displacements {[float(i) for i in node_3.round(3)]}")
print()

# print(
# 	"The resulting stress for part a: \n"
# 	f"e2:{_e2} sigma 2 {sigma2} N2:{N2}"
# )
# print(
# 	"The resulting stress for part a: \n"
# 	f"e1:{_e1} sigma 1{sigma1} N1:{N1}"
# )
