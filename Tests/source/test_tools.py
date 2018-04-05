# MIT License

# Copyright (c) 2018 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# @file   test_tools.py
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   August, 2017
# @brief  Python routines enabling computation of the analytical gravity acceleration created by an ellispoidal shape


import numpy as np
import scipy.integrate as integrate


'''
Computes the roots of the Phi function as defined in Scheeres' "Orbital Motion in Strongly Perturbed Environment",
page 66
Inputs:
-------
- x : x coordinate of query point in the principal axis frame on the ellispoid along largest axis
- y : y coordinate of query point in the principal axis frame on the ellispoid along intermediate axis
- z : z coordinate of query point in the principal axis frame on the ellispoid along shortest axis
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- roots: roots of the Phi function 
'''
def phi_roots(x,y,z,alpha,beta,gamma):

	coefs = np.array([
		-1, 
		- (alpha ** 2 + beta ** 2 + gamma ** 2) + x ** 2 + y ** 2 + z ** 2,
		x ** 2 * (beta ** 2 + gamma ** 2) 
		+ y ** 2 * (alpha ** 2 + gamma ** 2) 
		+ z ** 2 * (beta ** 2 + alpha ** 2)
		- (alpha ** 2 * beta ** 2 + gamma ** 2 * beta ** 2  + alpha ** 2 * gamma ** 2 ),
		(x * beta * gamma) ** 2 + (y * alpha * gamma) ** 2 + (z * alpha * beta) ** 2
		- (alpha * beta * gamma) ** 2
		])

	return np.roots(coefs)

'''
Proxy function used to computes the roots of Phi
Used for validation purposes only
Inputs:
-------
- u : input
- x : x coordinate of query point in the principal axis frame on the ellispoid along largest axis
- y : y coordinate of query point in the principal axis frame on the ellispoid along intermediate axis
- z : z coordinate of query point in the principal axis frame on the ellispoid along shortest axis
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- phi_proxy: evaluated phi_proxy function

'''
def phi_proxy(u,x,y,z,alpha,beta,gamma):

	coefs = np.array([
			-1, 
			- (alpha ** 2 + beta ** 2 + gamma ** 2) + x ** 2 + y ** 2 + z ** 2,
			x ** 2 * (beta ** 2 + gamma ** 2) 
			+ y ** 2 * (alpha ** 2 + gamma ** 2) 
			+ z ** 2 * (beta ** 2 + alpha ** 2)
			- (alpha ** 2 * beta ** 2 + gamma ** 2 * beta ** 2  + alpha ** 2 * gamma ** 2 ),
			(x * beta * gamma) ** 2 + (y * alpha * gamma) ** 2 + (z * alpha * beta) ** 2
			- (alpha * beta * gamma) ** 2
			])

	u_power = np.array([u ** 3,u ** 2, u ** 1, 1])
	return np.inner(coefs,u_power)





'''
Phi function as defined in Scheeres' "Orbital Motion in Strongly Perturbed Environment",
page 66
Inputs:
-------
- u : input
- x : x coordinate of query point in the principal axis frame on the ellispoid along largest axis
- y : y coordinate of query point in the principal axis frame on the ellispoid along intermediate axis
- z : z coordinate of query point in the principal axis frame on the ellispoid along shortest axis
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- phi: evaluated phi function
'''
def phi(u,x,y,z,alpha,beta,gamma):

	return x ** 2 / (alpha ** 2 + u) + y ** 2 / (beta ** 2 + u) + z ** 2 / (gamma ** 2 + u) - 1



'''
delta_u function as defined in Scheeres' "Orbital Motion in Strongly Perturbed Environment",
page 66
Inputs:
-------
- u : input
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- delta_u : evaluated delta_u function
'''
def delta_u(u,alpha,beta,gamma):
	return np.sqrt((alpha ** 2 + u) * (beta ** 2 + u) * (gamma ** 2 + u))


'''
Integrand of the x_axis acceleration
Inputs:
-------
- u : input
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- ax_integrand : evaluated ax_integrand function
'''
def ax_integrand(u,alpha,beta,gamma):
	return 1./((alpha ** 2 + u) * delta_u(u,alpha,beta,gamma))

'''
Integrand of the y_axis acceleration
Inputs:
-------
- u : input
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- ay_integrand : evaluated ay_integrand function
'''
def ay_integrand(u,alpha,beta,gamma):
	return 1./((beta ** 2 + u) * delta_u(u,alpha,beta,gamma))

'''
Integrand of the z_axis acceleration
Inputs:
-------
- u : input
- alpha : largest sma 
- beta : intermediate sma 
- gamma : shortest sma
Outputs:
-------
- ay_integrand : evaluated az_integrand function
'''
def az_integrand(u,alpha,beta,gamma):
	return 1./((gamma ** 2 + u) * delta_u(u,alpha,beta,gamma))



'''
Computes the analytic value of the acceleration due to gravity about 
an ellipsoidal shape
Inputs:
-------
- x : longest-axis coordinate of query point
- y : intermediate-axis coordinate of query point
- z : shortest-axis coordinate of query point
- alpha : longest semi-major axis of the body
- beta : intermediate semi-major axis of the body
- gamma : shortest semi-major axis of the body
- mu : standard gravitational parameter of the body
Outputs:
-------
- acceleration : evaluated acceleration 
'''
def acceleration(x,y,z,alpha,beta,gamma,mu):

	max_axis = max(alpha,beta,gamma)
	alpha = alpha / max_axis
	beta = beta / max_axis
	gamma = gamma / max_axis

	x = x / max_axis
	y = y / max_axis
	z = z / max_axis

	T = np.sqrt(max_axis ** 3 / mu)

	roots = phi_roots(x,y,z,alpha,beta,gamma)
	
	root = np.max(roots)

	a_x = - 3./2. * x * integrate.quad(ax_integrand,root,np.inf,args = (alpha,beta,gamma))[0]
	a_y = - 3./2. * y * integrate.quad(ay_integrand,root,np.inf,args = (alpha,beta,gamma))[0]
	a_z = - 3./2. * z * integrate.quad(az_integrand,root,np.inf,args = (alpha,beta,gamma))[0]

	return np.array([a_x,a_y,a_z]) * (max_axis / T ** 2)


def gfun(x,args):
	return args["y_min"]


def hfun(x,args):
	return args["y_max"]



def qfun(x,y,args):
	return args["z_min"]


def rfun(x,y,args):
	return args["z_max"]

def integrand(z,y,x):
	return 1./np.linalg.norm(np.array([x,y,z]))


def compute_cube_potential():
	h = 0.5
	x_point = 1.
	y_point = 2.
	z_point = 3.

	z_min = z_point - h
	z_max = z_point + h

	y_min = y_point - h
	y_max = y_point + h

	x_min = x_point - h
	x_max = x_point + h

	gfun = lambda x : y_min
	hfun = lambda x : y_max
	qfun = lambda x,y : z_min
	rfun = lambda x,y : z_max

	integral = integrate.tplquad(integrand,x_min,x_max,gfun, hfun, qfun, rfun)

	print integral



compute_cube_potential()
