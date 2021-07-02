"""
Compare the output of the code to that for the analytical thin current loop which uses elliptical integrals of the 
first and second kind. The solutions to these are generated in this script using python and scipy and output into the
file current_loop_solutions.txt which are compared against the result produced in solenoid rust program.


B_z(r,z) = B_0*(t/(2*pi*sqrt(nu)))*[(k**2-nu)/(nu*(1-k^2))*E(t) + K(t)]
B_r(r,z) = B_0*(t*x)/(2*pi*nu*sqrt(nu))[1/(1-k^2) * E(t) + K(t)]          (Eq. 40)

where K,E are complete elliptic integrals of the first and second kind
B_0 = (mu_0*current)/(2*radius)
x = (z-z0)/radius
nu = r/radius
k^2 = 2*nu/(1+nu^2+x^2)
t^2 = (2k^2)/(1+k^2)
W_0 = (1/2)*radius^2*B_0

These are defined in Eq. 41 of reference paper.

"""
import scipy
from scipy.constants import mu_0,pi
from scipy.special import ellipk,ellipe
from math import sqrt 

def get_x(z,z0,radius):
    return (z-z0)/radius

def get_b0(current,radius):
    return (mu_0*current)/(2*radius)

def get_nu(r,radius):
    return r/radius

def get_k_squared(nu,x):
    return (2*nu)/(1+nu**2+x**2)

def get_t(k_squared):
    return ((2*k_squared)/(1+k_squared))**0.5 

def get_w0(radius,b0):
    return 0.5*(radius**2)*b0

def magnetic_field_axial(r,z,current,radius,z0):
    b0 = get_b0(current,radius)
    x = get_x(z,z0,radius)

    if r != 0.0:
        nu = get_nu(r,radius)
        k_squared = get_k_squared(nu,x)
        t = get_t(k_squared)
        w0 = get_w0(radius,b0)

        front_term = (b0*t)/(2*pi*sqrt(nu))
        second_elliptical_factor = (k_squared -nu)/(nu*(1-k_squared))
        first_ellip = ellipk(t)
        second_ellip = ellipe(t)
        return front_term*(second_elliptical_factor*second_ellip + first_ellip)
    else:
        return b0 * (1/(1+x**2)**1.5)   
     
def magnetic_field_radial(r,z,current,radius,z0):
    b0 = get_b0(current,radius)
    x = get_x(z,z0,radius)

    if r != 0.0:
        nu = get_nu(r,radius)
        k_squared = get_k_squared(nu,x)
        t = get_t(k_squared)
        w0 = get_w0(radius,b0)

        front_term = (b0*t*x)/(2*pi*nu*sqrt(nu))
        second_elliptical_factor = 1/(1-k_squared)
        first_ellip = ellipk(t)
        second_ellip = ellipe(t)
        return front_term*(second_elliptical_factor*second_ellip + first_ellip)
    else:
        return 0.0 

radius = 1.0
z = 0.5  
r = 0.5
current = -1.0
z0 = 0.0

print(magnetic_field_radial(r,z,current,radius,z0))
print(magnetic_field_axial(r,z,current,radius,z0))
