# generate points in a cone

# import packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import useful_functions as uf

def generate_r(r_min, r_max, N):
    u = np.random.rand(N)
    R = (u*(r_max**3 - r_min**3) + r_min**3)**(1./3.)
    return R

def generate_phi(N):
    phi = np.random.rand(N)*2*np.pi
    return phi

def generate_theta(theta_max, N):
    u = np.random.rand(N)
    theta = np.arccos(1 - u*(1-np.cos(theta_max)))
    return theta

def xyz_Rphitheta(R, phi, theta):
    x = R*np.sin(theta)*np.cos(phi)
    y = R*np.sin(theta)*np.sin(phi)
    z = R*np.cos(theta)
    return x, y, z

def rotate_b(x, y, z, b):
    """
    rotation matrix about b galactic latitude
    note: differs from usual spherical coordinates theta by pi/2
    2018frankel
    """
    x1 = x*np.cos(np.pi/2 - b) - z*np.sin(np.pi/2 - b)
    y1 = y
    z1 = x*np.sin(np.pi/2 - b) + z*np.cos(np.pi /2 - b)
    return x1, y1, z1

def rotate_l(x, y, z, l):
    """
    rotation matrix about galactic latitude
    rotate of minus l because 
    Galactic coordinates go with the
    sense of rotation
    2018frankel
    """
    x1 = x*np.cos(l) - y*np.sin(l)
    y1 = x*np.sin(l) + y*np.cos(l)
    z1 = z
    return x1, y1, z1

def rotate_lb(x, y, z, l, b):
    """
    rotates a cone to be centered on l and b

    2018 frankel
    """
    x1, y1, z1 = rotate_b(x, y, z, b)
    x2, y2, z2 = rotate_l(x1, y1, z1, l)
    return x2, y2, z2

def make_cone_xyz(R, l_center, b_center, radius):
    """
    Samples points in a cone centered on l_center, b_center
    uniformly on the solid surface, and with the distribution
    of R in radius (R = np array)

    retuns Cartesian coordinates

    2018 frankel
    """
    N = len(R)
    # take l_center, b_center, radius in degrees
    # and convert them to work in radians
    [l, b, rad] = np.radians([l_center, b_center, radius])
    phi = generate_phi(N)
    theta = generate_theta(rad, N)
    x, y, z = xyz_Rphitheta(R, phi, theta)
    x1, y1, z1 = rotate_lb(x, y, z, l, b)
    return x1, y1, z1

def make_cone_lbd(R, l_center, b_center, radius):
    """
    Samples points in a cone centered on l_center, b_center
    uniformly on the solid surface, and with the distribution
    of R in radius (R = np array)

    retuns Galactic coordinates

    2018 frankel
    """
    x, y, z = make_cone_xyz(R, l_center, b_center, radius)
    # l and b come out in degrees l [0,360], b[-90,90]
    l, b, D = uf.convert_xyz_lbd(x, y, z)
    return l, b, D












