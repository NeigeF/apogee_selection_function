#import useful_functions as uf
import numpy as np
from astropy.table import Table
import astropy.units as un
import astropy.coordinates as coord
from vertical_structure_yst import *
import cone
import random

   ##########################################
   #                                        #
   #   Toy model for the MW disk            #
   #                                        #
   #                                        #  
   #  Highly inspired from the work         #
   #  of Bovy+12,14,16, Mackereth+17        #
   #  and Ting+19                           #
   #                                        #
   ##########################################

#
#
# Notes: this code is not a tool to use for publication
#        but to help understand and get inspried from
#
#        it might contain bugs and slow to run
#  
#        beware of the limitations
#


# put Milky Way model parameters together
class params(object):
    def __init__(self, rd, hz):
        self.rd = rd
        self.hz = hz
        self.R_flare=-8

def p_vertical_yst(z, R):
    # simplify without birth radius and for 2 Gyr old population
    return p_z(R, R, 2, z)
        
# Make a sort of flaring model for the vertical distribution
# (good if taking data slices in age)
def p_vertical_mackereth(R, z, hz, R_flare):
    """
    flared vertical profile for the disk
    Borrowed from Bovy+12 and Mackereth+17

    Arguments:
       R = Galactocentric radius (kpc)
       z = height above the plane (kpc)
       hz = scale height
       R_flare

    returns the pdf of z (same dimension as input z)
    (normalized to 1)

    hyper parameter: R_flare = 8 kpc

    2018-07-14 frankel
    """
    ln_p = -1./hz*np.exp((R-8)/R_flare)*np.abs(z)
    norm_ct = 2*hz*np.exp((8 - R)/R_flare)
    return np.exp(ln_p) / norm_ct

# 2D model for the Galactic disk
def nu_Rz(R, z, pm):
    """
    3D spatial density p(X, Y, z) as a function of R and z
    """
    y = np.exp(-R/pm.rd)/pm.rd
    y *= p_vertical_yst(z, R)#p_vertical_mackereth(R, z, pm.hz, pm.R_flare)
    return y


# For a standard candle of given absolute magnitude, give the fraction(Distance)
# This is a lazy way to do it - not the good one! Ok just for this quick example.
def interpolate_effective_selfunc(eff_frac, mag):
    """
    2019-08-21  frankel
    Interpolate the effective selection fraction 
    (or take the weighted mean between 2 values) as given in 
    the selection function table, for a standard candle of 
    absolute magnitude mag.

    Probably buggy

    Arguments
        eff_frac: 3D array: N_fields x N_distances x N_mag
                  as given in the selection function table grids

        mag: float: absolute magnitude of the standard candle
             we want to find the selection function for (not
             hidden behind dust)

    Returns:
        eff_frac(fields, distance | mag): 2D array 
        N_fields x N_distances
        of the selection dust fractions

    """
    # find the corresponding values in magnitude
    mag_c = np.linspace(-5, -1, 27)
    diff = mag_c - mag
    m1 = diff <= 0
    i_less = np.sum(m1) - 1
    i_great = i_less + 1
    
    # find the weight to take the weighted mean
    fact1 = (mag_c[i_great] - mag)/(mag_c[i_great] - mag_c[i_less])
    fact2 = (mag - mag_c[i_less])/(mag_c[i_great] - mag_c[i_less])
    
    # return the effective fractions: take the weighted mean
    frac = eff_frac[:,:,i_less]*fact1 + eff_frac[:,:,i_great]*fact2
    return frac
    

def f_cumul(x1, y1, axes=0):
    """
    cumulative of y1
    x1 needs to be a linspace
    """
    z = np.cumsum(y1, axis=axes)*(x1[1]-x1[0])
    return z

def draw_single_pointing_expo(n, l_p, b_p, D_linear, frac_dust,
               pm, radius_plate, dust=True):
    """
    Draws positions R, z, l, b, D in single pointings

    arguments
    n = initial number of points to sample from
    l_p = longitude pointing (deg)
    b_p = latitude pointing (deg)
    Dmin_p, D_maxp = min and max distances in pointing
    Rd = scale length (approximate)
    hz = scale height (approximate)
    
    frankel 2018
    """

    # Galactocentric dist & height along the pointing
    c = coord.Galactic(l=l_p*un.degree, b=b_p*un.degree, distance=D_linear*un.kpc)
    cc = c.transform_to(coord.Galactocentric)
    x, y, z = cc.x.value, cc.y.value, cc.z.value
    R_along = np.sqrt(x*x + y*y)
    z_along = z
    
    # distribution evaluation
    p = nu_Rz(R_along, z_along, pm)*D_linear*D_linear

    # dust
    if dust == True:
        p *= frac_dust

    # cumulative distribution:
    F_cumul = f_cumul(D_linear, p)

    u = np.random.rand(n)*np.max(F_cumul)
    D = np.interp(u, F_cumul, D_linear)

    l_s, b_s, D = cone.make_cone_lbd(D, l_p, b_p, radius_plate)

    cs = coord.Galactic(l=l_s*un.degree, b=b_s*un.degree, distance=D*un.kpc)
    ccs = cs.transform_to(coord.Galactocentric)
    xs, ys, zs = ccs.x.value, ccs.y.value, ccs.z.value
    R = np.sqrt(xs*xs + ys*ys)
    z = zs
    return R, z, l_s, b_s, D, np.max(F_cumul)
    

def draw_several_pointings_expo(n, l_p, b_p, D_linear, frac_dust,
                                 pm, radius_plate, area, selfrac,
                            locids, dust=True, return_maxF=False,
                                             return_Nkeep=False):

    """
    Draw data in each pointing provided
 
    return
       R = array of galactocentric radii  [kpc]
       z = array of heights above the plane [kpc]
       l = array of longitudes [deg]
       b = array of latitudes [deg]
       D = array of distances to the Sun [kpc]
       loc = list of locations of each field
       maxF = list normalization constant for each field 
              (including selfrac & area)
       Nkeep = array (same size as R) 
               number of points generated in field
               
    frankel 2018
    """
    # make empty lists for the data
    R, z, l, b, D, locs = ([] for i in range(6))

    # max of the f cumul to downsample
    maxF = np.zeros(len(l_p))

    for i in range(len(l_p)):
        Ri, zi, li, bi, Di, mF = draw_single_pointing_expo(n, l_p[i], b_p[i],
                             D_linear[i], frac_dust[i], pm, radius_plate[i])

        R.append(Ri)
        z.append(zi)
        l.append(li)
        b.append(bi)
        D.append(Di)
        maxF[i] = mF

    # find the max F to downsample, and scale also selfrac and area
    # Added 12/12/2018 
    max_Fcumul = np.max(maxF*np.array(selfrac)*np.array(area))
    F_cumul_scaled = np.array(selfrac)*maxF*np.array(area)/max_Fcumul

    # make empty arrays for the data
    R_arr = np.array([])
    z_arr = np.array([])
    l_arr = np.array([])
    b_arr = np.array([])
    D_arr = np.array([])
    Nkeep = np.array([])

    # now downsample each. Set random seed for smoothness.
    random.seed(1)
    for i in range(len(l_p)):
        n_keep = int(n*F_cumul_scaled[i])
        indices = random.sample(range(n), n_keep)

        # concatenate
        R_arr = np.concatenate((R_arr, R[i][indices]))
        z_arr = np.concatenate((z_arr, z[i][indices]))
        l_arr = np.concatenate((l_arr, l[i][indices]))
        b_arr = np.concatenate((b_arr, b[i][indices]))
        D_arr = np.concatenate((D_arr, D[i][indices]))
        locs = np.concatenate((locs, locids[i]*np.ones(len(indices))))
        Nkeep = np.concatenate((Nkeep, n_keep*np.ones(n_keep)))

    # return the data
    if return_maxF == False:
        return R_arr, z_arr, l_arr, b_arr, D_arr, locs
    else:
        if return_Nkeep == False:
            return R_arr, z_arr, l_arr, b_arr, D_arr, locs, maxF

        else:
            return R_arr, z_arr, l_arr, b_arr, D_arr, locs, maxF, Nkeep  
