

#
#
#  Borrowing the best fit vertical structure from Ting and Rix (2018)
#
#



# packages
import numpy as np
from galpy.potential import MWPotential2014 as pot
from galpy import potential
from galpy.potential import verticalfreq, epifreq, omegac

#=================================================================================
# vertical heating
#=================================================================================

# galpy scale units:
_REFR0 = 8.   #[kpc]  --> galpy length unit
_REFV0 = 220. #[km/s] --> galpy velocity unit

# useful to have: sech2 distribution
def sech2(x):
    """
    sech squared distribution

    Arguments:
        x = flot, int or array

    Returns:
        sech2(x), same format as x

    2019-01-19
    """
    return 1./np.cosh(x)**2



#----------------------------------------------------------------
# vertical heating in action
#----------------------------------------------------------------
def vertical_action(R, t):
    """
    Model for the mean vertical action for the Galactic disk,
    as a function of radius and time. From Ting and Rix (2018)

    Arguments:
        R = Galactocentric radius (mean between birth and present day)
        t = stellar age (Gyr)

    Returns:
        <J_z>(R, t)


    2019-01-19  frankel
    """
    a0 = 1.05
    a1 = 0.55*1e-1
    b0 = 1.79
    b1 = 0.5*1e-1
    c0 = 0.91
    c1 = 1.83*1e-1
    c2 = 8.7e-2
    c3 = 14e-3
    r = (R - 8)
    gamma = a0 + a1*r
    dJzdt = (b0 + b1*r)
    Jz0 = c0 + c1*r + c2*r*r + c3*r*r*r
    return Jz0 + dJzdt * t**(gamma)

#--------------------------------------------------------------
# isothermal distribution, epicycle approximation
#--------------------------------------------------------------
def scale_height(R, R0, tau):
    """
    Scale height in function of the mean vertical action
    and the vertical frequency as a function of birth radius,
    current radius and time. For isothermal sech2 distribution

    Arguments:
        R = present day Galctocentric radius [kpc]
        R0 = birth Galactocentric radius [kpc]
        t = stellar age [Gyr]

    Returns
        hz = the scale-height [kpc]

    2019-01-19 frankel
    """
    jz = vertical_action((R+R0)/2., tau)
    nu = verticalfreq(pot, R/_REFR0)*_REFV0/_REFR0
    hz = np.sqrt(2*jz/nu)
    return hz



# vertical distribution
def p_z(R, R0, tau, z):
    """
    Vertical isothermal distribution for the Galactic disk
    as a function of Galactocentric radius and age. Based on
    Ting and Rix analysis of the vertical heating history
    of the Galactic disk.

    Arguments
        Rm = (R+R0)/2.
        h = scale_height(R, R0, tau)
        return sech2(z/h)/(2*h)

    Returns
        p_z(z | R0, R, t) = proba of being at height z given
                            the R0, R and age of a star

    2019-01-22  frankel
    """
    h = scale_height(R, R0, tau)
    return sech2(z/h)/(2*h)


