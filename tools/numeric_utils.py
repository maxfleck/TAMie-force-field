
import numpy as np
from scipy.constants import Avogadro, epsilon_0, e

def calc_mie(r: np.ndarray, sigma: float, eps: float, n: int, m: int=6, key: str= "energy"):
    """
    This function computes the energy / force of the Mie potential for a given sigma, epsilon, distance, and repuslive, attrative exponent.

    Args:
        r (np.ndarray): Distance of the two atoms. Given in the same unit as sigma.
        sigma (float): Mixed sigma of the two atoms. Given in the same unit as the distance.
        eps (float): Mixed epsilon of the two atoms. Given in the energy unit.
        n (int): Repulsive exponent.
        m (int): Attractive exponent. Defaults to 6.
        key (str, optional): Key specifing if the energy (enery) or the force (force) should be computed. Defaults to "energy".

    Returns:
        mf (float): Either the energy or the force acting between the two atoms.
    """

    if eps > 0:

        # Mie constant
        c0 = n /  ( n - m ) * ( n / m )**( m / ( n - m ) )

        # energy U(r)
        if key == "energy":
            mf = c0 * eps * ( ( sigma / r ) ** n - ( sigma / r ) ** m) 
        
        # force = -dU(r)/dr
        elif key == "force":
            mf = -c0 * eps * ( m / r * ( sigma / r )**m - n / r * ( sigma/ r )**n )

    else:
        mf = np.zeros(len(r))

    return mf

def calc_coulomb(r: np.ndarray | np.ndarray, q1: float, q2: float, key: str= "energy"):
    """
    This function computes the energy / force of the Coulomb potential for given distance, and partial charges of the two atoms.

    Args:
        r (np.ndarray): Distance of the two atoms. Given in Angstrom.
        q1 (float): Partial charge of atom one. Given in portions of eletron charge.
        q2 (float): Partial charge of atom two. Given in portions of eletron charge.
        key (str, optional): Key specifing if the energy (enery) or the force (force) should be computed. Defaults to "energy".

    Returns:
        qf (float): Either the energy or the force acting between the two atoms. Unit is either kcal/(mol*AA) or kcal/(mol*AA^2)
    """
    # epsilon_0: As/Vm = C^2/(Nm^2)
    # distance: m = 10^-10AA
    # electron charge: C
    # q1*q2/(4*pi*eps_0*r) = const * q1*q2/r --> Unit: C^2 / ( C^2/(Nm^2) * m) --> Nm = J
    # * Na (1/mol) / 4184 kcal/J --> kcal/mol

    if q1 != 0.0 and q2!= 0.0:
        
        const = e**2 / ( 4 * np.pi * epsilon_0 * 1e-10 )  *  Avogadro / 4184

        # Energy U(r) = const*q1q2/r
        if key == "energy":
            qf = q1 * q2 * const / r 
        
        # Force = -dU(r)/dr = -d/dr (const*q1*q2*r^(-1)) = - (-const*q1*q2/r^2) = const*q1*q2/r^2
        elif key == "force":
            qf = q1 * q2 * const / r**2
    else:
        qf = np.zeros(len(r))

    return qf

def calc_bond(r: np.ndarray, r0: float, K: float, key = "energy"):
    """
    This function computes the energy / force of the harmonic bond potential for given distance, spring constant and equilibrium length.

    Args:
        r (np.ndarray): Distance of the two atoms. Given in the same unit as sigma.
        r0 (float): Equilibrium length of the bond.
        K (float): Harmonic spring constant of the bond.
        key (str, optional): Key specifing if the energy (enery) or the force (force) should be computed. Defaults to "energy".

    Returns:
        bf (float): Either the energy or the force acting between the two atoms. 
    """
    # Energy U(r) = K * (r-r0)^2
    if key == "energy":
        bf = K * ( r - r0 )**2
    
    # Force -dU(r)/dr = -K * 2 * (r-r0)
    elif key == "force":
        bf = -K * 2 * ( r - r0 )
        
    return bf