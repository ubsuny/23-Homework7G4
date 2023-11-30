"""
nacl.py - Computational Model for Na-Cl Ionic Cluster

This module defines a computational model for a cluster of Na (sodium) and Cl (chlorine) ions.
The model uses a simplified potential energy model to describe the interactions between the ions
in the cluster. The primary purpose is to simulate and optimize the configuration of the ions
to find the minimum potential energy, representing a stable or equilibrium state.

Constants:
- KE2: Coulomb force charge in electron-volts times nanometers.
- ALPHA, RHO, B, C: Parameters of the model.

Functions:
- cp(l): Generate combinations of pairs from a list.

Class:
- Cluster: Represents a cluster of Na and Cl ions.
    Attributes:
        - positions: Concatenated positions of Na and Cl ions.
        - charges: Concatenated charges of Na and Cl ions.
        - combs: Combinations of pairs of ion indices.
        - chargeprods: Product of charges for each ion pair.
        - rij: Distances between ion pairs.
    Methods:
        - calculate_potentials(): Calculate potentials between ion pairs.
        - calculate_total_potential(): Calculate the total potential energy.
        - calculate_binding_energy(): Calculate the binding energy.
        - get_flattened_positions(): Get flattened positions.
        - set_positions_from_flattened(vals): Set positions from flattened input.
        - __call__(vals): Function that scipy.optimize.minimize will call.

Usage:
- The module is used to create an instance of the Cluster class, define the initial positions
  of Na and Cl ions, and optimize the configuration to minimize the potential energy.

Example:
- The module includes an example where the initial configuration of Na and Cl ions is visualized
  in a 3D plot. 
  The scipy.optimize.minimize function is then used to find the optimized configuration,
  and the final positions and potential energy are printed.

This library serves as a computational tool for exploring the behavior of a simplified 
ionic cluster system, particularly focusing on the optimization of ion positions to achieve 
a stable configuration with minimal potential energy.
"""
# Have adjusted Uij to uij to follow sbakecase sugeested by pylint
import itertools
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import scipy.optimize



# Constants
KE2 = 197 / 137  # eV-nm   Coulomb force charge
ALPHA = 1.09e3  # eV      parameter of model
RHO = 0.0321    # nm      parameter of model
B = 1.0         # eV      regular
C = 0.01        # nm


def cp(l):
    """
    Generate combinations of pairs from a list.

    Parameters:
    - l (list): List of values.

    Returns:
    - numpy.ndarray: Array containing all combinations of pairs from the list.
    """
    return np.fromiter(itertools.chain(*itertools.combinations(l, 2)), dtype=int).reshape(-1, 2)


class Cluster:
    """
    Class representing a cluster of Na and Cl ions.
    """

    def __init__(self, r_na, r_cl):
        '''
        Initialize the Cluster.

        Parameters:
        - r_na (numpy.ndarray): Positions of Na ions.
        - r_cl (numpy.ndarray): Positions of Cl ions.
        '''
        self.positions = np.concatenate((r_na, r_cl))
        self.charges = np.concatenate([np.ones(r_na.shape[0]), np.full(r_cl.shape[0], -1)])
        self.combs = cp(np.arange(self.charges.size))
        self.chargeprods = self.charges[self.combs][:, 0] * self.charges[self.combs][:, 1]
        self.rij = np.linalg.norm(
    self.positions[self.combs][:, 0] - self.positions[self.combs][:, 1],
    axis=1
)

        self.uij_values = None  # Initialize Uij values to None

    def calculate_potentials(self):
        '''
        Calculate potentials between ion pairs.

        Returns:
        - numpy.ndarray: Vector of potentials for all combinations.
        '''
        uij = np.zeros_like(self.rij)
        pos = self.chargeprods > 0
        neg = ~pos
        uij[pos] = KE2 / self.rij[pos] + B * (C / self.rij[pos]) ** 12
        uij[neg] = (
    -KE2 / self.rij[neg]
    + ALPHA * np.exp(-self.rij[neg] / RHO)
    + B * (C / self.rij[neg]) ** 12
)

        self.uij_values = uij  # store Uij values
        return uij

    def calculate_total_potential(self):
        '''
        Calculate the total potential energy.

        Returns:
        - float: Total potential energy.
        '''
        if self.uij_values is None:
            self.calculate_potentials()
        return np.sum(self.uij_values)

    def calculate_binding_energy(self):
        '''
        Calculate the binding energy.

        Returns:
        - float: Binding energy.
        '''
        if self.uij_values is None:
            self.calculate_potentials()
        binding_energy = 0.5 * np.sum(self.uij_values)
        return binding_energy

    def get_flattened_positions(self):
        '''
        Get flattened positions.

        Returns:
        - numpy.ndarray: Flattened positions.
        '''
        return np.reshape(self.positions, -1)

    def set_positions_from_flattened(self, vals):
        '''
        Set positions from flattened input.

        Parameters:
        - vals (numpy.ndarray): Flattened positions.
        '''
        self.positions = vals.reshape(self.positions.shape)
        self.rij = np.linalg.norm(
    self.positions[self.combs][:, 0] - self.positions[self.combs][:, 1],
    axis=1
)

        self.uij_values = None  # Reset Uij values

    def __call__(self, vals):
        '''
        Function that scipy.optimize.minimize will call.

        Parameters:
        - vals (numpy.ndarray): Flattened positions.

        Returns:
        - float: Potential energy.
        '''
        self.set_positions_from_flattened(vals)
        return self.calculate_total_potential()
