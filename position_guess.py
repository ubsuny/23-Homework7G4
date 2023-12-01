'''
This is the "brute force" python file I was using to find/test
different initial positions, which I've added here for completeness.
It's sloppy but (I think) it was able to converge on most of them given the
proper filters
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
import scipy.optimize


ke2 = 197 / 137  # eV-nm   Coulomb force charge
alpha = 1.09e3  # eV      parameter of model
rho = 0.0321  # nm      parameter of model
b = 1.0  # eV      regular
c = 0.01  # nm

# Helpful solution to convert itertools combinations to numpy arrays here:
## https://stackoverflow.com/questions/33282369/convert-itertools-array-into-numpy-array
def cp(l):
    return np.fromiter(
        itertools.chain(*itertools.combinations(l, 2)), dtype=int
    ).reshape(-1, 2)


class Cluster:
    def __init__(self, r_na, r_cl):
        """
        Inputs the list of Na and Cl positions. Na has charge +1, Cl has -1.
        The array of ions itself does not change throughout the calculation, and
        neither do the charges. As such, we can just compute the combinations one time
        and refer to it throughout the calculation.
        """
        self.positions = np.concatenate((r_na, r_cl))
        self.charges = np.concatenate(
            [np.ones(r_na.shape[0]), np.full(r_cl.shape[0], -1)]
        )
        self.combs = cp(np.arange(self.charges.size))
        self.chargeprods = (
            self.charges[self.combs][:, 0] * self.charges[self.combs][:, 1]
        )
        self.rij = np.linalg.norm(
            self.positions[self.combs][:, 0] - self.positions[self.combs][:, 1], axis=1
        )

    def Vij(self):
        """Calculate a numpy vector of all of the potentials of the combinations"""
        self.Vij_ = np.zeros_like(self.rij)
        pos = self.chargeprods > 0
        neg = ~pos
        self.Vij_[pos] = ke2 / self.rij[pos] + b * (c / self.rij[pos]) ** 12
        self.Vij_[neg] = (
            -ke2 / self.rij[neg]
            + alpha * np.exp(-self.rij[neg] / rho)
            + b * (c / self.rij[neg]) ** 12
        )
        return self.Vij_

    def V(self):
        """Total potential, which is a sum of the Vij vector"""
        return np.sum(self.Vij())

    def get_vals(self):
        """Positions interpreted as a flat shape"""
        return np.reshape(self.positions, -1)

    def set_vals(self, vals):
        """Inputs flat shape of positions, used by __call__"""
        self.positions = vals.reshape(self.positions.shape)
        self.rij = np.linalg.norm(
            self.positions[self.combs][:, 0] - self.positions[self.combs][:, 1], axis=1
        )

    def __call__(self, vals):
        """Function that  scipy.optimize.minimize will call"""
        self.set_vals(vals)
        return self.V()


a = 0.2

# r_na = np.array([[0, 0, 0], [0, 0, 2*a],[0,a,a],[0,2*a,-0.5*a]])
# r_cl = np.array([[0, 0, a], [0, a, 2*a],[a,a,0],[a,- a,0]])

while True:

    r_na = np.random.rand(4, 3)*.75
    r_cl = np.random.rand(4, 3)*.75

    cluster = Cluster(r_na, r_cl)
    vals_init = cluster.get_vals()
    # print("initial Na positions:\n", r_na)
    # print("initial Cl positions:\n", r_cl)
    # print("initial positions flattened shape:\n", vals_init)
    # print("initial V  :", cluster.V())

    res = scipy.optimize.minimize(fun=cluster, x0=vals_init, tol=1e-3, method="BFGS")
    cluster.set_vals(
        res.x
    )  # For some reason, "minimize" is not updating the class at the last iteration

    if (
        res.fun < -26.9999
        and not np.isclose(res.fun, -27.7298, 0.0005)  # rectangular potential
        and not np.isclose(res.fun, -28.2358, 0.0005)  # cubic potential
        and not np.isclose(res.fun, -26.1369, 0.0005)  # just a straight line...
        and not np.isclose(res.fun, -27.7824, 0.0005)  # octagon
        and not np.isclose(res.fun, -27.3422, 0.0005)
    ):
        break
print("Final optimized cluster positions")
print(cluster.positions)
print("Final potential:", res.fun)



fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

charges = cluster.charges
x, y, z = cluster.positions[:, 0], cluster.positions[:, 1], cluster.positions[:, 2]
ax.scatter(x, y, z, c=charges, cmap="coolwarm", s = 200)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
