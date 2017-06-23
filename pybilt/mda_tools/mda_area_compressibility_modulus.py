""" Implementation of area compressibilty modulus
The are compressibilty modulus is an elastic property that can be computed from the total area (in the lateral
dimension) and its fluctuations. The quantity is given by
    K_A = kT<A>[ <(A - <A>)^2>]^-1,
where A is the area per lipid, k is Boltzmann's constant
and T is the temperature (in Kelvin).
See references:
Coarse Grained Model for Semiquantitative Lipid Simulations
Siewert J. Marrink, Alex H. de Vries, and Alan E. Mark
J. Phys. Chem. B, 2004, 108 (2), pp 750-760
DOI: 10.1021/jp036508g

Molecular Dynamics Simulations of Cardiolipin Bilayers
Martin Dahlberg and Arnold Maliniak
J. Phys. Chem. B, 2008, 112 (37), pp 11655-11663
DOI: 10.1021/jp803414g

"""

# we are going to use the MDAnalysis to read in topo and traj
# numpy
import numpy as np
# import my running stats class
from pybilt.common.running_stats import *
import scipy.constants as scicon


def area_compressibility_modulus(mda_trajectory, bilayer_selection, temperature, normal='z', first=0, last=-1,
                                 interval=1):
    # set the indices for indexing the lateral dimensions
    lateral_index = [0, 1]
    if normal == 'z':
        pass
    elif normal == 'x':
        lateral_index = [1, 2]
    elif normal == 'y':
        lateral_index = [0, 2]
    nlipids = len(bilayer_selection.residues)
    per_leaflet = nlipids / 2
    # get the lateral area at each frame
    area_run = RunningStats()
    times = []
    K_a = []
    for frame in mda_trajectory[first:last:interval]:
        lateral_edges = frame.dimensions[lateral_index]
        area = lateral_edges.prod()
        area_run.push(area)
        # compute the modulus
        Ka = (scicon.k * temperature * area_run.mean())/area_run.variance()
        #conversion factor for Joules/Angstrom^2 to milliNewtons/meter
        Ka*=10.0**23
        times.append(frame.time)
        K_a.append(Ka)
    # approximate area per lipid
    times = np.array(times)
    K_a = np.array(K_a)
    K_a_run = gen_running_average(K_a)
    time_series = np.zeros((len(K_a_run), 2))
    time_series[:, 0] = times[:]
    time_series[:, 1] = K_a[:]
    return time_series[1:]
