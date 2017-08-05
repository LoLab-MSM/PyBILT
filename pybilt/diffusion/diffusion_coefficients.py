"""Estimate self diffusion coefficients
This module provides a set of functions to estimate self diffusion coefficeints from mean squared displacement (MSD)
time trajectories.


"""


import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

def diffusion_coefficient_Einstein(times, msd_vals, dim=2, time_range=None):
    """Estimate diffusion coefficient from MSD
    A function to estimate the diffusion constant from a mean squared displacement time series.
    This function uses the long time mean squared displacement approximation (Einstein relation),
        MSD = lim_(t->inf) <||r_i(t) - r_i(0)||**2>_(nsels) = 2*dim*D*t
    where D is the diffusion coefficient.

    Returns

    """
    t = times
    msd = msd_vals
    if time_range is not None:
        t, msd = _time_block(times, msd_vals, time_range[0], time_range[1])
    nvals = len(msd)
    dt = t - t[0]
    D = msd[nvals-1]/(2.0*dim*dt[nvals-1])
    return D

def diffusion_coefficient_linear_fit(times, msd_vals, dim=2, time_range=None):
    t = times
    msd = msd_vals
    if time_range is not None:
        t, msd = _time_block(times, msd_vals, time_range[0], time_range[1])
    slope, unused_intercept, unused_r_value, unused_p_value, \
        std_err = stats.linregress(t,msd)
    return (slope/(2.0*dim), std_err/(2.0*dim))


#Anomalous diffusion


def _msd_anom_1d(time, D_alpha, alpha):
    return 2.0*D_alpha*time**alpha

def _msd_anom_2d(time, D_alpha, alpha):
    return 4.0*D_alpha*time**alpha

def _msd_anom_3d(time, D_alpha, alpha):
    return 6.0*D_alpha*time**alpha

def diffusion_coefficient_anomalous_fit(times, msd_vals, dim=2, time_range=None):

    # MSD(t) = 2*D_alpha*t^alpha
    # 0 < alpha < 1 - subdiffusion
    # 1 < alpha < 2 - superdiffusion
    #From
    # "Consistent picture of lateral subdiffusion in lipid bilayers: Molecular dynamics simulation and
    # exact results"
    # Gerald R. Kneller, Krzysztof Baczynski, and Marta Pasenkiewicz-Gierula
    # The Journal of Chemical Physics 135, 141105 (2011); doi: http://dx.doi.org/10.1063/1.3651800
    t = times
    msd = msd_vals

    func = _msd_anom_2d
    if dim == 3:
        func = _msd_anom_3d
    elif dim == 1:
        func = _msd_anom_1d
    if time_range is not None:
        t, msd = _time_block(times, msd_vals, time_range[0], time_range[1])
    popt, unused_pcov = curve_fit(func, t, msd)
    return popt


#helper function
def _time_block(times, values, t_lower, t_upper):
    npoints = len(times)
    t = []
    val = []
    for i in range(npoints):
        time = times[i]
        value = values[i]
        if time >= t_lower and time <= t_upper:
            t.append(time)
            val.append(value)
    return (np.array(t), np.array(val))
