"""Estimate self diffusion coefficients from MSD(t) data.

This module provides a set of functions to estimate self diffusion coefficients
from mean squared displacement (MSD) time trajectories.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from six.moves import range


def diffusion_coefficient_Einstein(times, msd_vals, dim=2, time_range=None):
    """Estimate diffusion coefficient from MSD(t) via Einstein relation.

    A function to estimate the diffusion constant from a mean squared
    displacement time series. This function uses the long time mean squared
    displacement approximation (Einstein relation):
        MSD = lim_(t->inf) <||r_i(t) - r_i(0)||**2>_(nsels) = 2*dim*D*t
    where D is the diffusion coefficient.

    Args:
        times (numpy.array): Array of the simulation times of each MSD point.
        msd_vals (numpy.array): Array of the MSD values.
        dim (Optional[int]): Set the dimension of the data and fit function: 1
            for 1-dimensional , 2 for 2-dimensional, or 3 for 3-dimensional.
            Defaults to 2.
        time_range (Optional[list]): Specify a time range via a list with
            format [time_start, time_end]; the range should be given in the
            same units as the time data in the Arg times. Defaults to None,
            which uses the whole time range in Args times.

    Returns:
        float: The value of the diffusion coefficient.

    References:
        1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical
            Properties of a Hydrated Lipid Bilayer from a Multinanosecond
            Molecular Dynamics Simulation, Biophysical Journal, Volume 81,
            Issue 5, 2001, Pages 2484-2494, ISSN 0006-3495,
            http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
            (http://www.sciencedirect.com/science/article/pii/S0006349501758948)
        2. Section 8.7,
            http://manual.gromacs.org/documentation/5.1.4/manual-5.1.4.pdf
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
    """Estimate the diffusion coefficent via linear least squares fit of MSD(t) data.

    Assumes the MSD data has the linear form,
        MSD(t) = 2*dim*D*t ,
    and uses a linear least squares fit to estimate the diffusion coefficent
    D.

    Args:
        times (numpy.array): Array of the simulation times of each MSD point.
        msd_vals (numpy.array): Array of the MSD values.
        dim (Optional[int]): Set the dimension of the data and fit function: 1
            for 1-dimensional , 2 for 2-dimensional, or 3 for 3-dimensional.
            Defaults to 2.
        time_range (Optional[list]): Specify a time range via a list with
            format [time_start, time_end]; the range should be given in the
            same units as the time data in the Arg times. Defaults to None,
            which uses the whole time range in Args times.

    Returns:
        tuple: Returns a tuple with format (D, Error_D), where D is the
        estimated diffusion coefficent and Error_D is the estimated error for
        D.

    References:
        1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical
            Properties of a Hydrated Lipid Bilayer from a Multinanosecond
            Molecular Dynamics Simulation, Biophysical Journal, Volume 81,
            Issue 5, 2001, Pages 2484-2494, ISSN 0006-3495;
            http://dx.doi.org/10.1016/S0006-3495(01)75894-8
            (http://www.sciencedirect.com/science/article/pii/S0006349501758948)

    """
    t = times
    msd = msd_vals
    if time_range is not None:
        t, msd = _time_block(times, msd_vals, time_range[0], time_range[1])
    #print(t, msd)
    dt = t - t[0]
    print("in diffco:", t[0], t[-1])
    slope, dummy_intercept, dummy_r_value, dummy_p_value, \
        std_err = stats.linregress(dt,msd)
    return (slope/(2.0*dim), std_err/(2.0*dim))


# Anomalous diffusion

def _msd_anom_1d(time, D_alpha, alpha):
    """1d anomalous diffusion function."""
    return 2.0*D_alpha*time**alpha


def _msd_anom_2d(time, D_alpha, alpha):
    """2d anomalous diffusion function."""
    return 4.0*D_alpha*time**alpha


def _msd_anom_3d(time, D_alpha, alpha):
    """3d anomalous diffusion function."""
    return 6.0*D_alpha*time**alpha


def diffusion_coefficient_anomalous_fit(times, msd_vals, dim=2,
                                        time_range=None):
    """Fit MSD time series data to an anomalous diffusion model.

    The anomalous diffusion function has the form:
        MSD(t) = 2 * dim * D_alpha * t**alpha
    For alpha values 0 < alpha < 1 corresponds to subdiffusion, while
    1 < alpha < 2 corresponds to superdiffusion.

    Args:
        times (numpy.array): Array of the simulation times of each MSD point.
        msd_vals (numpy.array): Array of the MSD values.
        dim (Optional[int]): Set the dimension of the data and fit function: 1
            for 1-dimensional , 2 for 2-dimensional, or 3 for 3-dimensional.
            Defaults to 2.
        time_range (Optional[list]): Specify a time range via a list with
            format [time_start, time_end]; the range should be given in the
            same units as the time data in the Arg times. Defaults to None,
            which uses the whole time range in Args times.

    Returns:
        tuple: The return is a tuple with values (D_alpha, alpha) where D_alpha is the
        anomalous diffusion coefficent and alpha is the anomalous diffusion
        power.

    References:
        1. Gerald R. Kneller, Krzysztof Baczynski, and Marta
            Pasenkiewicz-Gierula, Consistent picture of lateral subdiffusion in
            lipid bilayers: Molecular dynamics simulation and exact results,
            The Journal of Chemical Physics 135, 141105 (2011);
            doi: http://dx.doi.org/10.1063/1.3651800

    """
    t = times
    msd = msd_vals

    func = _msd_anom_2d
    if dim == 3:
        func = _msd_anom_3d
    elif dim == 1:
        func = _msd_anom_1d
    if time_range is not None:
        t, msd = _time_block(times, msd_vals, time_range[0], time_range[1])
    t -= t[0]
    popt, dummy_pcov = curve_fit(func, t, msd)
    return popt


# helper function
def _time_block(times, values, t_lower, t_upper):
    """Create a new sub block of data for the specified time range."""
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
