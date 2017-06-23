import numpy as np
from pybilt.mda_tools import diffusion_coefficients as dc

def test_diffusion_coefficient_estimation():
    #set known diffusion coefficient for test data
    diff_coeff = 2.0
    #dimensions
    dim = 2
    #times for the time series
    times = np.linspace(0,100.0,50, endpoint=True)
    #msd values
    #noise
    delta = delta = np.random.uniform(-10,10, size=(50,))
    msd_vals = 2.0*dim*diff_coeff*times + delta

    #Simple application of Einstein relation
    D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, dim=dim)
    #Use linear fit
    D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, dim=dim)
    #Use anomalous diffusion fit
    D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, dim=dim)
    print("Test data has diffusion coefficient: {}".format(diff_coeff))
    print("Values from estimators:")
    print("  Basic Einstein relation:")
    print("    Diffusion coefficient: {}".format(D_e))
    print("  Linear fit:")
    print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
    print("  Anomalous diffusion fit:")
    print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))

    return

if __name__ == '__main__':
    test_diffusion_coefficient_estimation()


