from __future__ import print_function
from pybilt.bilayer_analyzer import BilayerAnalyzer
import numpy as np
from pybilt.plot_generation import plot_generation_functions as pgf
from pybilt.plot_generation.plot_generation_functions import _color_list
from scipy.ndimage.filters import gaussian_filter
import scipy.constants as scicon

def test_lipid_grid_surface_area():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")
    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    nbins = 100
    analyzer.rep_settings['lipid_grid']['n_xbins'] = nbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = nbins
    analyzer.set_frame_range(0,1,1)
    i =0
    zgrid = None
    for _frame in analyzer:
        sa_upper, sa_lower = analyzer.reps['lipid_grid'].surface_area(filter_sigma=5.0)
        print(("surface for upper leaflet: {} lower leaflet {}".format(sa_upper, sa_lower)))
        box_area = analyzer.reps['current_mda_frame'].dimensions[0:2].prod()
        print(("box area: {}".format(box_area)))
        area_diff = (sa_upper+sa_lower)-(2.0*box_area)
        print(area_diff)
        area_ratio = (sa_upper + sa_lower)/(2.0*box_area)
        print(area_ratio)
        kc = (scicon.k*303.15 * np.log(300.0))/(8.0*np.pi*(area_ratio-1.0))
        print(kc)
    return


if __name__ == '__main__':
    test_lipid_grid_surface_area()
