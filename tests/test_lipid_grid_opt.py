from __future__ import print_function
from pybilt.bilayer_analyzer import BilayerAnalyzer
import numpy as np
from pybilt.plot_generation import plot_generation_functions as pgf
from pybilt.plot_generation.plot_generation_functions import _color_list
import pybilt.lipid_grid.lipid_grid as lg
import pybilt.lipid_grid.lipid_grid_opt as lgo


def test_lipid_grid_opt():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")
    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    #analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    nbins = 100
    #analyzer.rep_settings['lipid_grid']['n_xbins'] = nbins
    #analyzer.rep_settings['lipid_grid']['n_ybins'] = nbins
    analyzer.set_frame_range(0,1,1)
    i =0
    zgrid = None
    for _frame in analyzer:
        lipid_grid = lg.LipidGrids(analyzer.reps['com_frame'], analyzer.reps['leaflets'], [0,1], nxbins=nbins, nybins=nbins)
        lipid_grid_opt = lgo.LipidGrids(analyzer.reps['com_frame'], analyzer.reps['leaflets'], [0,1], nxbins=nbins, nybins=nbins)

        lg_u = lipid_grid.leaf_grid['upper'].lipid_grid
        lgo_u = lipid_grid_opt.leaf_grid['upper'].lipid_grid
        lg_l = lipid_grid.leaf_grid['lower'].lipid_grid
        lgo_l = lipid_grid_opt.leaf_grid['lower'].lipid_grid
        print((np.array_equal(lg_u, lgo_u)))
        print((np.array_equal(lg_l, lgo_l)))
        comp_l = lg_l == lgo_l
        wh = np.where(comp_l == False)
        print(wh)
        print((lg_l[wh]))
        print((lgo_l[wh]))
        # break
    return


if __name__ == '__main__':
    test_lipid_grid_opt()
