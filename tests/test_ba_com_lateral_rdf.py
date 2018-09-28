from __future__ import print_function
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
import matplotlib.pyplot as plt
import numpy as np
#define tests
def test_ba_com_lateral_rdf():
    print("testing the com_lateral_rdf protocol of the BilayerAnalyzer...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    #ba = BilayerAnalyzer(
    #    structure='../pybilt/sample_bilayer/sample_bilayer.psf',
    #    trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
    #    selection=sel_string,
    #)
    trajectory =['/home/blake/cl_bilayers/step7_prod_c0_p199.trr','/home/blake/cl_bilayers/step7_prod_c0_p201.trr']
    ba = BilayerAnalyzer(
        structure='/home/blake/cl_bilayers/step6.6_equilibration.gro',
        trajectory=trajectory,
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('com_lateral_rdf com_lateral_rdf resname_1 all resname_2 all n_bins 80 range_outer 40.0 leaflet upper')
    ba.set_frame_range(first=0, last=-1, interval=100)
    ba.rep_settings['com_frame']['name_dict'] = {'DOPE':['P'],'POPC':['P'],'TLCL2':['C2']}
    ba.settings['print_interval'] = 1
    print "frame_range: ",ba.settings['frame_range']
    ba.run_analysis()
    print('com_lateral_rdf: ')
    rdf, bins = ba.get_analysis_data('com_lateral_rdf')
    plt.plot(bins, rdf)

    hl_x = np.array([0.0, 40.0])
    hl_y = np.array([1.0, 1.0])
    plt.plot(hl_x, hl_y)

    return

if __name__ == '__main__':
    test_ba_com_lateral_rdf()
