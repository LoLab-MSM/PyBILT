from __future__ import print_function
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
import matplotlib.pyplot as plt

#define tests
def test_ba_com_lateral_rdf():
    print("testing the com_lateral_rdf protocol of the BilayerAnalyzer...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    # ba = BilayerAnalyzer(
    #    structure='../pybilt/sample_bilayer/sample_bilayer.psf',
    #    trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
    #    selection=sel_string,
    # )
    ba = BilayerAnalyzer(
        structure='/home/blake/cl_bilayers/step6.6_equilibration.gro',
        trajectory='/home/blake/cl_bilayers/step7_prod_c0_p199.trr',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('com_lateral_rdf com_lateral_rdf resname_1 DOPE resname_2 POPC n_bins 100 range_outer 50.0 leaflet lower')
    ba.set_frame_range(interval=10)
    ba.rep_settings['com_frame']['name_dict'] = {'DOPE':['P'],'POPC':['P'],'TLCL2':['C2']}
    ba.run_analysis()
    print('com_lateral_rdf: ')
    rdf, bins = ba.get_analysis_data('com_lateral_rdf')
    plt.plot(bins, rdf)
    #plt.show()
    return

if __name__ == '__main__':
    test_ba_com_lateral_rdf()
