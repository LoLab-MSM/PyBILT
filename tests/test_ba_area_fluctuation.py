from __future__ import print_function
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer

#define tests
def test_ba_area_fluctuation():
    print("testing the area fluctuation protocol of the BilayerAnalyzer...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('area_fluctuation area_fluctuation')
    ba.run_analysis()
    print('area_fluctuation: ')
    print((ba.get_analysis_data('area_fluctuation')))
    return

if __name__ == '__main__':
    test_ba_area_fluctuation()
