from __future__ import print_function
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from pybilt.bilayer_analyzer.analysis_protocols import valid_analysis

#define tests

#input options
def test_ba_acm():
    print("testing various input options...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('acm acm temperature 303.15')
    ba.run_analysis()
    print('ACM: ')
    print((ba.get_analysis_data('acm')))
    return

if __name__ == '__main__':
    test_ba_acm()
