from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
import matplotlib.pyplot as plt

#define tests
def test_ba_com_lateral_rdf():
    print("testing the com_lateral_rdf protocol of the BilayerAnalyzer...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('com_lateral_rdf com_lateral_rdf resname_1 POPC resname_2 POPC n_bins 50 range_outer 50.0')
    ba.run_analysis()
    print('com_lateral_rdf: ')
    rdf, bins = ba.get_analysis_data('com_lateral_rdf')
    plt.plot(bins, rdf)
    plt.show()
    return

if __name__ == '__main__':
    test_ba_com_lateral_rdf()
