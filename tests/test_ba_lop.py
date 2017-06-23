from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from pybilt.bilayer_analyzer.analysis_protocols import valid_analysis

#define tests

#input options
def test_ba_lop():
    print("testing various input options...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('lop lop_u resname POPC')
    ba.add_analysis('lop lop_l resname POPC leaflet lower')
    ba.run_analysis()
    print('Lateral Orientation Parameter (lop) upper: ')
    print(ba.get_analysis_data('lop_u'))
    print('Lateral Orientation Parameter (lop) lower: ')
    print(ba.get_analysis_data('lop_l'))
    return

if __name__ == '__main__':
    test_ba_lop()
