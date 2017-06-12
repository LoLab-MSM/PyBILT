from vorbilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from vorbilt.bilayer_analyzer.analysis_protocols import valid_analysis

#define tests

#input options
def test_ba_loa():
    print("testing various input options...")
    #initialize analyzer with keyword options--and default analysis
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    ba = BilayerAnalyzer(
        structure='../vorbilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.remove_analysis('msd_1')
    ba.add_analysis('loa loa_u resname POPC')
    ba.add_analysis('loa loa_l resname POPC leaflet lower')
    ba.run_analysis()
    print('Lateral Orientation Angle (loa) upper: ')
    print(ba.get_analysis_data('loa_u'))
    print('Lateral Orientation Angle (loa) lower: ')
    print(ba.get_analysis_data('loa_l'))
    return

if __name__ == '__main__':
    test_ba_loa()
