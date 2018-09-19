from __future__ import print_function
import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_leaflet_update_interval():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")


    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('nnf nnf_a resname_1 DOPE resname_2 POPC leaflet upper n_neighbors 6')
    analyzer.add_analysis('nnf nnf_b resname_1 DOPE resname_2 DOPE leaflet upper n_neighbors 4')
    analyzer.add_analysis('nnf nnf_c resname_1 DOPE resname_2 TLCL2 leaflet upper n_neighbors 8')
    print("Default Rep Settings:")
    analyzer.print_rep_settings()
    analyzer.adjust_rep_setting('leaflets', 'update_interval', 2)
    print("\nAdjusted Rep Settings:")
    analyzer.print_rep_settings()
    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print(analyzer.get_analysis_data('nnf_a'))
    print(" ")
    print(analyzer.get_analysis_data('nnf_b'))
    print(" ")
    print(analyzer.get_analysis_data('nnf_c'))

if __name__ == '__main__':
    test_leaflet_update_interval()
