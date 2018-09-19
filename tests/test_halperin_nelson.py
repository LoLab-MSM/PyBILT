from __future__ import print_function
import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_halperin_nelson():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")


    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('halperin_nelson hn leaflet upper')

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print(analyzer.get_analysis_data('hn'))

if __name__ == '__main__':
    test_halperin_nelson()
