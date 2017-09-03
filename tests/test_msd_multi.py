import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_msd_multi():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    #analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('msd_multi msd_multi_1 resname all leaflet both n_tau 10 n_sigma 10')
    analyzer.add_analysis('msd_multi msd_multi_2 resname all leaflet both n_tau 3 n_sigma 3')
    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print analyzer.get_analysis_data('msd_multi_1')
    print analyzer.get_analysis_data('msd_1')
    print analyzer.get_analysis_data('msd_multi_2')

    return

if __name__ == '__main__':
    test_msd_multi()
