import vorbilt.bilayer_analyzer.bilayer_analyzer as ba


def test_ndcorr():
    analyzer = ba.BilayerAnalyzer(structure='../vorbilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")
    analyzer.add_analysis('ndcorr ndcorr_a')
    analyzer.remove_analysis('msd_1')
    analyzer.print_analysis_protocol()
    #use the phosphorous atoms instead of full lipid center of mass
    analyzer.run_analysis()
    print(analyzer.get_analysis_data('ndcorr_a'))

test_ndcorr()
