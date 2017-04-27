import vorbilt.bilayer_analyzer.bilayer_analyzer as ba
def test_nnf():
    analyzer = ba.BilayerAnalyzer(psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('nnf nnf_a resname_1 DOPE resname_2 POPC leaflet upper n_neighbors 6')
    analyzer.add_analysis('nnf nnf_b resname_1 DOPE resname_2 DOPE leaflet upper n_neighbors 4')
    analyzer.add_analysis('nnf nnf_c resname_1 DOPE resname_2 TLCL2 leaflet upper n_neighbors 8')


    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print analyzer.get_analysis_data('nnf_a')
    print " "
    print analyzer.get_analysis_data('nnf_b')
    print " "
    print analyzer.get_analysis_data('nnf_c')

test_nnf()


