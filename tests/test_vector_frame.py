import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_vector_frame():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")

    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis("lipid_length ll leaflet upper resname TLCL2")
    analyzer.adjust_rep_setting('vector_frame', 'ref_atoms', {'DOPE':{'start':
                                                ['C218','C318'], 'end':'P'},
                                                'POPC':{'start':['C218', 'C316'],
                                                'end':'P'}, 'TLCL2':
                                                {'start':['CA18','CB18','CC18',
                                                 'CD18'], 'end':['P1', 'P3']}})
    analyzer.run_analysis()
    ll_dat = analyzer.get_analysis_data('ll')
    print('Lipid Lengths (vs time):')
    print(ll_dat)

if __name__ == '__main__':
    test_vector_frame()
