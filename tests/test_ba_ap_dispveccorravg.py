import pybilt.bilayer_analyzer.bilayer_analyzer as ba


def test_ba_ap_dispveccorravg():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")

    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('disp_vec_corr_avg disp_vec_corr_avg leaflet upper interval 2')

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print(analyzer.get_analysis_data('disp_vec_corr_avg'))
    return


if __name__ == '__main__':
    test_ba_ap_dispveccorravg()
