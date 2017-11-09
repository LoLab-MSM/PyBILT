import pybilt.bilayer_analyzer.bilayer_analyzer as ba
import matplotlib.pyplot as plt

def test_ba_ap_dispveccorravg():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")

    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('spatial_velocity_corr spatial_velocity_corr leaflet lower interval 9 n_bins 120 range_outer 60.0 resname_1 POPC resname_2 POPC')

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    bins, averages = analyzer.get_analysis_data('spatial_velocity_corr')
    plt.plot(bins, averages)
    plt.show()
    return


if __name__ == '__main__':
    test_ba_ap_dispveccorravg()
