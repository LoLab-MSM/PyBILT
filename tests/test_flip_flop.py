import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_flip_flop():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames_flipflop.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('flip_flop flip_flop')
    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    print analyzer.get_analysis_data('flip_flop')

    return

if __name__ == '__main__':
    test_flip_flop()
