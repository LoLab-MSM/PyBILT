
from pybilt.bilayer_analyzer.canned_analysis_protocols import dispvector_correlation


def test_analysis_module_dispvector_correlation():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    dispvector_correlation(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string, frame_interval=9)

    return

if __name__ == '__main__':
    test_analysis_module_dispvector_correlation()
