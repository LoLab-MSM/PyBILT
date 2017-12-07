
from pybilt.bilayer_analyzer.prefab_analysis_protocols import halperin_nelson


def test_analysis_module_halperin_nelson():
    sel_string = "resname POPC DOPE TLCL2"
    halperin_nelson(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string, frame_interval=1)

    return

if __name__ == '__main__':
    test_analysis_module_halperin_nelson()
