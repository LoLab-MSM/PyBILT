
from pybilt.bilayer_analyzer.prefab_analysis_protocols import msd_diffusion


def test_analysis_module_msd_diffusion():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    msd_diffusion(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  resnames=['POPC', 'DOPE'])

    return

if __name__ == '__main__':
    test_analysis_module_msd_diffusion()
