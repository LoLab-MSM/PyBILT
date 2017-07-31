from pybilt.mda_tools.mda_density_profile import mass_density_profile_multi_align
import MDAnalysis as mda
from pybilt.plot_generation.plot_generation_functions import plot_density_profile

def test_mass_density_profile_multi_align():

    u = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    ref = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    sel_popc = u.select_atoms('resname POPC')
    sel_bilayer = u.select_atoms('resname POPC DOPE TLCL2')
    centers, counts = mass_density_profile_multi_align(u, {'POPC':sel_popc, 'BILAYER':sel_bilayer}, ref, 'resid 1', refsel=sel_bilayer)

    print(sel_popc)
    print(sel_bilayer)

    print(centers)
    print(counts)

    plot_density_profile([(centers[1:-1], counts['POPC'][1:-1]), (centers[1:-1], counts['BILAYER'][1:-1])],
                         save=False, label_list=['POPC', 'BILAYER'],  show=True)
    return

if __name__ == '__main__':
    test_mass_density_profile_multi_align()