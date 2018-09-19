from __future__ import print_function
from pybilt.mda_tools.mda_density_map import position_density_map_2d_multi_align
import MDAnalysis as mda
from pybilt.plot_generation.plot_generation_functions import plot_position_density_map_2d

def test_position_density_map_2d_multi_align():
    u = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    ref = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    sel_popc = u.select_atoms('resname POPC')
    sel_dope = u.select_atoms('resname DOPE')
    sel_bilayer = u.select_atoms('resname POPC DOPE TLCL2')

    print(sel_popc)
    print(sel_dope)
    sel_bilayer = u.select_atoms('resname POPC POPS')
    print(sel_bilayer)
    x_centers, y_centers, counts = position_density_map_2d_multi_align(u, {'POPC':sel_popc, 'DOPE':sel_dope}, ref, 'resid 1', fstep=1, refsel=sel_bilayer, scale_to_max=True)
    #print(centers)
    #print(counts)
    #print(len(x_centers))
    #print(len(y_centers))
    plot_position_density_map_2d(x_centers, y_centers, counts['POPC'], save=False, show=False, scaled_to_max=True)
    plot_position_density_map_2d(x_centers, y_centers, counts['DOPE'], save=False, show=False, scaled_to_max=True)
    return

if __name__ == '__main__':
    test_position_density_map_2d_multi_align()
