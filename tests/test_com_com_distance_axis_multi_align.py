from pybilt.mda_tools.mda_distance import com_com_distance_axis_multi_align
import MDAnalysis as mda

def test_com_com_distance_axis_multi_align():

    u = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    ref = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    sel_1 = u.select_atoms('resid 1')
    sel_2 = u.select_atoms('resid 10')
    print(sel_1)
    print(sel_2)

    times, dists = com_com_distance_axis_multi_align(u, [[sel_1, sel_2]], ref, 'resid 1')

    print(times)
    print(dists)

    print("average distance: {} +- {}".format(dists[0].mean(), dists[0].std()))
    return

if __name__ == '__main__':
    test_com_com_distance_axis_multi_align()
