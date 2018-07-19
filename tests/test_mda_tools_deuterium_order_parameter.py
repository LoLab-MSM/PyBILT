from pybilt.mda_tools import deuterium_order_parameter
import MDAnalysis as mda

def test_mda_tools_deuterium_order_parameter():
    u = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf',
                     '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')

    lipid_segids = {'POPC':'MEMB'}
    lipid_acyl_carbons = {'POPC':{1:'C22', 2:'C23'}}
    lipid_acyl_hydrogens = {'POPC':{1:['H2R', 'H2S'], 2:['H3R', 'H3S']}}
    avgs, stds = deuterium_order_parameter(u, lipid_segids, lipid_acyl_carbons,
                                           lipid_acyl_hydrogens)

    print avgs['POPC']
    print stds['POPC']

if __name__ == '__main__':
    test_mda_tools_deuterium_order_parameter()
