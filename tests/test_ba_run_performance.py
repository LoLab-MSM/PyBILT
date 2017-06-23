import pybilt.bilayer_analyzer.bilayer_analyzer as ba

import timeit
def run_serial(analyzer):
    #make sure a fresh run
    analyzer.reset()
    analyzer.run_analysis()
    return

def test_ba_run_performance():
    setup = """\
import pybilt.bilayer_analyzer.bilayer_analyzer as ba
from __main__ import run_serial
analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                              trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                              selection="not resname CLA and not resname TIP3 and not resname POT")

analyzer.add_analysis('nnf nnf_a resname_1 DOPE resname_2 POPC leaflet upper n_neighbors 6')
analyzer.add_analysis('nnf nnf_b resname_1 POPC resname_2 POPC leaflet upper n_neighbors 10')
analyzer.add_analysis('bilayer_thickness bt')
analyzer.set_frame_range(interval=2)
"""
    print("timing the serial execution (10 repetitions):")
    time_s = timeit.timeit('run_serial(analyzer)', setup=setup, number=10)
    print("average serial execution time: {}".format(time_s/10.0))



if __name__ == '__main__':
    test_ba_run_performance()
