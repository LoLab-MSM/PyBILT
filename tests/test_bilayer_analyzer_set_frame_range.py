from __future__ import print_function
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer

def test_bilayer_analyzer_set_frame_range():
    #initialize analyzer with keyword options--and default analysis
    sel_string = "resname POPC DOPE TLCL2"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )

    #normal
    first = 1
    last = 8
    ba.set_frame_range(first=first, last=last)
    print("checking normal integer input:")
    print(("input: first {} and last {}".format(first, last)))
    print(("analyzer internal: first {} and last {}".format(ba.settings['frame_range'][0], ba.settings['frame_range'][1])))

    print("checking fraction inputs:")
    first = 0.5
    last = 0.9
    ba.set_frame_range(first=first, last=last)
    print(("input: first {} and last {}".format(first, last)))
    print(("analyzer internal: first {} and last {}".format(ba.settings['frame_range'][0], ba.settings['frame_range'][1])))

    return

if __name__ == '__main__':
    test_bilayer_analyzer_set_frame_range()