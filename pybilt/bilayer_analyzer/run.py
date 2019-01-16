from __future__ import absolute_import
from .bilayer_analyzer import BilayerAnalyzer
import sys
import os.path

if __name__ == '__main__':
    # get the input script from command line imputs
    input_script = os.path.abspath(sys.argv[1])
    # get the location to dump the output
    out_path = os.path.abspath(sys.argv[2])
    # initialize the analyzer with the script
    analyzer = BilayerAnalyzer(input_file=input_script)
    # run the analyzer
    analyzer.run_analysis()
    # dump the output
    print("dumping output to: " + out_path)
    analyzer.dump_data(path=out_path)
