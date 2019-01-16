from __future__ import absolute_import
from .bilayer_analyzer import BilayerAnalyzer
import sys
import os.path

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_script', metavar='input_script', type=str, help='The input script with its file path.')
    parser.add_argument('dump_path', metavar='dump_path', type=str, help='The file path where the analysis outputs should be dumped after the analysis is complete.')
    args = parser.parse_args()
    # get the input script from command line imputs
    input_script = os.path.abspath(args.input_script)
    # get the location to dump the output
    out_path = os.path.abspath(args.dump_path)
    # initialize the analyzer with the script
    analyzer = BilayerAnalyzer(input_file=input_script)
    # run the analyzer
    analyzer.run_analysis()
    # dump the output
    print("dumping output to: " + out_path)
    analyzer.dump_data(path=out_path)
