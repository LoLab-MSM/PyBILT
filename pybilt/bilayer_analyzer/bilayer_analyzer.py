"""The bilayer_analyzer module.

This module defines the BilayerAnalyzer class that can be used to build protocols for analyzing and generating plots
for the analysis of a lipid bilayer trajectory.

Example:
    >> import pybilt.bilayer_analyzer as ba
    >> analyzer = ba.BilayerAnalyzer(input_file="setup.in")
    >> analyzer.run_analysis()
    >> analyzer.dump_data()
    >> analyzer.save_all_plots()
"""

# imports
import MDAnalysis as mda
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle

import os
import multiprocessing

# PyBILT imports
import com_frame as cf
import leaflet as lf
import pybilt.lipid_grid.lipid_grid as lg
import analysis_protocols as ap
import plot_protocols as pp
from pybilt.common.running_stats import *
import mda_data as md

# import the coordinate wrapping function--for unwrapping
from pybilt.mda_tools.mda_unwrap import wrap_coordinates, \
    wrap_coordinates_parallel

default_analysis_commands = ['msd msd_1']


def word_list_to_string(word_list, delimeter=" "):
    string = ""
    for word in word_list:
        string+=word+delimeter
    nchar = len(string)
    return str(string[0:nchar-1])


def print_valid_analyses(show_settings=True):
    analyses = ap.Analyses([])
    for key in ap.valid_analysis:
        a_id = key+"_t"
        analyses.add_analysis(key+" "+a_id)
    for a_id in analyses.command_protocol.keys():
        a_key = analyses.command_protocol[a_id].analysis_key
        descript = analyses.command_protocol[a_id]._short_description
        print("analysis_key: {} ---> {}".format(a_key,descript))
        if show_settings:
            print("  with settings:")
            for setting in analyses.command_protocol[a_id].settings.keys():
                value = analyses.command_protocol[a_id].settings[setting]
                print("    {} --> {}".format(setting,type(value)))
    return

def print_analysis_settings(analysis_key):
    analyses = ap.Analyses([])
    if analysis_key in ap.valid_analysis:
        a_id = analysis_key+"_t"
        analyses.add_analysis(analysis_key+" "+a_id)
        for a_id in analyses.command_protocol.keys():
            a_key = analyses.command_protocol[a_id].analysis_key
            descript = analyses.command_protocol[a_id]._short_description
            print("analysis_key: {} ---> {}".format(a_key,descript))
            print("  with settings:")
            for setting in analyses.command_protocol[a_id].settings.keys():
                value = analyses.command_protocol[a_id].settings[setting]
                print("    {} --> {}".format(setting,type(value)))
    else:
        print("{} is not a valid analysis".format(analysis_key))
    return

def print_available_plots():
    print(pp.valid_plots)
    return


def _run_analysis_alias(protocol_analyzer):
    protocol = protocol_analyzer[0]
    protocol.run_analysis(protocol_analyzer[1])
    return protocol

# the main analyzer class

class BilayerAnalyzer:
    """An analyzer class to facilitate building analyses of the bilayers

    Attributes:
        valid_commands (list of str): A list of acceptable commands that can parsed from an input script.
        required_commands (list of str): A list of the commands that are required when parsing an input script.
        required_command_error_strings (dict): A dictionary of error messages to print when required commands are
            missing in an input script.
        input_script_name (str): The name of the input script if one was supplied.
        commands (dict): A dictionary of lists keyed to the setup commands.
        analysis_protocol (obj:Analyses): An instance of Analyses used to setup and store the defined
            analysis operations to perform during analysis.
        plot_protocol (obj:PlotProtocol): An instance of PlotProtocol used to setup and store the defined plotting
            operations to perform after the analysis is done.
        mda_data (obj:MDAData): An object used to store all the MDAnalysis object data for the input structure and
            trajectory.
        norm (int): An integer representing the index of 3 element xyz coordinate arrays for the bilayer normal
            dimension.
            Default: 2
        lateral (list): A two element list of integers representing the indices of 3 element xyz coordinate arrays for
            the bilayer lateral dimensions.
            Default: [0, 1]
        normal_dimension (str): A string representing the bilayer's normal dimension (i.e. 'x', 'y', or 'z').
            Default: 'z'
        lateral_dimension (str): A string representing the bilayer's lateral dimensions (i.e. 'xy', 'yz', 'xz', etc.).
            Default: 'xy'
        current_mda_frame (obj:MDAnalsysis-->Timestep): This variable is used during the analysis loop to store a copy
            of current frame in the MDAnalysis trajectory.
        frame_range (list): A three element list containing the range of frames and the skipping interval for the
            analysis (i.e. which frames to include in the loop over the trajectory).
            Structure: [first_frame, last_frame, interval]
            Default: [0, -1, 1]
        frame_index (int): An integer used to store the index of the current frame in the trajectory during the analysis
            loop.
        com_frame (obj:COMFrame): This variable is used to store an instance of the COMFrame object for the current
            frame during the analysis loop.
        dump_com_frame (bool): A boolean that is used to determine if the COMFrame objects should be dumped ater each
            iteraction of the analysis loop. The objects are dumped as pickle files.
            Default: False
        dump_com_frame_path (str): A string file path for where COMFrame objects are dumped.
            Default: './'
        leaflets: (dict): A dictionary of the Leaflet objects extracted from the bilayer at each frame of the analysis.
            The keys are 'upper' and 'lower' corresponding to the upper and lower leaflets of the bilayer.
        dump_leaflet (bool): Determines whether the leaflets are dumped (as pickles) after each frame in the analysis
            loop.
            Default: False
        dump_leaflet_path (str): A string file path for where the leaflet object dictionaries are dumped.
            Default: './'
        lipid_grid (obj:LipidGrids): An instance of the LipidGrids object for the current frame of the trajectory
            (analysis) loop.
        lg_nxbins (int): An integer defining the number of bins to use in the 'x' dimension of the LipidGrid
            objects.
            Default: 10
        lg_nybins (int): An integer defining the number of bins to use in the 'y' dimension of the LipidGrid
            objects.
            Default: 10
        dump_lipid_grid (bool): Determines wheter the LipidGrids objects built during interations of the analysis loop
            are dumped as pickle files.
            Default: False
        dump_lipid_grid_path (str): A string containing the path where dumped LipidGrids should be dumped.
            Default: './'
        first_com_frame (obj:COMFrame): Stores an instance of the COMFrame object from the first frame of the analysis
            loop.
    """
    # for the input script parser
    valid_commands = ["structure", "trajectory", "analysis", "selection", "frames",
                      "plot", "lipid_grid"]
    required_commands = ['structure', 'trajectory', 'selection']
    required_command_error_strings = {'structure': "the structure file needs to specified with command: \"structure path/psf_file_name\""}
    required_command_error_strings[
        'trajectory'] = "the trajectory file needs to specified with command: \"trajectory path/trajectory_file_name\""
    required_command_error_strings[
        'selection'] = "an MDAnalysis syntax selection needs to specified with command: \"selction \'selection string\'\""

    def __init__(self, structure=None, trajectory=None, selection=None,
                 input_file=None, input_dict=None):
        """BilayerAnalyzer initialization.
        The initialization initialially parses all the inputs and sets some the attribute values. It also calls the
        initialization of the Analyses and PlotProtocol objects and builds the MDAData object.
        Args:
            structure (str): Optional, the path and filename of the structure file.
            trajectory (str;list): Optional, the path and filename of the trajectory file. Also accepts a list of
                filenames.
            selection (str): Optional, the MDAnalysis compatible string to select the bilayer components.
            input_file (str): Optional, the path and filename of input setup file.
            input_dict (dict): Optional, a dictionary of keyed by valid commands and their input values.
        """
        #adjustable settings -- used internally in analysis loops and
        #   are passed to the individual analyses.
        self.settings = dict()
        #set default settings
        self.settings['frame_range'] = [0, -1, 1]
        self.settings['print_interval'] = 5
        self.settings['norm'] = 2
        self.settings['lateral'] = [0, 1]
        self.settings['lateral_dimension'] = "xy"
        self.settings['normal_dimension'] = "z"
        self.settings['frame_index'] = 0

        #buildable representation objects and their adjustable settings
        self.reps = dict()
        self.rep_settings = dict()
        #mda frame
        self.reps['current_mda_frame'] = None
        #com_frame
        self.reps['com_frame'] = None
        self.rep_settings['com_frame'] = {'dump' : False,
                                         'dump_path' : "./",
                                         'name_dict' : None}
        #leaflets
        self.reps['leaflets'] = None
        self.rep_settings['leaflets'] = {'dump': False,
                                         'dump_path' : "./"}
        #lipid_grid
        self.reps['lipid_grid'] = None
        self.rep_settings['lipid_grid'] = {'dump' : False,
                                           'dump_path' : "./",
                                           'n_xbins' : 10,
                                           'n_ybins' : 10}
        #first_com_frame
        self.reps['first_com_frame'] = None

        #non-adjustable internal settings and data
        self._input_script_name = None
        if input_file is not None:
            self._input_script_name = os.path.abspath(input_file)
        self._frame_loop_count = 0

        #input file parsing
        if input_file is not None:
            print ("parsing input file \'" + input_file + "\'...")
            self._commands = self.parse_input_script(self._input_script_name)
            # parse 'frames' input key--
            # for setting the frame range of the analysis
            if 'frames' in self._commands.keys():
                f_args = self._commands['frames']
                for i in range(0, len(f_args), 2):
                    arg_key = f_args[i]
                    arg_value = f_args[i + 1]
                    if arg_key == 'first':
                        self.settings['frame_range'][0] = arg_value
                    elif arg_key == 'last':
                        self.settings['frame_range'][1] = arg_value
                    elif arg_key == 'interval':
                        self.settings['frame_range'][2] = arg_value

            # parse inputs for lipid_grid settings
            if 'lipid_grid' in self._commands.keys():
                lg_args = self._commands['lipid_grid']
                for i in range(0, len(lg_args), 2):
                    arg_key = lg_args[i]
                    arg_value = lg_args[i + 1]
                    if arg_key == 'n_xbins':
                        self.rep_settings['lipid_grid']['n_xbins'] = arg_value
                    elif arg_key == 'n_ybins':
                        self.rep_settings['lipid_grid']['n_ybins'] = arg_value

        elif (structure is not None) and ((trajectory is not None) and
                                             (selection is not None)):
            print ("parsing inputs...")
            self._commands = dict()
            self._commands['structure'] = structure
            self._commands['trajectory'] = trajectory
            self._commands['selection'] = selection

        elif input_dict is not None and isinstance(input_dict, dict):
            self._commands = dict()
            id_keys = input_dict.keys()
            for key in self.required_commands:
                if key not in id_keys:
                    raise RuntimeError("key \'{}\' needs to be included in the input dictionary.".format(key))
            for key in id_keys:
                if key in self.valid_commands:
                    self._commands[key] = input_dict[key]
            # parse 'frames' input key--
            # for setting the frame range of the analysis
            if 'frames' in self._commands.keys():
                f_args = self._commands['frames']
                for i in range(0, len(f_args), 2):
                    arg_key = f_args[i]
                    arg_value = f_args[i + 1]
                    if arg_key == 'first':
                        self.settings['frame_range'][0] = int(arg_value)
                    elif arg_key == 'last':
                        self.settings['frame_range'][1] = int(arg_value)
                    elif arg_key == 'interval':
                        self.settings['frame_range'][2] = int(arg_value)

                        # parse inputs for lipid_grid settings
            if 'lipid_grid' in self._commands.keys():
                lg_args = self._commands['lipid_grid']
                for i in range(0, len(lg_args), 2):
                    arg_key = lg_args[i]
                    arg_value = lg_args[i + 1]
                    if arg_key == 'n_xbins':
                        self.rep_settings['lipid_grid']['n_xbins'] = int(arg_value)
                    elif arg_key == 'n_ybins':
                        self.rep_settings['lipid_grid']['n_ybins'] = int(arg_value)
        else:
            error = 'Must provide input_file or all three options for' \
                    '"structure", "trajectory", and "selection"'
            raise RuntimeError(error)

            # set up the analysis protocol
        print ("setting up analysis protocol:")
        if 'analysis' in self._commands.keys():
            self.analysis_protocol = ap.Analyses(
                self._commands['analysis'])
        else:
            self.analysis_protocol = ap.Analyses(
                default_analysis_commands)
        self.print_analysis_protocol()
        # set up the plot protocol
        print ("setting up plot protocol")
        if "plot" in self._commands.keys():
            self.plot_protocol = pp.PlotProtocol(self._commands['plot'],
                                                 self.analysis_protocol)
        else:
            self.plot_protocol = pp.PlotProtocol(None, self.analysis_protocol)
        for i in self._commands:
            print(i, self._commands[i])
        # build selection string for the MDAData object
        sel_string = self._commands['selection']

        # print (item)
        # sel_string = "not resname CLA and not resname TIP3 and not resname POT"
        # print "bilayer selection string:"
        # print sel_string
        # build the MDAData object
        print ('building the MDAnalysis objects...')

        self.mda_data = md.MDAData(self._commands['structure'],
                                   self._commands['trajectory'],
                                   sel_string)



        #more non-adjustable internal settings and data

        self._first_frame = True
        self._first_com = True
        #   for the iterator version of analysis loop
        self._current_frame = self.settings['frame_range'][0]
        self._last_frame = self.settings['frame_range'][1]
        self._oldcoords = None
        if self._last_frame < 0:
            self._last_frame += len(self.mda_data.mda_trajectory)

        return

    def __repr__(self):
        return 'BilayerAnalyzer'

    def parse_input_script(self, input_script_name):
        """Parses input setup scripts.
        Args:
            input_script_name (str): A string with the file path and name of
            the input script to be parsed.

        Returns:
            (dict): A dictionary containing lists of commands for each type of
             command key (e.g. analysis, plot
                etc.).
        Raises:
            RuntimeError: A runtime error is given if there is an invalid
            command key in the input script. A runtime error is also raised if
            the required commands are not provided in the input script.
        """
        commands = {}
        # args = {}
        with open(input_script_name) as ifile:
            for line in ifile:
                line = line[:-1]
                words = line.split()

                skip_line = (len(words) <= 1) or ((words[0] == '#') or
                                                  (words[0][0] == '#'))
                if skip_line:
                    continue
                key = words[0]
                if key in ['structure', 'trajectory', 'selection']:
                    second_term = line.replace(key, '', 1).strip()
                    if key == 'selection':
                        second_term = eval(second_term)
                    commands[key] = second_term
                else:
                    if key in self.valid_commands:
                        if key in commands.keys():
                            input_string = word_list_to_string(words[1:])
                            print (input_string)
                            commands[key].append(input_string)
                        else:
                            commands[key] = []
                            input_string = word_list_to_string(words[1:])
                            print(input_string)
                            commands[key].append(input_string)
                    else:
                        print("input command {} is not a valid"
                              " command".format(key))
                        error_str = "invalid input command '{}'".format(key)
                        raise RuntimeError(error_str)
        for required in self.required_commands:
            if required not in commands.keys():
                error_string = self.required_command_error_strings[required]
                raise RuntimeError(error_string)
                exit

        return commands

    def print_analysis_protocol(self):
        """Print the analysis protocol."""
        self.analysis_protocol.print_protocol()
        return

    def add_analysis(self, analysis_in):
        """Add a analysis to the analysis protocol.
        Args:
            analysis_in (str, list/tuple, or dict): The inpute defining the analysis key, analysis
             id, and settings for the new analysis.

        """
        self.analysis_protocol.add_analysis(analysis_in)
        return

    def remove_analysis(self, analysis_id):
        """Remove a specified analysis from the comptue protocol.
        Args:
            analysis_id (str): The string analysis id of the analysis to be
             removed from the protocol.

        """
        self.analysis_protocol.remove_analysis(analysis_id)
        return

    def remove_all_analyses(self):
        """Remove all analyses from the analysis protocol.
        """
        self.analysis_protocol.remove_all()
        return

    def get_analysis_ids(self):
        """Return the ids of currently initialized analysiss in the analysis
        protocol.
        Returns:
            (list): A list string analysis ids.

        """
        return self.analysis_protocol.analysis_ids

    def get_analysis_data(self, analysis_id):
        """Return the analysis output for the specified analysis.
        Args:
            analysis_id (str): A string analysis id.

        Returns:
            Variable: The analysis output of the specified analysis.
            The type/structure will depend on the analysis.

        """
        return self.analysis_protocol.command_protocol[analysis_id].get_data()

    def dump_data(self, path=None):
        """Dump all the anlysis outputs from the analysiss as pickle files."""
        self.analysis_protocol.dump_data(path=path)
        return

        ## plot data/access

    def print_plot_protocol(self):
        """Print the protocol for the plots that have initialized."""
        self.plot_protocol.print_protocol()
        return

    def add_plot(self, plot_string):
        """Add a new plot to the plot protocol.
        Args:
            plot_string (str): A string with the plot key, id, and arguments.

        """
        self.plot_protocol.add_plot(plot_string, self.analysis_protocol)
        return

    def remove_plot(self, plot_id):
        """Remove the specified plot from the protocol.
        Args:
            plot_id (str): Plot id of the plot to be removed.

        """
        self.plot_protocol.remove_plot(plot_id)
        return

    def print_plot_ids(self):
        """Print the ids of initialized plots in the plot protocol."""
        print(self.plot_protocol.plot_ids)

    def get_plot_ids(self):
        """Return the ids of the plots in the plot protocol.

        Returns:
            (list): A list of string plot ids that have been initialized in
             the plot protocol.
        """
        return self.plot_protocol.plot_ids

    def show_plot(self, plot_id):
        """Show the specified plot via matplotlibs .show() function.
        Args:
            plot_id (str): The string plot id to generate and show the plot for.

        """
        self.plot_protocol.command_protocol[plot_id].show_plot(
            self.analysis_protocol)
        return

    def generate_plot(self, plot_id):
        """Generate the plot image file (.eps) for the specified plot.
        Args:
            plot_id (str): The string plot id for the plot generate.

        """
        self.plot_protocol.command_protocol[plot_id].generate_plot(
            self.analysis_protocol)
        return

    def save_all_plots(self):
        """Generates the image files (.eps) for all plots in the plot protocol.
        """
        self.plot_protocol.save_plots(self.analysis_protocol)
        return

    # mda_trajectory data/access
    def update_mda_trajectory(self, new_trajectory):
        """Set a new trajectory file to read frames from.
        Args:
            new_trajectory (str): A file path and file name for the new
            trajectory.

        """
        print ("updating mda trajectory to:", new_trajectory)
        self._commands['trajectory'] = [new_trajectory]
        self.mda_data.update_trajectory(new_trajectory)
        return

    # buildable objects functions
    def set_dump_com_frame(self, on=True, path="./"):
        """Turn on (or off) the dumping of COMFrame objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
             dumped.
                Default: "./"

        """
        self.rep_settings['com_frame']['dump'] = on
        if path != self.rep_settings['com_frame']['dump_path']:
            self.rep_settings['com_frame']['dump_path'] = path
        return

    def set_dump_leaflet(self, on=True, path="./"):
        """Turn on (or off) the dumping of leaflet objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
             dumped.
                Default: "./"

        """
        self.rep_settings['leaflets']['dump'] = on
        if path != self.rep_settings['leaflets']['dump_path']:
            self.rep_settings['leaflets']['dump_path'] = path
        return

    def set_dump_lipid_grid(self, on=True, path="./"):
        """Turn on (or off) the dumping of LipidGrid objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
            dumped.
                Default: "./"

        """
        self.rep_settings['lipid_grid']['dump'] = on
        if path != self.rep_settings['lipid_grid']['dump_path']:
            self.rep_settings['lipid_grid']['dump_path'] = path
        return

    def set_frame_range(self, first=0, last=-1, interval=1):
        """Set the frame range and interval to use in the analysis.
        Args:
            first (int): The first frame.
            last (int): The last frame.
            interval (int): The interval to skip between frames in the analysis
             loop.

        """
        if first != self.settings['frame_range'][0]:
            self.settings['frame_range'][0] = first
        if last != self.settings['frame_range'][1]:
            self.settings['frame_range'][1] = last
            self._last_frame = last
            if self._last_frame < 0:
                self._last_frame+=len(self.mda_data.mda_trajectory)
        if interval != self.settings['frame_range'][2]:
            self.settings['frame_range'][2] = interval
        return

    def reset(self):
        """ Clears the analysis output stored in the analysiss.
        """
        self.settings['frame_index'] = 0
        # for iterator
        self._first_frame = True
        self._first_com = True
        self._current_frame = self.settings['frame_range'][0]
        self._frame_loop_count = 0
        #analyses
        self.analysis_protocol.reset()
        return

    # analysis
    def run_analysis(self, nprocs=1):
        """ Runs the analsysis.
        The function performs the loop over the trajectory. At each frame it
        builds the necessary objects (e.g. COMFrame) and then executes the
        analysis of each analysis that was initialized in the setup.

        Args:
            nprocs (int): An integer specifying the number of cores to use in
            multithreaded parallelelization.

        """

        parallel = False
        if nprocs > 1:
            parallel = True

        # now we need to unwrap the coordinates
        natoms = self.mda_data.natoms
        oldcoord = np.zeros((natoms, 3))
        currcoord = np.zeros((natoms, 3))
        wrapcoord = np.zeros((natoms, 3))
        first_frame_coord = np.zeros((natoms, 3))
        index = self.mda_data.bilayer_sel.indices
        firstframe = True
        first_com = True
        self.settings['frame_index'] = self.settings['frame_range'][0]
        n_mda_frames = len(self.mda_data.mda_trajectory)
        if self.settings['frame_range'][1] < 0:
            self.settings['frame_range'][1]+=n_mda_frames+1
        elif self.settings['frame_range'][1] > 0:
            self.settings['frame_range'][1]+=1
        for frame in self.mda_data.mda_trajectory[
                     self.settings['frame_range'][0]:self.settings['frame_range'][1]:self.settings['frame_range'][
                         2]]:
            self.reps['current_mda_frame'] = frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                first_frame_coord = np.copy(oldcoord)
                firstframe = False
                wrapcoord = np.copy(currcoord)
            else:
                abc = frame.dimensions[0:3]
                #oldcoord = self._oldcoords
                if parallel:
                    wrapcoord = wrap_coordinates_parallel(abc, currcoord,
                                                          oldcoord,
                                                          nprocs=nprocs)
                else:
                    #oldcoord = self._oldcoords.copy()
                    wrapcoord = wrap_coordinates(abc, currcoord, oldcoord)
                # frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
            #self._oldcoords = oldcoord
                # print ("wrapped coords:")
            # print (wrapcoord)
            #if self.analysis_protocol.use_objects['com_frame']:
            # now build the COMFrame
            self.reps['com_frame'] = cf.COMFrame(frame, self.mda_data.bilayer_sel,
                                         wrapcoord, name_dict=self.rep_settings['com_frame']['name_dict'])
            if first_com:
                self.reps['first_com_frame'] = self.reps['com_frame']
                first_com = False
            # now we can assign the lipids to the leaflets
            self.reps['leaflets'] = {'upper': lf.Leaflet('upper'),
                             'lower': lf.Leaflet('lower')}
            if self.rep_settings['com_frame']['dump']:
                ofname = self.rep_settings['com_frame']['dump_path'] + "com_frame_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['com_frame'], ofile)
            if self.rep_settings['leaflets']['dump']:
                ofname = self.rep_settings['leaflets']['dump_path'] + "leaflets_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['leaflets'], ofile)

            # first compute the average position along the normal direction
            zstat = RunningStats()
            for lipcom in self.reps['com_frame'].lipidcom:
                zstat.push(lipcom.com_unwrap[self.settings['norm']])
            zavg = zstat.mean()
            # now loop over the lipids
            l = 0
            for lipcom in self.reps['com_frame'].lipidcom:
                pos = ""
                # decide which leaflet
                #    print (lipcom.com_unwrap)
                #    print (lipcom.com)
                if lipcom.com_unwrap[self.settings['norm']] > zavg:
                    pos = 'upper'
                elif lipcom.com_unwrap[self.settings['norm']] < zavg:
                    pos = 'lower'
                # add to the chosen leaflet
                self.reps['com_frame'].lipidcom[l].leaflet = pos
                self.reps['leaflets'][pos].add_member(l, lipcom.type, lipcom.resid)
                l += 1

            if self.analysis_protocol.use_objects['lipid_grid']:
                self.reps['lipid_grid'] = lg.LipidGrids(self.reps['com_frame'], self.reps['leaflets'],
                                                self.settings['lateral'],
                                                nxbins=self.rep_settings['lipid_grid']['n_xbins'],
                                                nybins=self.rep_settings['lipid_grid']['n_ybins'])
                if self.rep_settings['lipid_grid']['dump']:
                    ofname = self.rep_settings['lipid_grid']['dump_path'] + "lipid_grid_" + str(
                        frame.frame) + ".pickle"
                    with open(ofname, 'wb') as ofile:
                        pickle.dump(self.reps['lipid_grid'], ofile)

                        # lipid_grid = None
            # now do analyses
            if self._frame_loop_count % self.settings['print_interval'] == 0:
                print("Frame", frame.frame)
            i = 0
            for analysis_id in self.analysis_protocol.analysis_ids:
                if self._frame_loop_count % self.settings['print_interval'] == 0:
                    print ("analysis " + analysis_id)
                self.analysis_protocol.command_protocol[analysis_id].run_analysis(
                    self.settings, self.reps, self.mda_data)
                # comp_out = analysis.run_analysis(self)
                # print (comp_out)
                #   analysis_out[i].append(comp_out)
                i += 1
            if self._frame_loop_count % self.settings['print_interval'] == 0:
                print(" ")
            self.settings['frame_index'] += self.settings['frame_range'][2]
            self._frame_loop_count+=1
            # print ('analysis_out:')
            # print (analysis_out)

    # analysis
    def run_analysis_mp(self, nprocs=1):
        """ Runs the analsysis.
        The function performs the loop over the trajectory. At each frame it
        builds the necessary objects (e.g. COMFrame) and then executes the
        analysis of each analysis that was initialized in the setup.

        Args:
            nprocs (int): An integer specifying the number of cores to use in
            multithreaded parallelelization.

        """

        parallel = False
        if nprocs > 1:
            parallel = True

        # now we need to unwrap the coordinates
        natoms = self.mda_data.natoms
        oldcoord = np.zeros((natoms, 3))
        currcoord = np.zeros((natoms, 3))
        wrapcoord = np.zeros((natoms, 3))
        first_frame_coord = np.zeros((natoms, 3))
        index = self.mda_data.bilayer_sel.indices
        firstframe = True
        first_com = True
        self.settings['frame_index'] = self.settings['frame_range'][0]
        n_mda_frames = len(self.mda_data.mda_trajectory)
        if self.settings['frame_range'][1] < 0:
            self.settings['frame_range'][1] += n_mda_frames + 1
        #build lists of which analysis protocols can be implemented
        # in parallel with multiprocessing (i.e. which can be properly pickled
        # and used without the MDAnalysis objects in mda_data
        run_serial = []
        run_parallel = []
        for a_id in self.analysis_protocol.command_protocol.keys():
        #    print a_id, " ", self.analysis_protocol.command_protocol[a_id]._pickleable
            if self.analysis_protocol.command_protocol[a_id]._pickleable:
        #        print(a_id)
                run_parallel.append(a_id)
            else:
                run_serial.append(a_id)

        #print(run_serial)
        #print(run_parallel)
        # create process pool
        pool = multiprocessing.Pool(processes=nprocs)
        #def in_func(input):
        #    ap._run_analysis_alias(input)

        in_func = _run_analysis_alias
        for frame in self.mda_data.mda_trajectory[
                     self.settings['frame_range'][0]:self.settings['frame_range'][1]:self.settings['frame_range'][
                         2]]:
            self.reps['current_mda_frame'] = frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                first_frame_coord = np.copy(oldcoord)
                firstframe = False
                wrapcoord = np.copy(currcoord)
            else:
                abc = frame.dimensions[0:3]
                # oldcoord = self._oldcoords
                if parallel:
                    wrapcoord = wrap_coordinates_parallel(abc, currcoord,
                                                          oldcoord,
                                                          nprocs=nprocs)
                else:
                    # oldcoord = self._oldcoords.copy()
                    wrapcoord = wrap_coordinates(abc, currcoord, oldcoord)
                # frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
                # self._oldcoords = oldcoord
                # print ("wrapped coords:")
            # print (wrapcoord)
            # if self.analysis_protocol.use_objects['com_frame']:
            # now build the COMFrame
            self.reps['com_frame'] = cf.COMFrame(frame, self.mda_data.bilayer_sel,
                                         wrapcoord, name_dict=self.rep_settings['com_frame']['name_dict'])
            if first_com:
                self.reps['first_com_frame'] = self.reps['com_frame']
                first_com = False
            # now we can assign the lipids to the leaflets
            self.reps['leaflets'] = {'upper': lf.Leaflet('upper'),
                             'lower': lf.Leaflet('lower')}
            if self.rep_settings['com_frame']['dump']:
                ofname = self.rep_settings['com_frame']['dump_path'] + "com_frame_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['com_frame'], ofile)
            if self.rep_settings['leaflets']['dump']:
                ofname = self.rep_settings['leaflets']['dump_path'] + "leaflets_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['leaflets'], ofile)

            # first compute the average position along the normal direction
            zstat = RunningStats()
            for lipcom in self.reps['com_frame'].lipidcom:
                zstat.push(lipcom.com_unwrap[self.settings['norm']])
            zavg = zstat.mean()
            # now loop over the lipids
            l = 0
            for lipcom in self.reps['com_frame'].lipidcom:
                pos = ""
                # decide which leaflet
                #    print (lipcom.com_unwrap)
                #    print (lipcom.com)
                if lipcom.com_unwrap[self.settings['norm']] > zavg:
                    pos = 'upper'
                elif lipcom.com_unwrap[self.settings['norm']] < zavg:
                    pos = 'lower'
                # add to the chosen leaflet
                self.reps['com_frame'].lipidcom[l].leaflet = pos
                self.reps['leaflets'][pos].add_member(l, lipcom.type, lipcom.resid)
                l += 1

            if self.analysis_protocol.use_objects['lipid_grid']:
                self.reps['lipid_grid'] = lg.LipidGrids(self.reps['com_frame'], self.reps['leaflets'],
                                                self.settings['lateral'],
                                                nxbins=self.rep_settings['lipid_grid']['n_xbins'],
                                                nybins=self.rep_settings['lipid_grid']['n_ybins'])
                if self.rep_settings['lipid_grid']['dump']:
                    ofname = self.rep_settings['lipid_grid']['dump_path'] + "lipid_grid_" + str(
                        frame.frame) + ".pickle"
                    with open(ofname, 'wb') as ofile:
                        pickle.dump(self.reps['lipid_grid'], ofile)

                        # lipid_grid = None
            # now do analyses
            if self._frame_loop_count % self.settings['print_interval'] == 0:
                print("Frame", frame.frame)
            r_analyses = []
            for analysis_id in self.analysis_protocol.analysis_ids:
                if self._frame_loop_count % self.settings['print_interval'] == 0:
                    print ("analysis " + analysis_id)

            for a_id in run_parallel:
            #    print(a_id)
                r_analyses.append([self.analysis_protocol.command_protocol[a_id], self])
                #r_analyses.append('hello')
            #print(r_analyses)

            #execute the analyses using pool.map
            #print("in_func")
            #print(in_func)
            #print("r_analyses")
            #print(r_analyses)
            results = pool.map(in_func, r_analyses)

            #print(results)
            #extract the results
            for item in results:
                analysis_id = item.analysis_id
                self.analysis_protocol.command_protocol[analysis_id] = item
            del results
            for a_id in run_serial:
                self.analysis_protocol.command_protocol[a_id].run_analysis(
                    self.settings,
                    self.reps,
                    self.mda_data)

            if self._frame_loop_count % self.settings['print_interval'] == 0:
                    print(" ")
            self.settings['frame_index'] += self.settings['frame_range'][2]
            self._frame_loop_count += 1
            # print ('analysis_out:')
            # print (analysis_out)
        #close and join the pool
        pool.close()
        pool.join()
        return

    @staticmethod
    def print_available_analysis():
        """Prints the keys of analysiss that can initialized.

        """
        print(ap.command_protocols.keys())
        return

    @staticmethod
    def print_available_plots():
        """Prints a list of the keys of the available plot types.        """
        print (pp.command_protocols.keys())
        return

    @staticmethod
    def available_analysis():
        """Returns the list of the keys of the available plot types.

        Returns:
            (list): A list of string keys for available plot types.

        """
        return pp.command_protocols.keys()

    @staticmethod
    def available_analysis():
        """ Returns the available analysiss.
        Returns:
            (list): A list of string keys corresponding to the available
            analysiss.

        """
        return ap.command_protocols.keys()
    #iterator functions --
    def __iter__(self):
        return self
    def next(self):
        """ Runs the analsysis as an iterator.
        The function performs the loop over the trajectory. At each frame it
        builds the necessary objects (e.g. COMFrame) and then executes the
        analysis of each analysis that was initialized in the setup.

        Args:
            nprocs (int): An integer specifying the number of cores to use in
            multithreaded parallelelization.

        """
        if self._current_frame > self._last_frame:
            raise StopIteration()
        nprocs = 1
        parallel = False
        if nprocs > 1:
            parallel = True

        # now we need to unwrap the coordinates
        natoms = self.mda_data.natoms
        oldcoord = np.zeros((natoms, 3))
        currcoord = np.zeros((natoms, 3))
        wrapcoord = np.zeros((natoms, 3))
        first_frame_coord = np.zeros((natoms, 3))
        index = self.mda_data.bilayer_sel.indices
        firstframe = self._first_frame
        first_com = self._first_com
        #print("first com: ",first_com)
        #print("first frame: ", firstframe)
        self.settings['frame_index'] = self.settings['frame_range'][0]
        #for frame in self.mda_data.mda_trajectory[
        #             self._current_frame:self._current_frame+self.frame_range[2]:self.frame_range[
        #                 2]]:
        #with self.mda_data.mda_trajectory[self._current_frame] as frame:
        #print(self._current_frame)
        #frame = self.mda_data.mda_trajectory[self._current_frame]
        for frame in self.mda_data.mda_trajectory[self._current_frame:self._current_frame+1]:
            self.reps['current_mda_frame'] = frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                first_frame_coord = np.copy(oldcoord)
                self._first_frame = False
                wrapcoord = np.copy(currcoord)
            else:
                abc = frame.dimensions[0:3]
                oldcoord = self._oldcoords
                if parallel:
                    wrapcoord = wrap_coordinates_parallel(abc, currcoord,
                                                          oldcoord,
                                                          nprocs=nprocs)
                else:
                    wrapcoord = wrap_coordinates(abc, currcoord, oldcoord)
                # frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
                # print ("wrapped coords:")
            self._oldcoords = oldcoord
            # print (wrapcoord)
            #if self.analysis_protocol.use_objects['com_frame']:
            # now build the COMFrame
            self.reps['com_frame'] = cf.COMFrame(frame, self.mda_data.bilayer_sel,
                                         wrapcoord,
                                         name_dict =
                                    self.rep_settings['com_frame']['name_dict'])
            if first_com:
                self.reps['first_com_frame'] = self.reps['com_frame']
                self._first_com = False
            # now we can assign the lipids to the leaflets
            self.reps['leaflets'] = {'upper': lf.Leaflet('upper'),
                             'lower': lf.Leaflet('lower')}
            if self.rep_settings['com_frame']['dump']:
                ofname = self.rep_settings['com_frame']['dump_path'] + "com_frame_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['com_frame'], ofile)
            if self.rep_settings['leaflets']['dump']:
                ofname = self.rep_settings['leaflets']['dump_path'] + "leaflets_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.reps['leaflets'], ofile)

            # first compute the average position along the normal direction
            zstat = RunningStats()
            for lipcom in self.reps['com_frame'].lipidcom:
                zstat.push(lipcom.com_unwrap[self.settings['norm']])
            zavg = zstat.mean()
            # now loop over the lipids
            l = 0
            for lipcom in self.reps['com_frame'].lipidcom:
                pos = ""
                # decide which leaflet
                #    print (lipcom.com_unwrap)
                #    print (lipcom.com)
                if lipcom.com_unwrap[self.settings['norm']] > zavg:
                    pos = 'upper'
                elif lipcom.com_unwrap[self.settings['norm']] < zavg:
                    pos = 'lower'
                # add to the chosen leaflet
                self.reps['com_frame'].lipidcom[l].leaflet = pos
                self.reps['leaflets'][pos].add_member(l, lipcom.type, lipcom.resid)
                l += 1

            if self.analysis_protocol.use_objects['lipid_grid']:
                self.reps['lipid_grid'] = lg.LipidGrids(self.reps['com_frame'], self.reps['leaflets'],
                                                self.settings['lateral'],
                                                nxbins=self.rep_settings['lipid_grid']['n_xbins'],
                                                nybins=self.rep_settings['lipid_grid']['n_ybins'])
                if self.rep_settings['lipid_grid']['dump']:
                    ofname = self.rep_settings['lipid_grid']['dump_path'] + "lipid_grid_" + str(
                        frame.frame) + ".pickle"
                    with open(ofname, 'wb') as ofile:
                        pickle.dump(self.reps['lipid_grid'], ofile)

                        # lipid_grid = None
            # now do analyses
            if self._frame_loop_count % self.settings['print_interval'] == 0:
                print("Frame: {}".format(frame.frame))
            i = 0
            for analysis_id in self.analysis_protocol.analysis_ids:
                if self._frame_loop_count % self.settings['print_interval'] == 0:
                    print ("  analysis: {}".format(analysis_id))
                self.analysis_protocol.command_protocol[analysis_id].run_analysis(
                    self.settings, self.reps, self.mda_data)
                # comp_out = analysis.run_analysis(self)
                # print (comp_out)
                #   analysis_out[i].append(comp_out)
                i += 1
        if self._frame_loop_count % self.settings['print_interval'] == 0:
            print(" ")
        self.settings['frame_index'] += self.settings['frame_range'][2]
        self._frame_loop_count+=1
        self._current_frame+=self.settings['frame_range'][2]
        if self._current_frame <= self._last_frame:

        #    return self._current_frame-self.frame_range[2]
            return
        #else:
         #   raise StopIteration()
