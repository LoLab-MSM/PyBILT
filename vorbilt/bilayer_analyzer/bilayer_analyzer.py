"""The bilayer_analyzer module.

This module defines the BilayerAnalyzer class that can be used to build protocols for analyzing and generating plots
for the analysis of a lipid bilayer trajectory.

Example:
    >> import vorbilt.bilayer_analyzer as ba
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
# VORBILT imports
import com_frame as cf
import leaflet as lf
import vorbilt.lipid_grid.lipid_grid as lg
import analysis_protocols as ap
import plot_protocols as pp
from vorbilt.common.running_stats import *
import mda_data as md

# import the coordinate wrapping function--for unwrapping
from vorbilt.mda_tools.mda_unwrap import wrap_coordinates, \
    wrap_coordinates_parallel

default_analysis_commands = [['msd', 'msd_1']]


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
    valid_commands = ["psf", "trajectory", "analysis", "selection", "frames",
                      "plot"]
    required_commands = ['psf', 'trajectory', 'selection']
    required_command_error_strings = {'psf': "the psf file needs to specified with command: \"psf path/psf_file_name\""}
    required_command_error_strings[
        'trajectory'] = "the trajectory file needs to specified with command: \"trajectory path/trajectory_file_name\""
    required_command_error_strings[
        'selection'] = "an MDAnalysis syntax selection needs to specified with command: \"selction \'selection string\'\""

    def __init__(self, psf_file=None, trajectory=None, selection=None,
                 input_file=None):
        """BilayerAnalyzer initialization.
        The initialization initialially parses all the inputs and sets some the attribute values. It also calls the
        initialization of the Analyses and PlotProtocol objects and builds the MDAData object.
        Args:
            psf_file (str): Optional, the path and filename of the structure file.
            trajectory (str;list): Optional, the path and filename of the trajectory file. Also accepts a list of
                filenames.
            selection (str): Optional, the MDAnalysis compatible string to select the bilayer components.
            input_file (str): Optional, the path and filename of input setup file.
        """
        self.input_script_name = input_file
        if input_file is not None:
            print ("parsing input file \'" + input_file + "\'...")
            self.commands = self.parse_input_script(input_file)
        elif (psf_file is not None) and ((trajectory is not None) and
                                             (selection is not None)):
            print ("parsing inputs...")
            self.commands = dict()
            self.commands['psf'] = psf_file
            self.commands['trajectory'] = trajectory
            self.commands['selection'] = selection
        else:
            error = 'Must provide input_file or all three options for' \
                    '"psf_file", "trajectory", and "selection"'
            raise RuntimeError(error)

            # set up the analysis protocol
        print ("setting up analysis protocol:")
        if 'analysis' in self.commands.keys():
            self.analysis_protocol = ap.Analyses(
                self.commands['analysis'])
        else:
            self.analysis_protocol = ap.Analyses(
                default_analysis_commands)
        self.print_analysis_protocol()
        # set up the plot protocol
        print ("setting up plot protocol")
        if "plot" in self.commands.keys():
            self.plot_protocol = pp.PlotProtocol(self.commands['plot'],
                                                 self.analysis_protocol)
        else:
            self.plot_protocol = pp.PlotProtocol(None, self.analysis_protocol)
        for i in self.commands:
            print(i, self.commands[i])
        # build selection string for the MDAData object
        sel_string = self.commands['selection']

        # print (item)
        # sel_string = "not resname CLA and not resname TIP3 and not resname POT"
        # print "bilayer selection string:"
        # print sel_string
        # build the MDAData object
        print ('building the MDAnalysis objects...')

        self.mda_data = md.MDAData(self.commands['psf'],
                                   self.commands['trajectory'],
                                   sel_string)
        # analysis defaults
        self.norm = 2
        self.lateral = [0, 1]
        self.lateral_dimension = "xy"
        self.normal_dimension = "z"
        self.current_mda_frame = None

        # settings for frame loop
        # frame_range[0]=first,frame_range[1]=last,frame_range[2]=interval
        self.frame_range = [0, -1, 1]
        self.frame_index = 0
        # buildable objects
        # com frame
        self.com_frame = None
        self.dump_com_frame = False
        self.dump_com_frame_path = "./"
        # leaflet
        self.leaflets = None
        self.dump_leaflet = False
        self.dump_leaflet_path = "./"
        # lipid grid- with default settings
        self.lipid_grid = None
        self.lg_nxbins = 10
        self.lg_nybins = 10
        self.dump_lipid_grid = False
        self.dump_lipid_grid_path = "./"

        # parse 'frames' input key--
        # for setting the frame range of the analysis
        if 'frames' in self.commands.keys():
            f_args = self.commands['frames']
            for i in range(0, len(f_args), 2):
                arg_key = f_args[i]
                arg_value = f_args[i + 1]
                if arg_key == 'first':
                    self.frame_range[0] = arg_value
                elif arg_key == 'last':
                    self.frame_range[1] = arg_value
                elif arg_key == 'interval':
                    self.frame_range[2] = arg_value

                    # parse inputs for lipid_grid settings
        if 'lipid_grid' in self.commands.keys():
            lg_args = self.commands['lipid_grid']
            for i in range(0, len(lg_args), 2):
                arg_key = lg_args[i]
                arg_value = lg_args[i + 1]
                if arg_key == 'n_xbins':
                    self.lg_nxbins = arg_value
                elif arg_key == 'n_ybins':
                    self.lg_nybins = arg_value

        self.first_com_frame = None
        return

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
                if key in ['psf', 'trajectory', 'selection']:
                    second_term = line.replace(key, '', 1).strip()
                    if key == 'selection':
                        second_term = eval(second_term)
                    commands[key] = second_term
                else:
                    if key in self.valid_commands:
                        if key in commands.keys():
                            commands[key].append(words[1:])
                        else:
                            commands[key] = []
                            commands[key].append(words[1:])
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

    def add_analysis(self, analysis_string):
        """Add a analysis to the analysis protocol.
        Args:
            analysis_string (str): A string defining the analysis key, analysis
             id, and arguments for the new analysis.

        """
        self.analysis_protocol.add_analysis(analysis_string)
        return

    def remove_analysis(self, analysis_id):
        """Remove a specified analysis from the comptue protocol.
        Args:
            analysis_id (str): The string analysis id of the analysis to be
             removed from the protocol.

        """
        self.analysis_protocol.remove_analysis(analysis_id)
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

    def dump_data(self):
        """Dump all the anlysis outputs from the analysiss as pickle files."""
        self.analysis_protocol.dump_data()
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
        self.commands['trajectory'] = [new_trajectory]
        self.mda_data.update_trajectory(new_trajectory)
        return

    # buildable objects functions
    def dump_com_frame(self, on=True, path="./"):
        """Turn on (or off) the dumping of COMFrame objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
             dumped.
                Default: "./"

        """
        self.dump_com_frame = on
        if path != self.dump_com_frame_path:
            self.dump_com_frame_path = path
        return

    def dump_leaflet(self, on=True, path="./"):
        """Turn on (or off) the dumping of leaflet objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
             dumped.
                Default: "./"

        """
        self.dump_leaflet = on
        if path != self.dump_leaflet_path:
            self.dump_leaflet_path = path
        return

    def dump_lipid_grid(self, on=True, path="./"):
        """Turn on (or off) the dumping of LipidGrid objects built during analysis.
        Args:
            on (bool, optional): Defines whether the COMFrame objects should be
             dumped (as pickle files) after each frame in the analysis loop.
                Default: True
            path (str, optional): Sets the path where the COMFrame objects are
            dumped.
                Default: "./"

        """
        self.dump_lipid_grid = on
        if path != self.dump_lipid_grid_path:
            self.dump_lipid_grid_path = path
        return

    def set_frame_range(self, first=0, last=-1, interval=1):
        """Set the frame range and interval to use in the analysis.
        Args:
            first (int): The first frame.
            last (int): The last frame.
            interval (int): The interval to skip between frames in the analysis
             loop.

        """
        if first != self.frame_range[0]:
            self.frame_range[0] = first
        if last != self.frame_range[1]:
            self.frame_range[1] = last
        if interval != self.frame_range[2]:
            self.frame_range[2] = interval
        return

    def reset(self):
        """ Clears the analysis output stored in the analysiss.
        """
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
        self.frame_index = self.frame_range[0]
        for frame in self.mda_data.mda_trajectory[
                     self.frame_range[0]:self.frame_range[1]:self.frame_range[
                         2]]:
            self.current_mda_frame = frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                first_frame_coord = np.copy(oldcoord)
                firstframe = False
                wrapcoord = np.copy(currcoord)
            else:
                abc = frame.dimensions[0:3]
                if parallel:
                    wrapcoord = wrap_coordinates_parallel(abc, currcoord,
                                                          oldcoord,
                                                          nprocs=nprocs)
                else:
                    wrapcoord = wrap_coordinates(abc, currcoord, oldcoord)
                # frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
                # print ("wrapped coords:")
            # print (wrapcoord)
            #if self.analysis_protocol.use_objects['com_frame']:
            # now build the COMFrame
            self.com_frame = cf.COMFrame(frame, self.mda_data.bilayer_sel,
                                         wrapcoord)
            if first_com:
                self.first_com_frame = self.com_frame
                first_com = False
            # now we can assign the lipids to the leaflets
            self.leaflets = {'upper': lf.Leaflet('upper'),
                             'lower': lf.Leaflet('lower')}
            if self.dump_com_frame:
                ofname = self.dump_com_frame_path + "com_frame_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.com_frame, ofile)
            if self.dump_leaflet:
                ofname = self.dump_leaflet_path + "leaflets_" + str(
                    frame.frame) + ".pickle"
                with open(ofname, 'wb') as ofile:
                    pickle.dump(self.leaflets, ofile)

            # first compute the average position along the normal direction
            zstat = RunningStats()
            for lipcom in self.com_frame.lipidcom:
                zstat.push(lipcom.com_unwrap[self.norm])
            zavg = zstat.mean()
            # now loop over the lipids
            l = 0
            for lipcom in self.com_frame.lipidcom:
                pos = ""
                # decide which leaflet
                #    print (lipcom.com_unwrap)
                #    print (lipcom.com)
                if lipcom.com_unwrap[self.norm] > zavg:
                    pos = 'upper'
                elif lipcom.com_unwrap[self.norm] < zavg:
                    pos = 'lower'
                # add to the chosen leaflet
                self.com_frame.lipidcom[l].leaflet = pos
                self.leaflets[pos].add_member(l, lipcom.type, lipcom.resid)
                l += 1

            if self.analysis_protocol.use_objects['lipid_grid']:
                self.lipid_grid = lg.LipidGrids(self.com_frame, self.leaflets,
                                                self.lateral,
                                                nxbins=self.lg_nxbins,
                                                nybins=self.lg_nybins)
                if self.dump_lipid_grid:
                    ofname = self.dump_lipid_grid_path + "lipid_grid_" + str(
                        frame.frame) + ".pickle"
                    with open(ofname, 'wb') as ofile:
                        pickle.dump(self.lipid_grid, ofile)

                        # lipid_grid = None
            # now do analyses

            print("Frame", frame.frame)
            i = 0
            for analysis_id in self.analysis_protocol.analysis_ids:
                print ("analysis " + analysis_id)
                self.analysis_protocol.command_protocol[analysis_id].run_analysis(
                    self)
                # comp_out = analysis.run_analysis(self)
                # print (comp_out)
                #   analysis_out[i].append(comp_out)
                i += 1
            print(" ")
            self.frame_index += self.frame_range[2]
            # print ('analysis_out:')
            # print (analysis_out)

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