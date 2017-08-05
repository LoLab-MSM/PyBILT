import pybilt.plot_generation.plot_generation_functions as pgf


command_protocols = {}
valid_plots = []
plot_func_dict = {}
plot_analysis_dict = {}


# protocol for the  plots to run after the frame loop
class PlotProtocol(object):
    def __init__(self, plot_commands, analysis_protocol):
        self.in_commands = plot_commands
        self.analysis_protocol = analysis_protocol
        # check analysis
        arguments = []
        plot_keys = []
        command_protocol = {}
        plot_ids = []
        if plot_commands is not None:
            for command in plot_commands:
                command = command.split()
                if len(command) < 2:
                    raise RuntimeError(
                        "wrong number of arguments for plot " + str(command))
                plot_key = command[0]
                plot_id = command[1]
                comp_args = command[1:]
                if plot_key in valid_plots:
                    if len(comp_args) >= 1:
                        if plot_id not in plot_ids:
                            arguments.append(comp_args)
                            plot_ids.append(comp_args[0])
                            plot_keys.append(plot_key)
                            command_protocol[comp_args[0]] = command_protocols[
                                plot_key](comp_args,
                                          analysis_protocol.analysis_keys,
                                          analysis_protocol.analysis_ids)
                        else:
                            raise RuntimeError("plot id '{}' has already "
                                               "been used".format(plot_id))
                    else:
                        raise RuntimeError("wrong number of arguments for plot"
                                           " {}".format(command))
                else:
                    raise RuntimeError("invalid plot type '{}' : "
                                       "{}".format(plot_key, command))
        self.arguments = arguments
        self.plot_keys = plot_keys
        self.plot_ids = plot_ids
        self.command_protocol = command_protocol
        self.n_commands = len(command_protocol)

        return

    def add_plot(self, plot_string, analysis_protocol):
        command = plot_string.split()
        plot_key = command[0]
        plot_id = command[1]
        comp_args = command[1:]
        if plot_key in valid_plots:
            if len(comp_args) >= 1:
                if plot_id not in self.plot_ids:
                    # comp_object = analysis_obj_name_dict[plot_key]
                    # self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.plot_ids.append(comp_args[0])
                    self.plot_keys.append(plot_key)
                    self.command_protocol[comp_args[0]] = command_protocols[
                        plot_key](comp_args, analysis_protocol.analysis_keys,
                                  analysis_protocol.analysis_ids)
                else:
                    raise RuntimeError("plot id '{}' has already "
                                       "been used".format(plot_id))
            else:
                raise RuntimeError("wrong number of arguments for plot"
                                   " {}".format(command))
        else:
            raise RuntimeError("invalid plot type '{}' : "
                               "{}".format(plot_key, command))
        self.n_commands += 1
        return

    def remove_plot(self, plot_id):
        if plot_id in self.plot_ids:
            del self.command_protocol[plot_id]
            index = self.plot_ids.index(plot_id)
            del self.arguments[index]
            del self.plot_keys[index]
            del self.plot_ids[index]
            self.n_commands -= 1
        else:
            raise RuntimeWarning("no plot with id '{}'".format(plot_id))
        return

    def print_protocol(self):
        print ("with plots:")
        if len(self.plot_ids) > 0:
            for plot_id in self.plot_ids:
                self.command_protocol[plot_id].print_protocol()
        else:
            print('None')
        return

    def save_plots(self, analysis_protocol):
        print ('dumping plot data to pickle files...')
        for plot_id in self.plot_ids:
            print ("plot id: " + plot_id + " ---> " + self.command_protocol[
                plot_id].save_file_name)
            self.command_protocol[plot_id].generate_plot(analysis_protocol)


class PlotFunctionProtocol(object):
    def __init__(self, args):

        # required
        self.valid_args = ['none']
        self.return_length = 1
        self.plot_key = 'none'
        self.plot_id = args[0]
        # default function settings
        self.save_file_name = self.plot_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.my_plot = None
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            dummy_arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                pass
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for plot " + self.plot_key)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("Parent protocol class for analysis.")
        return

    def generate_plot(self, bilayer_analyzer):
        # generate and save plot
        print(bilayer_analyzer)
        return

    def show_plot(self, analysis_protocol):
        print(analysis_protocol)
        # generate and show/display plot
        return


# define a new plot type 'msd'
valid_plots.append('msd')
plot_analysis_dict['msd'] = 'msd'


class MSDPlotProtocol(PlotFunctionProtocol):
    def __init__(self, args, analysis_keys, analysis_ids):
        ''' Initialization for Plotting protocol for MSD analysis
            args are the input arguments of format "plot_id analysis_id_1 legend_key_1..."
        '''
        # required
        self.plot_key = 'msd'
        self.plot_id = args[0]
        # pick up valid arguments- analysis_id of the matching analysis types
        self.valid_args = []
        i = 0
        for analysis_key in analysis_keys:
            # print (analysis_key)
            if analysis_key == plot_analysis_dict[self.plot_key]:
                self.valid_args.append(analysis_ids[i])

            i += 1

        # default function settings
        self.save_file_name = self.plot_id + ".eps"
        self.interval = 1
        # parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.my_plot = None
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for plot " + self.plot_key)
        return

        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
        "Plot " + self.plot_id + " for MSD analysis:" + str(self.include))
        return

    def generate_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(analysis_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])

        pgf.plot_msd(data, name_list=names, filename=self.save_file_name,
                     time_in='ps', time_out='ns', interval=self.interval)
        return

    def show_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(analysis_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])

        pgf.plot_msd(data, name_list=names, filename=self.save_file_name,
                     time_in='ps', time_out='ns', save=False, show=True,
                     interval=self.interval)
        return


command_protocols['msd'] = MSDPlotProtocol

# define a new plot type 'msd'
valid_plots.append('apl')
plot_analysis_dict['apl'] = ['apl_box', 'apl_grid']


class APLPlotProtocol(PlotFunctionProtocol):
    def __init__(self, args, analysis_keys, analysis_ids):
        ''' Initialization for Plotting protocol for MSD analysis
            args are the input arguments of format "plot_id analysis_id_1 legend_key_1..."
        '''
        # required
        self.plot_key = 'apl'
        self.plot_id = args[0]
        # pick up valid arguments- analysis_id of the matching analysis types
        self.valid_args = []
        i = 0
        for analysis_key in analysis_keys:
            # print (analysis_key)
            if analysis_key in plot_analysis_dict[self.plot_key]:
                self.valid_args.append(analysis_ids[i])

            i += 1

        # default function settings
        self.save_file_name = self.plot_id + ".eps"
        self.interval = 1
        # parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.my_plot = None
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for plot " + self.plot_key)
        return

        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
        "Plot " + self.plot_id + " for area per lipid (APL) analysis:" + str(
            self.include))
        return

    def generate_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            if analysis_protocol.command_protocol[c_id].analysis_key == 'apl_grid':
                get_data = analysis_protocol.command_protocol[c_id].get_data()
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:
                data.append(analysis_protocol.command_protocol[c_id].get_data())
                names.append(self.names[c_id])

        pgf.plot_area_per_lipid(data, name_list=names,
                                filename=self.save_file_name, time_in='ps',
                                time_out='ns', interval=self.interval)
        return

    def show_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            get_data = analysis_protocol.command_protocol[c_id].get_data()
            if analysis_protocol.command_protocol[c_id].analysis_key == 'apl_grid':
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:
                data.append(get_data)
                names.append(self.names[c_id])
        pgf.plot_area_per_lipid(data, name_list=names,
                                filename=self.save_file_name, time_in='ps',
                                time_out='ns', save=False, show=True,
                                interval=self.interval)
        return


command_protocols['apl'] = APLPlotProtocol


# define a new plot type 'disp_vec'
#  unfinished, so leave out of valid plots for now
# valid_plots.append('disp_vec')
# plot_analysis_dict['disp_vec'] = ['disp_vec']

class DispVecPlotProtocol(PlotFunctionProtocol):
    def __init__(self, args, analysis_keys, analysis_ids):
        ''' Initialization for Plotting protocol for MSD analysis
            args are the input arguments of format "plot_id analysis_id_1 legend_key_1..."
        '''
        # required
        self.plot_key = 'disp_vec'
        self.plot_id = args[0]
        # pick up valid arguments- analysis_id of the matching analysis types
        self.valid_args = []
        i = 0
        for analysis_key in analysis_keys:
            # print (analysis_key)
            if analysis_key in plot_analysis_dict[self.plot_key]:
                self.valid_args.append(analysis_ids[i])

            i += 1

        # default function settings
        self.save_file_name = self.plot_id + ".eps"
        self.interval = 1
        # parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.my_plot = None
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 1):
            arg_key = args[i]
            # arg_arg = args[i+1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
            # self.names[arg_key] = arg_arg
            else:
                raise RuntimeError("ignoring invalid argument key {} for "
                                   "plot {}".format(arg_key, self.plot_key))
        return

        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("Plot {} for displacement vector "
               "analysis: {}".format(self.plot_id, self.include))
        return

    def generate_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            if analysis_protocol.command_protocol[c_id].analysis_key == 'apl_grid':
                get_data = analysis_protocol.command_protocol[c_id].get_data()
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:
                data.append(analysis_protocol.command_protocol[c_id].get_data())
                names.append(self.names[c_id])

        pgf.plot_area_per_lipid(data, name_list=names,
                                filename=self.save_file_name, time_in='ps',
                                time_out='ns', interval=self.interval)
        return

    def show_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            get_data = analysis_protocol.command_protocol[c_id].get_data()
            if analysis_protocol.command_protocol[
                c_id].analysis_key == 'apl_grid':
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:
                data.append(get_data)
                names.append(self.names[c_id])
        pgf.plot_area_per_lipid(data, name_list=names,
                                filename=self.save_file_name, time_in='ps',
                                time_out='ns', save=False, show=True,
                                interval=self.interval)
        return


command_protocols['disp_vec'] = DispVecPlotProtocol

# define a new plot type 'bt'
valid_plots.append('bilayer_thickness')
plot_analysis_dict['bilayer_thickness'] = ['bilayer_thickness']


class BTPlotProtocol(PlotFunctionProtocol):
    def __init__(self, args, analysis_keys, analysis_ids):
        ''' Initialization for Plotting protocol for MSD analysis
            args are the input arguments of format "plot_id analysis_id_1 legend_key_1..."
        '''
        # required
        self.plot_key = 'bilayer_thickness'
        self.plot_id = args[0]
        # pick up valid arguments- analysis_id of the matching analysis types
        self.valid_args = []
        i = 0
        for analysis_key in analysis_keys:
            # print (analysis_key)
            if analysis_key in plot_analysis_dict[self.plot_key]:
                self.valid_args.append(analysis_ids[i])

            i += 1

        # default function settings
        self.save_file_name = self.plot_id + ".eps"
        self.interval = 1
        # parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.my_plot = None
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs-1, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for plot " + self.plot_key)
        return

        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("Plot {} for bilayer thickness "
               "analysis: {}".format(self.plot_id, self.include))
        return

    def generate_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(analysis_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])
        # don't need legend for single plots
        if len(names) < 2:
            names = None
        pgf.plot_bilayer_thickness(data, name_list=names,
                                   filename=self.save_file_name, time_in='ps',
                                   time_out='ns', interval=self.interval)
        return

    def show_plot(self, analysis_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(analysis_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])
        # don't need legend for single plots
        if len(names) < 2:
            names = None
        pgf.plot_bilayer_thickness(data, name_list=names,
                                   filename=self.save_file_name, time_in='ps',
                                   time_out='ns', save=False, show=True,
                                   interval=self.interval)
        return


command_protocols['bilayer_thickness'] = BTPlotProtocol
