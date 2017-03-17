import numpy as np
import vorbilt.plot_generation.plot_generation_functions as pgf 
import compute_protocols as cp

command_protocols = {}
valid_plots = []
plot_func_dict = {}
plot_compute_dict = {}


#protocol for the computes to run during the frame loop
class PlotProtocol:

    def __init__(self, plot_commands, compute_protocol):
        self.in_commands = plot_commands
        self.compute_protocol = compute_protocol
        #check computes
        arguments = []
        plot_keys = []
        command_protocol = {}
        plot_ids = []
        if plot_commands is not None:
            for command in plot_commands:
                if len(command)<2:
                    raise RuntimeError("wrong number of arguments for plot "+str(command))
                plot_key = command[0]
                plot_id = command[1]
                comp_args = command[1:]
                if (plot_key in valid_plots):
                    if (len(comp_args) >= 1):
                        if plot_id not in plot_ids:
                            arguments.append(comp_args)
                            plot_ids.append(comp_args[0])
                            plot_keys.append(plot_key)
                            command_protocol[comp_args[0]]= command_protocols[plot_key](comp_args, compute_protocol.compute_keys, compute_protocol.compute_ids)
                        else:
                            raise RuntimeError("plot id \'"+plot_id+"\' has already been used")
                    else:
                        raise RuntimeError("wrong number of arguments for plot "+str(command))
                else:
                    raise RuntimeError("invalid plot \'"+plot_key+"\' :"+str(command))   
        self.arguments = arguments
        self.plot_keys = plot_keys
        self.plot_ids = plot_ids
        self.command_protocol = command_protocol
        self.n_commands = len(command_protocol)

        return

    def add_plot(self, plot_string, compute_protocol):
        command = plot_string.split()
        plot_key = command[0]
        plot_id = command[1]
        comp_args = command[1:]
        if (plot_key in valid_plots):
            if (len(comp_args) >= 1):
                if plot_id not in self.plot_ids:
                    #comp_object = compute_obj_name_dict[plot_key]
                    #self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.plot_ids.append(comp_args[0])
                    self.plot_keys.append(plot_key)
                    self.command_protocol[comp_args[0]]= command_protocols[plot_key](comp_args, compute_protocol.compute_keys, compute_protocol.compute_ids)
                else:
                    raise RuntimeError("plot id \'"+plot_id+"\' has already been used")
            else:
                raise RuntimeError("wrong number of arguments for plot "+str(command))
        else:
            raise RuntimeError("invalid plot type \'"+plot_key+"\' :"+str(command))
        self.n_commands+=1    
        return
 
    def remove_plot(self, plot_id):
        if plot_id in self.plot_ids:
            del self.command_protocol[plot_id]
            index = self.plot_ids.index(plot_id)
            del self.arguments[index]
            del self.plot_keys[index]
            del self.plot_ids[index]
            self.n_commands-=1
        else:
            raise RuntimeWarning("no plot with id \'"+plot_id+"\'")
        return

    def print_protocol(self):
        print ("with plots:")
        if len(self.plot_ids)>0:        
            for plot_id in self.plot_ids:
                self.command_protocol[plot_id].print_protocol()
        else:
            print('None')
        return  


    def save_plots(compute_protocol):
        print ('dumping plot data to pickle files...')
        for plot_id in self.plot_ids:
            print ("plot id: "+plot_id+" ---> "+self.command_protocol[plot_id].save_file_name)
            self.command_protocol[plot_id].generate_plot(compute_protocol)


class PlotFunctionProtocol:

    def __init__(self, args):

        #required
        self.valid_args = ['none']
        self.return_length = 1
        self.plot_key = 'none'
        self.plot_id = args[0]
        #default function settings
        self.save_file_name = self.plot_id+".pickle"
        #parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)
        #storage for output
        self.my_plot = None  
        return

    #required- function to parse the input arguments    
    def parse_args(self, args):
        #print args
        nargs = len(args)
        for i in range(1,nargs,2):
            arg_key = args[i]
            arg_arg = args[i+1]
            if arg_key in self.valid_args:
                pass
            else:
                raise RuntimeError("ignoring invalid argument key "+arg_key+" for plot "+self.plot_key)
        return 
    #required - a check protocol function which reports relevant settings    
    def print_protocol(self):
        print ("Parent protocol class for computes.")
        return

    def generate_plot(self, bilayer_analyzer):
        #generate and save plot
        return

    def show_plot(self):
        #generate and show/display plot
        return



#define a new plot type 'msd'
valid_plots.append('msd')
plot_compute_dict['msd'] = 'msd'

class MSDPlotProtocol(PlotFunctionProtocol):

    def __init__(self, args, compute_keys, compute_ids):
    	''' Initialization for Plotting protocol for MSD computes
			args are the input arguments of format "plot_id compute_id_1 legend_key_1..."
    	'''
        #required
        self.plot_key = 'msd'
        self.plot_id = args[0]
        #pick up valid arguments- compute_id of the matching compute types
        self.valid_args = []
        i = 0
        for compute_key in compute_keys:
        	#print (compute_key)
        	if compute_key == plot_compute_dict[self.plot_key]:
        		self.valid_args.append(compute_ids[i])

    		i+=1

        #default function settings
        self.save_file_name = self.plot_id+".eps"
        self.interval = 1
        #parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        #storage for output
        self.my_plot = None   
        return

    #required- function to parse the input arguments    
    def parse_args(self, args):
        #print args
        nargs = len(args)
        for i in range(1,nargs,2):
            arg_key = args[i]
            arg_arg = args[i+1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError("ignoring invalid argument key "+arg_key+" for plot "+self.plot_key)
        return 

    #required - a check protocol function which reports relevant settings    
    def print_protocol(self):
        print ("Plot "+self.plot_id+" for MSD computes:"+str(self.include))
        return

    def generate_plot(self, compute_protocol):
    	data = []
    	names = []
    	for c_id in self.include:
    		data.append(compute_protocol.command_protocol[c_id].get_data())
    		names.append(self.names[c_id])

    	pgf.plot_msd(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns', interval=self.interval)	
    	return

    def show_plot(self, compute_protocol):
    	data = []
    	names = []
    	for c_id in self.include:
    		data.append(compute_protocol.command_protocol[c_id].get_data())
    		names.append(self.names[c_id])

    	pgf.plot_msd(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns',save=False,show=True, interval=self.interval)
    	return
command_protocols['msd'] = MSDPlotProtocol

#define a new plot type 'msd'
valid_plots.append('apl')
plot_compute_dict['apl'] = ['apl_box','apl_grid']

class APLPlotProtocol(PlotFunctionProtocol):

    def __init__(self, args, compute_keys, compute_ids):
        ''' Initialization for Plotting protocol for MSD computes
            args are the input arguments of format "plot_id compute_id_1 legend_key_1..."
        '''
        #required
        self.plot_key = 'apl'
        self.plot_id = args[0]
        #pick up valid arguments- compute_id of the matching compute types
        self.valid_args = []
        i = 0
        for compute_key in compute_keys:
            #print (compute_key)
            if compute_key in plot_compute_dict[self.plot_key]:
                self.valid_args.append(compute_ids[i])

            i+=1

        #default function settings
        self.save_file_name = self.plot_id+".eps"
        self.interval = 1
        #parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        #storage for output
        self.my_plot = None   
        return

    #required- function to parse the input arguments    
    def parse_args(self, args):
        #print args
        nargs = len(args)
        for i in range(1,nargs,2):
            arg_key = args[i]
            arg_arg = args[i+1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError("ignoring invalid argument key "+arg_key+" for plot "+self.plot_key)
        return 

    #required - a check protocol function which reports relevant settings    
    def print_protocol(self):
        print ("Plot "+self.plot_id+" for area per lipid (APL) computes:"+str(self.include))
        return

    def generate_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            if compute_protocol.command_protocol[c_id].compute_key == 'apl_grid':
                get_data = compute_protocol.command_protocol[c_id].get_data()
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:        
                data.append(compute_protocol.command_protocol[c_id].get_data())
                names.append(self.names[c_id])

        pgf.plot_area_per_lipid(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns', interval=self.interval)  
        return

    def show_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            get_data = compute_protocol.command_protocol[c_id].get_data()
            if compute_protocol.command_protocol[c_id].compute_key == 'apl_grid':
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:        
                data.append(get_data)
                names.append(self.names[c_id])
        pgf.plot_area_per_lipid(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns',save=False,show=True, interval=self.interval)
        return
command_protocols['apl'] = APLPlotProtocol

#define a new plot type 'disp_vec'
#  unfinished, so leave out of valid plots for now
#valid_plots.append('disp_vec')
#plot_compute_dict['disp_vec'] = ['disp_vec']

class DispVecPlotProtocol(PlotFunctionProtocol):

    def __init__(self, args, compute_keys, compute_ids):
        ''' Initialization for Plotting protocol for MSD computes
            args are the input arguments of format "plot_id compute_id_1 legend_key_1..."
        '''
        #required
        self.plot_key = 'disp_vec'
        self.plot_id = args[0]
        #pick up valid arguments- compute_id of the matching compute types
        self.valid_args = []
        i = 0
        for compute_key in compute_keys:
            #print (compute_key)
            if compute_key in plot_compute_dict[self.plot_key]:
                self.valid_args.append(compute_ids[i])

            i+=1

        #default function settings
        self.save_file_name = self.plot_id+".eps"
        self.interval = 1
        #parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        #storage for output
        self.my_plot = None   
        return

    #required- function to parse the input arguments    
    def parse_args(self, args):
        #print args
        nargs = len(args)
        for i in range(1,nargs,1):
            arg_key = args[i]
            #arg_arg = args[i+1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
            #   self.names[arg_key] = arg_arg
            else:
                raise RuntimeError("ignoring invalid argument key "+arg_key+" for plot "+self.plot_key)
        return 

    #required - a check protocol function which reports relevant settings    
    def print_protocol(self):
        print ("Plot "+self.plot_id+" for displacement vector computes:"+str(self.include))
        return

    def generate_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            if compute_protocol.command_protocol[c_id].compute_key == 'apl_grid':
                get_data = compute_protocol.command_protocol[c_id].get_data()
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:        
                data.append(compute_protocol.command_protocol[c_id].get_data())
                names.append(self.names[c_id])

        pgf.plot_area_per_lipid(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns', interval=self.interval)  
        return

    def show_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            get_data = compute_protocol.command_protocol[c_id].get_data()
            if compute_protocol.command_protocol[c_id].compute_key == 'apl_grid':
                for key in get_data.keys():
                    data.append(get_data[key])
                    names.append(key)
            else:        
                data.append(get_data)
                names.append(self.names[c_id])
        pgf.plot_area_per_lipid(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns',save=False,show=True, interval=self.interval)
        return
command_protocols['disp_vec'] = DispVecPlotProtocol

#define a new plot type 'bt'
valid_plots.append('bilayer_thickness')
plot_compute_dict['bilayer_thickness'] = ['bilayer_thickness']

class BTPlotProtocol(PlotFunctionProtocol):

    def __init__(self, args, compute_keys, compute_ids):
        ''' Initialization for Plotting protocol for MSD computes
            args are the input arguments of format "plot_id compute_id_1 legend_key_1..."
        '''
        #required
        self.plot_key = 'bilayer_thickness'
        self.plot_id = args[0]
        #pick up valid arguments- compute_id of the matching compute types
        self.valid_args = []
        i = 0
        for compute_key in compute_keys:
            #print (compute_key)
            if compute_key in plot_compute_dict[self.plot_key]:
                self.valid_args.append(compute_ids[i])

            i+=1

        #default function settings
        self.save_file_name = self.plot_id+".eps"
        self.interval = 1
        #parse input arguments if given
        self.include = []
        self.names = {}
        if len(args) > 1:
            self.parse_args(args)
        #storage for output
        self.my_plot = None   
        return

    #required- function to parse the input arguments    
    def parse_args(self, args):
        #print args
        nargs = len(args)
        for i in range(1,nargs,2):
            arg_key = args[i]
            arg_arg = args[i+1]
            if arg_key in self.valid_args:
                self.include.append(arg_key)
                self.names[arg_key] = arg_arg
            else:
                raise RuntimeError("ignoring invalid argument key "+arg_key+" for plot "+self.plot_key)
        return 

    #required - a check protocol function which reports relevant settings    
    def print_protocol(self):
        print ("Plot "+self.plot_id+" for bilayer thickness computes:"+str(self.include))
        return

    def generate_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(compute_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])
        #don't need legend for single plots    
        if len(names)<2:
            names = None
        pgf.plot_bilayer_thickness(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns', interval=self.interval)  
        return

    def show_plot(self, compute_protocol):
        data = []
        names = []
        for c_id in self.include:
            data.append(compute_protocol.command_protocol[c_id].get_data())
            names.append(self.names[c_id])
        #don't need legend for single plots    
        if len(names)<2:
            names = None
        pgf.plot_bilayer_thickness(data,name_list=names,filename=self.save_file_name,time_in='ps',time_out='ns',save=False,show=True, interval=self.interval)
        return

command_protocols['bilayer_thickness'] = BTPlotProtocol