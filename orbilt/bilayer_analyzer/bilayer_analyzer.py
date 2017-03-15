


#imports
import MDAnalysis as mda
import numpy as np
import pickle
#ORBILT imports
import com_frame as cf
import leaflet as lf 
import orbilt.lipid_grid.lipid_grid as lg
import compute_protocols as cp
import plot_protocols as pp
from orbilt.common.running_stats import *
import mda_data as md

# import the coordinate wrapping function--for unwrapping
from orbilt.mda_tools.mda_unwrap import wrap_coordinates,wrap_coordinates_parallel



#the main analyzer class

class BilayerAnalyzer:
    #for the input script parser
    valid_commands = ["psf", "trajectory", "compute", "selection", "frames", "com_frame", "lipid_grid", "plot"]
    required_commands = ['psf','trajectory', 'selection']
    required_command_error_strings = {'psf':"the psf file needs to specified with command: \"psf path/psf_file_name\""}
    required_command_error_strings['trajectory']="the trajectory file needs to specified with command: \"trajectory path/trajectory_file_name\""
    required_command_error_strings['selection']="an MDAnalysis syntax selection needs to specified with command: \"selction \'selection string\'\""


    def __init__(self, psf_file = None, trajectory=None, selection=None, input_file=None):

        self.input_script_name = input_file
        if input_file is not None:
            print ("parsing input file \'"+input_file+"\'...")
            self.commands = self.parse_input_script(input_file)
        elif (psf_file is not None) and ((trajectory is not None) and (selection is not None)):
            print ("parsing inputs...")
            self.commands = {}
            self.commands['psf'] = [psf_file]
            self.commands['trajectory'] = [trajectory]  
            self.commands['selection'] = [selection]
        print ("setting up compute protocol:")
        if 'compute' in self.commands.keys():
            self.compute_protocol = cp.ComputeProtocol(self.commands['compute'])
        else:
            self.compute_protocol = cp.ComputeProtocol(default_compute_commands)
        self.print_compute_protocol()   
        print ("setting up plot protocol")
        if "plot" in self.commands.keys():
            self.plot_protocol = pp.PlotProtocol(self.commands['plot'], self.compute_protocol)
        else:
            self.plot_protocol = pp.PlotProtocol(None, self.compute_protocol)     
        #build selection string
        sel_string = ""
        for item in self.commands['selection'][0]:
            word = ""
            for char in item:
                if char != "\"":
                    word+=char
            sel_string+=" "+word
        #   print (item)
        #sel_string = "not resname CLA and not resname TIP3 and not resname POT"    
        #print "bilayer selection string:"
        #print sel_string
        print ('building the MDAnalysis objects...')
        self.mda_data = md.MDAData(self.commands['psf'][0][0],self.commands['trajectory'][0][0],sel_string)
        self.norm = 2
        self.lateral = [0,1]
        self.lateral_dimension = "xy"
        self.normal_dimension = "z"
        self.current_mda_frame = None
        #buildable objects
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

        # parse inputs for lipid_grid settings
        if 'lipid_grid' in self.commands.keys():
            lg_args = self.commands['lipid_grid']
            for i in range(0,len(lg_args),2):
                arg_key = lg_args[i]
                arg_value = lg_args[i+1]
                if arg_key == 'n_xbins':
                    self.lg_nxbins=arg_value
                elif arg_key == 'n_ybins':
                    self.lg_nybins = arg_value

        self.first_com_frame = None   
        return  

    def parse_input_script(self, input_script_name):
        
        commands = {}
        #args = {}
        with open(input_script_name) as ifile:
            for line in ifile:
                words = line.split()
                skip_line = (len(words)<=1) or ((words[0]=='#') or (words[0][0]=='#'))
                if not skip_line:
                    if words[0] in self.valid_commands:
                        if words[0] in commands.keys():
                            commands[words[0]].append(words[1:])
                        else:
                            commands[words[0]] = []
                            commands[words[0]].append(words[1:])
                    else:
                        print ("input command ",words[0]," is not a valid command")
                        error_string = "invalid input command \""+words[0]+"\""
                        raise RuntimeError(error_string)
        for required in self.required_commands:
            if required not in commands.keys():
                error_string = self.required_command_error_strings[required]
                raise RuntimeError(error_string)
                exit

        return commands

    ### compute data/access
    def print_available_computes():
        print (cp.command_protocols.keys())
        return
    def available_computes():
        return cp.command_protocols.keys()

    def print_compute_protocol(self):
        self.compute_protocol.print_protocol()
        return
    def add_compute(self, compute_string):
        self.compute_protocol.add_compute(compute_string)
        return
    def remove_compute(self, compute_id):
        self.compute_protocol.remove_compute(compute_id)
        return

    def get_compute_ids(self):
        return self.compute_protocol.compute_ids

    def get_compute_data(self, compute_id):
        return self.compute_protocol.command_protocol[compute_id].get_data()
    
    def dump_data():
        self.compute_protocol.dump_data()
        return    
    
    ## plot data/access   
    def print_available_plots():
        print (pp.command_protocols.keys())
        return

    def available_computes():
        return pp.command_protocols.keys()

    def print_plot_protocol(self):
        self.plot_protocol.print_protocol()
        return

    def add_plot(self, plot_string):
        self.plot_protocol.add_plot(compute_string)
        return

    def remove_plot(self, plot_id):
        self.plot_protocol.remove_plot(plot_id)
        return

    def print_plot_ids(self):
        print (self.plot_protocol.plot_ids)


    def get_plot_ids(self):
         return self.plot_protocol.plot_ids

    def show_plot(self, plot_id):
        self.plot_protocol.command_protocol[plot_id].show_plot(self.compute_protocol)
        return     

    def generate_plot(self, plot_id):
        self.plot_protocol.command_protocol[plot_id].generate_plot(self.compute_protocol)
        return

    def save_all_plots():
        self.plot_protocol.save_plots(self.compute_protocol)
        return

    #mda_trajectory data/access
    def update_mda_trajectory(self, new_trajectory):
        print ("updating mda trajectory to:",new_trajectory)
        self.commands['trajectory'] = [new_trajectory]
        self.mda_data.update_trajectory(new_trajectory)
        return

    #buildable objects functions
    def dump_com_frame(on=True, path="./"):
        self.dump_com_frame=on
        if path != self.dump_com_frame_path:
            self.dump_com_frame_path = path
        return

    def dump_leaflet(on=True, path="./"):
        self.dump_leaflet=on 
        if path != self.dump_leaflet_path:
            self.dump_leaflet_path = path
        return
    def dump_lipid_grid(on=True, path="./"):
        self.dump_lipid_grid=on
        if path != self.dump_lipid_grid_path:
            self.dump_lipid_grid_path=path
        return
    #analysis     
    def run_analysis(self, nprocs=1):
        parallel = False
        if nprocs>1:
            parallel=True
        #now we need to unwrap the coordinates
        natoms = self.mda_data.natoms
        oldcoord = np.zeros((natoms,3))
        currcoord = np.zeros((natoms,3))
        wrapcoord = np.zeros((natoms,3))
        first_frame_coord = np.zeros((natoms,3))
        index = self.mda_data.bilayer_sel.indices
        firstframe = True
        first_com = True
        for frame in self.mda_data.mda_trajectory:
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
                    wrapcoord = wrap_coordinates_parallel(abc, currcoord, oldcoord,nprocs=nprocs)
                else:
                    wrapcoord = wrap_coordinates(abc, currcoord, oldcoord)
                #frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
           # print ("wrapped coords:")
            #print (wrapcoord)
            if self.compute_protocol.use_objects['com_frame']:    
                #now build the COMFrame
                self.com_frame = cf.COMFrame(frame, self.mda_data.bilayer_sel, wrapcoord)
                if first_com:
                    self.first_com_frame = self.com_frame
                    first_com = False
                # now we can assign the lipids to the leaflets 
                self.leaflets = {'upper':lf.Leaflet('upper'),'lower':lf.Leaflet('lower')}
                if self.dump_com_frame:
                    ofname = self.dump_com_frame_path+"com_frame_"+str(frame.frame)+".pickle"
                    with open(ofname,'wb') as ofile:
                        pickle.dump(self.com_frame,ofile)
                if self.dump_leaflet:
                    ofname = self.dump_leaflet_path+"leaflets_"+str(frame.frame)+".pickle"
                    with open(ofname,'wb') as ofile:
                        pickle.dump(self.leaflets,ofile)

                
                #first- compute the average position along the normal direction
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
                    if lipcom.com_unwrap[self.norm]>zavg:
                        pos = 'upper'                
                    elif lipcom.com_unwrap[self.norm]<zavg:
                        pos = 'lower'
                    #add to the chosen leaflet
                    self.com_frame.lipidcom[l].leaflet = pos
                    self.leaflets[pos].add_member(l, lipcom.type)    
                    l+=1
            if self.compute_protocol.use_objects['lipid_grid']:
                self.lipid_grid = lg.LipidGrids(self.com_frame,self.leaflets,self.lateral,nxbins=self.lg_nxbins,nybins=self.lg_nybins)
                if self.dump_lipid_grid:
                    ofname = self.dump_lipid_grid_path+"lipid_grid_"+str(frame.frame)+".pickle"
                    with open(ofname,'wb') as ofile:
                        pickle.dump(self.lipid_grid,ofile)

           #lipid_grid = None        
            #now do analyses - computes
            #compute_out = []
            #for i in range(self.compute_protocol.n_commands):
            #     compute_out.append([])
            print("Frame",frame.frame)
            i = 0
            for compute_id in self.compute_protocol.compute_ids:
                print ("compute " + compute_id)
                self.compute_protocol.command_protocol[compute_id].run_compute(self) 
               # comp_out = compute.run_compute(self)
               # print (comp_out)
             #   compute_out[i].append(comp_out)
                i+=1
            print(" ")    
            #print ('compute_out:')     
            #print (compute_out)    