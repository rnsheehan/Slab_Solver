# Import libraries
# You should try an import the bare minimum of modules
import sys # access system routines
import os
import glob
import re

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

# add path to our file
sys.path.append('c:/Users/robertsheehan/Programming/Python/Common/')
sys.path.append('c:/Users/robertsheehan/Programming/Python/Plotting/')

import Common
import Plotting

MOD_NAME_STR = "Plots" # use this in exception handling messages

def mode_plot():
    # make a plot of the slab waveguide mode profiles
    # R. Sheehan 3 - 9 - 2020

    FUNC_NAME = ".mode_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #filename = "TE_Mode_Profiles.txt"
        #filename = "Coupled_Mode_Profiles.txt"
        #filename = "integrate_KAA_field_values.txt"
        #filename = "integrate_KBB_field_values.txt"
        #filename = "integrate_KAB_field_values.txt"
        #filename = "integrate_KBA_field_values.txt"
        filename = "TE_Mode_Profiles_FL_A.txt"
        #filename = "TM_Disp_Eqns_FL_A.txt"
            
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            
            if len(data) > 2:      
                # multi-curve plot required
                hv_data = []; labels = []; marks = [];
                for i in range(1, len(data), 1):
                    hv_data.append([data[0], data[i]]); 
                    labels.append('M$_{%(v1)d}$'%{"v1":i}); 
                    marks.append(Plotting.labs_lins[(i-1)%len(Plotting.labs_lins)]); 
           
                # make the plot of the data set
                args = Plotting.plot_arg_multiple()

                args.loud = True
                args.crv_lab_list = labels
                args.mrk_list = marks
                args.x_label = 'Position ($\mu$m)'
                args.y_label = 'Field Value'
                args.fig_name = filename.replace('.txt','')

                Plotting.plot_multiple_curves(hv_data, args)

                del hv_data; del labels; del marks; 
            else:
                # single curve plot required
                
                args = Plotting.plot_arg_single()
                
                args.loud = True
                args.curve_label = 'M$_{1}$'
                args.marker = Plotting.labs_lins[0]
                args.x_label = 'Position ($\mu$m)'
                args.y_label = 'Field Value'
                args.fig_name = filename.replace('.txt','')
                
                Plotting.plot_single_curve(data[0], data[1], args)
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)
        
def disp_eqn_plot():
    # make a plot of the slab waveguide dispersion equations
    # R. Sheehan 3 - 9 - 2020

    FUNC_NAME = ".disp_eqn_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #filename = "TE_Dispersion_Eqn.txt"
        filename = "TE_Disp_Eqns_FL_AA.txt"
            
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            
            # single curve plot required
                
            args = Plotting.plot_arg_single()
            
            args.loud = True
            #args.curve_label = '$\beta$'
            args.marker = Plotting.labs_lins[0]
            args.x_label = 'Search Space'
            args.y_label = 'Disp. Eqn.'
            args.fig_name = filename.replace('.txt','')
            
            Plotting.plot_single_curve(data[0], data[1], args)
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)
        
def chuang_plots():
    # make plots for the computed data from the Chuang benchmark calculation
    # R. Sheehan 16 - 10 - 2020
    
    FUNC_NAME = ".chuang_plots()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
    
    try:
        filename = "Chuang_Benchmark.txt"
        
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            
            # plot the overlap integral versus RI difference
            args = Plotting.plot_arg_single()
            
            args.loud = True
            #args.curve_label = '$\beta$'
            args.marker = Plotting.labs_pts[0]
            args.x_label = 'RI Difference $\Delta n$'
            args.y_label = 'Overlap Integral $C_{AB}$'
            args.fig_name = filename.replace('.txt','') + '_OL'
            
            Plotting.plot_single_curve(data[0], data[1], args)
            
            # plot the asynchronism versus RI difference
            args.loud = True
            #args.curve_label = '$\beta$'
            args.marker = Plotting.labs_pts[1]
            args.x_label = 'RI Difference $\Delta n$'
            args.y_label = 'Asynchronism'
            args.fig_name = filename.replace('.txt','') + '_Async'
            
            Plotting.plot_single_curve(data[0], data[4], args)
            
            # multi-curve plot required
            hv_data = []; labels = []; marks = [];
            
            hv_data.append([data[0], data[2]]); labels.append('$k_{ab}$'); marks.append(Plotting.labs_pts[2]); 
            hv_data.append([data[0], data[3]]); labels.append('$k_{ba}$'); marks.append(Plotting.labs_pts[3]); 
       
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'RI Difference $\Delta n$'
            args.y_label = 'Coupling Coefficients'
            args.fig_name = filename.replace('.txt','') + '_CC'

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)
        
def coupling_ampl():
    # make of the coupling amplitudes for each waveguide
    # R. Sheehan 17 - 10 - 2020
    
    FUNC_NAME = ".coupling_ampl()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
    
    try:
        filename = "Coupled_Mode_Coefficients.txt"
        
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            
            # multi-curve plot required
            hv_data = []; labels = []; marks = [];
            
            hv_data.append([data[0], data[1]]); labels.append('|a(z)|'); marks.append(Plotting.labs_lins[2]); 
            hv_data.append([data[0], data[2]]); labels.append('|b(z)|'); marks.append(Plotting.labs_lins[3]); 
            
            print("max(a(z)) occurs at z = ",data[0][np.argmax(data[1])])
            print("min(b(z)) occurs at z = ",data[0][np.argmin(data[2])])
       
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Waveguide Length um'
            args.y_label = 'Coupling Amplitudes'
            args.fig_name = filename.replace('.txt','') + '_CA'

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)
        
def coupled_fields(stepsize, stepnum):
    # make a plot of the coupled fields as they propagate
    # R. Sheehan 18 - 10 - 2020
    
    FUNC_NAME = ".coupled_fields()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
    
    try:
        filename = "Coupled_Field_dz_%(v1)d_step_%(v2)s.txt"%{"v1":stepsize, "v2":stepnum}
        
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            
            # plot the overlap integral versus RI difference
            args = Plotting.plot_arg_single()
            
            args.loud = True
            #args.curve_label = '$\beta$'
            args.marker = Plotting.labs_lins[1]
            args.x_label = 'Position ($\mu$m)'
            args.y_label = 'Field Value'
            args.plt_title = 'z = %(v1)d $\mu$m'%{"v1":(stepsize*stepnum)}
            args.fig_name = filename.replace('.txt','')
            
            Plotting.plot_single_curve(data[0], data[1], args)
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory
    
    print(pwd)
    
    disp_eqn_plot()
    
    mode_plot()
    
    #chuang_plots()
    
    #coupling_ampl()
    
    #for i in range(0, 120, 10):
    #    coupled_fields(10, i)
