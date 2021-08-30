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

def wfn_plot():
    # make a plot of the computed wavefunctions
    # R. Sheehan 3 - 9 - 2020

    FUNC_NAME = ".wfn_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "Step_Solution.txt"
            
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
                args.y_label = '$\psi$(x)'
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
                args.y_label = '$\psi$(x)'
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
    
    wfn_plot()