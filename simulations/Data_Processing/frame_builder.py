'''
This script visualizes the data in the form of a video. It reads the data from the csv file and creates a video of the data.

The script is divided into the following sections:
1. Importing the required libraries
2. Reading the data from the csv file
   The data is stored in a directory structure as follows:
   {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
   where terms in {} are variables that are provided as input to the script.
   In these csvs, the columns are as follows:
   a_c	  x	  P(x; t)	 G(x; t)	 Pr(x; t)	 W(x; t)	 O(x; t) 
   where the 2nd column onwards (indexed from 0) are the data columns and represent the species concentration (for each species) at each grid point.
   The columns are of length {g*g} (excluding the header row) where g is the grid size.

3. To visualise the data, the script reads the data from the csv files and stores them in a Pandas DataFrame.
   The data for each frame is stored in a separate DataFrame, and is visualised by reshaping the data into a 2D grid.
   The data is then plotted using Seaborn's heatmap function.
   Two types of plots are created:
    a. A plot of the concentration of each species.
    b. Subplots of the concentration of each species.
    These plots are saved as images in the {out_dir} and {outdir/Combined} directory as png files.

4. The images of the individual plots and the subplots are then seperately combined to create two videos using the cv2.videoWriter. 
   The videos are saved in the {out_dir/Videos} directory.
5. Unless otherwise specified, the images are then deleted from the {out_dir} and {outdir/Combined} directories.

'''

# Importing the required libraries
import os
import math
import pandas as pan
import numpy as np
import seaborn as sea
import matplotlib.pyplot as plt
import cv2
import shutil
import argparse
import glob
import regex as re
from pathlib import Path

import matplotlib.ticker as mtick
import matplotlib.colors as mcolours
from matplotlib.collections import LineCollection

# User inputs.

in_dir = "StdParam_20_100_Test/"
out_dir = "../../Images/AdvDiffTest/"
PREFIX = "DiC-NEW"
out_dir = "../../Images/AdvDiffTest/" + PREFIX + "/"
g = 256;  dP = 11000; Geq = 4.802; R_max= 1;
T_vals =[]
a_vals = [0.043]    

#Making the output directories if they don't exist
Path(out_dir + 'Videos/').mkdir(parents=True, exist_ok=True)
Path(out_dir + 'Combined/').mkdir(parents=True, exist_ok=True)
Path(out_dir + 'Singular/').mkdir(parents=True, exist_ok=True)

hex_list= [['D8F3DC', 'B7E4C7', '95D5B2', '74C69D', '57CC99', '38A3A5', '006466', '0B525B', '1B4332'], 
           ['B7094C', 'A01A58', '892B64', '723C70', '5C4D7D', '455E89', '2E6F95', '1780A1', '0091AD'], 
           ['B80068', "BB0588", "CC00AA", "CC00CC", "BC00DD", "B100E8", "A100F2", "8900F2", '6A00F4'],
           ["c4fff9", "9ceaef", "68d8d6", "3dccc7", "07beb8"], 
           ['CAF0F8', '0096C7']]

float_list = [0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1]
vmax =[500, 35, 40]

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolours.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def home_video(dir, savepath, a_vals, T_vals, Rmax, pngformat= "Rho_s_{s}_a_{a}_T_{T}_n_{R}.png", delpng=True):
    '''
    This function creates a video of the data in the directory dir.
    The data is stored in the csv files in the directory structure as follows:

    '''
    return



def frame_visualiser(delpng = False):
    # If a_vals and T_vals are not provided, the script reads the relevant text files and extracts the a_vals and T_vals from them.
    # List of a_vals provided {in_dir}/{PREFIX}/a_vals.txt
    # List of T_vals provided {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/T_vals.txt

    global a_vals, T_vals
    if a_vals == []:
        try:
            with open(in_dir + PREFIX + "/a_vals.txt", 'r') as f:
                a_vals = [float(line.strip()) for line in f]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory and no a_vals specified. Terminating the program.")
            return
    print(f"List of a values: {a_vals}")
    for a in a_vals:
        # Convert a and T to int format if they are integers.
        a = int(a) if a.is_integer() else a
        if T_vals == []:
            try:
                with open(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt", 'r') as f:
                    T_vals = [float(line.strip()) for line in f]
            except FileNotFoundError:
                print(f"T_vals.txt not found in the input directory and no T_vals specified for a = {a}. Skipping this a-value.")
                continue
        print(f"List of T values: {T_vals}")
        for T in T_vals:
            # Convert a and T to int format if they are integers. 
            T = int(T) if T.is_integer() else T
            print(f"Processing a = {a}, T = {T} at dP = {dP}, Geq = {Geq}")
            for R in range(0, R_max):
            # Reading the data from the csv file
                # Convert a and T to int format if they are integers.

                data = pan.read_csv(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a}_R_{R}.csv", header= 0)
                # Extracting the data columns
                data = data.iloc[:, 2:]
                fig, axs = plt.subplots(1, 3, figsize = (20, 6))
                '''
                for s in range(3):
                    data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                    # Plot data first individually and then add to subplot.
                    fig, ax =plt.subplots(figsize=(9, 6))
                    ax = sea.heatmap(data_Sp, ax=ax, cmap=get_continuous_cmap(hex_list[s]), cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    plt.savefig(out_dir + f"Singular/Rho_s_{s}_a_{a}_T_{T}_n_{R}.png")
                    plt.close()
                '''
                # Creating subplots
                fig, axs = plt.subplots(1, 3, figsize = (20, 6))
                for s in range(3):
                    ax = axs[s]
                    data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                    ax = sea.heatmap(data_Sp , ax=ax, cmap=get_continuous_cmap(hex_list[s], float_list),  cbar=True, vmax = vmax[s])
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                plt.savefig(out_dir + f"Combined/Concentration_a_{a}_T_{T}_n_{R}.png")
                plt.close()


                data_gamma = pan.read_csv(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/GAMMA_T_{T}_a_{a}_R_{R}.csv", header= 0)
                data_gamma = data_gamma.iloc[:, 2:]
                fig, axs = plt.subplots(2, 3, figsize = (20, 12))

                for s in range(3):
                    ax1 = axs[0, s]
                    ax2 = axs[1, s]
                    data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                    data_GAM = np.array(data_gamma.iloc[:, s]).reshape(g, g)
                    ax1 = sea.heatmap(data_Sp, ax=ax1, cmap=get_continuous_cmap(hex_list[s], float_list), cbar=True, vmax = vmax[s])
                    ax1.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax1.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax1.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax2 = sea.heatmap(data_GAM, ax=ax2, cmap=get_continuous_cmap(hex_list[4]), cbar=True, vmin=0, vmax = 1)
                    ax2.set_title(r'$\gamma_%g(t = %g)$,' % (s, T))
                    ax2.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax2.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))

                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                plt.savefig(out_dir + f"Combined/Gamma_a_{a}_T_{T}_n_{R}.png")
                plt.close()
            
            # End of R loop
        # End of T loop

        # Combining images to create video for each a value.
        # Images are arranged to create video in increasing order of T, followed by increasing order of R.
        img_array = []
        '''
        for s in range(3):
            for R in range(R_max):
            
                for T in T_vals:
                    T = int(T) if T.is_integer() else T
                
                    try:
                        img = cv2.imread(out_dir + f"Singular/Rho_s_{s}_a_{a}_T_{T}_n_{R}.png")
                        img_array.append(img)
                    except:
                        print(f"Image not found for a = {a}, T = {T}, R = {R}, s = {s}, Skipping....")
                        continue

        height, width, layers = img_array[0].shape; size = (width, height)
        video_name = out_dir + f'Videos/IndRho{s}_a_{a}.mp4'
        out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 1, size)

        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()
        '''
        # Combining images to create video for each a value.
        # Images are arranged to create video in increasing order of T, followed by increasing order of R.
        img_array = []; img_array2 = []
        for R in range(R_max):
            for T in T_vals:
                T = int(T) if T.is_integer() else T
                try:
                    img = cv2.imread(out_dir + f"Combined/Concentration_a_{a}_T_{T}_n_{R}.png")
                    img_array.append(img)

                except:
                    print(f"Image not found for a = {a}, T = {T}, R = {R}, Skipping....")
                    continue
                try:
                    img2 = cv2.imread(out_dir + f"Combined/Gamma_a_{a}_T_{T}_n_{R}.png")
                    img_array2.append(img2)
                except:
                    print(f"Image not found for a = {a}, T = {T}, R = {R}, Skipping....")
                    continue
                
        height, width, layers = img_array[0].shape; size = (width, height)
        video_name = out_dir + f'Videos/CombinedRho_a_{a}.mp4'
        out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 1, size)

        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()

        if(len(img_array2) > 0):
            height, width, layers = img_array2[0].shape; size = (width, height)
            video_name = out_dir + f'Videos/CombinedGamma_a_{a}.mp4'
            out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 1, size)

            for i in range(len(img_array2)):
                out.write(img_array2[i])
            out.release()
        # Deleting the images
        if delpng:
            for R in range(R_max):
                for T in T_vals:
                    for s in range(3):
                        try:
                            os.remove(out_dir + f"Singular/Rho_s_{s}_a_{a}_T_{T}_n_{R}.png")
                            os.remove(out_dir + f"Combined/Concentration_a_{a}_T_{T}_n_{R}.png")
                        except:
                            print(f"Image not found for a = {a}, T = {T}, R = {R}, s = {s}, Skipping....")
                            continue



frame_visualiser()
            

        

