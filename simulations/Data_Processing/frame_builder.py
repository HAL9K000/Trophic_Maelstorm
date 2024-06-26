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
import fnmatch
from pathlib import Path

import matplotlib.ticker as mtick
import matplotlib.colors as mcolours
from matplotlib.collections import LineCollection


# User inputs.

in_dir = "../Data/Remote/Rietkerk/Reorg_Frames/3Sp/StdParam_20_100_MFTEQ/"
out_dir = "../../Images/3Sp/StdParam_20_100_MFTEQ/"
Path(out_dir).mkdir(parents=True, exist_ok=True)
prefixes = ["DiC-REF-1.1HI", "DiC-REF-0.5LI", "DiC-NREF-1.1HI", "DiC-NREF-0.5LI"]

g = 128;  dP = 10000; Geq = 4.802; R_max= -1;
T_vals =[]
a_vals = []    

#Making the output directories if they don't exist
#Path(out_dir + 'Videos/').mkdir(parents=True, exist_ok=True)
#Path(out_dir + 'Combined/').mkdir(parents=True, exist_ok=True)
#Path(out_dir + 'Singular/').mkdir(parents=True, exist_ok=True)

hex_list= [['#D8F3DC', '#B7E4C7', '#95D5B2', '#74C69D', '#57CC99', '#38A3A5', '#006466', '#0B525B', '#1B4332'], 
           ['#B7094C', '#A01A58', '#892B64', '#723C70', '#5C4D7D', '#455E89', '#2E6F95', '#1780A1', '#0091AD'], 
           ['#B80068', "#BB0588", "#CC00AA", "#CC00CC", "#BC00DD", "#B100E8", "#A100F2", "#8900F2", '#6A00F4'],
           ["#c4fff9", "#9ceaef", "#68d8d6", "#3dccc7", "#07beb8"], 
           ['#CAF0F8', '#0096C7']]

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


''' # Summary of auto_guess_avals_tvals(...)
    NOTE: This only works when the indir structure is of the form: 
    indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/
    This function attempts to recreate a_vals (and DEPENDENT t_vals nested within) from crawling relevant txt files in the above directory structure.
    Default behaviour is to guess a common set of  a_vals and t_vals for all sub-directories in indir + PREFIX.
    When a_vals is provided, it returns a_vals (unmodified) and a guessed set of t_vals 
    that exists in indir + PREFIX + f"/L_{g}_a_{a}/ for all such a vals
    When t_vals is passed as None, it only returns a_vals.
'''
def auto_guess_avals_tvals(indir, PREFIX, a_vals = [], t_vals= []):

    if len(a_vals) == 0:
        try:
            with open(indir + PREFIX + "/a_vals.txt", 'r') as f:
                a_vals = [float(line.strip()) for line in f]
                a_vals = [int(a) if a.is_integer() else a for a in a_vals]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory and no a_vals specified. Terminating the program.")
            
            #a_vals = None; t_vals = None;
            return (None, None)
    #print(t_vals)        
    if len(t_vals) == 0:
        for a in a_vals:
            a = int(a) if a.is_integer() else a
            try:
                with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt", 'r') as f:
                    t_vals.extend([float(line.strip()) for line in f])
            except FileNotFoundError:
                print(f"T_vals.txt not found in the input directory and no T_vals specified for a = {a}. Skipping this a-value.")
                continue
        t_vals= [int(T) if T.is_integer() else T for T in t_vals]
        # Sort and remove duplicates
        t_vals = sorted(list(set(t_vals)))
        # End of T loop
    return (a_vals, t_vals)
    
''' # Summary of home_video(...)
This function creates a video of the png data in the directory dir which has the structure:
{filedir}/{prefixes}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/{pngformat}
The data is stored in the png files with the format specified in the pngformat parameter.
The video is saved in the out_dir directory with the name specified in the videoname parameter.
The function reads the data from the png files and stores them in a list of images.
Optionally, one can provide the relative path (from parendir) to the text files delineating a_vals, T_vals and maxR (if these values are not specified/ set to 0)
by the relative paths given by avals_txtdir, Tvals_txtdir and maxR_txtdir. 
The same applies for pngdir (relative path to png files for given a, T, dP, Geq etc values)
'''
def home_video(file_parendir, out_dir, prefixes, a_vals, T_vals, maxR, pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir =  "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
               videoname="guess", video_relpath="guess", avals_txtdir = "{Pre}/a_vals.txt", Tvals_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt",
                maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt"):

    min_T = 0; max_T = 0; inputmax_R = maxR;
    if len(T_vals) > 0:
        min_T= min(T_vals); max_T = max(T_vals);# These are used to determine savepath and video name.

    for Pre in prefixes:
        img_array = []

        if len(a_vals) == 0:
            try:
                with open(file_parendir + avals_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq), 'r') as f:
                    a_vals = [float(line.strip()) for line in f]
            except FileNotFoundError:
                print(f"a_vals.txt not found in the input directory and no a_vals specified for {Pre}. skipping this prefix.")
                continue
        for a in a_vals:
            a = int(a) if float(a).is_integer() else a
            if len(T_vals) == 0:
                try:
                    with open(file_parendir + Tvals_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq, a=a), 'r') as f:
                        T_vals = [float(line.strip()) for line in f]
                    if(max(T_vals) > max_T):
                        max_T = max(T_vals)
                    if(min(T_vals) < min_T):
                        min_T = min(T_vals)
                except FileNotFoundError:
                    print(f"T_vals.txt not found in the input directory and no T_vals specified for a = {a}, {Pre}. Skipping this a-value.")
                    continue
            print(f"List of T values: {T_vals}")
            
            for T in T_vals:
                T = int(T) if float(T).is_integer() else T
                print(f"Processing a = {a}, T = {T} at dP = {dP}, Geq = {Geq}")
                if inputmax_R <= 0:
                    # This indicates a case where the number of frames is not known apriori and must be read from the data.
                    # This data is stored in the first line of the following text file:
                    # in_dir/PREFIX/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt
                    try:
                        with open(file_parendir + maxR_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq, a=a, T=T), 'r') as f:
                            maxR = int(f.readline().strip())

                    except FileNotFoundError:
                        print(f"maxR.txt not found in the input directory for a = {a}, T = {T}. Skipping this a-value.")
                        continue
                
                for R in range(0, maxR):
                    try:
                        img = cv2.imread(file_parendir +pngdir.format(s=0, Pre=Pre, g=g, dP=dP, Geq=Geq, a=a, T=T, R=R) 
                                         + pngformat.format(s=0, L=g, dP=dP, Geq=Geq,  a=a, T=T, R=R))
                        img_array.append(img)
                    except:
                        print(f"Image not found for a = {a}, T = {T}, R = {R}, Pref = {Pre}, in Path: \n"  + pngdir.format(s=0, Pre=Pre, L=g, dP=dP, Geq=Geq, a=a, T=T, R=R) 
                                         + pngformat.format(s=0, L=g, dP=dP, Geq=Geq,  a=a, T=T, R=R) + " \n Skipping....")
                        continue
                # End of R loop
            # End of T loop
        # End of a loop
        print(f"Number of images: {len(img_array)}")
        if(len(img_array) > 0):
            height, width, layers = img_array[0].shape; size = (width, height)
            if videoname == "guess":
                #Extract string preceding the first underscore in the pngformat string.
                Prefix = re.sub("/", "", Pre)
                videoname = f"{Prefix}_Tmax_{max_T}.mp4"
            
            if video_relpath == "guess":
                # Modify video save path depending on number of a_vals and T_vals encountered.
                if(len(a_vals) > 1):
                    video_relpath = out_dir + f'{Pre}/Videos/'+ pngformat.split("_")[0] + f'/{a_vals[0]}_{a_vals[-1]}/' 
                else:
                    video_relpath = out_dir + f'{Pre}/Videos/'+ pngformat.split("_")[0] + f'/{a_vals[0]}/'
                if(len(T_vals) > 1):
                    video_relpath += f'Tmin_{min_T}_Tmax_{max_T}/'
                else:
                    video_relpath += f'Tmax_{max_T}/'
            else:
                video_relpath = out_dir + video_relpath.format(Pre=Pre, g=g, dP=dP, Geq=Geq, a=a_vals[-1], T= max_T, R= maxR,  
                                                               amin= a_vals[0], amax = a_vals[-1], Tmin = min_T, Tmax = max_T)
            #Path(video_path).mkdir(parents=True, exist_ok=True) # Create the directory if it doesn't exist.
            video_name = video_relpath + videoname
            os.makedirs(video_relpath, exist_ok=True)
            out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 1, size)

            print(f"Saving {videoname} to Relative Path :\n {video_relpath}")

            for i in range(len(img_array)):
                out.write(img_array[i])
            out.release()
        else:
            print(f"No images found for {Pre}. Skipping....")
        
    # End of Prefix loop


def frame_visualiser(in_dir, out_dir, PREFIX, a_vals, T_vals, maxR, plt_gamma= False, delpng = False):
    # If a_vals and T_vals are not provided, the script reads the relevant text files and extracts the a_vals and T_vals from them.
    # List of a_vals provided {in_dir}/{PREFIX}/a_vals.txt
    # List of T_vals provided {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/T_vals.txt

    input_avals = a_vals; input_tvals = T_vals; input_Rmax = maxR
    print(f"Processing data for {PREFIX}...")
    #global a_vals, T_vals
    
    if a_vals == []:
        try:
            with open(in_dir + PREFIX + "/a_vals.txt", 'r') as f:
                a_vals = [float(line.strip()) for line in f]
                a_vals = [int(a) if a.is_integer() else a for a in a_vals]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory and no a_vals specified. Terminating the program.")
            return
    print(f"List of a values: {a_vals}")
    for a in a_vals:
        # Convert a and T to int format if they are integers.
        a = int(a) if float(a).is_integer() else a
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
            T = int(T) if float(T).is_integer() else T
            print(f"Processing a = {a}, T = {T} at dP = {dP}, Geq = {Geq}")
            if input_Rmax <= 0:
                # This indicates a case where the number of frames is not known apriori and must be read from the data.
                # This data is stored in the first line of the following text file:
                # in_dir/PREFIX/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt
                try:
                    #print(f"Reading maxR.txt for {PREFIX}, L= {g}, dP= {dP}, a = {a}, T = {T}...")
                    with open(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", 'r') as f:
                        maxR = int(f.readline().strip())
                except FileNotFoundError:
                    print(f"maxR.txt not found in the input directory for a = {a}, T = {T}. Skipping this a-value.")
                    continue
            
            savepath = out_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/"
            Path(savepath + "CombinedImg/").mkdir(parents=True, exist_ok=True)

            print(f"Creating images for a = {a}, T = {T} at dP = {dP}, Geq = {Geq}, maxR = {maxR}...")

            for R in range(0, maxR):
            # Reading the data from the csv file
                # Convert a and T to int format if they are integers.
                try:
                    data = pan.read_csv(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a}_R_{R}.csv", header= 0)
                except FileNotFoundError:
                    print(f"Data not found for a = {a}, T = {T}, R = {R}. Skipping....")
                    continue

                # Extracting the data columns
                data = data.iloc[:, 2:]
                fig, axs = plt.subplots(1, 3, figsize = (20, 6))
                ''' # Individual plots as opposed to combined ones.
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
                    ax = sea.heatmap(data_Sp , ax=ax, cmap=get_continuous_cmap(hex_list[s], float_list),  cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                
                plt.savefig(savepath + f"CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png")
                plt.close()

                # Creating subplots
                fig, axs = plt.subplots(1, 2, figsize = (12, 6))
                for s in range(2):
                    ax = axs[s]
                    data_Sp = np.array(data.iloc[:, s+3]).reshape(g, g)
                    ax = sea.heatmap(data_Sp , ax=ax, cmap=get_continuous_cmap(hex_list[s+3]),  cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s+3, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                
                plt.savefig(savepath + f"CombinedImg/WatConc_a_{a}_T_{T}_n_{R}.png")
                plt.close()


                if plt_gamma:
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
                    plt.savefig(savepath + f"CombinedImg/Gamma_a_{a}_T_{T}_n_{R}.png")
                    plt.close()

            # End of R loop
            plt.close()
        # End of T loop

        # Combining images to create video for each a value.
        # Images are arranged to create video in increasing order of T, followed by increasing order of R.
        img_array = []
        
        ''' PREVIOUS VERSION OF HOME_VIDEO FUNCTION
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
        '''
        """
        home_video(out_dir, out_dir, [PREFIX], input_avals, input_tvals, maxR=  input_Rmax, 
               pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
               videoname = "guess")
        """
    

''' # Summary Description of get_timeseriesData(...) (called by analyse_timeseriesData(...))
This function creates a multi-Indexed DataFrame with the timeseries data for each species.
The indices are as follows:
Prefix, a, T
The columns are determined by reading the header of tab-delineated txt file MEAN_INDV_REPLICATES.txt
which is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/
The columns of each such txt file follow the following structure:
R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
where {var} is one of the species names.
The second line of the file contains the values for the columns.

Each such file is read and the data is stored in the DataFrame and returned. 
Note that since R_max can vary for each a and T value, the DataFrame will have NaN values for the missing columns.

def get_timeseriesData(indir, PREFIX, a_vals, T_vals, filename = "MEAN_INDV_REPLICATES.txt", construct_means = False):

    auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    # Guess all a_vals and T_vals from the directory structure if not provided.
    print(f"List of a values: {a_vals}")
    print(f"List of T values: {T_vals}")
    # Create a multi-indexed DataFrame
    data = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            try:
                # Read tab-delineated txt file
                df = pan.read_table(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                # Use the header row as the column names and the second row as the data
                df.columns = df.columns.astype(str)
                #print(f"Data for a = {a}, T = {T}:")
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])
                #Assign the data df to the appropriate index in the data DataFrame
                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
            except FileNotFoundError:
                print(f"Data not found for a = {a}, T = {T}. Skipping....")
                continue
    return data
'''


'''# Summary Description of get_EQdata(...) (called by analyse_EQdata(...))
This function creates a multi-Indexed DataFrame with the equilibrium data for each species.
The indices are as follows:
Prefix, a, T
The columns are determined by reading the header of tab-delineated txt file MEAN_STD_Surviving_runs.txt
which is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/
These column names should be the same for all such files in all such directories.
The second line of the file contains the values for the columns.
All the data is stored in the DataFrame and returned.
# If a_vals and T_vals are not provided, the script reads the relevant text files and extracts the a_vals and T_vals from them.
# List of a_vals provided {in_dir}/{PREFIX}/a_vals.txt
# List of T_vals provided {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/T_vals.txt
'''
def get_EQdata(indir, PREFIX, a_vals, T_vals, filename = "MEAN_STD_Surviving_Runs.txt"):

    a_vals, T_vals = auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    if a_vals == None:
        print("Exiting..."); return;
    print(f"List of a values: {a_vals}")
    print(f"List of T values: {T_vals}")
    
    # Create a multi-indexed DataFrame
    #data = pan.MultiIndex.from_product([[PREFIX], a_vals, T_vals], names = ['Prefix', 'a', 'T'])
    data = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            try:
                # Read tab-delineated txt file
                df = pan.read_table(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                #data = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                # Use the header row as the column names and the second row as the data

                df.columns = df.columns.astype(str)
                #print(f"Data for a = {a}, T = {T}:")

                df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])

                '''#print(df.info())
                #print(df.head())
                #print(df.columns)'''

                #Assign the data df to the appropriate index in the data DataFrame
                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
                #print(data.columns)
                '''
                df = df.iloc[0, :]
                print(df.info()); print(data.shape);
                print(df.columns)
                df = df.to_frame().transpose()
                # Check and set column names
                if not hasattr(df, 'columns'):
                    df.columns = df.columns.astype(str)
                # Set column names to the first row of the data
                df.columns = df.iloc[0, :]
                
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])
                '''
            except FileNotFoundError:
                print(f"Data not found for a = {a}, T = {T}. Skipping....")
                continue
    return data

'''# Summary Description of analyse_timeseriesData(...)
This function analyses the timeseries data for each species and creates plots for each species for each Prefix and a value, with T as the x-axis.
NOTE: This function assumes that the data is stored in the FILENAME the following format:
R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
where {var} is one of the species names.
'''
def analyse_timeseriesData(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], filename = "MEAN_REPLICATES.txt"):

    savefilename = re.sub(r'\.txt$', '', filename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    init_a_vals = a_vals; init_T_vals = T_vals
    

    for Pre in prefixes:
        
        data = get_EQdata(indir, Pre, init_a_vals, [], filename)
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())

        # Find Tmax as the maximum value of T in the data.
        Tmax = data.index.get_level_values('T').max()
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax
        savecsvdir = out_dir + f"{Pre}/TimeSeries/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)

        # Save the data to a csv file.
        data.to_csv(savecsvdir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        # Use the data to create plots.
        # The plots are created for each species for each Prefix and a value, with T as the x-axis.
        # For each species, create a plot of T vs all the data for that species as lines, with the mean across surviving and all replicates as darker lines.
        # The means for individual replicates are lighter lines.
        
        # Try determining number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in data.columns if "_SURV" in col])
        except IndexError:
            print(f"No columns with _SURV found in the data. Not plotting data for Prefix {Pre}. Skipping....")
            num_species = None
            continue

        # Find column index of first column with _R_0 in the name. If not found, print a warning and set to 0.
        try:
            R0_index = data.columns.get_loc([col for col in data.columns if "_R_0" in col][0])
        except IndexError:
            print(f"No column with _R_0 found in the data. Not plotting individual replicate data for Prefix {Pre}.")
            R0_index = None
        print(f"R0_index = {R0_index}, num_species = {num_species}")
        # Get the list of a values and T values from the data.
        a_vals = data.index.get_level_values('a').unique().to_list()
        print(a_vals, T_vals)
        for a in a_vals:
            a = int(a) if float(a).is_integer() else a
            savepngdir = out_dir + f"{Pre}/TimeSeries/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{Tmax}/"
            Path(savepngdir).mkdir(parents=True, exist_ok=True)
            fig, axs = plt.subplots(1, 3, figsize = (20, 6))
            for s in range(3):
                data_A = data.loc[(slice(None), a, slice(None)), :]
                # Plotting for a fixed a value, species and Prefix for all T values.
                ax = axs[s]
                # RECALL, First column is R_max, then the AVG columns for surviving replicates, then AVG columns for all replicates, then AVG columns for each individual replicate.
                # Plotting mean of surviving replicates
                ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, 1+ s], label = r"$\mu_{surv}(\rho_{%g})$" %(s),
                        color = colours[s], linestyle = 'solid', linewidth = 2, alpha = 0.9)
                # Plotting mean of all replicates
                ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, 1 + num_species +  s], label = r"$\mu_{all}(\rho_{%g})$" %(s),
                        color = colours[s], linestyle = 'dashed', linewidth = 1.35, alpha = 0.75)
                # Plotting mean of individual replicates (suppress labels that represent columns with solely NaN values).
                if R0_index is not None and num_species is not None:
                    for r in range(R0_index, len(data_A.columns), num_species):
                        if not data_A.iloc[:, r].isnull().all():
                            ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, r + s],
                                    color ="grey", linestyle = 'solid', alpha = 0.45)
                ax.set_title(r"$\langle \rho_{%g}(t) \rangle_{x}$ vs $t$" %(s)  + f" at a = {a}")
                ax.set_xlabel(r'Time $(d)$')
                ax.set_ylabel(r'$\langle \rho_{%g}(t) \rangle_{x}$' % (s))
                ax.legend()
            #Set ymax for axs[3] to be 500.
            axs[2].set_ylim(-5, 500)
            fig.suptitle(r"$\mu(\rho(x, t))$" + f" For {Pre} at a = {a}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savepngdir + f"TimeSeries_a_{a}_dP_{dP}_Geq_{Geq}.png")
            #plt.show()
            plt.close()
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/TimeSeries"], a_vals= a_vals, T_vals=[Tmax], maxR= 1, 
                   pngformat= "TimeSeries_a_{a}_dP_{dP}_Geq_{Geq}.png", pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/", 
                   videoname = f"TimeSeries_{Pre}_Tmax_{Tmax}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")
        # Note that there is only one png per a value, so the video will be a simple slideshow (and R_max = 1).
        input("Press F to Continue...")
    # End of Prefix loop
    '''
    # Finally if prefixes > 1, use the combined data to create plots of T vs. the mean of the data for each species, with Prefix as hue.
    if len(prefixes) > 1:
        # Estimate number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in combined_data.columns if "_SURV" in col])
        except IndexError:
            print(f"No columns with _SURV found in the data. Skipping....")
            num_species = None
            return
        # Estimate R0_index by finding the first column with _R_0 in the name. If not found, print a warning and set to 0.  
        try:
            R0_index = combined_data.columns.get_loc([col for col in combined_data.columns if "_R_0" in col][0])
        except IndexError:
            print(f"No column with _R_0 found in the data. Not plotting individual replicate data for multiple prefixes.")
            R0_index = None
        for a in a_vals:
            fig, axs = plt.subplots(1, 3, figsize = (20, 6))
            for s in range(3):
                data_A = combined_data.loc[(slice(None), a, slice(None)), :]
                ax = axs[s]
                sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:, 1+ s], hue = 'Prefix', ax = ax, alpha = 1)
                sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:,1 + num_species + s], hue = 'Prefix', ax = ax, alpha = 0.75)
                if R0_index is not None and num_species is not None:
                    for r in range(R0_index, len(data_A.columns), num_species):
                        if not data_A.iloc[:, r].isnull().all():
                            sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:, r + s], hue = 'Prefix', ax = ax, alpha = 0.35)
                ax.set_title(r"$\langle \rho_{%g}(t) \rangle_{x}$ vs $t$" %(s) + f" at a = {a}")
                ax.set_xlabel(r'Time $(d)$')
                ax.set_ylabel(r'$\langle \rho_{%g}(t) \rangle_{x}$' % (s))
                ax.legend()
            fig.suptitle(r"$\mu$ of" + f" Species For All Prefixes at a = {a}, dP = {dP}, Geq = {Geq}")
            Path(out_dir + "TimeSeries/Combined/").mkdir(parents=True, exist_ok=True)
            plt.savefig(out_dir + f"TimeSeries/Combined/TimeSeries_AllPrefixes_a_{a}_dP_{dP}_Geq_{Geq}.png")
            plt.show()
            plt.close()
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, ["TimeSeries/Combined"], a_vals, [Tmax], 1, pngformat= "TimeSeries_AllPrefixes_a_{a}_dP_{dP}_Geq_{Geq}.png", videoname = f"TimeSeries_All_T_{Tmax}.mp4")
    '''

def analyse_EQdata(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], Tavg_window_index = [-30000, 0], filename = "MEAN_STD_Surviving_Runs.txt"):

    savefilename = re.sub(r'\.txt$', '', filename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    for Pre in prefixes:

        savedir = out_dir + f"{Pre}/PhaseDiagrams/"
        Path(savedir).mkdir(parents=True, exist_ok=True)
        
        data = get_EQdata(indir, Pre, a_vals, [], filename)
        # Data has columns for each species, of the form:
        # AVG[P(x; t)]	VAR[P(x; t)]	COUNT[P(x; t)]	AVG[G(x; t)]	VAR[G(x; t)]	COUNT[G(x; t)] ...
        # We are interested in the AVG and VAR columns for each species.
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print(f"Data for {Pre}:")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())
        # Save the data to a csv file
        data.to_csv(savedir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        Avals_list = sorted(data.index.get_level_values('a').unique().to_list())
        Tavg = pan.DataFrame()
        # Iterating over the Pre, a indices, average over the last Tavg_window_index values of T and assign to new DataFrame.
        for a in Avals_list:
            # Get Tmax for this a value
            Tmax = data.loc[(slice(None), a, slice(None)), :].index.get_level_values('T').max()
            Twin_min = Tmax + Tavg_window_index[0]; Twin_max = Tmax + Tavg_window_index[1]
            print(f"Averaging data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}:")
            df = data.loc[(slice(None), a, slice(Twin_min, Twin_max)), :].groupby(['Prefix', 'a']).mean()
            df.index = pan.MultiIndex.from_tuples([(Pre, a, Twin_max)], names = ['Prefix', 'a', 'Twindow'])
            Tavg = pan.concat([Tavg, df], axis = 0).sort_index(axis = 0)
        # End of a loop

        print(f"Data for {Pre} after averaging over last {Tavg_window_index[1]} values of T:")
        print(Tavg.head())
        print(Tavg.tail())
        print(Tavg.columns)
        print(Tavg.index)

        #Save the data to a csv file
        Tavg.to_csv(savedir + f"DEBUG_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.csv")

        # Use Tavg to create plots of a vs. the mean of the data for each species, with STD as infill.
        fig, axs = plt.subplots(1, 3, figsize = (20, 6))
        for s in range(3):
            ax = axs[s]
            ax.scatter(Tavg.index.get_level_values('a'), Tavg.iloc[:, 3*s], label = r"$\mu_{%g}$" %(s), color = colours[s], s = 15, alpha = 0.75)
            err = 2 * np.sqrt(Tavg.iloc[:, 3*s + 1]) # Get 95% confidence interval
            ax.fill_between(Tavg.index.get_level_values('a'), Tavg.iloc[:, 3*s] - err, Tavg.iloc[:, 3*s] + err, alpha = 0.35, color = colours[s])
            # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
            for i, txt in enumerate(Tavg.iloc[:, 3*s + 2]):
                ax.annotate(int(txt), (Tavg.index.get_level_values('a')[i], Tavg.iloc[i, 3*s]))
            ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r} $" % (s, Twin_max))
            ax.set_xlabel(r'R $(mm/hr)$')
            ax.set_ylabel(str(data.columns[3*s]))
            ax.legend()
        fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For {Pre}, dP = {dP}, Geq = {Geq} B/W {Twin_min} -- {Twin_max}")
        plt.savefig(savedir + f"MeanSTD_{Pre}_amin_{Avals_list[0]}_amax_{Avals_list[-1]}_dP_{dP}_Geq_{Geq}.png")
        plt.show()
        plt.close()
        

        ''' # Previous version of the code, where we iterate over T values and create plots for each T value.
        # Use the data to create plots.
        # For each T-value, create a plot of a vs. the mean of the data for each species, with STD as infill.
        for T in T_vals:
            T = int(T) if float(T).is_integer() else T
            data_T = data.loc[(slice(None), slice(None), T), :]
            fig, axs = plt.subplots(1, 3, figsize = (20, 6))
            for s in range(3):
                ax = axs[s]
                ax.scatter(data_T.index.get_level_values('a'), data_T.iloc[:, 3*s], label = r"$\rho_{%g}$" %(s), color = colours[s], s = 15, alpha = 0.75)
                err = 2 * np.sqrt(data_T.iloc[:, 3*s + 1]) # Get 95% confidence interval
                ax.fill_between(data_T.index.get_level_values('a'), data_T.iloc[:, 3*s] - err, data_T.iloc[:, 3*s] + err, alpha = 0.35, color = colours[s])
                # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
                for i, txt in enumerate(data_T.iloc[:, 3*s + 2]):
                    ax.annotate(txt, (data_T.index.get_level_values('a')[i], data_T.iloc[i, 3*s]))
                ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r} $" % (s, T))
                ax.set_xlabel(r'R $(mm/hr)$')
                ax.set_ylabel(str(data.columns[3*s]))
                ax.legend()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For {Pre} at T = {T}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savedir + f"MeanSTD_{Pre}_T_{T}_dP_{dP}_Geq_{Geq}.png")
            plt.show()
            plt.close()
        '''

    savedir = out_dir + "PhaseDiagrams/"
    Path(savedir).mkdir(parents=True, exist_ok=True)
    # Finally if prefixes > 1, use the combined data to create plots of a vs. the mean of the data for each species, with STD as infill, and Prefix as hue.
    if len(prefixes) > 1:
        for T in T_vals:
            T = int(T) if float(T).is_integer() else T
            combined_data_T = combined_data.loc[(slice(None), slice(None), T), :]
            fig, axs = plt.subplots(1, 3, figsize = (20, 6))
            for s in range(3):
                ax = axs[s]
                sea.lineplot(data = combined_data_T, x = 'a', y = combined_data_T.iloc[:, 3*s], hue = 'Prefix', ax = ax)
                err = 2 * np.sqrt(combined_data_T.iloc[:, 3*s + 1])
                ax.fill_between(combined_data_T.index.get_level_values('a'), combined_data_T.iloc[:, 3*s] - err, combined_data_T.iloc[:, 3*s] + err, alpha = 0.35)
                # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
                #for i, txt in enumerate(combined_data.iloc[:, 3*s + 2]):
                #    ax.annotate(txt, (combined_data.index.get_level_values('a')[i], combined_data.iloc[i, 3*s]))
                ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r}$" % (s, T))
                ax.set_xlabel('R $(mm/hr)$')
                ax.set_ylabel(str(data.columns[3*s]))
                ax.legend()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" For All Prefixes at T = {T}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savedir + f"MeanSTD_T_{T}_AllPrefixes.png")
            plt.show()
            plt.close()


''' # Summary Description of recursive_copydir(...)
This function recursively copies files from the source directory to the destination directory.
src: source directory
dst: destination directory
include_filetypes: list of filetypes to be copied
symlinks: If True, symlinks are copied as symlinks. If False, the files are copied as hardlinks.
'''    
def recursive_copydir(src, dst, include_filetypes = ["*.txt"], 
                      exclude_filetypes =["*.csv", "*.png", "*.jpg", "*.jpeg", "*.mp4"], symlinks=False):
    for root, dirs, files in os.walk(src):
        for file in files:
            if any(fnmatch.fnmatch(file, pattern) for pattern in include_filetypes):
                if not any(fnmatch.fnmatch(file, pattern) for pattern in exclude_filetypes):
                    src_file = os.path.join(root, file)
                    dst_file = os.path.join(dst, os.path.relpath(src_file, src))
                    os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                    if symlinks:
                        if os.path.exists(dst_file):
                            os.remove(dst_file)
                        os.symlink(src_file, dst_file)
                    else:
                        shutil.copy2(src_file, dst_file)

    print("\n\n__________________________________________________________________________________________\n\n")
    print(f"Done copying files from {src} to {dst}.")
    print("\n\n__________________________________________________________________________________________\n\n")


#recursive_copydir(in_dir, out_dir, include_filetypes = ["*.txt"], exclude_filetypes =["*.png", "*.jpg", "*.jpeg", "*.mp4"], symlinks=False)
#a_vals = [0.041, 0.0415, 0.042, 0.0425, 0.043, 0.044, 0.045, 0.047, 0.05, 0.053]; 
#T_vals= [229087]
'''
for Pre in prefixes:
    frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= 0, plt_gamma= False, delpng = False)
    print(f"Done with making frames for {Pre}")
    for a in a_vals:
        print(f"Making video for {Pre} at a = {a} \n\n")
        home_video(out_dir, out_dir, [Pre], [a], T_vals=[], maxR= -1, 
                   pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
                     maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", videoname = "guess",
                     video_relpath= "{Pre}/Videos/Conc/{a}/{Tmin}-{Tmax}/")
'''
#frame_visualiser(in_dir, out_dir, prefixes[0], a_vals, T_vals, maxR= -1, plt_gamma= False, delpng = False)
#home_video(out_dir, out_dir, prefixes, a_vals, T_vals, R_max, pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", videoname = "guess")

prefixes = ["DIC-REF-1.1HI", "DIC-REF-0.5LI"]
analyse_EQdata(in_dir, out_dir, prefixes, a_vals, T_vals, filename = "MEAN_STD_Surviving_Runs.txt")
#analyse_timeseriesData(in_dir, out_dir, prefixes, a_vals, T_vals, filename = "MEAN_REPLICATES.txt")