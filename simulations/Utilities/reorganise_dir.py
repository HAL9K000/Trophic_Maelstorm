import numpy as np
import os
import pandas as pan
import regex as re
from pathlib import Path
import glob
import shutil
import sys
import argparse
import copy
import warnings

import scipy.stats as stats
from scipy.interpolate import CubicSpline

from glow_up import *

# Show the first FutureWarning that occurs in the script, then ignore all future FutureWarnings.
warnings.simplefilter(action='once', category=FutureWarning)
warnings.simplefilter(action='once', category=pan.errors.PerformanceWarning)


'''
This script reorganises the directory structure of the data files in the Rietkerk model.
The original directory structure (filepaths) are as follows if jID (unique ID assigned to each replicate) is present in the filename:
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_jID_{jID}_*_R_{R}.csv OR
root_dir/{PREFIX}*_dP_{dP}_Veq_{Veq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_jID_{jID}_*_R_{R}.csv

The original directory structure (filepaths) are as follows if jID is not present in the filename:
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv OR
root_dir/{PREFIX}*_dP_{dP}_Veq_{Veq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv

The new directory structure (filepaths) will be as follows:
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv OR
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Veq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv

The script will:
1. Find all immediate subdirectories in root_dir.
2. For each subdirectory, find all unique values of a_val in the files in the subdirectory.
3. For each value of a_val, find all unique values of T and jID (if present) in the files in the subdirectory.
4. Looping over an ascending order of jID, for each value of T, copy all files in the subdirectory to the new directory structure in out_dir.
5. If a file with the same name already exists in the new directory, rename the file to avoid conflicts.

This script ASSUMES that the files in the original directory structure have the following format,
and WILL NOT WORK if the files are named differently:
I. Contains "_G_{L}_T_{Time}_" in the filename, which is followed by
II. "_a_{a_val}_" in the filename, which is followed by
III. "_R_{R}.csv" in the filename.
IV. "_jID_{jID}_" MAY OR MAY NOT be present in the filename.

The script accepts the following arguments:
1. prefixes: A list of prefixes to be used in the subdirectory names in out_dir.
2. root_dir: The root directory containing the original data files.
3. out_dir_noprefix: The output directory where the reorganised data files will be stored.
4. dP: The value of dP in the subdirectory names in root_dir.
5. Geq: The value of Geq in the subdirectory names in root_dir (if "NA", Geq is not used in the subdirectory name).
6. L: A list of values of L to be used in the subdirectory names in out_dir.
7. indx_vals_t: Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.

'''

#prefixes =["DiC-NREF-1.1HI", "DiC-NREF-0.5LI", "DiC-NREF-0.1LI"]#, "DiC-NEW"]
prefixes =["HX2001-UA125A125-5E2UNI", "HX1005-UA125A0-5E2UNI"]
#prefixes =["", "DiC", "BURNIN", "DiC-BURNIN", "DDM-DiC", "DDM-DiC-BURNIN"]
#root_dir = "../Data/Amarel/Rietkerk/Prelims/Stochastic/3Sp/"
root_dir = "../Data/Remote/Rietkerk/Frames/Stochastic/3Sp/"
#out_dir_noprefix = "../Data/Remote/Rietkerk/Reorg_Frames/1Sp/StdParam_MFT/"
out_dir_noprefix = "../Data/Remote/Rietkerk/Reorg_Frames/3Sp/ASCALE_20_100_BRNIN/"

#out_dir_noprefix = "../Data/Remote/Rietkerk/Reorg_Frames/3Sp/StdParam_20_100_MFTNu/"

dP = 10000
Geq = 0.19208 # Optional. If Geq is not used in the subdirectory name, set Geq = "NA".
Veq = "NA" # Optional. If Veq is not used in the subdirectory name, set Veq = "NA".
L= [128]
indx_vals_t = -25
#Extract n largest values of T if indx_vals_t = -n, 
# n smallest values of T if indx_vals_t = n.
tmin = None; tmax = None; 
dynamic_inspect = False;    # Set to True to stop script periodically and inspect values of parameters and outputs.

def set_frames_input():
    parser = argparse.ArgumentParser(description='Reorganise directory structure of data files in Rietkerk model.')
    # Boolean flag for dynamic input values.
    parser.add_argument("--dynamic", action="store_true", help="Use dynamic input values for the script.")
    parser.add_argument("--prefixes", nargs="+", help="List of prefixes to be used in the subdirectory names in out_dir.")
    parser.add_argument("--indir", help="Root directory containing the original data files.")
    parser.add_argument("--outdir", help="Output directory where the reorganised data files will be stored.")
    parser.add_argument("--dP", type=int, help="Value of dP in the subdirectory names in root_dir.")
    parser.add_argument("--Geq", help="Value of Geq in the subdirectory names in root_dir (if 'NA', Geq is not used in the subdirectory name).")
    parser.add_argument("--Veq", help="Value of Veq in the subdirectory names in root_dir if Geq is unknown (if 'NA', Veq is not used in the subdirectory name).")
    parser.add_argument("--L", nargs="+", type=int, help="List of values of L to be used in the subdirectory names in out_dir.")
    parser.add_argument("--indx_vals_t", type=int, help="Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.")
    parser.add_argument("--tmin", help="Minimum value of T to be extracted, must be integer or None.")
    parser.add_argument("--tmax", help="Maximum value of T to be extracted, must be integer or None.")
    args = parser.parse_args()

    global prefixes, root_dir, out_dir_noprefix, dP, Geq, Veq, L, indx_vals_t, tmin, tmax, dynamic_inspect

    if args.dynamic:
        # First list all immediate subdirectories in root_dir as information to the user.
        subdir = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(root_dir, '*'))]
        print("For given root_dir: " + root_dir)
        print("Found subdirectories: \n" )
        for sub in subdir:
            print(sub)

        print("Do you wish to provide new values for the following parameters? (y/n)")
        print("out_dir_noprefix: " + out_dir_noprefix)
        print("dP: " + str(dP))
        print("Geq: " + str(Geq))
        print("Veq: " + str(Veq))
        print("L: " + str(L))
        print("indx_vals_t: " + str(indx_vals_t))
        print(f"Tmin:   {tmin} \t Tmax:   {tmax}")

        dynamic_inspect = True

        choice = input("Enter choice: ")
        if choice == "y":
            out_dir_noprefix = input("Enter new out_dir_noprefix: ")
            dP = int(input("Enter new dP: "))
            Geq = input("Enter new Geq ('NA' where it isn't applicable): ")
            Geq = float(Geq) if Geq.isfloat() else Geq
            Veq = input("Enter new Veq ('NA' where it isn't applicable): ")
            Veq = float(Veq) if Veq.isfloat() else Veq
            L = [int(x) for x in input("Enter new L: ").split()]
            indx_vals_t = int(input("Enter new indx_vals_t: "))
            tmin = input("Enter new tmin: "); tmin = int(tmin) if tmin.isnumeric() else None
            tmax = input("Enter new tmax: "); tmax = int(tmax) if tmax.isnumeric() else None

        return

    else:
        if args.prefixes:
            prefixes = args.prefixes
            #Remove newline, tab and trailing or leading whitespaces from prefixes.
            prefixes = [pre.strip("\n\t ").lstrip().rstrip() for pre in prefixes]
        if args.indir:
            root_dir = args.indir
        if args.outdir:
            out_dir_noprefix = args.outdir
        if args.dP:
            try:
                dP = int(args.dP)
            except ValueError:
                print("Need to provide an integer value for dP. Terminating script..."); sys.exit(1)
        if args.Geq:
            try:
                Geq = float(args.Geq)
                # Next, check if Geq is an integer or a float. If it is an integer, convert it to an integer.
                if Geq.is_integer():
                    Geq = int(Geq)
            except ValueError:
                print("Setting Geq to NA..."); Geq = "NA"
        if args.Veq:
            try:
                Veq = float(args.Veq)
                # Next, check if Veq is an integer or a float. If it is an integer, convert it to an integer.
                if Veq.is_integer():
                    Veq = int(Veq)
            except ValueError:
                print("Setting Veq to NA..."); Veq = "NA"
        if args.L:
            L = args.L; L = [int(x) for x in L]
        if args.indx_vals_t:
            try:
                indx_vals_t = int(args.indx_vals_t)
            except ValueError:
                print("Setting indx vals to -10..."); indx_vals_t = -10;
        if args.tmin:
            tmin = int(args.tmin) if args.tmin.isnumeric() else None
        if args.tmax:
            tmax = int(args.tmax) if args.tmax.isnumeric() else None

        print("prefixes: " + str(prefixes))
        print("root_dir: " + root_dir)
        print("out_dir_noprefix: " + out_dir_noprefix)
        print("dP: " + str(dP))
        print("Geq: " + str(Geq))
        print("Veq: " + str(Veq))
        print("L: " + str(L))
        print("indx_vals_t: " + str(indx_vals_t))
        
        return



""" # Summary of the function post_process(...)
# A wrapper function that post-processes the data files in out_dir_noprefix, after the main() function has been executed/
# Out_dir_noprefix is the output directory where the reorganised data files are stored.
# It is structured as follows:
# out_dir_noprefix/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv OR
# out_dir_noprefix/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Veq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
# The function will:
# 1. Find all immediate subdirectories in out_dir_noprefix if prefixes is empty, or in out_dir_noprefix/{PREFIX} for each PREFIX in prefixes.
# 2. Recursively navigate to each sub-directory that contains files of the form "FRAME*.csv".
# 3. Next, for each such non-empty nested sub-directory, use glob to find all files of the form "FRAME*.csv".
# 4. Next, pass this list of files to other functions for further processing.
"""
def post_process(prefixes= []):
    #Assigning prefixes to all subdirectories in out_dir_noprefix if prefixes is empty.
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(out_dir_noprefix, '*/'))]
    
    for pre in prefixes:
        # Recursively navigate to each sub-directory within out_dir_noprefix/{pre} that contains files of the form "FRAME*.csv".
        # These subdirectories may be immediate or nested.
        # For each such non-empty nested sub-directory, use glob to find all files of the form "FRAME*.csv".
        # Pass this list of files to other functions for further processing.
        workdir = os.path.join(out_dir_noprefix, pre)
        all_valid_subdirpaths = [dir for dir, _, files in os.walk(workdir) if len(glob.glob(dir + "/FRAME*.csv")) > 0]
        for subdirpath in all_valid_subdirpaths:
            files = glob.glob(subdirpath + "/FRAME*.csv")
            # Now pass this list of files to other functions for further processing.
            # For example, one can calculate the mean and standard deviation of the files in the subdirectory.
            # Or one can calculate the autocorrelation of the files in the subdirectory.
            # Or one can calculate the power spectrum of the files in the subdirectory.
            # Or one can calculate the mean squared displacement of the files in the subdirectory.
            # Or one can calculate the velocity autocorrelation of the files in the subdirectory.
            # Or one can calculate the spatial correlation of the files in the subdirectory.
            # Or one can calculate the spatial power spectrum of the files in the subdirectory.

            '''# Find mean and standard deviation of each column in each file in files if the column values are non-zero.
            #Mean and Std Rho Density & Max Replicates
            df_surviving = gen_MEAN_SD_COLSfiledata(files, pathtodir=subdirpath, ext="csv", nonzero=True, add_counts= True)
            df_all = gen_MEAN_SD_COLSfiledata(files, pathtodir=subdirpath, ext="csv", nonzero=False, add_counts= True)
            
            # Save this df (excluding first two columns) to a txt file in subdirpath.
            try:
                if df_surviving is not None:
                    df_surviving.to_csv(subdirpath + "/MEAN_STD_Surviving_Runs.txt", sep="\t", index=False, header=True)
                if df_all is not None:
                    df_all.to_csv(subdirpath + "/MEAN_STD_All_Runs.txt", sep="\t", index=False, header=True)
            except Exception as e:
                print("Error: Could not write MEAN_STD to " + subdirpath + "/MEAN_STD.txt with error message: \n" + str(e))
            
            #'''

            '''# Find max R in files.
            Rvals=[]; find_vals(files, ["R_[\d]+"], Rvals) 
            # Assumes R is an integer and files are of the form "FRAME_T_{T}_a_{a_val}_R_{R}.csv"
            maxR = max([int(re.sub("R_", "", r)) for r in Rvals])
            # Save maxR+1 to a text file "maxR.txt" in subdirpath, overwriting the file if it already exists.
            try:
                with open(subdirpath + "/maxR.txt", "w") as f:
                    f.write(str(maxR + 1))
                f.close()
            except Exception as e:
                print("Error: Could not write maxR to */" + os.path.basename(subdirpath) + "/maxR.txt with error message: \n" + str(e))
            #'''

            '''# Find Mean for each column in each file in files and the mean across all files for each column using gen_MEAN_INDVL_Colsfiledata.
            # Save this df to a txt file in subdirpath.
            df_replicates = gen_MEAN_INDVL_Colsfiledata(files, pathtodir=subdirpath, ext="csv")
            try:
                if df_replicates is not None:
                    df_replicates.to_csv(subdirpath + "/MEAN_REPLICATES.txt", sep="\t", index=False, header=True)
            except Exception as e:
                print("Error: Could not write MEAN_REPLICATES to " + subdirpath + "/MEAN_REPLICATES.txt with error message: \n" + str(e))
            #'''

            # Finds the power spectrum of each column in each file using gen_FFT_PowerSpectra.
            # Bin_mask : "auto" (use GMMs to find thresholds), "read" (read thresholds from files in savedir/BIN_MASKS/{files}.txt), 
            # or None (no thresholding).
            df_fft_power = gen_FFT_PowerSpectra(files, pathtodir=subdirpath, ext="csv", Tmin = 150000, Tmax = 240000, bin_mask= "auto", 
                                                exclude_col_labels= ["a_c", "x", "L", "W(x; t)", "O(x; t)"], binwidth=0.5)
            Path(subdirpath + "/FFT_PowerSpect").mkdir(parents=True, exist_ok=True)

            print("Found FFT Power Spectra:")
            #print(df_fft_power)

            try:
                if df_fft_power is not None:
                    df_fft_power.to_csv(subdirpath + "/FFT_PowerSpect/FFT_POWERspectra.csv", sep=",", index=False, header=True)
            except Exception as e:
                print("Error: Could not write FFT_PowerSpectra to " + subdirpath + "/FFT_PowerSpect/FFT_PowerSpectra.csv with error message: \n" + str(e))
            
            '''# Find potential well data in files using gen_potential_well_data.
            # Save the potential well data to a csv file in subdirpath/Pot_Well.
            # If evaluate_local_minima is True, also save the local minima data to a csv file in subdirpath/Pot_Well.
            df_kde, df_local_minima = gen_potential_well_data(files, pathtodir=subdirpath, ext="csv", 
                    exclude_col_labels= ["a_c", "x", "L"], Tmin = 150000, evaluate_local_minima= True, bins =300)
            Path(subdirpath + "/Pot_Well").mkdir(parents=True, exist_ok=True)
            try:
                if df_kde is not None:
                    df_kde.to_csv(subdirpath + "/Pot_Well/Pot_Well.csv", sep=",", index=False, header=True)
                if df_local_minima is not None:
                    df_local_minima.to_csv(subdirpath + "/Pot_Well/LOCAL_MINIMA.csv", sep=",", index=False, header=True)
            except Exception as e:
                print("Error: Could not write potential well data to " + subdirpath + "/Pot_Well/Pot_Well.csv with error message: \n" + str(e))
            #'''  
        # Done with all subdirectories in out_dir_noprefix/{pre} that contain files of the form "FRAME*.csv".
        print("\n=====================================================================================================\n")
        print(f"Done post-processing for prefix: {pre}...")
        print("\n=====================================================================================================\n")
            

def main():
    
    #Find a list of immediate subdirectories in root_dir
    subdir = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(root_dir, '*/'))]
    
    for pre in prefixes:

        out_dir = out_dir_noprefix;
        if pre == "":
            out_dir = out_dir_noprefix + "NoPrefix/"
        else:
            out_dir = out_dir_noprefix + pre + "/"
        
        #Subdirectories are of the form "PREFIX*_dP_{dP}_Geq_{val}" OR "PREFIX*_dP_{dP}_Veq_{val}".
        #Using regex and glob, create a sorted list of subdir (in ascending order of val (can be float or integer)) with prefix pre.
        
        if (Geq == "NA" and Veq == "NA"):
            subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP))), 
                            key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
        else:

            if(Geq != "NA"):
                subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP) + "_Geq_" + str(Geq))), 
                        key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
            elif (Veq != "NA"):
                subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP) + "_Veq_" + str(Veq))), 
                        key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
            
        
        subdir_test = [os.path.basename(s) for s in subdir]
        print("Found subdirectories: " + str(subdir_test))
        #print(os.path.basename(subdir[0])[0])
        
        if pre == "":
            #If pre is empty, valid subdirectories start with a digit.
            subdir = [sub for sub in subdir if str(os.path.basename(sub))[0].isdigit()]

        print("Found subdirectories: " + str(subdir))
        if dynamic_inspect:
            input("Press Enter to continue...")
    
        for g in L:

            # In each selected subdir, all files are  of the form 
            # "FRAME_RAND*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*jiD_{jID}_*_R_{R}_*.csv" if the file contains jiD data.
            # "FRAME_RAND*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv" if the file does not contain jiD data.

            # First make a sorted list of a_val in each subdir.

            for sub in subdir:
                #Get all files in parent directory that begin with the above (but not subdirectories).
                files = glob.glob(sub + "/*_G_" + str(g) + "_T_*.csv")
                #Find all unique values of a_val in files (in ascending order) (using regex with finding "a_{a_val}" and then removing a_).
                a_vals = [re.findall(r'a_[\d]*[.][\d]+' , f)[0]
                            if re.findall(r'a_[\d]*[.][\d]+' , f) 
                            else re.findall(r'a_[\d]+' , f)[0] if re.findall(r'a_[\d]+' , f) else "a_0" for f in files]
                #print("For subdir: " + sub)
                #print("Found a_vals: " + str(a_vals))
                print(len(a_vals))
                # Remove a_ from a_vals and sort them in ascending order.
                a_vals = sorted([float(re.findall(r'[\d]*[.][\d]+', a)[0]) if re.findall(r'[\d]*[.][\d]+', a) 
                                 else int(re.findall(r'[\d]+', a)[0]) for a in a_vals])
                
                # Only consider unique values of a.
                a_vals =np.unique(a_vals)

                print("For subdir: " + sub)
                print("Found a_vals: " + str(a_vals))

                # Check if out_dir exists in savedir.txt in sub.
                # If it does, provide a warning followed by a prompt to skip ("skip"/1) or continue ("continue"/2).

                
                if os.path.exists(sub + "/savedir.txt"):
                    with open(sub + "/savedir.txt", "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            line.strip("\n")
                            if out_dir + "/L_" +str(g) + "_*" in line:
                                print("WARNING!!!: " + out_dir + " already exists in " + sub + "/savedir.txt")
                                if dynamic_inspect:
                                    print("Do you want to skip this directory ('skip'/1) or continue overwriting ('continue'/2)?")
                                    choice = input("Enter choice: ")
                                    if choice == "skip" or choice == "1":
                                        continue
                                    elif choice == "continue" or choice == "2":
                                        break
                                    else:
                                        print("Invalid choice. Skipping directory.")
                                        continue
                    f.close()
                                
                for a in a_vals:
                    a = int(a) if float(a).is_integer() else a
                    savedir = "L_" + str(g) + "_a_" + str(a) + "/dP_" + str(dP) + "/Geq_" 
                    if Geq == "NA" and Veq != "NA":
                        savedir += str(Veq) + "/"
                    else:
                        savedir += str(Geq) + "/"
                    outdir = out_dir + savedir
                    # Recursively create outdir (and parent directories) if it doesn't exist.
                    Path(outdir).mkdir(parents=True, exist_ok=True)

                    t_range = find_timerange(files, a, indx_vals_t, tmin, tmax)

                    if (0 not in t_range):
                        print("Warning: 0 not in t_range. Adding 0 to t_range.")
                        t_range = [0] + t_range

                    print(" For a = " + str(a) + " Found t_range: " + str(t_range))

                    # Finding jID values in files which contain a_{a} in their name.
                    jiD_range= [ re.findall(r'jID_[\d]+' , f)[0] if re.findall(r'jID_[\d]+' , f) else "jID_-1" for f in files
                                if (re.findall(r'a_[\d]*[.][\d]+' , f) and  re.findall(r'a_[\d]*[.][\d]+' , f)[0] == "a_" + str(a))
                                or (re.findall(r'a_[\d]+' , f) and  re.findall(r'a_[\d]+' , f)[0] == "a_" + str(a))]
                    # Extract jiD values from jiD_range.
                    jiD_vals = sorted([int(re.sub("jID_", "", j)) for j in jiD_range]); 
                    # Remove duplicates from jiD_vals, and append -1 at the end to indicate possible frames with no jiD data.
                    jiD_vals = np.unique(jiD_vals); jiD_vals = np.append(jiD_vals, -1)
                    
                    print("For a = " + str(a) + " Found jID_vals: " + str(jiD_vals))

                    print("Processing: " + savedir)
                    
                    for jiD in jiD_vals:    # Loop over jiD values first.
                        for t in t_range:   # Loop over t values next.
                            #Get all files in parent directory that begin with the above (but not subdirectories).
                            if jiD != -1:
                                selectedfiles = glob.glob(sub + "/*_G_" + str(g) + "_T_" + str(t) +"*_a_" + str(a) + "*_jID_" + str(jiD) + "_*.csv")
                            else:
                                # RECALL- A file with jiD = -1 indicates that the file does not contain unique jiD data.
                                selectedfiles = [f for f in glob.glob(sub + "/*_G_" + str(g) + "_T_" + str(t) +"*_a_" + str(a) + "_*.csv") 
                                                 if "jID" not in f]
                            if len(selectedfiles) == 0:
                                continue
                            t_subdir = outdir + "/T_" + str(t) + "/"
                            # If selectedfiles is empty, skip the loop over t.
                            # Recursively create t_subdir (and parent directories) if it doesn't exist.
                            Path(t_subdir).mkdir(parents=True, exist_ok=True)

                            print(f"Processing jID: {jiD} at t= {t}")
                            #print("Processing t= " + str(t))
                            
                            # Copy all files in files to t_subdir, with new names and avoiding conflicts.
                            # If a file with the same name exists in t_subdir, rename the file to avoid conflicts.
                            # This renaming will be done by sending the conflicting file name to a function that will rename it.

                            for source_file in selectedfiles:
                                source_filename = os.path.basename(source_file)
                                # Change source_filename as follows:
                                # Source file name is of the form "FRAME_RAND*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv"
                                # Change it to "FRAME_T_{T}_a_{a_val}_R_{R}.csv"

                                # First extract T, a_val, R from source_filename.
                                T = t; a_val = a; R = int(re.findall(r'[\d]+' , source_filename)[-1])
                                # Find filetype of source_filename, which are the characters before the first "_"
                                filetype = source_filename.split("_")[0]
                                # Now create new filename.
                                out_filename = filetype + "_T_" + str(T) + "_a_" + str(a_val) + "_R_" + str(R) + ".csv"
                                '''
                                if re.search("GAMMA", source_filename):
                                    out_filename = "GAMMA_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
                                else:
                                    out_filename = "FRAME_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
                                '''
                                outfile = t_subdir + out_filename
                                #Check if outfile already exists.
                                #fd_source = os.open(source_file, os.O_RDONLY);
                                #print("Status of source file: " + str(os.fstat(fd_source)))
                                if os.path.exists(outfile):
                                    print("Warning: File " + out_filename + " already exists, and conflicts with source file: "
                                        + source_filename + ". Renaming file.")
                                    # First check if outfile is the same as source_file (this is done by looking at the attributes of the files (filesize etc)).
                                    # If they are the same, skip the file.
                                    for out_f in glob.glob(t_subdir + "*.csv"):

                                        if os.path.samefile(source_file, out_f):
                                            print("WARNING!!!: File " + outfile + " is the same as " + source_file + ". Skipping file.")
                                            continue

                                        fd_source = os.open(source_file, os.O_RDONLY);
                                        fd_outfile = os.open(out_f, os.O_RDONLY);
                                        if os.fstat(fd_source).st_ino == os.fstat(fd_outfile).st_ino:
                                            print("WARNING!!!: File " + outfile + " is the same as " + source_file + ". Skipping file.")
                                            continue

                                        os.close(fd_source)
                                        os.close(fd_outfile)
                                        
                                    # If outfile already exists, rename it to avoid conflicts.
                                    max_R = max([int(re.findall(r'[\d]+' , f)[-1]) for f in glob.glob(t_subdir + filetype + "*.csv")])
                                    out_filename = filetype + "_T_" + str(T) + "_a_" + str(a_val) + "_R_" + str(max_R + 1) + ".csv"
                                    outfile = t_subdir + out_filename

                                # Copy source_file to outfile.
                                try:
                                    shutil.copy2(source_file, outfile)
                                except shutil.Error as e:
                                    print("Error: Could not copy file " + source_file + "\tto:\n" + outfile + "with error message: \n" + str(e))
                                    print("Skipping file.")
                                    continue
                                except FileNotFoundError as e:
                                    print("Error: Could not find file " + source_file + "  with error message: \n" + str(e))
                                    print("Skipping file.")
                                    continue
                            # This is the end of the loop over source_file in selectedfiles.

                        # This is the end of the loop over t in t_range.
                    # This is the end of the loop over jiD in jiD_vals.
                    #Create a t-vals.txt file in outdir.
                    test_txt_RW(outdir, "*/", "T", ["[\d]*[.][\d]+", "[\d]+"])

                # This is the end of the loop over a in a_vals.
                # Now appending list of savedir in a txt file in sub.
                
                #Open a file in sub (in append mode), and append savedir to it.
                try:
                    with open(sub + "/savedir.txt", "a") as f:
                        f.write(out_dir + "L_" + str(g) + "_*" +"\n")
                    f.close()
                except Exception as e:
                    print("Error: Could not write savedir to " + sub + "/savedir.txt with error message: \n" + str(e))
                
        print("Finished processing prefix: " + pre)
        

        # Check if out_dir exists.
        if os.path.exists(out_dir):
            print("Adding a-val.txt to " + out_dir)
            test_txt_RW(out_dir, "*/", "a", ["[\d]*[.][\d]+", "[\d]+"])

        # Read a_vals.txt in out_dir and print the values.
        if os.path.exists(out_dir + "/a_vals.txt"):
            print("Reading a_vals.txt in " + out_dir)
            with open(out_dir + "/a_vals.txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    print(line.strip())
            f.close()
    
    # This is the end of the loop over pre in prefixes.
    # Now post-process the data files in out_dir_noprefix.
    '''
    print("\n\n__________________________________________________________________________________________\n\n")
    print("Now, Post-processing data files in " + out_dir_noprefix)
    print("\n\n__________________________________________________________________________________________\n\n")
    post_process(prefixes)
    '''


set_frames_input()
#main()
post_process(prefixes)
#dir = "\\\\?\\D:\\cygwin64\\home\\koust\\Code\\Trophic_Maelstorm\\simulations\\Data\\Remote\\Rietkerk\\Frames\\Stochastic\\3Sp\\DDM_DiC_BURNIN_0.025-0.065_dP_30000_Geq_4.802";
#rename(dir, "_RANDDDMDiCBURNIN_ThreeSp_P_c_DP", "")

#dir = "../Data/Remote/Test_Rietkerk_Frames/Stochastic/3Sp/Reorganised/StandardParam_20_100/DiC_BURNIN_/"
#test_txt_RW(dir)