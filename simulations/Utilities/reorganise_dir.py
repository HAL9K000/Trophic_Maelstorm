import numpy as np
import os
import pandas as pan
import regex as re
from pathlib import Path
import glob
import shutil
import sys
import copy

import scipy.stats as stats
from scipy.interpolate import CubicSpline

'''
This script reorganises the directory structure of the data files in the Rietkerk model.
The original directory structure (filepaths) are as follows:
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv

The new directory structure (filepaths) will be as follows:
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv

The script will:
1. Find all immediate subdirectories in root_dir.
2. For each subdirectory, find all unique values of a_val in the files in the subdirectory.
3. For each value of a_val, find all unique values of T in the files in the subdirectory.
4. For each value of T, copy all files in the subdirectory to the new directory structure in out_dir.
5. If a file with the same name already exists in the new directory, rename the file to avoid conflicts.

This script ASSUMES that the files in the original directory structure have the following format,
and WILL NOT WORK if the files are named differently:
I. Contains "_G_{L}_T_{Time}_" in the filename, which is followed by
II. "_a_{a_val}_" in the filename, which is followed by
III. "_R_{R}.csv" in the filename.

The script accepts the following arguments:
1. prefixes: A list of prefixes to be used in the subdirectory names in out_dir.
2. root_dir: The root directory containing the original data files.
3. out_dir_noprefix: The output directory where the reorganised data files will be stored.
4. dP: The value of dP in the subdirectory names in root_dir.
5. Geq: The value of Geq in the subdirectory names in root_dir (if "NA", Geq is not used in the subdirectory name).
6. L: A list of values of L to be used in the subdirectory names in out_dir.
7. indx_vals_t: Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.

'''

prefixes =["DiC-REF-HI", "DiC-REF-LI", "DiC-NREF-HI", "DiC-NREF-LI"]#, "DiC-NEW"]
#prefixes =["", "DiC", "BURNIN", "DiC-BURNIN", "DDM-DiC", "DDM-DiC-BURNIN"]
root_dir = "../Data/Remote/Rietkerk/Frames/Stochastic/3Sp/Test/"
#out_dir_noprefix = "../Data/Rietkerk/Reorganised_Frames/Stoc/3Sp/StdParam_20_100_Test/"

out_dir_noprefix = "../Data/Remote/Rietkerk/Reorg_Frames/3Sp/StdParam_20_100_EQ"

dP = 10000
Geq = 4.802 # Optional. If Geq is not used in the subdirectory name, set Geq = "NA".
L= [128]
tmin = 60000; tmax = None; 
# Optional. If tmin and tmax are not provided, set tmin = None and tmax = None.
indx_vals_t = 26
#Extract n largest values of T (within the subset [tmin, tmax] if provided) if indx_vals_t = -n, 
# n smallest values of T if indx_vals_t = n.


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
print("L: " + str(L))
print("indx_vals_t: " + str(indx_vals_t))
print(f"Tmin:   {tmin} \t Tmax:   {tmax}")

choice = input("Enter choice: ")
if choice == "y":
    out_dir_noprefix = input("Enter new out_dir_noprefix: ")
    dP = int(input("Enter new dP: "))
    Geq = input("Enter new Geq ('NA' where it isn't applicable): ")
    L = [int(x) for x in input("Enter new L: ").split()]
    indx_vals_t = int(input("Enter new indx_vals_t: "))
    tmin = float(input("Enter new tmin: "))
    tmax = float(input("Enter new tmax: "))


def find_timerange(files, a, indx= indx_vals_t, min_t = None, max_t = None):
    # Files are a list of files in a directory.
    # of the form "FRAME_RAND*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv"
    # a is the value of a in the above.
    # indx provides the number of t-values to extract from the file name.
    # Negative values of indx (say -n) indicate the n largest values of T are to be extracted.
    # Positive values of indx (say n) indicate the n smallest values of T are to be extracted.
    # Return a list of T values in files.

    # First, extract all T values from files that have a_val = a. (t-values can be float or integer).
    T_vals = [re.findall(r'T_[\d]*[.][\d]+' , f)[0] if re.findall(r'T_[\d]*[.][\d]+' , f) 
              else re.findall(r'T_[\d]+' , f)[0] for f in files 
              if re.findall(r'a_[\d]*[.][\d]+' , f)[0] == "a_" + str(a) or re.findall(r'a_[\d]+' , f)[0] == "a_" + str(a)]
    # Remove T_ from T_vals and sort them in ascending order.
    T_vals = sorted([float(re.findall(r'[\d]*[.][\d]+', t)[0]) 
                     if re.findall(r'[\d]*[.][\d]+', t) 
                     else int(re.findall(r'[\d]+', t)[0]) for t in T_vals])
    # Only consider unique values of T.
    T_vals = sorted(list(set(T_vals)))
    # If min_t and max_t are provided, extract T values in the range [min_t, max_t].

    if np.abs(indx) > len(T_vals):
        print("Warning: indx greater than number of T values. Returning all T values.")
        return T_vals
    
    if min_t is not None:
        T_vals = [t for t in T_vals if t >= min_t]
    if max_t is not None:
        T_vals = [t for t in T_vals if t <= max_t]
    
    if np.abs(indx) > len(T_vals):
        print("Warning: indx greater than number of T values (as per constraints). Returning all T values (as per constraints).")
        return T_vals
    return T_vals[indx:] if indx < 0 else T_vals[:indx]

def find_vals(files, reg_search_patterns, vals, dir=""):
    # Recursively iterate over reg_search_patterns and find all values in files that match the pattern.
    # Append the values to vals.
    files_copy = copy.deepcopy(files)
    if len(reg_search_patterns) == 0:
        return
    else:
        for f in files:
            vals += [re.findall(reg_search_patterns[0], f)[0] for f in files_copy if re.findall(reg_search_patterns[0], f)]
        # Remove matching values from files.
        files_copy = [f for f in files_copy if not re.findall(reg_search_patterns[0], f)]
        find_vals(files_copy, reg_search_patterns[1:], vals, dir)



def test_txt_RW(parendir, type_id ="", reg_search_patterns = ["[\d]*[.][\d]+", "[\d]+"], savefile = "default", verbose= False):

    # Test reading and writing text files.
    # Dir has directories of the form "L_{G}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/"
    # Extract all vals in dir that match regex_search_patterns, with type_id as identifying prefix (if not empty) (such as "a" or "L" or "T" etc.)
    # and save them to a text file in the same directory (tab delineated), with the name "savefile"

    # First make a list of all immediate subdirectories in dir.

    files = glob.glob(os.path.join(parendir, "*/"))
    #print("Found files: " + str(files))
    #Find all unique values of vals in files (in ascending order) (using regex with finding "a_{a_val}" and then removing a_).
    if type_id != "":
        reg_search_type_patterns = [type_id + "_" + reg for reg in reg_search_patterns]
    else:
        reg_search_type_patterns = reg_search_patterns
    # Recursively iterate over regex_search_patterns and find all values in files that match the pattern.
    vals=[]
    find_vals(files, reg_search_type_patterns, vals, parendir)
    #print(f"Found {type_id} vals: " + str(vals))
    # Remove type_id from vals and sort them in ascending order if type_id is not empty.
    vals = sorted([float(re.sub(type_id + "_", "", v)) for v in vals])
    vals = np.unique(vals)
    print(f"Found {type_id} vals: " + str(vals))
    # Save vals to a text file in parendir (tab delineated), with the name savefile.
    # Int values are saved as integers, float values are saved as floats.

    if savefile == "default":
        savefile = type_id + "_vals.txt"
    try:
        np.savetxt(parendir + "/" + savefile, vals, delimiter="\t", fmt="%g")
        # Test reading the text file line by line.
        if(verbose):
            with open(parendir + "/" + savefile, "r") as f:
                lines = f.readlines()
                for line in lines:
                    print(line.strip())
                # Close the file.
                f.close()
    except Exception as e:
        print("Error: Could not write vals to " + parendir + "/" + savefile + " with error message: \n" + str(e))
        return

''' # Summary of the function gen_MEAN_SD_COLSfiledata(...)
# This function generates the mean and standard deviation of each of the columns in the files.
# excluding exclude_col_labels. Files are standardised to have the same names of columns (but not necessarily in the same order).
# If the files have non-standardised columns (different names other than exclude_col_labels), the function will return an error.
#If nonzero is True, averages and means are calculated over only non-zero columns across all files, else over all columns.
# If add_counts is True, the function will also add a column to the dataframe next to the mean and standard deviation columns for each column-labels
# that contains the number of non-zero columns for that column-label across all files. If nonzero is False, the added columns will simply contain the number of columns
# for that column-label across all files (which will be the same as the number of files).
'''

def gen_MEAN_SD_COLSfiledata(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], nonzero=True, add_counts= True):

    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    if(add_counts):
        col_names = [ y + x + "]" for x in col_labels  for y in ["AVG[", "VAR[", "COUNT["]]
    else:
        col_names = [ y + x + "]" for x in col_labels  for y in ["AVG[", "VAR["] ]
    
    # Create a df with col_names as column names.
    df_out = pan.DataFrame(columns= col_names)
    # Pooled Data Stats will store the pooled average, standard deviation, sample size and counts for each column in col_labels across all files.
    # It will be iteratively updated as the function iterates over all files in files.
    pooled_data_stats = pan.DataFrame(columns = [y + x for x in col_labels 
            for y in ["POOLED_AVG ", "POOLED_VAR ", "POOLED_SAMPLE_SIZE ", "POOLED_COUNTS "]])
    # Initialise the pooled_data_stats to 0.
    pooled_data_stats.loc[0] = 0
    # Iterate over all files in files.
    for file in files:

        # Read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]
        
        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        for col in col_labels:

            # Find mean, variance and sample size of each column in df.
            mean_col = df[col].mean()
            if nonzero:
                if mean_col == 0:
                    continue

            variance_col = df[col].var()
            # Store sample size in ss
            ss= len(df[col])

            # Update Pooled Data Stats iteratively.
            pooled_data_stats.loc[0, "POOLED_AVG " + col] = (pooled_data_stats.loc[0, "POOLED_AVG " + col] * pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] 
                    + mean_col * ss)/(pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] + ss)
            pooled_data_stats.loc[0, "POOLED_VAR " + col] =  (pooled_data_stats.loc[0, "POOLED_VAR " + col] * (pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] -1)
                    + variance_col * (ss -1))/(pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] + ss -2)
            pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] += ss
            pooled_data_stats.loc[0, "POOLED_COUNTS " + col] += 1

    # Done with all files in files.

    # Store the pooled_data_stats in df_out.
    for col in col_labels:
        df_out.loc[0, "AVG[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_AVG " + col]
        df_out.loc[0, "VAR[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_VAR " + col]
        if add_counts:
            df_out.loc[0, "COUNT[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_COUNTS " + col]
    
    return df_out


'''# Summary of the function gen_MEAN_INDVL_Colsfiledata(...)
# This function generates the mean for each column for each file in files, as well as the mean across all files for each column.
# It will exclude exclude_col_labels from the mean calculations. If the files have non-standardised columns 
# (different names other than exclude_col_labels), the function will return an error.
# It will create a dataframe with the following column structure:
#R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
# where {var} is one of the species names NOT in exclude_col_labels.
'''
def gen_MEAN_INDVL_Colsfiledata(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"]):
    
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    # Create a df with col_labels as column names. This will store the mean data, fot each column in col_labels for each file in files.
    # Create another df which will concatenate the mean data across all files.
    df_file = pan.DataFrame(columns= ["AVG[" + x + "]" for x in col_labels])
    df_pooled = pan.DataFrame()
    # Initialise the pooled_data_stats to 0.
    # Iterate over all files in files, sorted in ascending order of R.
    for file in sorted(files, key= lambda x: int(re.sub("R_", "", re.findall(r'R_[\d]+', x)[0]))):
        # Read the file.
        Rstr = re.findall(r'R_[\d]+', file)[0]
        #R = int(re.search(r'R_[\d]+', file).group(0).split("_")[1])
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None

        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        # Find mean of each column in df.
        for col in col_labels:
            mean_col = df[col].mean()
            # Store the mean in df_file.
            df_file.loc[0, "AVG[" + col + "]"] = mean_col
            # Store the mean in df_pooled.
            df_pooled.loc[0, "AVG[" + col + "]_" + Rstr] = mean_col
    
    # Done with all files in files.
    # Finally insert mean across all files for each column in col_labels in df_pooled at the start.
    for i in range(len(col_labels)):
        df_pooled.insert(i, "AVG[" + col_labels[i] + "]_ALL", df_pooled.filter(like= "AVG[" + col_labels[i] + "]_").mean(axis=1))

    # Also insert mean across all files for each column in col_labels in df_pooled that are non-zero at the start.
    for i in range(len(col_labels)):
        df_pooled.insert(i, "AVG[" + col_labels[i] + "]_SURV", df_pooled.filter(like= "AVG[" + col_labels[i] + "]_R").replace(0, np.nan).mean(axis=1))
    
    # Insert max R at the start of the dataframe.
    df_pooled.insert(0, "R", len(files))
    # Replace all NaNs with 0.
    df_pooled = df_pooled.fillna(0)
    return df_pooled

""" # Summary of the function gen_potential_well_data(...)
The following function gen_potential_well_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], evaluate_local_minima =True)
will generate the potential well data from the files in files.
It will read all files in files, and will concatenate matching columns in the files to generate the potential well data.
After concatenating the columns across all files, for each column, it will use stats.gaussian_kde to generate a kernel density estimate of the data.
It will then find the potential well for each column as potential well = -np.log10(kde(x)), where x is the data in the column.
If evaluate_local_minima is True, it will evaluate the local minima of the potential well data for each column.
It will do so by using the CubicSpline function in scipy.interpolate to find all the local minima of the potential well data.
In case of an entirely flat potential well, the function will return a warning and return local minima for that column as None.
It will then store the local minima in a different dataframe, along with the corresponding x-values.
The function will then return the potential well dataframe and (optionally) the local minima dataframes.
"""
def gen_potential_well_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], evaluate_local_minima= True):
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    kde_cols = [y + x + "]" for x in col_labels for y in ["KDE[", "POTENTIAL_WELL["]]
    # Create a df with col_labels as column names. This will store the KDE and potential well data.
    df_kde = pan.DataFrame(columns= kde_cols)
    if(evaluate_local_minima):
        min_cols = ["LOCAL_MINIMA[" + x + "]" for x in col_labels]
        # Create a df with col_labels as column names. This will store the local minima data.
        df_local_minima = pan.DataFrame(columns= min_cols)
    
    concat_df = pan.DataFrame()
    # Iterate over all files in files. For each file, concatenate the columns in col_labels across all files.
    # Store the concatenated data in concat_df.

    for file in files:

        #First read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        concat_df = pan.concat([concat_df, df[col_labels]], axis=0)

    #print(concat_df.head())
    #print(concat_df.tail())
    print(concat_df.shape)
    # Done with all files in files.

    # Now generate the KDE and potential well data for each column in col_labels.
    for col in col_labels:
        print("Evaluating column: " + col)
        #First check if the column is entirely zero. If it is, skip the column.
        if all(concat_df[col] == 0):
            print("Warning: Entirely zero column across all replicates in " + col + f". Skipping column.")
            continue
        # Check if the column is entirely flat, i.e. std < 1e-6. If it is, skip the column.
        if concat_df[col].std() < 1e-6:
            print("Warning: Entirely flat column across all replicates in " + col + f". Skipping column.")
            continue

        # Evaluate this over a certain range of x-values for each column.
        eval_range = np.linspace(concat_df[col].min(), concat_df[col].max(), 300)

        # Find the kernel density estimate of the data in col.
        kde = stats.gaussian_kde(concat_df[col]); kernel = kde.evaluate(eval_range)
        # Report indices of kernel where kernel values are zero.
        if np.any(kernel == 0):
            print("Warning: Kernel is zero at indices: " + str(np.where(kernel == 0)))
        
        # Find the potential well for the data in col, accounting for zero values in kernel.
        potential_well = -np.log10(kernel + 1e-10)
        # Store the KDE and potential well data in df_kde.
        
        df_kde["KDE[" + col + "]"] = pan.Series(kernel)
        df_kde["POTENTIAL_WELL[" + col + "]"] = pan.Series(potential_well)

        #print(df_kde.head())
        #print(df_kde.tail())
        print(df_kde.shape)
        #print(df_kde.describe())

        if(evaluate_local_minima):
            # Find the local minima of the potential well data in col.
            # Use CubicSpline to find the local minima.
            local_minima = []
            cs = CubicSpline(eval_range, potential_well)
            # Find the local minima.
            cs_1st_derivative = cs.derivative()
            cs_2nd_derivative = cs.derivative(nu=2)
            extrema = cs_1st_derivative.roots()
            for x in extrema:
                if cs_2nd_derivative(x) > 0:
                    local_minima.append(x)
            # If the potential well is entirely flat, return a warning and set local minima to None.
            if len(local_minima) == 0:
                print("Warning: Entirely flat potential well for column " + col + ". Local minima set to NAN.")
                local_minima.append(np.nan)
            # Store the list of local minima in df_local_minima.
            df_local_minima["LOCAL_MINIMA[" + col + "]"] = pan.Series(local_minima)

    # Done with all columns in col_labels.     
    return (df_kde, df_local_minima) if evaluate_local_minima else (df_kde, None)


""" # Summary of the function post_process(...)
# A wrapper function that post-processes the data files in out_dir_noprefix, after the main() function has been executed/
# Out_dir_noprefix is the output directory where the reorganised data files are stored.
# It is structured as follows:
# out_dir_noprefix/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
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

            # Find mean and standard deviation of each column in each file in files if the column values are non-zero.
            ''' #Mean and Std Rho Density & Max Replicates
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
            
            # Find
            # Find max R in files.
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
            '''

            # Find Mean for each column in each file in files and the mean across all files for each column using gen_MEAN_INDVL_Colsfiledata.
            # Save this df to a txt file in subdirpath.
            df_replicates = gen_MEAN_INDVL_Colsfiledata(files, pathtodir=subdirpath, ext="csv")
            try:
                if df_replicates is not None:
                    df_replicates.to_csv(subdirpath + "/MEAN_REPLICATES.txt", sep="\t", index=False, header=True)
            except Exception as e:
                print("Error: Could not write MEAN_REPLICATES to " + subdirpath + "/MEAN_REPLICATES.txt with error message: \n" + str(e))

            '''
            # Find potential well data in files using gen_potential_well_data.
            # Save the potential well data to a csv file in subdirpath/Pot_Well.
            # If evaluate_local_minima is True, also save the local minima data to a txt file in subdirpath.
            df_kde, df_local_minima = gen_potential_well_data(files, pathtodir=subdirpath, ext="csv", exclude_col_labels= ["a_c", "x", "L"], evaluate_local_minima= True)
            Path(subdirpath + "/Pot_Well").mkdir(parents=True, exist_ok=True)
            try:
                df_kde.to_csv(subdirpath + "/Pot_Well/Pot_Well.csv", sep=",", index=False, header=True)
                if df_local_minima is not None:
                    df_local_minima.to_csv(subdirpath + "/LOCAL_MINIMA.txt", sep="\t", index=False, header=True)
            except Exception as e:
                print("Error: Could not write potential well data to " + subdirpath + "/Pot_Well/Pot_Well.csv with error message: \n" + str(e))
            ''' 
        # Done with all subdirectories in out_dir_noprefix/{pre} that contain files of the form "FRAME*.csv".
        #print("\n=====================================================================================================\n")
        print(f"Done post-processing for prefix: {pre}...")
        #print("\n=====================================================================================================\n")
            

def main():
    
    #Find a list of immediate subdirectories in root_dir
    subdir = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(root_dir, '*/'))]
    
    for pre in prefixes:

        out_dir = out_dir_noprefix;
        if pre == "":
            out_dir = out_dir_noprefix + "NoPrefix/"
        else:
            out_dir = out_dir_noprefix + pre + "/"
        
        #Subdirectories are of the form "PREFIX*_dP_{dP}_Geq_{val}"
        #Using regex and glob, create a sorted list of subdir (in ascending order of val (can be float or integer)) with prefix pre.
        
        if(Geq == "NA"):
            subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP))), 
                            key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
        else:
            subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP) + "_Geq_" + str(Geq))), 
                        key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
        
        subdir_test = [os.path.basename(s) for s in subdir]
        print("Found subdirectories: " + str(subdir_test))
        #print(os.path.basename(subdir[0])[0])
        
        if pre == "":
            #If pre is empty, valid subdirectories start with a digit.
            subdir = [sub for sub in subdir if str(os.path.basename(sub))[0].isdigit()]

        print("Found subdirectories: " + str(subdir))
        input("Press Enter to continue...")
    
        for g in L:

            # In each selected subdir, all files are  of the form
            # "FRAME_RAND*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv"

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
                    savedir = "L_" + str(g) + "_a_" + str(a) + "/dP_" + str(dP) + "/Geq_" + str(Geq) + "/"
                    outdir = out_dir + savedir
                    # Recursively create outdir (and parent directories) if it doesn't exist.
                    Path(outdir).mkdir(parents=True, exist_ok=True)

                    t_range = find_timerange(files, a, indx_vals_t, tmin, tmax)

                    if (0 not in t_range):
                        print("Warning: 0 not in t_range. Adding 0 to t_range.")
                        t_range = [0] + t_range

                    print("Processing: " + savedir)

                    for t in t_range:
                        #Get all files in parent directory that begin with the above (but not subdirectories).
                        selectedfiles = glob.glob(sub + "/*_G_" + str(g) + "_T_" + str(t) +"*_a_" + str(a) + "_*.csv")
                        t_subdir = outdir + "/T_" + str(t) + "/"
                        # Recursively create t_subdir (and parent directories) if it doesn't exist.
                        Path(t_subdir).mkdir(parents=True, exist_ok=True)
                        print("Processing t= " + str(t))
                        
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

                    #Create a t-vals.txt file in outdir.
                    test_txt_RW(outdir, "T", ["[\d]*[.][\d]+", "[\d]+"])

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
            test_txt_RW(out_dir, "a", ["[\d]*[.][\d]+", "[\d]+"])

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
        

def rename(dir, keyword, new_keyword):
    # Rename all files in dir (including the entire filepath) that contain keyword to new_keyword.
    # This function is used to rename files that conflict with other files in the directory.
    # This function is called by main() when a file with the same name already exists in the new directory.
    # The function will rename the file to avoid conflicts.
    # The function will rename the file by replacing keyword with new_keyword in the filename.
    # The function will rename the file by adding a number to the end of the filename if new_keyword is already present in the filename.

    files = glob.glob(dir + "\\*" + keyword + "*")
    for file in files:
        new_file = file.replace(keyword, new_keyword)
        if os.path.exists(new_file):
            new_file = new_file.replace(new_keyword, new_keyword + "_1")
        os.rename(file, new_file)



#main()
post_process(prefixes)
#dir = "\\\\?\\D:\\cygwin64\\home\\koust\\Code\\Trophic_Maelstorm\\simulations\\Data\\Remote\\Rietkerk\\Frames\\Stochastic\\3Sp\\DDM_DiC_BURNIN_0.025-0.065_dP_30000_Geq_4.802";
#rename(dir, "_RANDDDMDiCBURNIN_ThreeSp_P_c_DP", "")

#dir = "../Data/Remote/Test_Rietkerk_Frames/Stochastic/3Sp/Reorganised/StandardParam_20_100/DiC_BURNIN_/"
#test_txt_RW(dir)