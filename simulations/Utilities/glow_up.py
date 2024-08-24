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

import scipy.stats as stats
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline

'''
This script contains a set of utility functions that can be used to read and write text files, read and write csv files,
and generally manipulate/reorganise data in given directories. It is intended to be invoked from reorganise_dir.py and other scripts
in the Utilities directory.

A detailed description of each function is provided below.
'''


'''# Summary of the function find_timerange(...)
    This function finds all T values in provided "files" (a list of files that may or may not be from the same directory)
    that have a_val = a. (t-values can be float or integer). It then extracts the "indx" smallest or largest values of T 
    ( if indx is negative/ positive respectively) from the list of T values.
    It then returns the list of T values that are in the range [min_t, max_t] (if min_t and max_t are provided).
'''
def find_timerange(files, a, indx, min_t = None, max_t = None):
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

''' # Summary of the function savetxt_to_dir(...)
    This function saves a dataframe or an array-like-object to a directory (savedir), with name savefilename.
    If the directory does not exist, it creates the directory. If the file already exists in the directory, 
    it renames the file to avoid conflicts.
'''
def savetxt_to_dir(data, savedir, savefilename, header = None, index = False, sep = "\t", overwrite = True, verbose = False):
    
    if not os.path.exists(savedir):
        Path(savedir).mkdir(parents=True, exist_ok=True)
    # Check if the file already exists in the directory.
    if os.path.exists(savedir + "/" + savefilename):
        # Provide a warning before overwriting the file.
        if not overwrite:
            print("Warning: File " + savedir + "/" + savefilename + " already exists. Skipping file.")
            return
    try:
        # Check if data is a dataframe or an array (or list-like object).
        if isinstance(data, pan.DataFrame):
            data.to_csv(savedir + "/" + savefilename, header = header, index = index, sep = sep, fmt="%g")
        else:
            np.savetxt(savedir + "/" + savefilename, data, delimiter = sep, fmt="%g")

        if verbose:
            with open(savedir + "/" + savefilename, "r") as f:
                lines = f.readlines()
                for line in lines:
                    print(line.strip())
                # Close the file.
                f.close()
    except Exception as e:
        print("Error: Could not write vals to " + savedir + "/" + savefilename + " with error message: \n" + str(e))
        return

''' # Summary of the function test_txt_RW(...)
    The function test_txt_RW(...) are used to read and write text files from a directory (parendir) by iterating over 
    some set of files in parendir. BY DEFAULT  test_txt_RW() iterates over all immediate subdirectories in parendir, but 
    can be used to iterate over other sets of files as well in parendir (BY CHANGING file_type).
    test_txt_RW(...) is used to read and write text files from a directory (parendir).
    If  set_of_files is provided, the function will iterate over the files in set_of_files, else it will iterate over the files in parendir as specified by file_type.
    Both functions extract all vals in files that match regex_search_patterns, with type_id as identifying prefix (if not empty)
    (such as "a" or "L" or "T" etc.) and save them to a text file in the same directory (tab delineated), with the name "savefile".
'''

def test_txt_RW(parendir, file_type="*/", type_id ="", reg_search_patterns = ["[\d]*[.][\d]+", "[\d]+"], 
                set_of_files= None,  savefile = "default", verbose= False):

    # Extract all vals in dir that match regex_search_patterns, with type_id as identifying prefix (if not empty) (such as "a" or "L" or "T" etc.)
    # and save them to a text file in the same directory (tab delineated), with the name "savefile"

    # First make a list of all immediate subdirectories in dir.
    if set_of_files is not None:
        files = set_of_files
    else:
        files = glob.glob(os.path.join(parendir, file_type))
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

    if savefile == "default":
        savefile = type_id + "_vals.txt"
    savetxt_to_dir(vals, parendir, savefile, header = None, index = False, sep = "\t", overwrite = True, verbose = verbose)

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


''' # Summary of the function quadratic_spline_roots(...)
    The function quadratic_spline_roots(spl) finds the roots of a quadratic spline function spl.
    It does so by finding the roots of the quadratic spline function in each interval [a, b] of the spline.
    Stolen from: https://stackoverflow.com/questions/50371298/find-maximum-minimum-of-a-1d-interpolated-function
'''
def quadratic_spline_roots(spl):
    roots = []
    knots = spl.get_knots()
    for a, b in zip(knots[:-1], knots[1:]):
        u, v, w = spl(a), spl((a+b)/2), spl(b)
        t = np.roots([u+w-2*v, w-u, 2*v])
        t = t[np.isreal(t) & (np.abs(t) <= 1)]
        roots.extend(t*(b-a)/2 + (b+a)/2)
    return np.array(roots)

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
def gen_potential_well_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], Tmin= None, Tmax = None, evaluate_local_minima= True, bins=300):
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None, None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels];
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    kde_cols = [y + x + "]" for x in col_labels for y in ["BINS[", "KDE[", "POTENTIAL_WELL["]]
    # Create a df with col_labels as column names. This will store the KDE and potential well data.
    df_kde = pan.DataFrame(columns= kde_cols)
    if(evaluate_local_minima):
        min_cols = ["LOCAL_MINIMA[" + x + "]" for x in col_labels]
        # Create a df with col_labels as column names. This will store the local minima data.
        df_local_minima = pan.DataFrame(columns= min_cols)
    
    concat_df = pan.DataFrame()
    # Iterate over all files in files. For each file, concatenate the columns in col_labels across all files.
    # Store the concatenated data in concat_df.

    # Get T and a values from the first file.
    Tval = re.search(r'T_[\d]+', files[0]).group(0).split("_")[1] if re.search(r'T_[\d]+', files[0]) else re.search(r'T_[\d]*[.][\d]+', files[0]).group(0).split("_")[1]
    aval = re.search(r'a_[\d]*[.][\d]+' , files[0]).group(0).split("_")[1] if re.search(r'a_[\d]*[.][\d]+' , files[0]) else re.search(r'a_[\d]+' , files[0]).group(0).split("_")[1]
    # Check if Tmin and Tmax are provided. If they are, only calculate KDE and potential well data for Tval in the range [Tmin, Tmax].
    if Tmin is not None:
        if float(Tval) < Tmin:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None, None
    if Tmax is not None:
        if float(Tval) > Tmax:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None, None
    print(f"Generating potential well data for T = {Tval} and a = {aval}.")
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
            return None, None
        
        concat_df = pan.concat([concat_df, df[col_labels]], axis=0)

    #print(concat_df.head())
    #print(concat_df.tail())
    print(concat_df.shape)
    # Done with all files in files.

    # Now generate the KDE and potential well data for each column in col_labels.
    for col in col_labels:
        is_zero = False; is_flat = False; # Flags to check if the column is entirely zero or flat. 
        #print("Evaluating column: " + col)
        #First check if the column is entirely zero. If it is, skip the column.
        if all(concat_df[col] == 0):
            print("Warning: Entirely zero column across all replicates in " + col + f" for T = {Tval} and a = {aval} ... Adjusting KDE ...")
            # Define eval_range between 0 and 1, with 300 bins. Store the KDE as 1e-10 for all values in eval_range, except 0.
            eval_range = np.linspace(0, 1, bins); kernel = np.ones(bins) * 1e-10; kernel[0] = 1; is_zero = True;
        
        # Check if the column is entirely flat, i.e. std < 1e-6. If it is, skip the column.
        elif concat_df[col].std() < 1e-6:
            print("Warning: Entirely flat column across all replicates in " + col + f"  for T = {Tval} and a = {aval} ... Adjusting KDE ...")
            # Define eval_range b/w min - 1 and max + 1, with 300 bins. Store the KDE as 1e-10 for all values in eval_range except the mean of the column.
            eval_range = np.linspace(concat_df[col].min() - 1, concat_df[col].max() + 1, bins); is_flat = True;
            kernel = np.ones(bins) * 1e-10;
            # Change the kernel value at the index bin closest to the mean of the column to 1.
            kernel[np.argmin(np.abs(eval_range - concat_df[col].mean()))] = 1;
        else:
            # Evaluate this over a certain range of x-values for each column.
            eval_range = np.linspace(concat_df[col].min(), concat_df[col].max(), bins)

            # Find the kernel density estimate of the data in col.
            try:
                kde = stats.gaussian_kde(concat_df[col]); kernel = kde.evaluate(eval_range)
            except stats.NearConstantInputWarning  or stats.ConstantInputWarning as e:
                print(f"Warning: {e} for column " + col + f" for T = {Tval} and a = {aval}. Skipping column.")
                continue

            # Report indices of kernel where kernel values are zero.
            if np.any(kernel == 0):
                #print(f"Warning: Kernel is zero for T = {Tval} and a = {aval} at indices: " + str(np.where(kernel == 0)))
                # Replace zero values in kernel with 1e-10.
                kernel[kernel == 0] = 1e-10
        
        # Find the potential well for the data in col, accounting for zero values in kernel.
        potential_well = -np.log10(kernel)
        # Store the KDE and potential well data in df_kde.
        df_kde["BINS[" + col + "]"] = pan.Series(eval_range)
        df_kde["KDE[" + col + "]"] = pan.Series(kernel)
        df_kde["POTENTIAL_WELL[" + col + "]"] = pan.Series(potential_well)

        #print(df_kde.head())
        #print(df_kde.tail())
        #print(df_kde.shape)
        #print(df_kde.describe())

        if(evaluate_local_minima):
            # Find the local minima of the potential well data in col.
            # Use CubicSpline to find the local minima.
            local_minima = []; min_pot =[];
            # local_minima will store the x-values corresponding to the local minima of the potential well.
            # min_pot will store the potential well values at the local minima (y-values).
            # First check if the potential well is entirely flat.
            if is_zero or is_flat:
                # Append the eval_range corresponding to potential well = 0 to local_minima.
                local_minima.append(eval_range[np.argmin(np.abs(potential_well))])
                min_pot.append(0)
                df_local_minima["LOCAL_MINIMA[" + col + "]"] = pan.Series(local_minima);
                df_local_minima["MIN_POT[" + col + "]"] = pan.Series(min_pot);
                continue

            cs = InterpolatedUnivariateSpline(eval_range, potential_well, k=3)
            # Find the local minima.
            cs_1st_derivative = cs.derivative()
            cs_2nd_derivative = cs.derivative(n=2)
            extrema = quadratic_spline_roots(cs_1st_derivative)
            for x in extrema:
                if cs_2nd_derivative(x) > 0:
                    local_minima.append(x)
                    min_pot.append(cs(x))
            # If the potential well is entirely flat, return a warning and set local minima to None.
            if len(local_minima) == 0:
                print("Warning: Entirely flat potential well for column " + col + ". Local minima set to NAN.")
                local_minima.append(np.nan)
                min_pot.append(np.nan)
            # Store the list of local minima in df_local_minima.
            df_local_minima["LOCAL_MINIMA[" + col + "]"] = pan.Series(local_minima)
            df_local_minima["MIN_POT[" + col + "]"] = pan.Series(min_pot)

    # Done with all columns in col_labels.     
    return (df_kde, df_local_minima) if evaluate_local_minima else (df_kde, None)



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


# =============================== FUNCTIONS SPECIFIC TO PRELIM FILES PROCESSING ============================

# Computes and returns ln-distributed points from t=10^0 to log10(t_max)
    # (the latter rounded down to 1 decimal place)
def logarithm10_time_bins(t_max, dt, size=0):
    
    
    if t_max <= 1:
        return [0.0]  # No points recorded if t_max < e^0

    tend_one = np.log10(t_max*1.01)
    power = np.arange(0.0, tend_one, 0.04)  # 25 time-stamps per decade

    t_fin = np.floor( np.power(10, power) / dt) * dt

    t_fin[t_fin <= 1.0] += dt  # Ensure t_fin is strictly above 1.0

    logspaced = np.unique(t_fin).tolist()

    # Removing time-stamps that are too close to each other (within dt/2.0)
    
    n = len(logspaced)
    for i in range(n - 1, 0, -1):
        if logspaced[i] < logspaced[i - 1] + dt / 2.0:
            logspaced.pop(i)

    
    # Add 0.0 to the list of time-stamps if it is not already present at the start
    if logspaced[0] > 0.0:
        logspaced.insert(0, 0.0)
    
    print(f"Generated {len(logspaced)} logarithmically spaced time-stamps from 0 to {tend_one} (log10({t_max})) with dt = {dt}.")
    #print(f"Time-stamps: {logspaced}")
    if size > 0 and len(logspaced) != size:
        print(f"Warning! The number of time-stamps generated {len(logspaced)} is not equal to the specified size = {size}. Returning the first  {size} time-stamps.")
        return logspaced[:size]
    else:
        return logspaced


'''# Summary of the function gen_MEAN_INDVL_Colsfiledata(...)
# This function generates the mean for each column for each file in files, as well as the mean across all files for each column.
# It will only include include_col_labels from the mean calculations. If the files have non-standardised columns 
# (a label in include_col_labels is not present in the column names of the files), the function will return an error.
# It will create a dataframe with the following column structure:
#R_max  t   AVG[{var}]_SURV   ...   SD[{var}]_SURV  ... R_SURV[{var}]   ...   AVG[{var}]_ALL   ...   SD[{var}]_ALL
# where {var} is one of the species names in include_col_labels (excluding "t").
# The optional argument dt is used to generate a logarithmically spaced time range from 0 to t_max, with time-stamps spaced by dt,
# if the time values in the files are all zero.
'''
def gen_MEAN_INDVL_Prelimsfiledata(files, pathtodir="", ext="csv", tmax =None, dt =None, include_col_labels= ["t",  "<<P(x; t)>_x>_r" , 
                                        "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r", "<<W(x; t)>_x>_r", "<<O(x; t)>_x>_r"]):
    
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
    df.columns = [col.strip() for col in df.columns]
    col_labels = [col for col in df.columns if col in include_col_labels and col != "t"];

    # Create a df with col_labels as column names. This will store the mean data, fot each column in col_labels for each file in files.
    df_file = pan.DataFrame(columns= ["Rmax", "t"] + ["AVG[" + x + "]_ALL" for x in col_labels] + ["VAR[" + x + "]_ALL" for x in col_labels] 
                            + ["AVG[" + x + "]_SURV" for x in col_labels] + ["VAR[" + x + "]_SURV" for x in col_labels] + ["R_SURV[" + x + "]" for x in col_labels])
    
    pooled_df = pan.DataFrame()

    if not all(df["t"] == 0) or dt == None or tmax == None:
        pooled_df["t"] = df["t"] # Store the time values in pooled_df.
    else:
        # Use a custom time range if the time values are all zero.
        pooled_df["t"] = logarithm10_time_bins(tmax, dt, len(df))


    maxR = len(files)
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
        df = df[[col for col in df.columns if col in include_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        # Also check if the time values in the files are the same. If not, return an error.
        if not all(df["t"] == pooled_df["t"]) and not all(df["t"] == 0):
            print("Error: Time values in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        Rstr = re.findall(r'R_[\d]+', file)[0]
        
        for col in col_labels:

            # Store column values in pooled_df
            pooled_df[col + "_" + Rstr] = df[col]

    # Done with all files in files.
    # Now generate element-wise mean and standard deviation for each column in col_labels, across all files.
    # Store the mean and standard deviation in df_file.

    for col in col_labels:

        # Find mean and standard deviation of each column in pooled_df.
        mean_col = pooled_df.filter(like= col).mean(axis=1)
        variance_col = pooled_df.filter(like= col).var(axis=1)
        # Store the mean and standard deviation in df_file.
        df_file["AVG[" + col + "]_ALL"] = mean_col
        df_file["VAR[" + col + "]_ALL"] = variance_col
        
        # Surviving data is calculated as the element-wise mean and standard deviation of the column across all files for each non-zero element.
        # Store the mean and standard deviation in df_file.
        mean_col_surv = pooled_df.filter(like= col).replace(0, np.nan).mean(axis=1)
        variance_col_surv = pooled_df.filter(like= col).replace(0, np.nan).var(axis=1)
        df_file["AVG[" + col + "]_SURV"] = mean_col_surv
        df_file["VAR[" + col + "]_SURV"] = variance_col_surv

        # Store the number of non-zero elements in the column across all files in df_file.
        R_surv_col = pooled_df.filter(like= col).replace(0, np.nan).count(axis=1)
        df_file["R_SURV[" + col + "]"] = R_surv_col
        # Store the Rmax in df_file, as a column with the same value for all rows.
        df_file["Rmax"] = [maxR for i in range(len(df_file))]
        df_file["t"] = pooled_df["t"]
    
    return df_file
        





