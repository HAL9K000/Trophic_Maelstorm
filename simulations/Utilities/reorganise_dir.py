import numpy as np
import os
import pandas as pan
import regex as re
from pathlib import Path
import glob
import shutil
import sys
import copy

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

prefixes =["DiC-NEW-REF-IMP"]#, "DiC-NEW"]
#prefixes =["", "DiC", "BURNIN", "DiC-BURNIN", "DDM-DiC", "DDM-DiC-BURNIN"]
root_dir = "../Data/Remote/Rietkerk/Frames/Stochastic/3Sp/Test/"
#out_dir_noprefix = "../Data/Rietkerk/Reorganised_Frames/Stoc/3Sp/StdParam_20_100_Test/"

out_dir_noprefix = "../Data_Processing/StdParam_20_100_Test/"

dP = 11000
Geq = 4.802 # Optional. If Geq is not used in the subdirectory name, set Geq = "NA".
L= [256]
indx_vals_t = 26
#Extract n largest values of T if indx_vals_t = -n, 
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

choice = input("Enter choice: ")
if choice == "y":
    out_dir_noprefix = input("Enter new out_dir_noprefix: ")
    dP = int(input("Enter new dP: "))
    Geq = input("Enter new Geq ('NA' where it isn't applicable): ")
    L = [int(x) for x in input("Enter new L: ").split()]
    indx_vals_t = int(input("Enter new indx_vals_t: "))


def find_timerange(files, a, indx= indx_vals_t ):
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
    if np.abs(indx) > len(T_vals):
        print("Warning: indx greater than number of T values. Returning all T values.")
        return T_vals
    return T_vals[indx:] if indx < 0 else T_vals[:indx]

def find_vals(files, dir, reg_search_patterns, vals):
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
        find_vals(files_copy, dir, reg_search_patterns[1:], vals)



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
    find_vals(files, parendir, reg_search_type_patterns, vals)
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
    
    
    '''
    #vals = [re.findall(r'a_[\d]*[.][\d]+' , f)[0]
    #            if re.findall(r'a_[\d]*[.][\d]+' , f) 
    #            else re.findall(r'a_[\d]+' , f)[0] for f in files]
    
    # Remove a_ from a_vals and sort them in ascending order.
    a_vals = sorted([float(re.findall(r'[\d]*[.][\d]+', a)[0]) if re.findall(r'[\d]*[.][\d]+', a) 
                     else int(re.findall(r'[\d]+', a)[0]) for a in a_vals])
    a_vals = np.unique(a_vals)
    print("Found a_vals: " + str(a_vals))
    # Save a_vals to a text file in dir.
    try:

        np.savetxt(dir + "/a_vals.txt", a_vals, delimiter="\t", fmt="%g")
        # Test reading the text file line by line.
        with open(dir + "/a_vals.txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                print(line.strip())
            

        # Close the file.
        f.close()
    except Exception as e:
        print("Error: Could not write a_vals to " + dir + "/a_vals.txt with error message: \n" + str(e))
        return
    '''
        
    
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

                    t_range = find_timerange(files, a, indx_vals_t)

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
                            # Now create new filename.
                            if re.search("GAMMA", source_filename):
                                out_filename = "GAMMA_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
                            else:
                                out_filename = "FRAME_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
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
                                max_R = max([int(re.findall(r'[\d]+' , f)[-1]) for f in glob.glob(t_subdir + "*.csv")])
                                out_filename = "FRAME_T_" + str(T) + "_a_" + str(a_val) + "_R_" + str(max_R + 1) + ".csv"
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



    
    


    

main()
#dir = "\\\\?\\D:\\cygwin64\\home\\koust\\Code\\Trophic_Maelstorm\\simulations\\Data\\Remote\\Rietkerk\\Frames\\Stochastic\\3Sp\\DDM_DiC_BURNIN_0.025-0.065_dP_30000_Geq_4.802";
#rename(dir, "_RANDDDMDiCBURNIN_ThreeSp_P_c_DP", "")

#dir = "../Data/Remote/Test_Rietkerk_Frames/Stochastic/3Sp/Reorganised/StandardParam_20_100/DiC_BURNIN_/"
#test_txt_RW(dir)