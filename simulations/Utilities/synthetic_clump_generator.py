import numpy as np
import pandas as pan
import scipy.ndimage as ndimage
import regex as re
import matplotlib.pyplot as plt
import os
import glob
import random

SPB =3
prefix ="HEXBLADE";
L = 128; r =5;
output_dir = f"../Input/DP/{SPB}Sp/{prefix}/L_{L}_a_0/"

# Create the output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)


hex_length = {"  P(x; t)": 10}
gauss_2sd_radius = {"  P(x; t)": 4}
gauss_amp = {"  P(x; t)": 10}
min_val = {"  P(x; t)": 0.0}
percent_missing = {"  P(x; t)": 0.1}
retain_cluster = {"  P(x; t)": 0}

def create_gaussian_nonperiodic(L, center, amplitude, sigma, min_val):
    """Create a 2D Gaussian centered at `center` in a grid of size `L`."""
    y, x = np.indices((L, L))
    y0, x0 = center
    # Calculate the distance with periodic boundary conditions
    gaussian = (amplitude) * np.exp( -(np.power(x - x0, 2)+ np.power(y - y0,2)) / (2*np.power(sigma,2)))
    #gaussian = (amplitude)/(2*np.power(sigma,2)*np.pi) * np.exp( -(np.power(x - x0, 2)+ np.power(y - y0,2)) / (2*np.power(sigma,2)))
    gaussian[gaussian < min_val] = min_val
    return gaussian

def create_gaussian_periodic(L, center, amplitude, sigma, min_val):
    """Create a 2D Gaussian centered at `center` in a grid of size `L`."""
    y, x = np.indices((L, L))
    y0, x0 = center
    # Calculate the distance with periodic boundary conditions
    dx = np.minimum(np.abs(x - x0), L - np.abs(x - x0))
    dy = np.minimum(np.abs(y - y0), L - np.abs(y - y0))
    gaussian = (amplitude/2*np.power(sigma,2)*np.pi) * np.exp( -(np.power(dx, 2)+ np.power(dy,2)) / (2*np.power(sigma,2)))
    gaussian[gaussian < min_val] = min_val
    return gaussian

def tile_hex_grid(L, hex_length):
    """Create a hexagonal grid pattern."""
    grid = np.zeros((L, L))
    hex_height = np.sqrt(3) * hex_length
    hex_width = 2 * hex_length
    
    for i in range(int(hex_height), L, int(hex_height)):
        for j in range(int(hex_width / 2), L, int(hex_width)):
            if (i // int(hex_height)) % 2 == 1:
                center = (i, j)
            else:
                center = (i, j + int(hex_width / 2))
            grid[i % L, j % L] = 1  # Mark hex center or set a pattern
    
    return grid

def generate_grid(L, hex_length, gauss_2sd_radius, gauss_amp, retain_cluster=0, percent_missing=0, min_val=0):
    if gauss_2sd_radius >= hex_length:
        raise ValueError("gauss_2sd_radius must be less than hex_length")
    
    # Get the Gaussian SD.
    gauss_sd = gauss_2sd_radius / 2.0
    
    grid = np.zeros((L, L))
    hex_grid = tile_hex_grid(L, hex_length)

    
    # Identify hexagon centers (non-zero values in hex_grid)
    centers = np.argwhere(hex_grid > 0)
    num_hex = len(centers)
    # If retain_cluster is set, randomly select a cluster of hexagons to retain, equal to retain_cluster
    if retain_cluster > 0 and retain_cluster <= num_hex:
        selected_centers = random.sample(list(centers), retain_cluster)
    else:
        # Randomly select (1 - percent_missing) fraction of hexagon centers
        selected_centers = random.sample(list(centers), int((1 - percent_missing) * num_hex))
    
    for center in selected_centers:
        gaussian = create_gaussian_nonperiodic(L, center, gauss_amp, gauss_sd, min_val)
        grid += gaussian
    
    grid[grid < min_val] = min_val
    return grid

def plot_grid(grid, title=""):
    plt.imshow(grid, cmap='hot', interpolation='nearest')
    plt.title(title)
    plt.show()
    plt.close()

def save_grid(df, replicate, savedir, frame_subdir ="SEP_{sep}_WID_{wid}/MSF_{missing}/MAMP_{amp}_MINVAL_{minval}/", overwrite=False):
    """Save the grid to CSV, ensuring no file overwrite."""
    # Creating the directory structure for the frame.
    first_key = list(df.keys())[2]
    save_frame_dir = os.path.join(savedir, frame_subdir.format(sep = hex_length[first_key], missing = percent_missing[first_key],
            wid = gauss_2sd_radius[first_key], amp = gauss_amp[first_key], minval = min_val[first_key]))
    os.makedirs(save_frame_dir, exist_ok=True)
    
    existing_files = glob.glob(os.path.join(save_frame_dir, 'FRAME_T_0_a_0_R_*.csv'))
    existing_replicates = [int(re.findall(r'[\d]+' , f)[-1]) for f in existing_files]
    
    if replicate in existing_replicates and not overwrite:
        replicate = max(existing_replicates) + 1
    
    filename = os.path.join(save_frame_dir, f"FRAME_T_0_a_0_R_{replicate}.csv")
    df.to_csv(filename, index=False)

def main(hex_length, gauss_2sd_radius, L, r, gauss_amp, min_val, output_dir, percent_missing):
    for replicate in range(0, r):
        # Create a pandas dataframe with columns: "a_c", "x" followed by keys in hex_length
        # Each column should be a 1D array of length L x L
        df = pan.DataFrame(columns=["a_c", "x"] + list(hex_length.keys()))
        #Store a 0 array of length L x L in the "a_c" column
        df["a_c"] = np.zeros(L*L); df["x"] = np.arange(L*L); 

        for key in hex_length:
            grid = generate_grid(L, hex_length[key], gauss_2sd_radius[key], gauss_amp[key], percent_missing = percent_missing[key], min_val = min_val[key])
            #plot_grid(grid, key)
            #Flatten the grid and store it in the respective column
            df[key] = grid.flatten()
            print("Max value: ", np.max(grid))
            #grid = df[key].values.reshape(L, L)
            plot_grid(grid, key)
        save_grid(df, replicate, output_dir, overwrite=False)

if __name__ == "__main__":
    # Example usage
    main(hex_length, gauss_2sd_radius, L, r, gauss_amp, min_val, output_dir, percent_missing)
