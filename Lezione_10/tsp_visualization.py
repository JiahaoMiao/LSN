import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import imageio.v2 as imageio
from tqdm import tqdm

def load_data(path, generation):
    filename = os.path.join(path, 'config', f'BestPath{generation}.dat')
    try:
        X, Y = np.loadtxt(filename, usecols=(0, 1), delimiter=' ', unpack=True)
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")
        X, Y = np.array([]), np.array([])
    return X, Y

def plot_generation(X, Y, generation, best_length, output_path, df, BBox, ruh_m):
    fig, ax = plt.subplots(figsize=(8, 8))
    # Add the map image
    ax.imshow(ruh_m, zorder=0, extent=BBox, aspect='equal')

    # Scatter plot for the cities
    ax.scatter(df.Latitude, df.Longitude, zorder=1, alpha=1, c='b', s=10, label='Cities')
    
    # Highlight the starting point
    if len(X) > 0 and len(Y) > 0:
        ax.scatter(X[0], Y[0], marker='*', s=100, color='green', label='Start')
    
    # Plot the path
    if len(X) > 1 and len(Y) > 1:
        ax.plot(X, Y, marker=' ', color='red', label='Best Path')
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(BBox[0], BBox[1])
    ax.set_ylim(BBox[2], BBox[3])
    ax.set_title(f"Generation {generation}, Euclidean Length = {best_length}")
    
    
    ax.grid(True)
    ax.legend()
    
    output_filename = os.path.join(output_path, f"{generation}.png")
    plt.savefig(output_filename, format="png", dpi=150)
    plt.close(fig)
    return output_filename

def create_gif(filenames, output_filename):
    with imageio.get_writer(output_filename, mode='I') as writer:
        for filename in tqdm(filenames, desc="Creating GIF"):
            try:
                image = imageio.imread(filename)
                writer.append_data(image)
            except Exception as e:
                print(f"Error reading {filename}: {e}")

def main():
    path = "DATA/Migration/0/"
    best_len_file = os.path.join(path, "BestLength.dat")
    output_path = "imgs"
    gif_output = os.path.join(output_path, "Italy.gif")
    
    os.makedirs(output_path, exist_ok=True)  # Create output directory if it doesn't exist
    
    try:
        gen, best = np.loadtxt(best_len_file, usecols=(0, 1), delimiter=' ', unpack=True)
    except Exception as e:
        print(f"Error loading best length data: {e}")
        return
    
    best_gen = 350  # Total number of generations to plot
    filenames = []

    print(f"Generating {best_gen} images...")

    for j in tqdm(range(best_gen), desc="Generating images"):
        X, Y = load_data(path, 10*(j+1))  # Load data every 10 generations
        if len(X) == 0 or len(Y) == 0:
            continue
        filename = plot_generation(X, Y, 10*(j+1), best[10*(j+1)-1], output_path, df, BBox, ruh_m)  # Plot the generation and save the image
        filenames.append(filename)

    print("\nGenerating GIF...")
    create_gif(filenames, gif_output)
    
    # Cleanup generated image files
    for filename in filenames:
        os.remove(filename)
    
    print("\nEND")

if __name__ == "__main__":
    # Global variables
    df = pd.read_csv('INPUT/cap_prov_ita.dat', sep=' ', header=None, names=['Latitude', 'Longitude'])
    # Define the bounding box
    BBox = ((df.Latitude.min(),   df.Latitude.max(),      
         df.Longitude.min(), df.Longitude.max()))
    # Load the map image once to avoid loading it multiple times
    ruh_m = plt.imread('imgs/italy.png')
    main()