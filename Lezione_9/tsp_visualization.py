import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import imageio.v2 as imageio
from tqdm import tqdm

def load_data(path, problem, generation):
    filename = os.path.join(path, problem, 'config', f'BestPath{generation}.dat')
    try:
        X, Y = np.loadtxt(filename, usecols=(0, 1), delimiter=' ', unpack=True)
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")
        X, Y = np.array([]), np.array([])
    return X, Y

def plot_generation(X, Y, generation, best_length, output_path, problem):
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Scatter plot for the cities
    ax.scatter(X, Y, marker='o',s=34, color='blue', label='Cities')
    
    # Highlight the starting point
    if len(X) > 0 and len(Y) > 0:
        ax.scatter(X[0], Y[0], marker='*',s = 100, color='green', label='Start')
    
    # Plot the path
    if len(X) > 1 and len(Y) > 1:
        ax.plot(X, Y, marker=' ', color='red', alpha=1, linewidth=2)
    
    # Add the geometric shape based on the problem
    if problem == "Circumference":
        circle = plt.Circle((0, 0), 1, linewidth=2, edgecolor='black',fill = False)
        ax.add_artist(circle)
        ax.set_aspect(1)
    elif problem == "Square":
        square = patches.Rectangle((-1, -1), 2, 2, linewidth=2, edgecolor='black',fill = False)
        ax.add_patch(square)
        ax.set_aspect(1)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_ylim([-1.5, 1.5])
    ax.set_xlim([-1.5, 1.5])
    ax.set_title(f"Generation {generation}, Length = {best_length}")
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
    path = "DATA/"
    # problem = "Square"
    problem = "Circumference"
    best_len_file = os.path.join(path, problem, "BestLength.dat")
    output_path = os.path.join("imgs", problem)
    gif_output = os.path.join(output_path, problem + ".gif")
    
    os.makedirs(output_path, exist_ok=True)
    
    try:
        gen, best = np.loadtxt(best_len_file, usecols=(0, 1), delimiter=' ', unpack=True)
    except Exception as e:
        print(f"Error loading best length data: {e}")
        return
    
    best_gen = 200  # Total number of generations to plot
    filenames = []

    print(f"Generating {best_gen} images...")
    
    for j in tqdm(range(best_gen), desc="Generating images"):
        X, Y = load_data(path, problem, j)
        if len(X) == 0 or len(Y) == 0:
            continue
        filename = plot_generation(X, Y, j, best[j], output_path, problem)
        filenames.append(filename)

    print("\nGenerating GIF...")
    create_gif(filenames, gif_output)
    
    # Cleanup generated image files
    for filename in filenames:
        os.remove(filename)
    
    print("\nEND")

if __name__ == "__main__":
    main()
