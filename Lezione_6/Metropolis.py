import os
import subprocess
from tqdm import tqdm

def modify_input_file(input_file, new_temp, h=0.0):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines): # Loop over the lines, i->line index, line->line content
        if line.startswith('TEMP'):
            lines[i] = f'TEMP {new_temp}\n' # Modify the line content then add a newline
        elif line.startswith('SIMULATION_TYPE'):
            values = line.split() # Split the line into a list of strings
            values[-1] = str(h) # Modify the last element
            lines[i] = ' '.join(values) + '\n' # Join the list into a string and add a newline

    with open(input_file, 'w') as f:
        f.writelines(lines) # Write the modified lines back to the file

def extract_last_line(output_file, destination_file, new_temp): # Extract the last line of the output file
    with open(output_file, 'r') as f: # Open the output file in read mode
        lines = f.readlines() 
        if lines: # Check if the file is not empty
            last_line = lines[-1].split() # Split the last line into a list of strings
            last_line[0] = f'{new_temp}'  # Modify the first column 
            last_line.pop(1)  # Remove the second column
            last_line_str = ' '.join(last_line) # Join the list into a string
             # Check if the destination file exists
    if os.path.exists(destination_file):
        # If the file exists, open it in append mode ('a')
        with open(destination_file, 'a') as dest:
            dest.write(last_line_str+ '\n')
    else:
        # If the file doesn't exist, open it in write mode ('w')
        with open(destination_file, 'w') as dest:
            dest.write(last_line_str+ '\n')

def main():

    # Define the variables
    input_file = '../INPUT/input.ising'  # input file name
    cxx_program = './simulator.exe'  # program to execute, argument is not needed here

    initial_temp = 0.5
    num_iterations = 25  # Change this to the number of iterations you want

    output_dir = '../OUTPUT'
    data_dir = '../DATA/Ising/Metropolis'
    prop = ["total_energy", "specific_heat", "susceptibility"]
    
    #--------------------------------------------------------------
    
    
    # Run the simulation for properties with h=0.0
    print("Running simulation for properties with h=0.0")

    with tqdm(total=num_iterations) as pbar: # Create a progress bar
        for i in range(num_iterations+1):
            new_temp = "{:.1f}".format(initial_temp + i * 0.1)  # Change this to the modification you need
            
            # Modify the input file with the new temperature
            modify_input_file(input_file, new_temp, h=0.0)
            
            #if the program requires an argument: subprocess.run([cxx_program, 'arg1', 'arg2',...])
            subprocess.run(cxx_program)  # Run the simulation 
            # Write results of the simulation to the corresponding output file
            for name in prop:
                output_file = os.path.join(output_dir, f'{name}.dat')
                destination_file = os.path.join(data_dir, f'{name}.dat')
                extract_last_line(output_file, destination_file, new_temp)
            
            pbar.update(1)

    # Run the simulation for magnetization with h=0.02
    h_val = 0.02

    print(f"Running simulation for magnetization with h={h}")
    with tqdm(total=num_iterations) as pbar:
        for i in range(num_iterations+1):
            new_temp = "{:.1f}".format(initial_temp + i * 0.1)  # Change this to the modification you need
            modify_input_file(input_file, new_temp, h=h) # Modify the input file with the new temperature and h value
            
            subprocess.run(cxx_program) # Run the simulation

            output_file = os.path.join(output_dir, f'magnetization.dat')
            destination_file = os.path.join(data_dir, f'magnetization.dat')
            extract_last_line(output_file, destination_file, new_temp)

            pbar.update(1)

if __name__ == "__main__":
    main()
