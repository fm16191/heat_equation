#!/bin/env python3

import sys
import numpy as np
import os
import shutil

# Get data file
if len(sys.argv) < 2:
    print("Export 3d utput to CSV files that Paraview can read")
    print(f"Usage: {sys.argv[0]} <file.data>")
    exit(1)

data_file = sys.argv[1]

# Create folder
output_directory = f"csv_{data_file}"
if os.path.exists(output_directory):
    print(f"Folder exists already : {output_directory}")
    if input("Overwrite ? [y/N] ") == "y": shutil.rmtree(output_directory)
    else: exit(1)

os.makedirs(output_directory)

# Get simulation dimensions
fo = open(data_file, "r")
line = fo.readline()
SIZE_X_AXIS, SIZE_Y_AXIS, SIZE_Z_AXIS = map(int, line.split()[:3])
lines = sum(1 for line in fo) + 1 # first line is already read
print(f"Simulation dimensions are {SIZE_X_AXIS} x {SIZE_Y_AXIS} x {SIZE_Z_AXIS}")
print(f"Simulation has {lines} steps")
fo.close()
fo = open(data_file, "r")


# Create CSV files
for t, line in enumerate(fo.readlines()):
    data = np.array(list(map(float, line.split()[3:])))

    output_filename = f"{output_directory}/{data_file}_{t:03d}.csv"
    fw = open(output_filename, "w")

    fw.write("x, y, z, t\n")

    print(SIZE_X_AXIS, SIZE_Y_AXIS, SIZE_Z_AXIS, len(data), SIZE_X_AXIS * SIZE_Y_AXIS * SIZE_Z_AXIS)

    for i in range(len(data)):
        fw.write(f"{int(i / SIZE_Y_AXIS / SIZE_Z_AXIS)}, {int(i / SIZE_Z_AXIS) % SIZE_Y_AXIS}, {i % SIZE_Z_AXIS}, {data[i]}\n")

    fw.close()
fo.close()

print(f"CSV files written to \"{output_directory}/\"")
