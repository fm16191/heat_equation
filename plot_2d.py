import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

import shutil

import sys
import os
from tempfile import mkdtemp

# Get data file
if len(sys.argv) < 2:
    print("Export simulation x y T to picture (if one line) or gif (if multiple lines)")
    print(f"Usage: {sys.argv[0]} <file.data>")
    exit(1)

data_file = sys.argv[1]


def export_frame(line:str, output_filename:str, frame:str=None, min_t:float=None, max_t:float=None):
    SIZE_Y_AXIS, SIZE_X_AXIS = map(int, line.split()[:2])

    data = np.array(list(map(float, line.split()[2:])))
    data = data.reshape((SIZE_Y_AXIS, SIZE_X_AXIS))

    fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'heatmap'}]])

    heatmap = go.Heatmap(z=data, colorscale='plasma', zmin=min_t, zmax=max_t)
    fig.add_trace(heatmap, row=1, col=1)

    fig.update_layout(
        xaxis_title="x (dx)",
        yaxis_title="y (dy)",
        width=1200,
        height=1200,
        title=f"Heat Equation 2D - timestep {frame}" if frame else f"Heat Equation 2D",
    )

    pio.write_image(fig, output_filename, format='png')

# Read line count
with open(data_file, 'r') as file:
    print("Reading lines count ...", end="\r")
    num_lines = sum(1 for line in file)
    print("Number of lines :", num_lines)

if num_lines == 0:
    print("Empty file")
    exit(1)

# Export one frame
elif num_lines == 1:
    print("Export heatmap")
    out_filename = "heatmap_2d.png"

    with open(data_file, 'r') as file:
        line = file.readline()
        export_frame(line, out_filename)

    print(f"Heatmap exported as \"{out_filename}\"")

# Export multiple frames (gif)
else:

    output_directory = f"frames_{data_file}"
    if os.path.exists(output_directory):
        print(f"Folder exists already : {output_directory}")
        fdel = input("Overwrite ? [y/N] ")
        if fdel == "y":
            shutil.rmtree(output_directory)
        else:
            exit(1)

    os.makedirs(output_directory)

    print(f"Export heatmap frames to folder \"{output_directory}\"")

    with open(data_file, 'r') as file:
        for frame, line in enumerate(file):
            if frame == 0:
                data = np.array(list(map(float, line.split()[2:])))
                min_t = data.min() # Should be 0
                max_t = data.max() * 1.05

            print(f"Export frame {frame: 4d}/{num_lines}", end='\r')
            frame_filename = os.path.join(output_directory, f'frame_{frame:04d}.png')
            export_frame(line, frame_filename, frame=frame, min_t=min_t, max_t=max_t)

    # shutil.rmtree(output_directory)
    print("Exporting animation frames complete.")
