import plotly.graph_objects as go
import plotly.io as pio
import numpy as np

import sys

# Get data file
if len(sys.argv) < 2:
    print("Export simulation x y T to picture (if one line) or gif (if multiple lines)")
    print(f"Usage: {sys.argv[0]} <file.data>")
    exit(1)

data_file = sys.argv[1]

# Read file
fo = open(data_file, "r")
line = fo.readline()
fo.close()

SIZE_X_AXIS, SIZE_Y_AXIS, SIZE_Z_AXIS = map(int, line.split()[:3])

data = np.array(list(map(float, line.split()[3:])))
temperature_values = data.reshape((SIZE_X_AXIS, SIZE_Y_AXIS, SIZE_Z_AXIS))

# Meshgrid for plotting
x, y, z = np.meshgrid(range(SIZE_X_AXIS), range(SIZE_Y_AXIS), range(SIZE_Z_AXIS))

# Filter out points the lowest temperature points
min_temp = np.min(temperature_values)
max_temp = np.max(temperature_values)
temp_diff = max_temp - min_temp

indices = temperature_values != min_temp
x = x[indices]
y = y[indices]
z = z[indices]
temperature_values = temperature_values[indices]

# Create a 3D subplot with different opacity for different temperature ranges
def subplot(fig, t_range, opacity=None, showlegend=False):
    # Reduce dataset with temperature range
    subset_indices = (temperature_values >= t_range[0] * max_temp) & (
        temperature_values <= t_range[1] * max_temp
    )
    subset = temperature_values[subset_indices]

    if not opacity:
        opacity = np.mean(subset) ** 2 / max_temp**2 * 0.8
        print(opacity)

    fig.add_trace(
        go.Scatter3d(
            x=x[subset_indices].flatten(),
            y=y[subset_indices].flatten(),
            z=z[subset_indices].flatten(),
            mode="markers",
            marker=dict(
                size=5,
                color=subset.flatten(),
                colorscale="plasma",
                colorbar=dict(title="Temperature") if showlegend else None,
            ),
            opacity=opacity,
            showlegend=showlegend,
            hovertext=[f"Temperature: {temp:.2f}" for temp in subset.flatten()],
        )
    )


fig = go.Figure()
subplot(fig, (0.00, 1.0), opacity=0.07, showlegend=True)

t_ranges = [(0.3, 0.9), (0.9, 1.0)]
for t_range in t_ranges:
    subplot(fig, t_range)

frame = 0
fig.update_layout(
    showlegend=False,
    # xaxis_title="x",
    # yaxis_title="y",
    width=1000,
    height=1000,
    title=f"Heat Equation 3D - timestep {frame}" if frame else f"Heat Equation 3D",
)

output_filename = f"heatmap_3d.png"
pio.write_image(fig, output_filename, format="png")
print(f"Heatmap 3D saved at {output_filename}")

fig.show()
