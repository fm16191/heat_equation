import numpy as np
import matplotlib.pyplot as plt

import sys
if len(sys.argv) < 2:
    print(f"Usage : {sys.argv[0]} <file.data>")
    exit(1)

data = np.loadtxt(sys.argv[1])

x, y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

plt.pcolormesh(x, y, data, cmap='plasma', shading='auto')
plt.colorbar(label='Temperature', orientation='vertical')
plt.title('Heat Equation 1D')
plt.xlabel('x (dx)')
plt.ylabel('t (dt)')

out_filename = 'heatmap_1d.png'
print(f"Heatmap exported as \"{out_filename}\"")
plt.savefig(out_filename)

plt.show()
