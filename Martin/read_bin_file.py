import scipy.io
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


"""
I matlab:

% Save the entire array_3d to a .mat file
save('phantom_data.mat', 'array_3d');
"""

# Load the .mat file
mat_data = scipy.io.loadmat('phantom_data.mat')

# Extract the 3D array from the loaded dictionary
array_3d = mat_data['array_3d']

# Print the shape of the array
print(array_3d.shape)

# Define the shape of the array (replace with the correct dimensions)
x, y, z = 256, 256, 1200


# Choose a slice index along the z-axis (middle slice in this case)
slice_idx = z // 2

# Display a slice along the z-axis
plt.imshow(array_3d[:, :, slice_idx], cmap='gray')
plt.title(f'Slice {slice_idx}')
plt.colorbar()
plt.show()




# Create a figure and axis for plotting the 2D slice
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(bottom=0.25)
initial_slice_index = z // 2  # Start at the middle slice

# Display the initial slice
img = ax.imshow(array_3d[:, :, initial_slice_index], cmap='gray')
ax.set_title(f'Slice {initial_slice_index}')
plt.colorbar(img, ax=ax)

# Add a slider for navigating through slices
ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Slice Index', 0, z-1, valinit=initial_slice_index, valstep=1)

# Update function for the slider
def update(val):
    slice_index = int(slider.val)
    img.set_data(array_3d[:, :, slice_index])
    ax.set_title(f'Slice {slice_index}')
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()


print(array_3d[100:150,100:150,600])