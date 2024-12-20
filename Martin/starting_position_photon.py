from imports_file import *

#   ----------------------------------------------------------------------
#   READ IN DATA
#   ----------------------------------------------------------------------

# Load the .mat file
mat_file = 'phantom_data.mat'  # Replace with your .mat file path
data = scipy.io.loadmat(mat_file)


# Check what keys are in the .mat file
print("Keys in the .mat file:", data.keys())

# Access the 3D array from the .mat file (replace 'your_matrix_key' with the actual key from the file)
array_3d = data['array_3d']  # Replace 'phantom_matrix' with the actual key in your .mat file

#   ----------------------------------------------------------------------
#   CREATE SLICED ARRAYS
#   ----------------------------------------------------------------------

sliced_array_phantom = array_3d[30:230, 40:210, 570:900]
# sliced_array_phantom = array_3d[:,:,:]
x, y, z = sliced_array_phantom.shape

sliced_array_njure = np.zeros((x, y, z))

target_values = [17, 18, 19, 20, 21, 22, 23]
mask = np.isin(sliced_array_phantom, target_values)
sliced_array_njure[mask] = sliced_array_phantom[mask]

#   ----------------------------------------------------------------------
#   INPUT ARRAY FOR PLOTTING
#   ----------------------------------------------------------------------

# replace old array with new array, which will then be plotted
array_3d = sliced_array_njure

#   ----------------------------------------------------------------------
#   PLOTTING
#   ----------------------------------------------------------------------

# Create a larger figure with a GridSpec layout
fig = plt.figure(figsize=(18, 10))  # Increase figure size
gs = GridSpec(1, 3, width_ratios=[1, 1, 1], height_ratios=[1])  # Adjust ratios to make plots fill space

# Create subplots with custom sizing
ax1 = plt.subplot(gs[0])  # First larger plot
ax2 = plt.subplot(gs[1])  # Second larger plot
ax3 = plt.subplot(gs[2])  # Smaller plot

# Initialize slice indices for each axis
slice_x_index = x // 2  # Start in the middle of x
slice_y_index = y // 2  # Start in the middle of y
slice_z_index = z // 2  # Start in the middle of z

# Display the initial slices
vmin, vmax = 0, np.max(array_3d)
img1 = ax1.imshow(array_3d[slice_x_index, :, :], cmap='gray', vmin=vmin, vmax=vmax)  # x-slice
ax1.set_title(f'X Slice: {slice_x_index}')
img2 = ax2.imshow(array_3d[:, slice_y_index, :], cmap='gray', vmin=vmin, vmax=vmax)  # y-slice
ax2.set_title(f'Y Slice: {slice_y_index}')
img3 = ax3.imshow(array_3d[:, :, slice_z_index], cmap='gray', vmin=vmin, vmax=vmax)  # z-slice
ax3.set_title(f'Z Slice: {slice_z_index}')

# Define consistent slider dimensions
slider_width = 0.75  # Consistent width for all sliders
slider_height = 0.03  # Consistent height for all sliders
slider_spacing = 0.05  # Vertical space between sliders

# Calculate vertical positions for the sliders
slider_x_ypos = 0.01
slider_y_ypos = slider_x_ypos + slider_spacing
slider_z_ypos = slider_y_ypos + slider_spacing

# Create slider axes with consistent dimensions
ax_slider_x = plt.axes([0.1, slider_x_ypos, slider_width, slider_height], facecolor='lightgoldenrodyellow')
ax_slider_y = plt.axes([0.1, slider_y_ypos, slider_width, slider_height], facecolor='lightgoldenrodyellow')
ax_slider_z = plt.axes([0.1, slider_z_ypos, slider_width, slider_height], facecolor='lightgoldenrodyellow')

# Create the sliders
slider_x = Slider(ax_slider_x, 'X Slice', 0, x - 1, valinit=slice_x_index, valstep=1)
slider_y = Slider(ax_slider_y, 'Y Slice', 0, y - 1, valinit=slice_y_index, valstep=1)
slider_z = Slider(ax_slider_z, 'Z Slice', 0, z - 1, valinit=slice_z_index, valstep=1)

# Explicitly reset Z-slider bounds again after creation to ensure correctness
ax_slider_z.set_position([0.1, slider_z_ypos, slider_width, slider_height])

def update(val):
    slice_x_index = int(slider_x.val)
    slice_y_index = int(slider_y.val)
    slice_z_index = int(slider_z.val)

    img1.set_data(array_3d[slice_x_index, :, :])
    ax1.set_title(f'X Slice: {slice_x_index}')

    img2.set_data(array_3d[:, slice_y_index, :])
    ax2.set_title(f'Y Slice: {slice_y_index}')

    img3.set_data(array_3d[:, :, slice_z_index])
    ax3.set_title(f'Z Slice: {slice_z_index}')

    fig.canvas.draw_idle()


# Attach the update function to the sliders
slider_x.on_changed(update)
slider_y.on_changed(update)
slider_z.on_changed(update)

# plt.show()

#   ----------------------------------------------------------------------
#   PRINTING
#   ----------------------------------------------------------------------

# print(array_3d[slice_x_index, :, :])  # X slice
# print(array_3d[:, slice_y_index, :])  # Y slice
# print(array_3d[:, :, slice_z_index])  # Z slice