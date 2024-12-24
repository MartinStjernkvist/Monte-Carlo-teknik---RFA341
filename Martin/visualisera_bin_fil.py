from imports import *


"""
I matlab:

% Save the entire array_3d to a .mat file
save('phantom_data.mat', 'array_3d');
"""

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

# Ensure array_3d is in the shape (x, y, z)
x, y, z = array_3d.shape


#   ----------------------------------------------------------------------
#   INPUT
#   ----------------------------------------------------------------------

# Choose a slice index along the z-axis (middle slice in this case)
slice_idx = z // 2

array_3d = array_3d


#   ----------------------------------------------------------------------
#   PLOTTING 2D
#   ----------------------------------------------------------------------

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

# plt.show()



#   ----------------------------------------------------------------------
#   INPUT 3D
#   ----------------------------------------------------------------------

# from starting_position_photon import sliced_array_njure
# array_3d = sliced_array_njure

# array_3d = array_3d[50:-50,50:200,600:1100]
# array_3d = array_3d[:,25:215,:] # denna matris inkluderar kroppen, försöker skära bort luft runtomkring
array_3d = array_3d[:, 25:215, 500:1150] # försöker skära en rimlig matris där växelverkan kan ske och fortfarande leda till energideponering i benmärgen i ryggen
x, y, z = array_3d.shape

#   ----------------------------------------------------------------------
#   PLOTTING 3D
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
img1 = ax1.imshow(array_3d[slice_x_index, :, :], cmap='gray')  # x-slice
ax1.set_title(f'X Slice: {slice_x_index}')
img2 = ax2.imshow(array_3d[:, slice_y_index, :], cmap='gray')  # y-slice
ax2.set_title(f'Y Slice: {slice_y_index}')
img3 = ax3.imshow(array_3d[:, :, slice_z_index], cmap='gray')  # z-slice
ax3.set_title(f'Z Slice: {slice_z_index}')

# Create sliders for controlling the slices
ax_slider_x = plt.axes([0.1, 0.01, 0.75, 0.03], facecolor='lightgoldenrodyellow')
slider_x = Slider(ax_slider_x, 'X Slice', 0, x - 1, valinit=slice_x_index, valstep=1)

ax_slider_y = plt.axes([0.1, 0.06, 0.75, 0.03], facecolor='lightgoldenrodyellow')
slider_y = Slider(ax_slider_y, 'Y Slice', 0, y - 1, valinit=slice_y_index, valstep=1)

ax_slider_z = plt.axes([0.1, 0.11, 0.75, 0.03], facecolor='lightgoldenrodyellow')
slider_z = Slider(ax_slider_z, 'Z Slice', 0, z - 1, valinit=slice_z_index, valstep=1)


# Update function for the sliders
def update(val):
    slice_x_index = int(slider_x.val)
    slice_y_index = int(slider_y.val)
    slice_z_index = int(slider_z.val)

    # Update the images for each slice
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

# Show the plot
plt.show()