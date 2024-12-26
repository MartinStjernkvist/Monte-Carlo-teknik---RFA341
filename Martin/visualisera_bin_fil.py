from imports import *


"""
I matlab:

% Save the entire array_3d to a .mat file
save('phantom_data.mat', 'array_3d');
"""

#   ----------------------------------------------------------------------
#   READ IN DATA
#   ----------------------------------------------------------------------

mat_file = 'phantom_data.mat'
data = scipy.io.loadmat(mat_file)

if __name__ == "__main__":
    print("Keys in the .mat file:", data.keys())

array_phantom = data['array_3d']

"""
#   ----------------------------------------------------------------------
#   INPUT 2D
#   ----------------------------------------------------------------------

# Choose a slice index along the z-axis (middle slice in this case)
slice_idx = z // 2

array_phantom = array_phantom

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
"""



"""
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
"""


#   ----------------------------------------------------------------------
#   PLOTTING 3D
#   ----------------------------------------------------------------------

class visualisera_matris:
    def __init__(self, array_3d):
        self.array_3d = array_3d
        self.x, self.y, self.z = array_3d.shape

        # Initialize slice indices for each axis
        self.slice_x_index = self.x // 2
        self.slice_y_index = self.y // 2
        self.slice_z_index = self.z // 2

        # Create the figure and layout
        self.fig = plt.figure(figsize=(15, 7.5))
        self.gs = GridSpec(1, 3, width_ratios=[1, 1, 1], height_ratios=[1])

        # Create subplots
        self.ax1 = plt.subplot(self.gs[0])
        self.ax2 = plt.subplot(self.gs[1])
        self.ax3 = plt.subplot(self.gs[2])

        # Display the initial slices
        # self.img1 = self.ax1.imshow(self.array_3d[self.slice_x_index, :, :], cmap='gray')
        # self.ax1.set_title(f'X Slice: {self.slice_x_index}')
        # self.img2 = self.ax2.imshow(self.array_3d[:, self.slice_y_index, :], cmap='gray')
        # self.ax2.set_title(f'Y Slice: {self.slice_y_index}')
        # self.img3 = self.ax3.imshow(self.array_3d[:, :, self.slice_z_index], cmap='gray')
        # self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        self.vmin, self.vmax = 0, np.max(array_3d)
        self.img1 = self.ax1.imshow(array_3d[self.slice_x_index, :, :], cmap='gray', vmin=self.vmin, vmax=self.vmax)  # x-slice
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')
        self.img2 = self.ax2.imshow(array_3d[:, self.slice_y_index, :], cmap='gray', vmin=self.vmin, vmax=self.vmax)  # y-slice
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')
        self.img3 = self.ax3.imshow(array_3d[:, :, self.slice_z_index], cmap='gray', vmin=self.vmin, vmax=self.vmax)  # z-slice
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        # Create sliders
        self.create_sliders()

        # Attach update function to sliders
        self.slider_x.on_changed(self.update)
        self.slider_y.on_changed(self.update)
        self.slider_z.on_changed(self.update)

    def create_sliders(self):
        ax_slider_x = plt.axes([0.1, 0.01, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_x = Slider(ax_slider_x, 'X Slice', 0, self.x - 1, valinit=self.slice_x_index, valstep=1)

        ax_slider_y = plt.axes([0.1, 0.06, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_y = Slider(ax_slider_y, 'Y Slice', 0, self.y - 1, valinit=self.slice_y_index, valstep=1)

        ax_slider_z = plt.axes([0.1, 0.11, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_z = Slider(ax_slider_z, 'Z Slice', 0, self.z - 1, valinit=self.slice_z_index, valstep=1)

    def update(self, val):
        # Update slice indices
        self.slice_x_index = int(self.slider_x.val)
        self.slice_y_index = int(self.slider_y.val)
        self.slice_z_index = int(self.slider_z.val)

        # Update the images
        self.img1.set_data(self.array_3d[self.slice_x_index, :, :])
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')

        self.img2.set_data(self.array_3d[:, self.slice_y_index, :])
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')

        self.img3.set_data(self.array_3d[:, :, self.slice_z_index])
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        # Redraw the canvas
        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()

#   ----------------------------------------------------------------------
#   INPUT 3D
#   ----------------------------------------------------------------------

if __name__ == "__main__":

    # from starting_position_photon import sliced_array_njure
    # array_3d = sliced_array_njure

    # array_3d = array_3d[50:-50,50:200,600:1100]
    # array_3d = array_3d[:,25:215,:] # denna matris inkluderar kroppen, försöker skära bort luft runtomkring

    # försöker skära en rimlig matris där växelverkan kan ske och fortfarande leda till energideponering i benmärgen i ryggen
    array_3d = array_phantom[:,25:215, 500:1150]

    viewer = visualisera_matris(array_3d)
    viewer.show()