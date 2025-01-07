from imports import *


class visualisera_matris:
    """
    Kod som visualiserar fantomen i tre olika tvärsnitt, med sliders.
    """
    def __init__(self, array_3d, visa_något=False):
        self.array_3d = array_3d
        self.x, self.y, self.z = array_3d.shape

        # Initialisera mitt-tvärsnittet för respektive dimension.
        self.slice_x_index = self.x // 2
        self.slice_y_index = self.y // 2
        self.slice_z_index = self.z // 2

        # Skapa figur.
        self.fig = plt.figure(figsize=(15, 7.5))
        self.gs = GridSpec(1, 3, width_ratios=[1, 1, 1], height_ratios=[1])

        # Skapa delfigurer.
        self.ax1 = plt.subplot(self.gs[0])
        self.ax2 = plt.subplot(self.gs[1])
        self.ax3 = plt.subplot(self.gs[2])

        # visa_något bestämmer ifall det som plottas ska ha olika färger beroende på voxelvärde, eller
        # ifall allt som plottas ska vara vitt (så att det syns)
        if visa_något == False:
            self.vmin, self.vmax = 0, np.max(array_3d)
        else:
            self.vmin, self.vmax = 0, 0.0001*np.max(array_3d)

        self.img1 = self.ax1.imshow(array_3d[self.slice_x_index, :, :], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # x-tvärsnitt
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')
        self.img2 = self.ax2.imshow(array_3d[:, self.slice_y_index, :], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # y-tvärsnitt
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')
        self.img3 = self.ax3.imshow(array_3d[:, :, self.slice_z_index], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # z-tvärsnitt
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        # Skapa sliders.
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
        # Updatera sliders.
        self.slice_x_index = int(self.slider_x.val)
        self.slice_y_index = int(self.slider_y.val)
        self.slice_z_index = int(self.slider_z.val)

        # Updatera figurerna.
        self.img1.set_data(self.array_3d[self.slice_x_index, :, :])
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')

        self.img2.set_data(self.array_3d[:, self.slice_y_index, :])
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')

        self.img3.set_data(self.array_3d[:, :, self.slice_z_index])
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()

#   ----------------------------------------------------------------------
#   READ IN DATA
#   ----------------------------------------------------------------------

"""
I matlab:

% Save the entire array_3d to a .mat file
save('phantom_data.mat', 'array_3d');
"""

data = scipy.io.loadmat(mat_file)
fantom_matris = data['array_3d']

if __name__ == "__main__":
    print("Keys in the .mat file:", data.keys())

    #   ----------------------------------------------------------------------
    #   PLOTTING 3D
    #   ----------------------------------------------------------------------

    # from starting_position_photon import sliced_array_njure
    # array_3d = sliced_array_njure

    # array_3d = array_3d[50:-50,50:200,600:1100]
    # array_3d = array_3d[:,25:215,:] # denna matris inkluderar kroppen, försöker skära bort luft runtomkring

    '''
    försöker skära en rimlig matris där växelverkan kan ske och fortfarande leda till energideponering i benmärgen i ryggen
    '''
    array_3d = fantom_matris[:, 25:215, 500:1150]

    viewer = visualisera_matris(array_3d)
    viewer.show()

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
