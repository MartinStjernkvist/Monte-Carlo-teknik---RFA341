#Monte Carlo (fantomen)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec
# Hitta filen
file_path = "C:/Users/Admin/Documents/GitHub/Monte Carlo Linnea/FullPhantom.bin"
x, y, z = 256, 256, 1200  # Antal voxlar i x-dimensionen

voxel_size = 0.15  # Voxelstorlek i cm

# Läs in de binära data från filen
try:
    raw_data = np.fromfile(file_path, dtype=np.float32)
except FileNotFoundError:
    print("Fil inte hittad!")
    exit()

# Omforma den råa datan till en 3D-array
array_3d = raw_data.reshape((x, y, z))
""""
# Visualisera med sliceViewer (skärvisualisering) i 2D
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.imshow(array_3d[:, :, z//2], cmap='gray')  # Välj mitten-slice för visualisering
ax.set_title('Mitten-slice av 3D-array')
plt.show()
"""

# Skapar figuren i 2D slice 
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(bottom=0.25)

#Skapar slice i z-axeln
initial_slice_index = z // 2  # Start at the middle slice
img = ax.imshow(array_3d[:, :, initial_slice_index], cmap='gray')
ax.set_title(f'Slice i z={initial_slice_index}')
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


# Voxelvärde som står för båda njurarna i arrayen 
Njurar_index=[np.where(array_3d == 17), np.where(array_3d == 18) ,np.where(array_3d == 19),np.where(array_3d == 20),np.where(array_3d == 21) ,np.where(array_3d == 22)]

#Skapar en annan array med bara njurarna och sätter in alla voxlar som har njurarna
Array_njurar=np.full_like(array_3d,0)
for i in range(len(Njurar_index)):
    Array_njurar(i)==1

#Skapar slice i z-axeln för njurarna
initial_slice_index = z // 2  # Start at the middle slice
img = ax.imshow(Array_njurar[:, :, initial_slice_index], cmap='gray')
ax.set_title(f'Slice i z={initial_slice_index}')
plt.colorbar(img, ax=ax)

# Add a slider for navigating through slices
ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Slice Index', 0, z-1, valinit=initial_slice_index, valstep=1)

# Update function for the slider
def update(val):
    slice_index = int(slider.val)
    img.set_data(Array_njurar[:, :, slice_index])
    ax.set_title(f'Slice {slice_index}')
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()
# Create a larger figure with a GridSpec layout
fig = plt.figure(figsize=(18, 10))  # Increase figure size
gs = GridSpec(1, 3, width_ratios=[3, 3, 1], height_ratios=[2])  # Adjust ratios to make plots fill space

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

plt.show()

