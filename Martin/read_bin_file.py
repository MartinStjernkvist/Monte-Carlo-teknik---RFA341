
# with open('FullPhantom.bin', 'rb') as file:
#     binary_data = file.read()
#     print(binary_data)  # This will print the binary content of the file
#
# text = binary_data.decode('utf-8')  # Decoding binary data to a string using UTF-8 encoding
# print(text)



# import pickle
#
# with open('FullPhantom.bin', 'rb') as file:
#     data = pickle.load(file)
#     print(data)  # Assuming the data is a dictionary, list, etc.


# import struct
# # Open the BMP file in binary read mode
# with open("FullPhantom.bin", "rb") as file:
#     header = file.read(14)
#     # Extract image width and height from the header
#     width, height = struct.unpack("<II", header[10:18]) # Little-endian integers
#     # Print the extracted dimensions
#     print("Image width:", width)
#     print("Image height:", height)
#     # Read the remaining image data (pixel values)
#     image_data = file.read()
#     # Process the image data (e.g., display or manipulate)
#     # This example simply prints a portion of the pixel data
#     print("First 10 pixel values:", image_data[:10])


# import numpy as np
#
# file_path = "FullPhantom.bin"
# matrix_shape = (256, 256, 1200)  # Update if necessary
# data_type = np.int16  # Check the actual data type (e.g., np.float32, np.uint8)
#
# # Read the binary file
# data = np.fromfile(file_path, dtype=data_type)
#
# # Ensure the total size matches
# expected_size = np.prod(matrix_shape)
# if data.size != expected_size:
#     raise ValueError(f"File size mismatch: expected {expected_size} elements, got {data.size}.")
#
# # Reshape the data into the matrix format
# matrix = data.reshape(matrix_shape)
#
# print("Matrix shape:", matrix.shape)
# print("Matrix data type:", matrix.dtype)




# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# # Define file path and properties
# file_path = "FullPhantom.bin"
# x, y, z = 256, 256, 1200  # Dimensions
# voxel_size = 0.15  # cm (not used directly here)
#
# # Read the binary file (int16 format, little-endian)
# data_type = np.int32
# data = np.fromfile(file_path, dtype=data_type)
#
# # Ensure the data matches the expected size
# expected_size = x * y * z
# if data.size != expected_size:
#     raise ValueError(f"File size mismatch: expected {expected_size} elements, got {data.size}.")
#
# # Reshape the data into a 3D array
# array_3d = data.reshape((x, y, z))
#
# # Visualization
# # Slice in the middle of each dimension
# fig = plt.figure(figsize=(10, 8))
#
# # Plot a slice along the Z-axis
# ax1 = fig.add_subplot(1, 3, 1)
# ax1.imshow(array_3d[:, :, z // 2], cmap="gray")
# ax1.set_title(f"Z Slice at z={z//2}")
# ax1.axis("off")
#
# # Plot a slice along the Y-axis
# ax2 = fig.add_subplot(1, 3, 2)
# ax2.imshow(array_3d[:, y // 2, :], cmap="gray", aspect="auto")
# ax2.set_title(f"Y Slice at y={y//2}")
# ax2.axis("off")
#
# # Plot a slice along the X-axis
# ax3 = fig.add_subplot(1, 3, 3)
# ax3.imshow(array_3d[x // 2, :, :], cmap="gray", aspect="auto")
# ax3.set_title(f"X Slice at x={x//2}")
# ax3.axis("off")
#
# plt.tight_layout()
# plt.show()
#
# # 3D visualization using slices
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Use the slice function to visualize 3D
# x_slice, y_slice, z_slice = x // 2, y // 2, z // 2
# ax.contourf(array_3d[:, :, z_slice], zdir='z', offset=z_slice, cmap="gray")
# ax.contourf(array_3d[:, y_slice, :], zdir='y', offset=y_slice, cmap="gray")
# ax.contourf(array_3d[x_slice, :, :], zdir='x', offset=x_slice, cmap="gray")
#
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
#
# plt.title("3D Slice Visualization")
# plt.show()






# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.widgets import Slider
#
# # Define file path and properties
# file_path = "FullPhantom.bin"
# x, y, z = 256, 256, 1200  # Dimensions
# data_type = np.int32
#
# # Read the binary file
# data = np.fromfile(file_path, dtype=data_type)
# array_3d = data.reshape((x, y, z))
#
# # Normalize the data for better visualization (optional)
# array_3d = (array_3d - np.min(array_3d)) / (np.max(array_3d) - np.min(array_3d))
#
# # Initial slice index
# initial_slice = z // 2
#
# # Create a figure for interactive visualization
# fig, ax = plt.subplots()
# plt.subplots_adjust(bottom=0.25)
# img = ax.imshow(array_3d[:, :, initial_slice], cmap="gray")
# ax.set_title(f"Z Slice: {initial_slice}")
#
# # Add a slider for Z-axis
# ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
# slider = Slider(ax_slider, 'Z Slice', 0, z - 1, valinit=initial_slice, valstep=1)
#
# # Update function for the slider
# def update(val):
#     slice_idx = int(slider.val)
#     img.set_data(array_3d[:, :, slice_idx])
#     ax.set_title(f"Z Slice: {slice_idx}")
#     fig.canvas.draw_idle()
#
# slider.on_changed(update)
#
# plt.show()



# import numpy as np
#
# # File properties
# file_path = "FullPhantom.bin"
# x, y, z = 256, 256, 1200
# data_type = np.int32  # Assumed type
#
# # Check file size
# import os
# file_size = os.path.getsize(file_path)
# expected_size = x * y * z * np.dtype(data_type).itemsize
# print(f"File size: {file_size} bytes, Expected size: {expected_size} bytes")
#
# # Try reading the file
# data = np.fromfile(file_path, dtype=data_type)
# if data.size == 0:
#     print("File is empty or contains no readable data.")
#
# # Inspect subset of data
# print("Raw data preview:", data[:10])
#
# # Test offset
# offset = 1024  # Adjust based on file inspection
# with open(file_path, "rb") as f:
#     f.seek(offset)
#     data_offset = np.fromfile(f, dtype=data_type)
#     print("Data preview after offset:", data_offset[:10])
#
# # Test alternative data types
# data_uint8 = np.fromfile(file_path, dtype=np.uint8)
# print("Preview as uint8:", data_uint8[:10])
# data_float32 = np.fromfile(file_path, dtype=np.float32)
# print("Preview as float32:", data_float32[:10])
#
#
#
# import gzip
# with gzip.open(file_path, 'rb') as f:
#     data = np.fromfile(f, dtype=np.int32)
#     print(data[:10])
#
#
# with open(file_path, "rb") as f:
#     # Read the first 100 bytes for inspection
#     header = f.read(100)
#     print(header)



import scipy.io
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

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
