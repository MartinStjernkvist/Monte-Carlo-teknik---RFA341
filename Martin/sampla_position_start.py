from imports import *
from matriser import sliced_array_njure

def position_start():
    x_size, y_size, z_size = sliced_array_njure.shape
    # Step 1: Randomly select coordinates until a non-zero value is found
    while True:
        # Randomly select indices (x, y, z)
        x = np.random.randint(0, x_size)
        y = np.random.randint(0, y_size)
        z = np.random.randint(0, z_size)

        # Step 2: Check if the selected voxel is non-zero
        if sliced_array_njure[x, y, z] != 0:
            # print(f"Randomly selected voxel: (x={x}, y={y}, z={z}), Value: {sliced_array_njure[x, y, z]}")
            break
    return x,y,z

