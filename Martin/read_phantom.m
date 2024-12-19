% Define the file path and file properties
file_path = "C:\MATLAB\Monte-Carlo teknik - RFA341\FullPhantom.bin";
x = 256; % Replace with your desired dimensions
y = 256;
z = 1200;
voxelSize = 0.15; % cm
% Open the binary file for reading
fid = fopen(file_path, 'rb');
if fid == -1
error('Could not open the file.');
end
% Read the binary data
raw_data = fread(fid, x * y * z, 'single', 0, 'ieee-le');
% Close the file
fclose(fid);
% Reshape the raw data into a 3D array
array_3d = reshape(raw_data, [x, y, z]);
%% Reading data into image slices
figure;
orthosliceViewer(array_3d)
figure;
sliceViewer(array_3d)





