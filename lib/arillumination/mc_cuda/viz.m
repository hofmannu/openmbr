clear all;
close all;

nA = 601;
nR = 701;
fileID = fopen("outputFile.bin");
clear A;
A = fread(fileID, 1e9, 'single');

A = reshape(A, [nR, nA]);

r = (0:(nR - 1)) * 10e-6;
z = (0:(nA - 1)) * 10e-6;

imagesc(r, z, log(A)')
xlabel('radial distance [m]');
ylabel('axial distance [m]')
colormap(bone(1024));
colorbar;

