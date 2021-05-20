# openmbr
An open source and fully freely usable model based reconstruction procedure framework for optoacoustic microscopy.

## Installation

To use this toolbox I recommend running it on Archlinux. The interface is made with ImGUI and allows the user to define most of the settings. From the hardware perspective you need a CUDA capable GPU for now since the computational requirements are quite high for model based reconstruction. CPU version in the fourier domain is still under development.

The following libraries are required to compile the code
*  `hdf5` file import and export, saving and loading of large datasets
*  `cuda` and `nvidia` or `nvidia-lts` 
*  `cmake` and `make` for compiling
*  `glfw-x11`, `glew` used to display stuff
*  `nlohmann-json` to save settings of reconstruction to json file

Archlinux installation command

```
pacman -S nlohmann-json cuda nvidia hdf5 cmake make glfw-x11 glew
```

Ubuntu installation command
```
apt-get install libhdf5-dev nvidia-cuda-toolkit cmake make libglfw3-dev libglfw3-dev 
```

Building the code through:
```
git clone git@github.com:hofmannu/openmbr.git
cd yamct
git submodule init
git submodule update
mkdir Debug
cd Debug
cmake .. && make all && ./gui
```

## Model building

## Format of datasets required for import and export

### Scan dataset

### Experimentally determined optoacoustic waveform

### Format of exported model matrix
For reusing of the acoustic or acoustic-optical models, the model matrix can be exported to a `h5` file. Variables in this file are:

*  `ntModel` number of elements along time domain 
*  `nrModel` number of elements along radial domain
*  `nzModel` number of elements along vertical domain
*  `rRes` radial resolution
*  `zRes` axial resolution
*  `tRes` temporal resolution
*  `t0` zero timepoint
*  `z0` zero axial point
*  `data` actual data matrix
*  `sir` spatial impulse response for debugging mostly
*  `nzSir` number of element in spatial impulse response along z

The data is thereby stored in the following order: `[it, ir, iz]` and needs to be reshaped after loading into MATLAB. 
