# openmbr
An open source and fully freely usable model based reconstruction procedure framework for optoacoustic microscopy.


# Format of exported model matrix
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