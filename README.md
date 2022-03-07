# openmbr

An open source and open access reconstruction platform for optoacoustic mesoscopy based on iterative inversion through LSQR. Runs on CUDA capable GPUs and requires most likely a Linux environment. I will leave testing on Windows for other users.

## General procedure

*  load dataset from `h5` file (format specifications below)
*  load model matrix from `h5` fille (format specifications below)
*  run reconstruction and watch the absorbance converging towards a hopefully nice result

Please make sure for now that the spatial resolution of the model matrix and the dataset are the same. The reconstruction GUI does not support interpolation of the datasets to the same grid yet.



## Format specifications

### Input dataset

The input dataset contains the time domain signals measured by the hopefully spherically focused transducer  at each x any y grid position. The file should contain the following variables

*  `vol`, dataset in x, y, t order; type: `float`
*  `dim`: dimensions of the `vol` array along x, y, t; type: `uint64_t`
*  `dr`: 3 element array describing resolution along x, y, t; type: `float`
*  `origin`: 3 element array describing the origin of the first element along x, y, t; type: `float`

### Model matrix

The model matrix is a 4 dimensional array and describes the model of the transducer in its local sensitivity field. The dataset needs to contain the following variables:

*  `dim`
*  `dr`
*  `origin`
*  `data`