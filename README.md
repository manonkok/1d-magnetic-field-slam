# Online One-Dimensional Magnetic Field SLAM with Loop-Closure Detection

This repository is the official implementation of the methods in the publication:

* Manon Kok and Arno Solin (2024). **Online One-Dimensional Magnetic Field SLAM with Loop-Closure Detection**. In *IEEE International Conference on Multisensor Fusion and Integration*.

### Motivation

We present a lightweight magnetic field simultaneous localisation and mapping (SLAM) approach for drift correction in odometry paths, where the interest is purely in the odometry and not in map building. We represent the past magnetic field readings as a one-dimensional trajectory against which the current magnetic field observations are matched. This approach boils down to sequential loop-closure detection and decision-making, based on the current pose state estimate and the magnetic field. We combine this setup with a path estimation framework using an extended Kalman filter and smoother which fuse the odometry increments with the detected loop-closure timings. We demonstrate the practical applicability of the model with several different real-world examples from a handheld iPad moving in indoor scenes.

https://github.com/user-attachments/assets/b77933b3-a214-4bcf-8241-0d638c2b8654

## Dependencies

The codes in this repository have been tested with **Mathworks MATLAB R2024a (Update 1)**. The core functions (those in `src`) do not depend on any additional toolboxes, but helper files in `tools`. The main file runSLAM.m and the plotting file makePlots.m in tools use the procrustes function from the following built-in toolboxes:
* Statistics and Machine Learning Toolbox (tested with version 24.1)
to compute the RMSE.

## Structure of the codes

### Main file to run

```
  runSLAM - Main file to make all results
```

### Core functions (under `src`)

```
  magSLAMwithLoopClosures - Run 1D magnetic field SLAM with loop closures
  run_filter_from_scratch - Runs the EKF from the start of the data set
```

### Helper functions (under `tools`)

```
             dynamics - Dynamic model
            makePlots - Generates plots from the paper
	  prepareData - Prepares data to be used in EKF, also adds a random noise realization to the odometry
    	     quat2eul - Converts quaternion to Euler angles
                 rotx - Computes a rotation matrix from a rotation angle around the x-axis
		 roty - Computes a rotation matrix from a rotation angle around the y-axis
		 rotz - Computes a rotation matrix from a rotation angle around the z-axis
```

## Data availability and access

Data used in the examples are provided in the folder `data`.

## Running experiment examples

In the `runSLAM.m` file it is possible to choose which of the four data sets to run in the variable `indDataSet` on line 41.

## License

This software is provided under the [MIT License](LICENSE).
