# 5G-SLAM
This repo contains a 5G simultaneous localization and mapping (SLAM) system written in the Matlab language. This code computes the performances of the vehicle state estimation and the mapping of the environment in 5G mmWave vehicular networks.


## Summary
The create_measurement.m generates vehicle trajectories, three-types of objects (base station, virtual anchors, scattering points), and clutter.

The main.m is used to estimate the vehicle state and the objects' state, and compute the RMSE of the vehicle sate and GOSPA of the objects.


## Authors
The code was developed by Yu Ge, Ph.D. student at Chalmers University of Technology, Gothenburg, Sweden, contact email yuge@chalmers.se

For more information, please see
```
@article{ge20205g,
  title={5G SLAM using the clustering and assignment approach with diffuse multipath},
  author={Ge, Yu and Wen, Fuxi and Kim, Hyowon and Zhu, Meifang and Jiang, Fan and Kim, Sunwoo and Svensson, Lennart and Wymeersch, Henk},
  journal={Sensors},
  volume={20},
  number={16},
  pages={4656},
  year={2020},
  publisher={Multidisciplinary Digital Publishing Institute}
}

```


## License
This project is licensed under the project of Chalmers University of Technology.


## Configuration Instructions
A working and licensed installation of Matlab is required. The version should be no earlier than 2016b.


## Installation Instructions
Clone/download all included files, and there is no further installation required.


## Operating Instructions
Run the file main.m, or type main directly to the command line.
To run the test, type run(Testclass) in the command line.
