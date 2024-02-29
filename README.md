To initialaze the repo:
``` bash
git clone --recursive https://github.com/Niccolo-Marcucci/Circular_cavity_with_spiral_outcoupler
```
or
``` bash
git clone https://github.com/Niccolo-Marcucci/Circular_cavity_with_spiral_outcoupler
git submodule update --init --recursive
```

Before running any script, verify that the submodules (and subsubmodules) are properly loaded.

To create a simulation file, run the script `s00_init_project.lsf`. This will include all the possible elements in the simulation domain. In this file choose a proper **base** name for the folder and the simulation files. 

With `s01_set_sim_parameters.lsf` you will load the simulation base file, edit it disabling the undesired components and assigning the proper values for all the various elements properties. 
Then saving it with a similar name and run the simulations.

`s03_*` is used to export farfield data