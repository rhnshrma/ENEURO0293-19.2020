`insilico`
========

![GNU GPLv3 License](http://img.shields.io/badge/license-GPLv3-green.svg)
[![insilico Trello](https://img.shields.io/badge/Trello-insilico-blue.svg)](https://trello.com/b/lkXzPGqD/insilico)

`insilico` is a Computational Neuroscience simulation library written in C++. `insilico` encourages ready-to-start approach for quick setup of simulation environment, without hindering programmers time and focus from intended experiment.

Library homepage: http://www.iiserpune.ac.in/~collins/insilico/

Install
=======

Please refer to `doc/INSTALL` file.

Build
=====
Clone the source and run the following commands on terminal.
```
   g++ -Wall -std=c++11  -I ./include -o insilico.out TH_network_mAch.cpp
```
use flag `-Ofast` to speed up simulation (no inaccuracies noticed)

Execute
=======

Run the python script `gen_nsets_Ca_N.py` to generate the different neuron configuration files

Run the following command on terminal to execute the code.
```
  insilico.out -o <output_file>.csv -n <neuron_file.isf> -s <synapse_file.isf> -t <time>

  Options:
    -o   Output file
    -n   Neuron configuration file
    -s   Synapse configuration file (ssets.isf)
    -t   simulation time in milliseconds (set to 15000 for 15 second simulations)

```
Please read `doc/FILES` file for details about the input and output files and their formats.

License
=======

This simulator library is licensed under GNU GPLv3 which can be found in `LICENSE` file under home directory of this project.


Biophysical basis of alpha rhythm disruption in Alzheimer’s Disease (AD)
=======

# https://doi.org/10.1523/ENEURO.0293-19.2020


