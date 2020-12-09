
-----------------------------------
INTERFACE BETWEEN LAMMPS AND AELAS
-----------------------------------

In this example we show an interface between ```LAMMPS``` and ```AELAS``` to calculate the elastic constants. The programs ```LAMMPS``` and ```AELAS``` can be downloaded here

[https://lammps.sandia.gov](https://lammps.sandia.gov)

[https://sites.google.com/site/zrfbuaa/softwares/aelas](https://sites.google.com/site/zrfbuaa/softwares/aelas)

In the file ```workflow_lmp_aelas.pdf``` you can see a diagram of the workflow of this interface.
As you can see, it makes use of the program ```ATOMSK``` in order to convert some files, you can download it here 

[https://atomsk.univ-lille.fr](https://atomsk.univ-lille.fr)

To illustrate this example we apply it to BCC Fe. The folder ```./src``` contains the three general source files ```gen_pos.sh```, ```run_lmp.sh``` and ```OSZICAR```.
To run this example make a new folder

```mkdir ./BCC_Fe/output_files_test```

and copy all files in folders ```./src/``` and in ```./BCC_Fe/input_files/``` into the folder ```./BCC_Fe/output_files_test```, that is

```cp ./src/* ./BCC_Fe/output_files_test/```

```cp ./BCC_Fe/input_files/* ./BCC_Fe/output_files_test/```

Next, go to the folder ```./BCC_Fe/output_files_test```

```cd ./BCC_Fe/output_files_test/```

and type 

```./gen_pos.sh```

This will generate the required input ```POSCAR``` files for the ```AELAS``` code. Next, type

```./run_lmp.sh```

It will calculate the energy for each distorted cell with ```LAMMPS``` and create ```OSZICAR```-like files, so that ```AELAS``` can read this energies and compute the elastic constants. The calculated elastic constants are printed in the file called ```ELADAT```. Please, check that you get the same files and results that those given in the folder ```./BCC_Fe/output_files/```.

More details can be found in reference: 
P. Nieves et al. “Spin lattice model for cubic crystals” arXiv (2020)