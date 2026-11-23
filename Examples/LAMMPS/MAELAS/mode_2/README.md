
--------------------------------------------
INTERFACE BETWEEN LAMMPS AND MAELAS (-mode 2)
--------------------------------------------

In this example we show an interface between ```LAMMPS``` and ```MAELAS``` with ```-mode 2``` to calculate the magnetoelastic constants. The program ```LAMMPS``` can be downloaded here

[https://lammps.sandia.gov](https://lammps.sandia.gov)

As you can see, it makes use of the program ```ATOMSK``` in order to convert some files, you can download it here 

[https://atomsk.univ-lille.fr](https://atomsk.univ-lille.fr)

To illustrate this example we apply it to BCC Fe. The folder ```./input``` contains the source files.
To run this example, go to the folder ```./input``` and type

```chmod +x run.sh```

```./run.sh```

The calculated magnetoelastic constants are printed in the file called ```results_derive.out```. Please, check that you get the same files and results that those given in the folder ```./output/```.

