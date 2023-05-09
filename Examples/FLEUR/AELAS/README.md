
--------------------------------------------
INTERFACE BETWEEN FLEUR AND AELAS
--------------------------------------------

In this example we show an interface between ```FLEUR``` and ```AELAS``` to calculate the elastic constants. 
The program ```FLEUR``` can be downloaded here

[https://www.flapw.de/master/downloads/](https://www.flapw.de/master/downloads/)

This interface makes use of the programs ```vasp2cif``` and ```cif2cell``` in order to convert VASP format into FLEUR format, you can download them here 

[https://github.com/egplar/vasp2cif](https://github.com/egplar/vasp2cif)

[https://github.com/torbjornbjorkman/cif2cell](https://github.com/torbjornbjorkman/cif2cell)

To illustrate this example we apply it to bcc Fe. The folder ```./input``` contain 5 files:

1) ```run.sh``` bash script to generate the deformed unit cells with ```AELAS``` and run ```FLEUR``` 

2) ```jsub``` job script to run this example in supercomputer facilities

3) ```OSZICAR``` template for the OSZICAR-like file format. Calculated energy will be overwritten on the word ```energy```, which is the place where ```AELAS``` reads the energy of VASP OSZICAR files 

4) ```POSCAR``` relaxed structure in ```VASP``` format 

5) ```cp_osz.sh``` bash script to extract energies from output files of ```FLEUR``` and derive elastic constants with ```AELAS```

To run this example first type

```chmod +x run.sh```

```./run.sh```

Once all ```FLEUR``` calculations are correctly done, then type

```chmod +x cp_osz.sh```

```./cp_osz.sh```

The obtained results are shown in folder ```./output```. In this example, the default input file for ```FLEUR``` is modified
 (e.g, k-point mesh, plane-wave cut-off, ...) through the script ```run.sh```. The calculated elastic constants are printed in the file ```ELADAT```.