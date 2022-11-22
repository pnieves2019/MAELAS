
--------------------------------------------
INTERFACE BETWEEN QuantumEspresso AND MAELAS (-mode 3)
--------------------------------------------

In this example we show an interface between ```QuantumEspresso``` and ```MAELAS``` with ```-mode 3``` to calculate the isotropic magnetoelastic constants, isotropic magnetostrictive coefficients and spontaneous volume magnetostriction without spin-orbit coupling. 
The program ```QuantumEspresso``` can be downloaded here

[https://www.quantum-espresso.org/](https://www.quantum-espresso.org/)

This interface makes use of the program ```ATOMSK``` in order to convert some files, you can download it here 

[https://atomsk.univ-lille.fr](https://atomsk.univ-lille.fr)

To illustrate this example we apply it to bcc Fe. The folder ```./input``` contain 5 files:

1) ```run.sh``` the executable file

2) ```jsub``` job script to run this example in supercomputer facilities

3) ```OSZICAR``` template for the OSZICAR file format

4) ```structure.pw``` relaxed structure in ```QuantumEspresso``` format (in angstrom units!)

5) ```ELADAT``` elastic tensor calculated with program AELAS, see example [https://github.com/pnieves2019/MAELAS/tree/master/Examples/QuantumEspresso/AELAS](https://github.com/pnieves2019/MAELAS/tree/master/Examples/QuantumEspresso/AELAS)

To run this example in your local machine just type

```./run.sh```

or to run it in supercomputer facilities

```qsub jsub```

The obtained results are shown in folder ```./output```. Note that ```ATOMSK``` generates a default input file for ```QuantumEspresso```. 
In this example, this default input file is modified (e.g, pseudopotential folder, k-point mesh, cut-off, starting magnetization, ...) through the script ```run.sh```.
The calculated isotropic magnetoelastic constants, isotropic magnetostrictive coefficients and spontaneous volume magnetostriction are printed in the file ```output_derivation.dat```.