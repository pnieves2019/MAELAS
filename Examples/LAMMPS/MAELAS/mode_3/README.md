---------------------------------------
Calculation of volume magnetostriction in LAMMPS
---------------------------------------

In this folder you can find two ways to compute volume magnetostriction (w_s) in LAMMPS:

1) The first method is based on ```-mode 3``` of MAELAS. To run this calculation go to folder ```calculation_volume_magnetostriction_MAELAS```. Here, there are two folders called ```input``` and ```output```. Copy all files in the folder ```input``` into the folder ```output```, go to folder ```output``` and execute the file called ```run``` 

```bash
./run
```  

This will generate a file called ```output.dat```, where you can find the calculated volume magnetostriction. In this example we obtain:
 
```bash
w_s=-0.010028
```  

This calculation is very sensitive to the lattice parameter used the initial POSCAR file (INPOS.lmp) which is the reference state from which MAELAS generates the deformations. Here, one needs to use the lattice parameter that corresponds to the equilibrium paramagnetic state (a0=2.86642687316 A, v_PM=23.551718415 A^3), which is found in the second method by fitting the energy vs volume to the Murnaghan equation of state (EOS) without offset in the exchange energy. In this case, the energy of the deformed unit cells should be evaluated without offset in the exchange energy (see files ```in.aelas``` and ```in.maelas```)

```bash
pair_coeff 	* * spin/exchange exchange 3.5 0.02726 0.2171 1.841 offset no 
``` 
 

2) The second method derives the volume magnetostriction using the equation w_s = (v_FM - v_PM)/v_PM, where v_FM and v_PM are the equilibrium volume at the ferromagnetic and paramagnetic states, respectively. These equilibrium volumes are calcualted by fitting the energy vs volume to the Murnaghan equation of state (EOS), without the offset in the exchange energy 

```bash
pair_coeff 	* * spin/exchange exchange 3.5 0.02726 0.2171 1.841 offset no 
``` 

To run this calculation go to folder ```calculation_volume_magnetostriction_alternative_method```. Here, there are two folders called ```input``` and ```output```. Copy all files in the folder ```input``` into the folder ```output```, go to folder ```output``` and execute the file called ```run``` 

```bash
./run
```  

This will generate a file called ```output.dat```, where you can find the calculated volume magnetostriction. In this example we obtain:
 
```bash
w_s=-0.0098989
```  

which is very close to the value obtained with the first method. A similar value for w_s is obtained if the offset is included. Note that the volume magnetostriction depends on the parameterization of the exchange parameters through the Bethe-Slater function. More details about the influence of this parameterization on w_s can be found in the following reference: 

P. Nieves, et al., “Spin-lattice model for cubic crystals”, Phys. Rev. B 103, 094437 (2021)

[https://doi.org/10.1103/PhysRevB.103.094437](https://doi.org/10.1103/PhysRevB.103.094437)



