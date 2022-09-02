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
w_s=-0.0101021
```  

Note that to get consistent results with the second method one needs consider the following two things:

	- This calculation is very sensitive to the lattice parameter used the initial POSCAR file (INPOS.lmp). Here, one needs to use the lattice parameter that corresponds to the equilibrium ferromagnetic state (a0=2.8664962944 A, v_FM=23.55342963513 A^3), which is found in the second method by fitting the energy vs volume to the Murnaghan equation of state (EOS) including the offset in the exchange energy (the offset should be included to find the equilibrium volume due to a pathological behaviour of classical spin-lattice model where the pressure is not zero at the minimum energy without offset).

	- However, to compute w_s the offset in the exchange energy must not be included because it sets the exchange energy to zero for any deformation at the ferromagentic state, which modifies the strain dependence of the isotropic magnetoelastic energy artificially.

```bash
pair_coeff 	* * spin/exchange exchange 3.4 0.02726 0.2171 1.841 offset no 
``` 
 

2) The second method derives the volume magnetostriction using the equation w_s = (v_FM - v_PM)/v_PM, where v_FM and v_PM are the equilibrium volume at the ferromagnetic and paramagnetic states, respectively. These equilibrium volumes are calcualted by fitting the energy vs volume to the Murnaghan equation of state (EOS), including the the offset in the exchange energy 

```bash
pair_coeff 	* * spin/exchange exchange 3.4 0.02726 0.2171 1.841 offset yes 
``` 

To run this calculation go to folder ```calculation_volume_magnetostriction_alternative_method```. Here, there are two folders called ```input``` and ```output```. Copy all files in the folder ```input``` into the folder ```output```, go to folder ```output``` and execute the file called ```run``` 

```bash
./run
```  

This will generate a file called ```output.dat```, where you can find the calculated volume magnetostriction. In this example we obtain:
 
```bash
w_s=-0.010026486
```  

which is very close to the value obtained with the first method. A similar value for w_s is obtained if the offset is not included. Note that the volume magnetostriction depends on the parameterization of the exchange parameters through the Bethe-Slater function. More details about the influence of this parameterization on w_s can be found in the following reference: 

P. Nieves, et al., “Spin-lattice model for cubic crystals”, Phys. Rev. B 103, 094437 (2021)

[https://doi.org/10.1103/PhysRevB.103.094437](https://doi.org/10.1103/PhysRevB.103.094437)



