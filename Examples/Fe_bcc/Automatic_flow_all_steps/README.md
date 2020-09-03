--------------------------------------
EXAMPLE: Fe bcc (automatic flow)
--------------------------------------


In this folder there is a procedure to run all steps in an automatic way using the bash script called "automatic_flow.sh".

The POSCAR and POTCAR files are the only inputs that we need to calculate the anisotropic magnetostrictive coefficients. 
Here, we also include the file ELADAT that contains the elastic tensor that will be used to calculate the magnetoelastic constants in the last step. 
The POSCAR (called POSCAR_Fe_bcc) and ELADAT files should be in the same folder as the bash script automatic_flow.sh.
The POTCAR is copied from a folder containing VASP potentials which is set in the beginning of the script.
To run this example one only needs to type:
```bash
nohup ./automatic_flow.sh > report.out &
```
In the file report.out and folder ./results you can find the generated outputs (we manually removed POTCAR files after all jobs finished).

Note that in this example we used a small k-mesh (-k 60) in order to just check this script quickly.
Hence, the calculated MAE and magnetostriction coefficients are not reliable.
Much higher number of k-points are required to obtain reliable and well-converged results. 
For instance, in our tests for Fe-bcc we set -k 180. 
