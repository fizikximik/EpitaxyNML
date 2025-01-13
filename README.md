# Edited for new version of numpy and scipy! Last update: 13.01.2025
EpitaxyNML is an application to check the commensurability between two slabs based on the paper of Zur and McGill, [J. Appl. Phys. 55 (2)](https://aip.scitation.org/doi/10.1063/1.333084)

Only two POSCAR files are required as input. 

Results are saved in results.dat file which contains:

* Supercell sizes of the thinfilm, N_TF, and the substrate, N_Sub.
* Percentual mismatch between the surfaces (Err_a, Err_b, Err_alpha). 
* Generator matrices for the supercells (Mat_TF, MatSub).

You can contact J. J. Plata for more information or reporting bugs at: jplata@us.es.
