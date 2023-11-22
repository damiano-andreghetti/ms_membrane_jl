This julia software runs simulation of an elastic membrane on which particles attach and locally impose a rigidity different from that of bulk membrane.
The main function is run_sim(...) inside MCLangevin.jl
In the following surface is discretized into L by L lattice of spacing a.
Needed arguments are (in order):
h: Matrix L by L of Floats, which are the heights of the membrane at each lattice site

part: Matrix L by L of Bool, which is the occupation of each lattice site (1 if particle present 0 otherwise)

ld: when sampling new membrane configurations each of them is obtained displacing each height of the membrane within the interval of size 2*ld around the previous position. Changing this affects acceptance of MC moves for membrane

N: number of total MC sweeps for the simulation 
(1 sweep= go through each site in random order and attempt displacement of hegiht of memrbane at that site+ update through Langevin eq the displacement of particles+if dispalcement exceed threshold a/2 make particle jump)

Nav: number of sweeps over which take moving average of effective diffusion rate to obtain smooth values

Ninf: number of sweeps between two consecutive file IO (both saving configurations and measures of observables)

k: rigidity of membrane where particle is present (kBT units)

k0: rigidity of bulk membrane (kBT units)

W: direct interaction parameter (kBT units)

sigma: surface tension parameter (kBT/a^2 units)

phi_kd_ratio: ratio between phi=kI(1-rho) and kd that we want to have for the simulation. At each sweep insertion rate kI is updated based on measured values of rho and kd in roder to keep phi/kd fixed.

Ne: extraction size. Clusters (connected components) containing >=Ne particles are removed

folder: string with the name of the folder where results will be saved
