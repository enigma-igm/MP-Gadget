# First, convert the CICASS ICs to MP-Gadget format
# This needs to be done with one process and one thread.
# (note that cdm and xray use the same ICs)
export OMP_NUM_THREADS=1
mpirun -np 1 ../../genic/MP-GenIC-CICASS paramfile_cdm.genic || exit 1
mpirun -np 1 ../../genic/MP-GenIC-CICASS paramfile_wdm.genic || exit 1

# Now run the hydro
# Threads optimized for igm, but you might need NMPI=1 on Mac
export NMPI=32
export OMP_NUM_THREADS=4
# Each of these takes about 5 minutes to run to z = 5 on igm.
mpirun -np $NMPI --use-hwthread-cpus ../../gadget/MP-Gadget paramfile_cdm.gadget || exit 1
mpirun -np $NMPI --use-hwthread-cpus ../../gadget/MP-Gadget paramfile_wdm.gadget || exit 1
mpirun -np $NMPI --use-hwthread-cpus ../../gadget/MP-Gadget paramfile_xray.gadget || exit 1
