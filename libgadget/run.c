#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "utils.h"

#include "allvars.h"
#include "walltime.h"
#include "gravity.h"
#include "density.h"
#include "domain.h"
#include "run.h"
#include "init.h"
#include "cooling.h"
#include "checkpoint.h"
#include "petaio.h"
#include "petapm.h"
#include "timestep.h"
#include "drift.h"
#include "forcetree.h"
#include "blackhole.h"
#include "hydra.h"
#include "sfr_eff.h"
#include "metal_return.h"
#include "slotsmanager.h"
#include "hci.h"
#include "fof.h"
#include "cooling_qso_lightup.h"
#include "lightcone.h"
#include "timefac.h"

/* stats.c only used here */
void energy_statistics(FILE * FdEnergy, const double Time,  struct part_manager_type * PartManager);
/*!< file handle for energy.txt log-file. */
static FILE * FdEnergy;
static FILE  *FdCPU;    /*!< file handle for cpu.txt log-file. */
static FILE *FdSfr;     /*!< file handle for sfr.txt log-file. */
static FILE *FdBlackHoles;  /*!< file handle for blackholes.txt log-file. */
static FILE *FdBlackholeDetails;  /*!< file handle for BlackholeDetails binary file. */
static FILE *FdHelium; /* < file handle for the Helium reionization log file helium.txt */

static struct ClockTable Clocks;

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */
static void write_cpu_log(int NumCurrentTiStep, const double atime, FILE * FdCPU);

/* Updates the global storing the current random offset of the particles,
 * and stores the relative offset from the last random offset in rel_random_shift*/
static void update_random_offset(double * rel_random_shift);
static void set_units();
/* Construct a unit system struct*/
static struct UnitSystem get_unitsystem(void);

static void
open_outputfiles(int RestartSnapNum);

static void
close_outputfiles(void);

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file). Please avoid adding new variables in favour of things which are local to a module.
 */
struct global_data_all_processes All;

/*Set the global parameters*/
void
set_all_global_params(ParameterSet * ps)
{
    int ThisTask;
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    if(ThisTask == 0) {
        /* Start reading the values */
        param_get_string2(ps, "InitCondFile", All.InitCondFile, sizeof(All.InitCondFile));
        param_get_string2(ps, "OutputDir", All.OutputDir, sizeof(All.OutputDir));
        param_get_string2(ps, "SnapshotFileBase", All.SnapshotFileBase, sizeof(All.SnapshotFileBase));
        param_get_string2(ps, "FOFFileBase", All.FOFFileBase, sizeof(All.FOFFileBase));
        param_get_string2(ps, "EnergyFile", All.EnergyFile, sizeof(All.EnergyFile));
        All.OutputEnergyDebug = param_get_int(ps, "OutputEnergyDebug");
        param_get_string2(ps, "CpuFile", All.CpuFile, sizeof(All.CpuFile));

        All.CP.CMBTemperature = param_get_double(ps, "CMBTemperature");
        All.CP.RadiationOn = param_get_int(ps, "RadiationOn");
        All.CP.Omega0 = param_get_double(ps, "Omega0");
        All.CP.OmegaBaryon = param_get_double(ps, "OmegaBaryon");
        All.CP.OmegaLambda = param_get_double(ps, "OmegaLambda");
        All.CP.Omega_fld = param_get_double(ps, "Omega_fld");
        if(All.CP.OmegaLambda > 0 && All.CP.Omega_fld > 0)
            endrun(0, "Cannot have OmegaLambda and Omega_fld (evolving dark energy) at the same time!\n");
        All.CP.w0_fld = param_get_double(ps,"w0_fld");
        All.CP.wa_fld = param_get_double(ps,"wa_fld");
        All.CP.Omega_ur = param_get_double(ps, "Omega_ur");
        All.CP.HubbleParam = param_get_double(ps, "HubbleParam");

        All.OutputDebugFields = param_get_int(ps, "OutputDebugFields");

        All.TimeMax = param_get_double(ps, "TimeMax");
        All.Asmth = param_get_double(ps, "Asmth");
        All.ShortRangeForceWindowType = (enum ShortRangeForceWindowType) param_get_enum(ps, "ShortRangeForceWindowType");
        All.Nmesh = param_get_int(ps, "Nmesh");

        All.CoolingOn = param_get_int(ps, "CoolingOn");
        All.HydroOn = param_get_int(ps, "HydroOn");
        All.DensityOn = param_get_int(ps, "DensityOn");
        All.TreeGravOn = param_get_int(ps, "TreeGravOn");
        All.LightconeOn = param_get_int(ps, "LightconeOn");
        All.FastParticleType = param_get_int(ps, "FastParticleType");
        All.PairwiseActiveFraction = param_get_double(ps, "PairwiseActiveFraction");
        All.TimeLimitCPU = param_get_double(ps, "TimeLimitCPU");
        All.AutoSnapshotTime = param_get_double(ps, "AutoSnapshotTime");
        All.TimeBetweenSeedingSearch = param_get_double(ps, "TimeBetweenSeedingSearch");
        All.RandomParticleOffset = param_get_double(ps, "RandomParticleOffset");

        All.PartAllocFactor = param_get_double(ps, "PartAllocFactor");
        All.SlotsIncreaseFactor = param_get_double(ps, "SlotsIncreaseFactor");

        All.SnapshotWithFOF = param_get_int(ps, "SnapshotWithFOF");

        All.RandomSeed = param_get_int(ps, "RandomSeed");

        All.BlackHoleOn = param_get_int(ps, "BlackHoleOn");
        All.WriteBlackHoleDetails = param_get_int(ps,"WriteBlackHoleDetails");

        All.StarformationOn = param_get_int(ps, "StarformationOn");
        All.MetalReturnOn = param_get_int(ps, "MetalReturnOn");
        All.MaxDomainTimeBinDepth = param_get_int(ps, "MaxDomainTimeBinDepth");

        /*Massive neutrino parameters*/
        All.MassiveNuLinRespOn = param_get_int(ps, "MassiveNuLinRespOn");
        All.HybridNeutrinosOn = param_get_int(ps, "HybridNeutrinosOn");
        All.CP.MNu[0] = param_get_double(ps, "MNue");
        All.CP.MNu[1] = param_get_double(ps, "MNum");
        All.CP.MNu[2] = param_get_double(ps, "MNut");
        All.HybridVcrit = param_get_double(ps, "Vcrit");
        All.HybridNuPartTime = param_get_double(ps, "NuPartTime");
        if(All.MassiveNuLinRespOn && !All.CP.RadiationOn)
            endrun(2, "You have enabled (kspace) massive neutrinos without radiation, but this will give an inconsistent cosmology!\n");
        /*End massive neutrino parameters*/

        if(All.StarformationOn != 0 && All.CoolingOn == 0)
        {
                endrun(1, "You try to use the code with star formation enabled,\n"
                          "but you did not switch on cooling.\nThis mode is not supported.\n");
        }
    }
    MPI_Bcast(&All, sizeof(All), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
int begrun(int RestartFlag, int RestartSnapNum)
{
    if(RestartFlag == 1) {
        RestartSnapNum = find_last_snapnum(All.OutputDir);
        message(0, "Last Snapshot number is %d.\n", RestartSnapNum);
    }

    hci_init(HCI_DEFAULT_MANAGER, All.OutputDir, All.TimeLimitCPU, All.AutoSnapshotTime, All.SnapshotWithFOF);

    petapm_module_init(omp_get_max_threads());
    petaio_init();
    walltime_init(&Clocks);

    petaio_read_header(RestartSnapNum);

    slots_init(All.SlotsIncreaseFactor * PartManager->MaxPart, SlotsManager);
    /* Enable the slots: stars and BHs are allocated if there are some,
     * or if some will form*/
    if(All.NTotalInit[0] > 0)
        slots_set_enabled(0, sizeof(struct sph_particle_data), SlotsManager);
    if(All.StarformationOn || All.NTotalInit[4] > 0)
        slots_set_enabled(4, sizeof(struct star_particle_data), SlotsManager);
    if(All.BlackHoleOn || All.NTotalInit[5] > 0)
        slots_set_enabled(5, sizeof(struct bh_particle_data), SlotsManager);

    set_units();

#ifdef DEBUG
    char * pidfile = fastpm_strdup_printf("%s/%s", All.OutputDir, "PIDs.txt");

    MPIU_write_pids(pidfile);
    myfree(pidfile);
#endif

    init_forcetree_params(All.FastParticleType);

    const struct UnitSystem units = get_unitsystem();
    init_cooling_and_star_formation(All.CoolingOn, All.StarformationOn, &All.CP, All.MassTable[0], All.BoxSize, units);

    All.MinEgySpec = get_MinEgySpec();

    gravshort_fill_ntab(All.ShortRangeForceWindowType, All.Asmth);

    set_random_numbers(All.RandomSeed);

    if(All.LightconeOn)
        lightcone_init(&All.CP, All.TimeInit, All.UnitLength_in_cm, All.OutputDir);
    return RestartSnapNum;
}

/* Small function to decide - collectively - whether to use pairwise gravity this step*/
static int
use_pairwise_gravity(ActiveParticles * Act, struct part_manager_type * PartManager)
{
    /* Find total number of active particles*/
    int64_t total_active, total_particle;
    MPI_Allreduce(&Act->NumActiveParticle, &total_active, 1, MPI_INT64, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&PartManager->NumPart, &total_particle, 1, MPI_INT64, MPI_SUM, MPI_COMM_WORLD);

    /* Since the pairwise step is O(N^2) and tree is O(NlogN) we should scale the condition like O(N)*/
    return total_active < All.PairwiseActiveFraction * total_particle;
}

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void
run(int RestartSnapNum)
{
    /*Number of timesteps performed this run*/
    int NumCurrentTiStep = 0;
    /*Is gas physics enabled?*/
    int GasEnabled = All.NTotalInit[0] > 0;

    int SnapshotFileCount = RestartSnapNum;
    PetaPM pm = {0};
    gravpm_init_periodic(&pm, All.BoxSize, All.Asmth, All.Nmesh, All.CP.GravInternal);

    /* ... read initial model and initialise the times*/
    inttime_t ti_init = init(RestartSnapNum, All.TimeIC, All.TimeInit, All.TimeMax, &All.CP, All.SnapshotWithFOF, All.MassiveNuLinRespOn, All.MassTable, All.NTotalInit);

    DriftKickTimes times = init_driftkicktime(ti_init);

    DomainDecomp ddecomp[1] = {0};
    domain_decompose_full(ddecomp);	/* do initial domain decomposition (gives equal numbers of particles) */

    if(All.DensityOn) {
        double uu_in_cgs = All.UnitEnergy_in_cgs / All.UnitMass_in_g;
        setup_smoothinglengths(RestartSnapNum, ddecomp, &All.CP, All.BlackHoleOn, All.MinEgySpec, uu_in_cgs, ti_init, All.TimeInit, All.NTotalInit[0]);
    }

    /* Stored scale factor of the next black hole seeding check*/
    double TimeNextSeedingCheck = All.TimeInit;

    walltime_measure("/Misc");

    open_outputfiles(RestartSnapNum);

    write_cpu_log(NumCurrentTiStep, All.TimeInit, FdCPU); /* produce some CPU usage info */

    double atime = set_global_time(times.Ti_Current);

    while(1) /* main loop */
    {
        /* Find next synchronization point and the timebins active during this timestep.
         *
         * Note that on startup, P[i].TimeBin == 0 for all particles,
         * all bins except the zeroth are inactive and so we return 0 from this function.
         * This ensures we run the force calculation for the first timestep.
         */
        inttime_t Ti_Next = find_next_kick(times.Ti_Current, times.mintimebin);
        inttime_t Ti_Last = times.Ti_Current;

        times.Ti_Current = Ti_Next;

        /*Convert back to floating point time*/
        double newatime = set_global_time(times.Ti_Current);
        if(newatime < atime)
            endrun(1, "Negative timestep: %g New Time: %g Old time %g!\n", newatime - atime, newatime, atime);
        atime = newatime;

        /* Compute the list of particles that cross a lightcone and write it to disc.*/
        if(All.LightconeOn)
            lightcone_compute(atime, PartManager->BoxSize, &All.CP, Ti_Last, Ti_Next);

        int is_PM = is_PM_timestep(&times);

        SyncPoint * next_sync; /* if we are out of planned sync points, terminate */
        SyncPoint * planned_sync; /* NULL; if the step is not a planned sync point. */

        next_sync = find_next_sync_point(times.Ti_Current);
        planned_sync = find_current_sync_point(times.Ti_Current);

        HCIAction action[1];

        hci_action_init(action); /* init to no action */

        int stop = 0;

        if(is_PM) {
            /* query HCI requests only on PM step; where kick and drifts are synced */
            stop = hci_query(HCI_DEFAULT_MANAGER, action);

            if(action->type == HCI_TERMINATE) {
                endrun(0, "Human triggered termination.\n");
            }
        }

        double rel_random_shift[3] = {0};
        if(NumCurrentTiStep > 0 && is_PM  && All.RandomParticleOffset > 0) {
            update_random_offset(rel_random_shift);
        }

        int extradomain = is_timebin_active(times.mintimebin + All.MaxDomainTimeBinDepth, times.Ti_Current);
        /* drift and ddecomp decomposition */
        /* at first step this is a noop */
        if(extradomain || is_PM) {
            /* Sync positions of all particles */
            drift_all_particles(Ti_Last, times.Ti_Current, &All.CP, rel_random_shift);
            /* full decomposition rebuilds the domain, needs keys.*/
            domain_decompose_full(ddecomp);
        } else {
            /* FIXME: add a parameter for ddecomp_decompose_incremental */
            /* currently we drift all particles every step */
            /* If it is not a PM step, do a shorter version
             * of the ddecomp decomp which just exchanges particles.*/
            struct DriftData drift;
            drift.CP = &All.CP;
            drift.ti0 = Ti_Last;
            drift.ti1 = times.Ti_Current;
            domain_maintain(ddecomp, &drift);
        }
        update_lastactive_drift(&times);


        ActiveParticles Act = {0};
        rebuild_activelist(&Act, &times, NumCurrentTiStep, atime);

        set_random_numbers(All.RandomSeed + times.Ti_Current);

        /* Are the particle neutrinos gravitating this timestep?
         * If so we need to add them to the tree.*/
        int HybridNuGrav = All.HybridNeutrinosOn && atime <= All.HybridNuPartTime;

        /* Collective: total number of active particles must be small enough*/
        int pairwisestep = use_pairwise_gravity(&Act, PartManager);

        MyFloat * GradRho = NULL;
        if(sfr_need_to_compute_sph_grad_rho())
            GradRho = (MyFloat *) mymalloc2("SPH_GradRho", sizeof(MyFloat) * 3 * SlotsManager->info[0].size);

        /* Need to rebuild the force tree because all TopLeaves are out of date.*/
        ForceTree Tree = {0};
        force_tree_rebuild(&Tree, ddecomp, HybridNuGrav, !pairwisestep && All.TreeGravOn, All.OutputDir);

        /* density() happens before gravity because it also initializes the predicted variables.
        * This ensures that prediction consistently uses the grav and hydro accel from the
        * timestep before this one, which matches Gadget-2/3. It was tested to make a small difference,
        * since prediction is only really used for artificial viscosity.
        *
        * Doing it first also means the density is up to date for
        * adaptive gravitational softenings. */
        if(GasEnabled)
        {
            /*Allocate the memory for predicted SPH data.*/
            struct sph_pred_data sph_predicted = slots_allocate_sph_pred_data(SlotsManager->info[0].size);

            if(All.DensityOn)
                density(&Act, 1, DensityIndependentSphOn(), All.BlackHoleOn, All.MinEgySpec, times, &All.CP, &sph_predicted, GradRho, &Tree);  /* computes density, and pressure */

            /***** update smoothing lengths in tree *****/
            force_update_hmax(Act.ActiveParticle, Act.NumActiveParticle, &Tree, ddecomp);
            /***** hydro forces *****/
            MPIU_Barrier(MPI_COMM_WORLD);

            /* adds hydrodynamical accelerations  and computes du/dt  */
            if(All.HydroOn)
                hydro_force(&Act, atime, &sph_predicted, All.MinEgySpec, times, &All.CP, &Tree);

            /* Scratch data cannot be used checkpoint because FOF does an exchange.*/
            slots_free_sph_pred_data(&sph_predicted);
        }

        /* The opening criterion for the gravtree
        * uses the *total* gravitational acceleration
        * from the last timestep, GravPM+GravAccel.
        * So we must compute GravAccel for this timestep
        * before gravpm_force() writes the PM acc. for
        * this timestep to GravPM. Note initially both
        * are zero and so the tree is opened maximally
        * on the first timestep.*/
        const int NeutrinoTracer =  All.HybridNeutrinosOn && (atime <= All.HybridNuPartTime);
        const double rho0 = All.CP.Omega0 * 3 * All.CP.Hubble * All.CP.Hubble / (8 * M_PI * All.CP.GravInternal);

        if(All.TreeGravOn) {
            /* Do a short range pairwise only step if desired*/
            if(pairwisestep) {
                struct gravshort_tree_params gtp = get_gravshort_treepar();
                grav_short_pair(&Act, &pm, &Tree, gtp.Rcut, rho0, NeutrinoTracer, All.FastParticleType);
            }
            else
                grav_short_tree(&Act, &pm, &Tree, rho0, NeutrinoTracer, All.FastParticleType);
        }

        /* We use the total gravitational acc.
        * to open the tree and total acc for the timestep.
        * Note that any of (GravAccel, GravPM,
        * HydroAccel) may change much faster than
        * the total acc.
        * We do the same as Gadget-2, but one could
        * instead use short-range tree acc. only
        * for opening angle or short-range timesteps,
        * or include hydro in the opening angle.*/
        if(is_PM)
        {
            gravpm_force(&pm, &Tree, &All.CP, atime, All.UnitLength_in_cm, All.OutputDir, All.MassiveNuLinRespOn, All.TimeIC, All.HybridNeutrinosOn, All.FastParticleType, All.BlackHoleOn);

            /* compute and output energy statistics if desired. */
            if(All.OutputEnergyDebug)
                energy_statistics(FdEnergy, atime, PartManager);
        }

        MPIU_Barrier(MPI_COMM_WORLD);
        message(0, "Forces computed.\n");

        /* Update velocity to Ti_Current; this synchronizes TiKick and TiDrift for the active particles
         * and sets Ti_Kick in the times structure.*/
        if(is_PM) {
            apply_PM_half_kick(&All.CP, &times);
        }

        /* Need a scale factor for entropy and velocity limiters*/
        apply_half_kick(&Act, &All.CP, &times, atime, All.MinEgySpec);

        /* Cooling and extra physics show up as a source term in the evolution equations.
         * Formally you can write the structure of the partial differential equations:
           dU/dt +  div(F) = S
         * where the cooling, BH and SFR are the source term S.
         * The extra physics is done after the kick, using a Strang split operator.
         * Gadget3/Arepo tries to follow the general operator splitting ansatz (often called Strang splitting).
         * Here you alternate the evolution under the operator generating the time evolution of the homogenous system (ie, without S)
         * with the operator generating the time evolution under the source function alone.
         * This means to advance the full system by dt, you first evolve dU/dt = S
         * by dt, and then dU/dt +  div(F) = 0 by dt.
         * [Actually, for second-order convergence in time, you should rather evolve S for dt/2, then the homogenous part for dt, and then S again for dt/2.]

         * The operator-split approach offers a number of practical advantages when the source function is stiff.
         * You can, for example, solve dU/dt = S in a robust and stable fashion with an implict solver, whereas the Gadget2 approach is severely challenged
         * and either requires an artficial clipping of the maximum allowed cooling rate, or a severe reduction of the timestep, otherwise the
         * predicted entropy due to cooling someting during the timestep can become severely wrong. Also, the source term approach can be easily
         * used to treat effectively instantaneous injections of energy (like from BHs), which is again hard to properly incorporate in the
         * time-integration approach where you want to have a "full" dU/dt all times. (Volker Springel 2020).
         */
        if(GasEnabled)
        {
            /* Do this before sfr and bh so the gas hsml always contains DesNumNgb neighbours.*/
            if(All.MetalReturnOn) {
                double AvgGasMass = All.CP.OmegaBaryon * 3 * All.CP.Hubble * All.CP.Hubble / (8 * M_PI * All.CP.GravInternal) * pow(PartManager->BoxSize, 3) / All.NTotalInit[0];
                metal_return(&Act, ddecomp, &All.CP, atime, AvgGasMass);
            }

            /* this will find new black hole seed halos.
             * Note: the FOF code does not know about garbage particles,
             * so ensure we do not have garbage present when we call this.
             * Also a good idea to only run it on a PM step.
             * This does not break the tree because the new black holes do not move or change mass, just type.
             * It does not matter that the velocities are half a step off because they are not used in the FoF code.*/
            if (is_PM && ((All.BlackHoleOn && atime >= TimeNextSeedingCheck) ||
                (during_helium_reionization(1/atime - 1) && need_change_helium_ionization_fraction(atime)))) {

                /* Seeding: builds its own tree.*/
                FOFGroups fof = fof_fof(ddecomp, 0, MPI_COMM_WORLD);
                if(All.BlackHoleOn && atime >= TimeNextSeedingCheck) {
                    fof_seed(&fof, &Act, atime, MPI_COMM_WORLD);
                    TimeNextSeedingCheck = atime * All.TimeBetweenSeedingSearch;
                }

                if(during_helium_reionization(1/atime - 1)) {
                    /* Helium reionization by switching on quasar bubbles*/
                    do_heiii_reionization(atime, &fof, ddecomp, &All.CP, All.UnitEnergy_in_cgs / All.UnitMass_in_g, FdHelium);
                }
                fof_finish(&fof);
            }

            if(is_PM) {
                /*Rebuild the force tree we freed in gravpm to save memory. Means might be two trees during FOF.*/
                force_tree_rebuild(&Tree, ddecomp, HybridNuGrav, 0, All.OutputDir);
            }

            /* Black hole accretion and feedback */
            if(All.BlackHoleOn) {
                struct UnitSystem units = get_unitsystem();
                blackhole(&Act, atime, &All.CP, &Tree, units, FdBlackHoles, FdBlackholeDetails);
            }

            /**** radiative cooling and star formation *****/
            if(All.CoolingOn)
                cooling_and_starformation(&Act, atime, get_dloga_for_bin(times.mintimebin, times.Ti_Current), &Tree, &All.CP, GradRho, FdSfr);

        }
        /* We don't need this timestep's tree anymore.*/
        force_tree_free(&Tree);

        if(GradRho) {
            myfree(GradRho);
            GradRho = NULL;
        }

        /* If a snapshot is requested, write it.
         * write_checkpoint is responsible to maintain a valid ddecomp and tree after it is called.
         *
         * We only attempt to output on sync points. This is the only chance where all variables are
         * synchronized in a consistent state in a K(KDDK)^mK scheme.
         */

        int WriteSnapshot = 0;
        int WriteFOF = 0;

        if(planned_sync) {
            WriteSnapshot |= planned_sync->write_snapshot;
            WriteFOF |= planned_sync->write_fof;
        }

        if(is_PM) { /* the if here is unnecessary but to signify checkpointing occurs only at PM steps. */
            WriteSnapshot |= action->write_snapshot;
            WriteFOF |= action->write_fof;
        }
        if(WriteSnapshot || WriteFOF) {
            /* Get a new snapshot*/
            SnapshotFileCount++;
            /* The accel may have created garbage -- collect them before writing a snapshot.
             * If we do collect, reset active list size.*/
            int compact[6] = {0};
            if(slots_gc(compact, PartManager, SlotsManager))
                Act.NumActiveParticle = PartManager->NumPart;
        }
        FOFGroups fof = {0};
        if(WriteFOF) {
            /* Compute FOF and assign GrNr so it can be written in checkpoint.*/
            fof = fof_fof(ddecomp, 1, MPI_COMM_WORLD);
        }

        /* WriteFOF just reminds the checkpoint code to save GroupID*/
        write_checkpoint(SnapshotFileCount, WriteSnapshot, WriteFOF, atime, All.OutputDir, All.SnapshotFileBase, All.OutputDebugFields);

        /* Save FOF tables after checkpoint so that if there is a FOF save bug we have particle tables available to debug it*/
        if(WriteFOF) {
            fof_save_groups(&fof, All.OutputDir, All.FOFFileBase, SnapshotFileCount, All.PartAllocFactor, &All.CP, atime, All.MassTable, All.StarformationOn, All.BlackHoleOn, MPI_COMM_WORLD);
            fof_finish(&fof);
        }

        write_cpu_log(NumCurrentTiStep, atime, FdCPU);    /* produce some CPU usage info */

        report_memory_usage("RUN");

        if(!next_sync || stop) {
            /* out of sync points, or a requested stop, the run has finally finished! Yay.*/
            if(action->type == HCI_TIMEOUT)
                message(0, "Stopping: not enough time for another PM step before TimeLimitCPU is reached.\n");
            break;
        }

        /* more steps to go. */

        /* assign new timesteps to the active particles,
         * now that we know they have synched TiKick and TiDrift,
         * and advance the PM timestep.*/
        const double asmth = All.Asmth * PartManager->BoxSize / All.Nmesh;
        find_timesteps(&Act, &times, atime, All.FastParticleType, &All.CP, asmth, NumCurrentTiStep == 0, All.OutputDir);

        /* Update velocity and ti_kick to the new step, with the newly computed step size */
        apply_half_kick(&Act, &All.CP, &times, atime, All.MinEgySpec);

        if(is_PM) {
            apply_PM_half_kick(&All.CP, &times);
        }

        /* We can now free the active list: the new step have new active particles*/
        free_activelist(&Act);

        NumCurrentTiStep++;
    }

    close_outputfiles();
}

void write_cpu_log(int NumCurrentTiStep, const double atime, FILE * FdCPU)
{
    walltime_summary(0, MPI_COMM_WORLD);

    if(FdCPU)
    {
        int NTask;
        MPI_Comm_size(MPI_COMM_WORLD, &NTask);
        fprintf(FdCPU, "Step %d, Time: %g, MPIs: %d Threads: %d Elapsed: %g\n", NumCurrentTiStep, atime, NTask, omp_get_max_threads(), Clocks.ElapsedTime);
        walltime_report(FdCPU, 0, MPI_COMM_WORLD);
        fflush(FdCPU);
    }
}

/* We operate in a situation where the particles are in a coordinate frame
 * offset slightly from the ICs (to avoid correlated tree errors).
 * This function updates the global variable containing that offset, and
 * stores the relative shift from the last offset in the rel_random_shift output
 * array. */
static void
update_random_offset(double * rel_random_shift)
{
    int i;
    for (i = 0; i < 3; i++) {
        /* Note random number table is duplicated across processors*/
        double rr = get_random_number(i);
        /* Upstream Gadget uses a random fraction of the box, but since all we need
         * is to adjust the tree openings, and the tree force is zero anyway on the
         * scale of a few PM grid cells, this seems enough.*/
        rr *= All.RandomParticleOffset * PartManager->BoxSize / All.Nmesh;
        /* Subtract the old random shift first.*/
        rel_random_shift[i] = rr - PartManager->CurrentParticleOffset[i];
        PartManager->CurrentParticleOffset[i] = rr;
    }
    message(0, "Internal particle offset is now %g %g %g\n", PartManager->CurrentParticleOffset[0], PartManager->CurrentParticleOffset[1], PartManager->CurrentParticleOffset[2]);
#ifdef DEBUG
    /* Check explicitly that the vector is the same on all processors*/
    double test_random_shift[3] = {0};
    for (i = 0; i < 3; i++)
        test_random_shift[i] = PartManager->CurrentParticleOffset[i];
    MPI_Bcast(test_random_shift, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (i = 0; i < 3; i++)
        if(test_random_shift[i] != PartManager->CurrentParticleOffset[i])
            endrun(44, "Random shift %d is %g != %g on task 0!\n", i, test_random_shift[i], PartManager->CurrentParticleOffset[i]);
#endif
}

/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
static void
open_outputfiles(int RestartSnapNum)
{
    const char mode[3]="a+";
    char * buf;
    char * postfix;
    int ThisTask;
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    FdCPU = NULL;
    FdEnergy = NULL;
    FdBlackHoles = NULL;
    FdSfr = NULL;
    FdBlackholeDetails = NULL;

    if(RestartSnapNum != -1) {
        postfix = fastpm_strdup_printf("-R%03d", RestartSnapNum);
    } else {
        postfix = fastpm_strdup_printf("%s", "");
    }

    /* all the processors write to separate files*/
    if(All.BlackHoleOn && All.WriteBlackHoleDetails){
        buf = fastpm_strdup_printf("%s/%s%s/%06X", All.OutputDir,"BlackholeDetails",postfix,ThisTask);
        fastpm_path_ensure_dirname(buf);
        if(!(FdBlackholeDetails = fopen(buf,"a")))
            endrun(1, "Failed to open blackhole detail %s\n", buf);
        myfree(buf);
    }

    /* only the root processors writes to the log files */
    if(ThisTask != 0) {
        return;
    }

    buf = fastpm_strdup_printf("%s/%s%s", All.OutputDir, All.CpuFile, postfix);
    fastpm_path_ensure_dirname(buf);
    if(!(FdCPU = fopen(buf, mode)))
        endrun(1, "error in opening file '%s'\n", buf);
    myfree(buf);

    if(All.OutputEnergyDebug) {
        buf = fastpm_strdup_printf("%s/%s%s", All.OutputDir, All.EnergyFile, postfix);
        fastpm_path_ensure_dirname(buf);
        if(!(FdEnergy = fopen(buf, mode)))
            endrun(1, "error in opening file '%s'\n", buf);
        myfree(buf);
    }

    if(All.StarformationOn) {
        buf = fastpm_strdup_printf("%s/%s%s", All.OutputDir, "sfr.txt", postfix);
        fastpm_path_ensure_dirname(buf);
        if(!(FdSfr = fopen(buf, mode)))
            endrun(1, "error in opening file '%s'\n", buf);
        myfree(buf);
    }

    if(qso_lightup_on()) {
        buf = fastpm_strdup_printf("%s/%s%s", All.OutputDir, "helium.txt", postfix);
        fastpm_path_ensure_dirname(buf);
        if(!(FdHelium = fopen(buf, mode)))
            endrun(1, "error in opening file '%s'\n", buf);
        myfree(buf);
    }

    if(All.BlackHoleOn) {
        buf = fastpm_strdup_printf("%s/%s%s", All.OutputDir, "blackholes.txt", postfix);
        fastpm_path_ensure_dirname(buf);
        if(!(FdBlackHoles = fopen(buf, mode)))
            endrun(1, "error in opening file '%s'\n", buf);
        myfree(buf);
    }
}


/*!  This function closes the global log-files.
*/
static void
close_outputfiles(void)
{
    if(FdCPU)
        fclose(FdCPU);
    if(FdEnergy)
        fclose(FdEnergy);
    if(FdSfr)
        fclose(FdSfr);
    if(FdBlackHoles)
        fclose(FdBlackHoles);
    if(FdBlackholeDetails)
        fclose(FdBlackholeDetails);
}

/* Construct a unit system struct*/
struct UnitSystem
get_unitsystem(void)
{
    struct UnitSystem units;
    units.UnitMass_in_g = All.UnitMass_in_g;
    units.UnitVelocity_in_cm_per_s = All.UnitVelocity_in_cm_per_s;
    units.UnitLength_in_cm = All.UnitLength_in_cm;

    units.UnitTime_in_s = units.UnitLength_in_cm / units.UnitVelocity_in_cm_per_s;
    units.UnitDensity_in_cgs = units.UnitMass_in_g / pow(units.UnitLength_in_cm, 3);
    units.UnitEnergy_in_cgs = units.UnitMass_in_g * pow(units.UnitLength_in_cm, 2) / pow(units.UnitTime_in_s, 2);
    return units;
}

/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
static void
set_units(void)
{
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;

    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

    /* convert some physical input parameters to internal units */

    All.CP.Hubble = HUBBLE * All.UnitTime_in_s;
    All.CP.UnitTime_in_s = All.UnitTime_in_s;
    All.CP.GravInternal = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);

    init_cosmology(&All.CP, All.TimeIC);
    /* Detect cosmologies that are likely to be typos in the parameter files*/
    if(All.CP.HubbleParam < 0.1 || All.CP.HubbleParam > 10 ||
        All.CP.OmegaLambda < 0 || All.CP.OmegaBaryon < 0 || All.CP.OmegaG < 0 || All.CP.OmegaCDM < 0)
        endrun(5, "Bad cosmology: H0 = %g OL = %g Ob = %g Og = %g Ocdm = %g\n",
               All.CP.HubbleParam, All.CP.OmegaLambda, All.CP.OmegaBaryon, All.CP.OmegaCDM);

    /*Initialise the hybrid neutrinos, after Omega_nu*/
    if(All.HybridNeutrinosOn)
        init_hybrid_nu(&All.CP.ONu.hybnu, All.CP.MNu, All.HybridVcrit, LIGHTCGS/1e5, All.HybridNuPartTime, All.CP.ONu.kBtnu);

    message(0, "Hubble (internal units) = %g\n", All.CP.Hubble);
    message(0, "G (internal units) = %g\n", All.CP.GravInternal);
    message(0, "UnitLength_in_cm = %g \n", All.UnitLength_in_cm);
    message(0, "UnitMass_in_g = %g \n", All.UnitMass_in_g);
    message(0, "UnitTime_in_s = %g \n", All.UnitTime_in_s);
    message(0, "UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
    message(0, "UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
    message(0, "UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
    message(0, "Dark energy model: OmegaL = %g OmegaFLD = %g\n",All.CP.OmegaLambda, All.CP.Omega_fld);
    message(0, "Photon density OmegaG = %g\n",All.CP.OmegaG);
    if(!All.MassiveNuLinRespOn)
        message(0, "Massless Neutrino density OmegaNu0 = %g\n",get_omega_nu(&All.CP.ONu, 1));
    message(0, "Curvature density OmegaK = %g\n",All.CP.OmegaK);
    if(All.CP.RadiationOn) {
        /* note that this value is inaccurate if there is massive neutrino. */
        double OmegaTot = All.CP.OmegaG + All.CP.OmegaK + All.CP.Omega0 + All.CP.OmegaLambda;
        if(!All.MassiveNuLinRespOn)
            OmegaTot += get_omega_nu(&All.CP.ONu, 1);
        message(0, "Radiation is enabled in Hubble(a). "
               "Following CAMB convention: Omega_Tot - 1 = %g\n", OmegaTot - 1);
    }
    message(0, "\n");
}
