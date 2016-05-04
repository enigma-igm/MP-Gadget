#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>


#include "allvars.h"
#include "param.h"
#include "densitykernel.h"
#include "proto.h"
#include "cosmology.h"
#include "cooling.h"
#include "petaio.h"
#include "mymalloc.h"

#include "config.h"

/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
    if(ThisTask == 0)
    {
        /*    printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version()); */
        printf("\nThis is P-Gadget, version %s.\n", GADGETVERSION);
        printf("\nRunning on %d MPIs .\n", NTask);
        printf("\nRunning on %d Threads.\n", omp_get_max_threads());
        printf("\nCode was compiled with settings:\n %s\n", COMPILETIMESETTINGS);
        printf("\nSize of particle structure       %td  [bytes]\n",sizeof(struct particle_data));
        printf("\nSize of blackhole structure       %td  [bytes]\n",sizeof(struct bh_particle_data));
        printf("\nSize of sph particle structure   %td  [bytes]\n",sizeof(struct sph_particle_data));

        if(RestartFlag == 1) {
            fprintf(stderr, "Restarting from restart file is no longer supported. Use a snapshot instead.\n");
            abort();
        }
    }

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
    x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

    read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

    mymalloc_init();
    walltime_init(&All.CT);
    petaio_init();


#ifdef DEBUG
    write_pid_file();
    enable_core_dumps_and_fpu_exceptions();
#endif

    set_units();


#ifdef COOLING
    set_global_time(All.TimeBegin);
    InitCool();
#endif

#if defined(SFR)
    init_clouds();
#endif

#ifdef LIGHTCONE
    lightcone_init();
#endif

    boxSize = All.BoxSize;
    boxHalf = 0.5 * All.BoxSize;
    inverse_boxSize = 1. / boxSize;

    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, 42);	/* start-up seed */

    if(RestartFlag != 3 && RestartFlag != 4)
        long_range_init();

    All.TimeLastRestartFile = 0;


    set_random_numbers();

    init();			/* ... read in initial model */

    open_outputfiles();

    reconstruct_timebins();

#ifdef TWODIMS
    int i;

    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[2] = 0;
        P[i].Vel[2] = 0;

        P[i].GravAccel[2] = 0;

        if(P[i].Type == 0)
        {
            SPHP(i).VelPred[2] = 0;
            SPHP(i).a.HydroAccel[2] = 0;
        }
    }
#endif


    init_drift_table();

    if(RestartFlag == 2)
        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
    else
        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);


    All.TimeLastRestartFile = 0;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
    double meanweight;

    All.UnitVelocity_in_cm_per_s = 1e5; /* 1 km/sec */
    All.UnitLength_in_cm = 3.085678e21; /* 1.0 Kpc /h */
    All.UnitMass_in_g = 1.989e43;       /* 1e10 Msun/h*/

    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);

    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
    All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
    All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

    /* convert some physical input parameters to internal units */

    All.Hubble = HUBBLE * All.UnitTime_in_s;

    if(ThisTask == 0)
    {
        printf("\nHubble (internal units) = %g\n", All.Hubble);
        printf("G (internal units) = %g\n", All.G);
        printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
        printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
        printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
        printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
        printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
        printf("Radiation density Omega_R = %g\n",OMEGAR);

        printf("\n");
    }

    meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

    All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
    All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


#if defined(SFR)
    set_units_sfr();
#endif
}



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
    char mode[2], buf[200];
    char dumpdir[200];
    char postfix[128];

    if(RestartFlag == 0 || RestartFlag == 2)
        strcpy(mode, "w");
    else
        strcpy(mode, "a");

    if(RestartFlag == 2) {
        sprintf(postfix, "-R%03d", RestartSnapNum);
    } else {
        sprintf(postfix, "%s", "");
    }

    /* create spliced dirs */
    int chunk = 10;
    if (NTask > 100) chunk = 100;
    if (NTask > 1000) chunk = 1000;

    sprintf(dumpdir, "%sdumpdir-%d%s/", All.OutputDir, (int)(ThisTask / chunk), postfix);
    mkdir(dumpdir, 02755);

    if(ThisTask != 0)		/* only the root processors writes to the log files */
        return;

    sprintf(buf, "%s%s%s", All.OutputDir, All.CpuFile, postfix);
    if(!(FdCPU = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }

    sprintf(buf, "%s%s%s", All.OutputDir, All.InfoFile, postfix);
    if(!(FdInfo = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }

    sprintf(buf, "%s%s%s", All.OutputDir, All.EnergyFile, postfix);
    if(!(FdEnergy = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }

    sprintf(buf, "%s%s%s", All.OutputDir, All.TimingsFile, postfix);
    if(!(FdTimings = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }

#ifdef SFR
    sprintf(buf, "%s%s%s", All.OutputDir, "sfr.txt", postfix);
    if(!(FdSfr = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
#endif

#ifdef BLACK_HOLES
    sprintf(buf, "%s%s%s", All.OutputDir, "blackholes.txt", postfix);
    if(!(FdBlackHoles = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
#endif



}




/*!  This function closes the global log-files.
*/
void close_outputfiles(void)
{

    if(ThisTask != 0)		/* only the root processors writes to the log files */
        return;

    fclose(FdCPU);
    fclose(FdInfo);
    fclose(FdEnergy);
    fclose(FdTimings);

#ifdef SFR
    fclose(FdSfr);
#endif

#ifdef BLACK_HOLES
    fclose(FdBlackHoles);
#endif

}

#ifdef BLACK_HOLES
static void 
BlackHoleFeedbackMethodAction (ParameterSet * ps, char * name, void * data)
{
    int v = param_get_enum(ps, name);
    if(HAS(v, BH_FEEDBACK_TOPHAT) == HAS(v, BH_FEEDBACK_SPLINE)) {
        printf("error BlackHoleFeedbackMethod contains either tophat or spline, but both\n");
        endrun(0);
    }
    if(HAS(v, BH_FEEDBACK_MASS) ==  HAS(v, BH_FEEDBACK_VOLUME)) {
        printf("error BlackHoleFeedbackMethod contains either volume or mass, but both\n");
        endrun(0);
    }
}
#endif

static void
DensityKernelTypeAction(ParameterSet * ps, char * name, void * data)
{
    int v = param_get_int(ps, name);

    if(v >= density_kernel_type_end()) {
        printf("Error. DensityKernelType can be\n");
        int i;
        for(i = 0; i < density_kernel_type_end(); i++) {
            printf("%d %s\n", i, density_kernel_name(i));
        }
        endrun(111);
    }
}

#ifdef SFR
static void
StarformationCriterionAction(ParameterSet * ps, char * name, void * data)
{
    int v = param_get_enum(ps, name);
    if(!HAS(v, SFR_CRITERION_DENSITY)) {
        printf("error: At least use SFR_CRITERION_DENSITY\n");
        endrun(0);
    }
#if ! defined SPH_GRAD_RHO || ! defined METALS
    if(HAS(v, SFR_CRITERION_MOLECULAR_H2)) {
        printf("error: enable SPH_GRAD_RHO to use h2 criterion in sfr \n");
        endrun(0);
    }
    if(HAS(v, SFR_CRITERION_SELFGRAVITY)) {
        printf("error: enable SPH_GRAD_RHO to use selfgravity in sfr \n");
        endrun(0);
    }
#endif
}
#endif

/*! this fucking function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 *
 *  It needs to be rewriten.
 */

static void
OutputListFilenameAction(ParameterSet * ps, char * name, void * data)
{
    char * fname = param_get_string(ps, name);

    FILE *fd;
    int count, flag;
    char buf[512];

    if(!(fd = fopen(fname, "r")))
    {
        printf("can't read output list in file '%s'\n", fname);
        endrun(111);
    }

    All.OutputListLength = 0;

    while(1)
    {
        if(fgets(buf, 500, fd) != buf)
            break;

        count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

        if(count == 1)
            flag = 1;

        if(count == 1 || count == 2)
        {
            if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
            {
                if(ThisTask == 0)
                    printf("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n",
                            (int) MAXLEN_OUTPUTLIST);
                endrun(13);
            }

            All.OutputListFlag[All.OutputListLength] = flag;
            All.OutputListLength++;
        }
    }

    fclose(fd);

    printf("\nfound %d times in output-list.\n", All.OutputListLength);

    free(fname);
}

static char *
fread_all(char * filename)
{
    FILE * fp = fopen(filename, "r");
    if(!fp){
        printf("Could not open parameter file '%s' for reading\n",filename);
        endrun(1);
    }
    fseek(fp, 0, SEEK_END);
    int size = ftell(fp);
    char * r = malloc(size + 1);
    fseek(fp, 0, SEEK_SET);
    fread(r, 1, size, fp);
    r[size] = 0;
    fclose(fp);
    return r;
}

/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
    if(ThisTask == 0) {
        ParameterSet * ps = parameter_set_new();
#ifdef BLACK_HOLES
        param_set_action(ps, "BlackHoleFeedbackMethod", BlackHoleFeedbackMethodAction, NULL);
#endif
        param_set_action(ps, "DensityKernelType", DensityKernelTypeAction, NULL);
#ifdef SFR
        param_set_action(ps, "StarformationCriterion", StarformationCriterionAction, NULL);
#endif
        param_set_action(ps, "OutputListFilename", OutputListFilenameAction, NULL);

        char * content = fread_all(fname);
        param_parse(ps, content);
        free(content);

        All.NumThreads = omp_get_max_threads();
        All.ICFormat = 1;
        All.SnapFormat = 3;
        All.CompressionLevel = 4;

    /* Start reading the values */
        param_get_string2(ps, "InitCondFile", All.InitCondFile);
        param_get_string2(ps, "OutputDir", All.OutputDir);
        param_get_string2(ps, "TreeCoolFile", All.TreeCoolFile);
        param_get_string2(ps, "MetalCoolFile", All.MetalCoolFile);
        param_get_string2(ps, "UVFluctuationfile", All.UVFluctuationFile);
        param_get_string2(ps, "SnapshotFileBase", All.SnapshotFileBase);
        param_get_string2(ps, "EnergyFile", All.EnergyFile);
        param_get_string2(ps, "CpuFile", All.CpuFile);
        param_get_string2(ps, "InfoFile", All.InfoFile);
        param_get_string2(ps, "TimingsFile", All.TimingsFile);
        param_get_string2(ps, "RestartFile", All.RestartFile);
        param_get_string2(ps, "OutputListFilename", All.OutputListFilename);

        All.DensityKernelType = param_get_int(ps, "DensityKernelType");

        All.Omega0 = param_get_double(ps, "Omega0");
        All.OmegaBaryon = param_get_double(ps, "OmegaBaryon");
        All.OmegaLambda = param_get_double(ps, "OmegaLambda");
        All.HubbleParam = param_get_double(ps, "HubbleParam");
        All.BoxSize = param_get_double(ps, "BoxSize");

        All.MaxMemSizePerCore = param_get_int(ps, "MaxMemSizePerCore");
        All.TimeOfFirstSnapshot = param_get_double(ps, "TimeOfFirstSnapshot");
        All.CpuTimeBetRestartFile = param_get_double(ps, "CpuTimeBetRestartFile");
        All.TimeBetStatistics = param_get_double(ps, "TimeBetStatistics");
        All.TimeBegin = param_get_double(ps, "TimeBegin");
        All.TimeMax = param_get_double(ps, "TimeMax");
        All.TreeDomainUpdateFrequency = param_get_double(ps, "TreeDomainUpdateFrequency");
        All.ErrTolTheta = param_get_double(ps, "ErrTolTheta");
        All.ErrTolIntAccuracy = param_get_double(ps, "ErrTolIntAccuracy");
        All.ErrTolForceAcc = param_get_double(ps, "ErrTolForceAcc");
        All.Nmesh = param_get_int(ps, "Nmesh");

        All.MinGasHsmlFractional = param_get_double(ps, "MinGasHsmlFractional");
        All.MaxGasVel = param_get_double(ps, "MaxGasVel");
        All.MaxSizeTimestep = param_get_double(ps, "MaxSizeTimestep");

        All.MinSizeTimestep = param_get_double(ps, "MinSizeTimestep");
        All.MaxRMSDisplacementFac = param_get_double(ps, "MaxRMSDisplacementFac");
        All.ArtBulkViscConst = param_get_double(ps, "ArtBulkViscConst");
        All.CourantFac = param_get_double(ps, "CourantFac");
        All.DensityResolutionEta = param_get_double(ps, "DensityResolutionEta");

        All.DensityContrastLimit = param_get_double(ps, "DensityContrastLimit");
        All.MaxNumNgbDeviation = param_get_double(ps, "MaxNumNgbDeviation");

        All.NumFilesPerSnapshot = param_get_int(ps, "NumFilesPerSnapshot");
        All.NumWritersPerSnapshot = param_get_int(ps, "NumWritersPerSnapshot");
        All.NumFilesPerPIG = param_get_int(ps, "NumFilesPerPIG");
        All.NumWritersPerPIG = param_get_int(ps, "NumWritersPerPIG");

        All.CoolingOn = param_get_int(ps, "CoolingOn");
        All.StarformationOn = param_get_int(ps, "StarformationOn");
        All.TypeOfTimestepCriterion = param_get_int(ps, "TypeOfTimestepCriterion");
        All.TypeOfOpeningCriterion = param_get_int(ps, "TypeOfOpeningCriterion");
        All.TimeLimitCPU = param_get_double(ps, "TimeLimitCPU");
        All.SofteningHalo = param_get_double(ps, "SofteningHalo");
        All.SofteningDisk = param_get_double(ps, "SofteningDisk");
        All.SofteningBulge = param_get_double(ps, "SofteningBulge");
        All.SofteningGas = param_get_double(ps, "SofteningGas");
        All.SofteningStars = param_get_double(ps, "SofteningStars");
        All.SofteningBndry = param_get_double(ps, "SofteningBndry");
        All.SofteningHaloMaxPhys= param_get_double(ps, "SofteningHaloMaxPhys");
        All.SofteningDiskMaxPhys= param_get_double(ps, "SofteningDiskMaxPhys");
        All.SofteningBulgeMaxPhys= param_get_double(ps, "SofteningBulgeMaxPhys");
        All.SofteningGasMaxPhys= param_get_double(ps, "SofteningGasMaxPhys");
        All.SofteningStarsMaxPhys= param_get_double(ps, "SofteningStarsMaxPhys");
        All.SofteningBndryMaxPhys= param_get_double(ps, "SofteningBndryMaxPhys");

        All.BufferSize = param_get_double(ps, "BufferSize");
        All.PartAllocFactor = param_get_double(ps, "PartAllocFactor");
        All.TopNodeAllocFactor = param_get_double(ps, "TopNodeAllocFactor");

        All.InitGasTemp = param_get_double(ps, "InitGasTemp");
        All.MinGasTemp = param_get_double(ps, "MinGasTemp");

    #if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
        All.ReferenceGasMass = param_get_double(ps, "ReferenceGasMass");
    #endif

    #ifdef FOF
        All.FOFHaloLinkingLength = param_get_double(ps, "FOFHaloLinkingLength");
        All.FOFHaloMinLength = param_get_int(ps, "FOFHaloMinLength");
    #endif

    #ifdef BLACK_HOLES
        All.BlackHoleSoundSpeedFromPressure = 0;

        All.TimeBetBlackHoleSearch = param_get_double(ps, "TimeBetBlackHoleSearch");
        All.BlackHoleAccretionFactor = param_get_double(ps, "BlackHoleAccretionFactor");
        All.BlackHoleEddingtonFactor = param_get_double(ps, "BlackHoleEddingtonFactor");
        All.SeedBlackHoleMass = param_get_double(ps, "SeedBlackHoleMass");
        All.MinFoFMassForNewSeed = param_get_double(ps, "MinFoFMassForNewSeed");

        All.BlackHoleNgbFactor = param_get_double(ps, "BlackHoleNgbFactor");

        All.BlackHoleMaxAccretionRadius = param_get_double(ps, "BlackHoleMaxAccretionRadius");
        All.BlackHoleFeedbackFactor = param_get_double(ps, "BlackHoleFeedbackFactor");
        All.BlackHoleFeedbackRadius = param_get_double(ps, "BlackHoleFeedbackRadius");

        All.BlackHoleFeedbackRadiusMaxPhys = param_get_double(ps, "BlackHoleFeedbackRadiusMaxPhys");

        All.BlackHoleFeedbackMethod = param_get_enum(ps, "BlackHoleFeedbackMethod");

    #endif

    #ifdef SFR
        All.StarformationCriterion = param_get_enum(ps, "StarformationCriterion");
        All.CritOverDensity = param_get_double(ps, "CritOverDensity");
        All.CritPhysDensity = param_get_double(ps, "CritPhysDensity");

        All.FactorSN = param_get_double(ps, "FactorSN");
        All.FactorEVP = param_get_double(ps, "FactorEVP");
        All.TempSupernova = param_get_double(ps, "TempSupernova");
        All.TempClouds = param_get_double(ps, "TempClouds");
        All.MaxSfrTimescale = param_get_double(ps, "MaxSfrTimescale");
        All.WindModel = param_get_enum(ps, "WindModel");

        /* The following two are for VS08 and SH03*/
        All.WindEfficiency = param_get_double(ps, "WindEfficiency");
        All.WindEnergyFraction = param_get_double(ps, "WindEnergyFraction");

        /* The following two are for OFJT10*/
        All.WindSigma0 = param_get_double(ps, "WindSigma0");
        All.WindSpeedFactor = param_get_double(ps, "WindSpeedFactor");

        All.WindFreeTravelLength = param_get_double(ps, "WindFreeTravelLength");
        All.WindFreeTravelDensFac = param_get_double(ps, "WindFreeTravelDensFac");

        All.QuickLymanAlphaProbability = param_get_double(ps, "QuickLymanAlphaProbability");

    #endif

    #ifdef SOFTEREQS
        All.FactorForSofterEQS = param_get_double(ps, "FactorForSofterEQS");
    #endif

        printf("The Density Kernel type is %s\n", density_kernel_name(All.DensityKernelType));
        All.DesNumNgb = density_kernel_desnumngb(All.DensityKernelType,
                All.DensityResolutionEta);
        printf("The Density resolution is %g * mean separation, or %d neighbours\n",
                All.DensityResolutionEta, All.DesNumNgb);

        parameter_set_free(ps);
    }

    MPI_Bcast(&All, sizeof(All), MPI_BYTE, 0, MPI_COMM_WORLD);

    if(All.TypeOfTimestepCriterion >= 3)
    {
        if(ThisTask == 0)
        {
            printf("The specified timestep criterion\n");
            printf("is not valid\n");
        }
        endrun(0);
    }

#ifdef SFR

    if(All.StarformationOn == 0)
    {
        if(ThisTask == 0)
        {
            printf("StarformationOn is disabled!\n");
        }
    }
    if(All.CoolingOn == 0)
    {
        if(ThisTask == 0)
        {
            printf("You try to use the code with star formation enabled,\n");
            printf("but you did not switch on cooling.\nThis mode is not supported.\n");
        }
        endrun(0);
    }
#else
    if(All.StarformationOn == 1)
    {
        if(ThisTask == 0)
        {
            printf("Code was compiled with star formation switched off.\n");
            printf("You must set `StarformationOn=0', or recompile the code.\n");
        }
        All.StarformationOn = 0;
    }
#endif





    if(All.NumWritersPerSnapshot > NTask)
    {
       All.NumWritersPerSnapshot = NTask;
    }
    if(All.NumWritersPerPIG > NTask)
    {
       All.NumWritersPerPIG = NTask;
    }


#ifdef METALS
#ifndef SFR
    if(ThisTask == 0)
    {
        printf("Code was compiled with METALS, but not with SFR.\n");
        printf("This is not allowed.\n");
    }
    endrun(0);
#endif
#endif

}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
    int i;
    int64_t ti_end;

    if(sizeof(int64_t) != 8)
    {
        if(ThisTask == 0)
            printf("\nType 'int64_t' is not 64 bit on this platform\n\n");
        endrun(555);
    }

    if(ThisTask == 0)
    {
        printf("\nAll.TimeMax has been changed in the parameterfile\n");
        printf("Need to adjust integer timeline\n\n\n");
    }

    if(TimeMax_new < TimeMax_old)
    {
        if(ThisTask == 0)
            printf("\nIt is not allowed to reduce All.TimeMax\n\n");
        endrun(556);
    }

    ti_end = (int64_t) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);

    while(ti_end > TIMEBASE)
    {
        All.Timebase_interval *= 2.0;

        ti_end /= 2;
        All.Ti_Current /= 2;

        All.PM_Ti_begstep /= 2;
        All.PM_Ti_endstep /= 2;

        for(i = 0; i < NumPart; i++)
        {
            P[i].Ti_begstep /= 2;
            P[i].Ti_current /= 2;

            if(P[i].TimeBin > 0)
            {
                P[i].TimeBin--;
                if(P[i].TimeBin <= 0)
                {
                    printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
                    endrun(8765);
                }
            }
        }

        All.Ti_nextlineofsight /= 2;
    }

    All.TimeMax = TimeMax_new;
}
