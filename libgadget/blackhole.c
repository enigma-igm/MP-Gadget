#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <omp.h>

#include "physconst.h"
#include "cooling.h"
#include "gravity.h"
#include "densitykernel.h"
#include "treewalk.h"
#include "slotsmanager.h"
#include "blackhole.h"
#include "timestep.h"
#include "hydra.h"
#include "density.h"
#include "sfr_eff.h"
#include "winds.h"
#include "walltime.h"
/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

struct BlackholeParams
{
    double BlackHoleAccretionFactor;	/*!< Fraction of BH bondi accretion rate */
    double BlackHoleFeedbackFactor;	/*!< Fraction of the black luminosity feed into thermal feedback */
    enum BlackHoleFeedbackMethod BlackHoleFeedbackMethod;	/*!< method of the feedback*/
    double BlackHoleEddingtonFactor;	/*! Factor above Eddington */
    int BlackHoleRepositionEnabled; /* If true, enable repositioning the BH to the potential minimum*/
    
    int BlackHoleKineticOn; /*If 1, perform AGN kinetic feedback when the Eddington accretion rate is low */
    double BHKE_EddingtonThrFactor; /*Threshold of the Eddington rate for the kinetic feedback*/
    double BHKE_EddingtonMFactor; /* Factor for mbh-dependent Eddington threshold */
    double BHKE_EddingtonMPivot; /* Pivot MBH for mbh-dependent Eddington threshold */
    double BHKE_EddingtonMIndex; /* Powlaw index for mbh-dependent Eddington threshold */
    double BHKE_EffRhoFactor; /* Minimum kinetic feedback efficiency factor scales with BH density*/
    double BHKE_EffCap; /* Cap of the kinetic feedback efficiency factor */
    double BHKE_InjEnergyThr; /*Minimum injection of KineticFeedbackEnergy, controls the burstiness of kinetic feedback*/
    double BHKE_SfrCritOverDensity; /*for KE efficiency calculation, borrow from sfr.params */
    /**********************************************************************/
    int MergeGravBound; /*if 1, apply gravitational bound criteria for BH mergers */

    int BH_DynFrictionMethod;/*0 for off; 1 for Star Only; 2 for DM+Star; 3 for DM+Star+Gas */
    int BH_DFBoostFactor; /*Optional boost factor for DF */
    double BH_DFbmax; /* the maximum impact range, in physical unit of kpc. */
    int BH_DRAG; /*Hydro drag force*/

    double SeedBHDynMass; /* The initial dynamic mass of BH particle */

    double SeedBlackHoleMass;	/*!< (minimum) Seed black hole mass */
    double MaxSeedBlackHoleMass; /* Maximum black hole seed mass*/
    double SeedBlackHoleMassIndex; /* Power law index for BH seed mass*/
    /************************************************************************/
} blackhole_params;

int
BHGetRepositionEnabled(void)
{
    return blackhole_params.BlackHoleRepositionEnabled;
}

typedef struct {
    TreeWalkQueryBase base;
    MyFloat Density;
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat BH_Mass;
    MyFloat Vel[3];
    MyFloat Accel[3];
    MyIDType ID;
    MyFloat Mtrack;
} TreeWalkQueryBHAccretion;

typedef struct {
    TreeWalkResultBase base;
    MyFloat BH_MinPotPos[3];
    MyFloat BH_MinPotVel[3];
    MyFloat BH_MinPot;
    int BH_minTimeBin;
    int encounter;
    MyFloat FeedbackWeightSum;

    MyFloat SmoothedEntropy;
    MyFloat GasVel[3];
    /* used for AGN kinetic feedback */
    MyFloat V2sumDM;
    MyFloat V1sumDM[3];
    MyFloat NumDM;
    MyFloat MgasEnc;
} TreeWalkResultBHAccretion;

typedef struct {
    TreeWalkNgbIterBase base;
    DensityKernel accretion_kernel;
    DensityKernel feedback_kernel;
} TreeWalkNgbIterBHAccretion;


/*****************************************************************************/
typedef struct {
    TreeWalkQueryBase base;
    MyFloat Hsml;
} TreeWalkQueryBHDynfric;

typedef struct {
    TreeWalkResultBase base;

    MyFloat SurroundingVel[3];
    MyFloat SurroundingDensity;
    int SurroundingParticles;
    MyFloat SurroundingRmsVel;

} TreeWalkResultBHDynfric;

typedef struct {
    TreeWalkNgbIterBase base;
    DensityKernel dynfric_kernel;
} TreeWalkNgbIterBHDynfric;

/*****************************************************************************/


typedef struct {
    TreeWalkQueryBase base;
    MyFloat Hsml;
    MyFloat Mtrack;
    MyFloat BH_Mass;
    MyIDType ID;
    MyFloat Density;
    MyFloat FeedbackEnergy;
    MyFloat FeedbackWeightSum;
    MyFloat KEFeedbackEnergy;
    int FdbkChannel; /* 0 thermal, 1 kinetic */
} TreeWalkQueryBHFeedback;

typedef struct {
    TreeWalkResultBase base;
    MyFloat Mass; /* the accreted Mdyn */
    MyFloat AccretedMomentum[3];
    MyFloat BH_Mass;
    int BH_CountProgs;
    MyFloat acMtrack; /* the accreted Mtrack */
    int alignment; /* Ensure alignment*/
} TreeWalkResultBHFeedback;

typedef struct {
    TreeWalkNgbIterBase base;
    DensityKernel feedback_kernel;
} TreeWalkNgbIterBHFeedback;

struct BHPriv {
    /* Temporary array to store the IDs of the swallowing black hole for gas.
     * We store ID + 1 so that SwallowID == 0 can correspond to the unswallowed case. */
    MyIDType * SPH_SwallowID;
    /* These are temporaries used in the accretion treewalk*/
    MyFloat * MinPot;
    MyFloat * BH_Entropy;
    MyFloat (*BH_SurroundingGasVel)[3];

    /*************************************************************************/
    /* used in the dynamic friction treewalk*/
    MyFloat * BH_SurroundingDensity;
    int * BH_SurroundingParticles;
    MyFloat (*BH_SurroundingVel)[3];
    MyFloat * BH_SurroundingRmsVel;

    /*************************************************************************/

    MyFloat (*BH_accreted_momentum)[3];

    /* These are temporaries used in the feedback treewalk.*/
    MyFloat * BH_accreted_Mass;
    MyFloat * BH_accreted_BHMass;
    MyFloat * BH_accreted_Mtrack;

    /* This is a temporary computed in the accretion treewalk and used
     * in the feedback treewalk*/
    MyFloat * BH_FeedbackWeightSum;
    
    /* temporary computed for kinetic feedback energy threshold*/
    MyFloat * NumDM;
    MyFloat (*V1sumDM)[3];
    MyFloat * V2sumDM;
    MyFloat * MgasEnc;
    /* mark the state of AGN kinetic feedback, 1 accumulate, 2 release */
    int * KEflag;

    /* Time factors*/
    double atime;
    double a3inv;
    double hubble;
    struct UnitSystem units;
    Cosmology * CP;
    /* Counters*/
    int64_t * N_sph_swallowed;
    int64_t * N_BH_swallowed;
};
#define BH_GET_PRIV(tw) ((struct BHPriv *) (tw->priv))

/* Structure needs to be packed to ensure disc write is the same on all architectures and the record size is correct. */
struct __attribute__((__packed__)) BHinfo{
    /* Stores sizeof(struct BHinfo) - 2 * sizeof(size1) . Allows size of record to be stored in the struct.*/
    int size1;
    MyIDType ID;
    MyFloat Mass;
    MyFloat Mdot;
    MyFloat Density;
    int minTimeBin;
    int encounter;

    double  MinPotPos[3];
    MyFloat MinPot;
    MyFloat BH_Entropy;
    MyFloat BH_SurroundingGasVel[3];
    MyFloat BH_accreted_momentum[3];

    MyFloat BH_accreted_Mass;
    MyFloat BH_accreted_BHMass;
    MyFloat BH_FeedbackWeightSum;

    MyIDType SPH_SwallowID;
    MyIDType SwallowID;

    int CountProgs;
    int Swallowed;

    /****************************************/
    double Pos[3];
    MyFloat BH_SurroundingDensity;
    MyFloat BH_SurroundingParticles;
    MyFloat BH_SurroundingVel[3];
    MyFloat BH_SurroundingRmsVel;

    double BH_DFAccel[3];
    double BH_DragAccel[3];
    double BH_GravAccel[3];
    double Velocity[3];
    double Mtrack;
    double Mdyn;
    
    double KineticFdbkEnergy;
    double NumDM;
    double V1sumDM[3];
    double V2sumDM;
    double MgasEnc;
    int KEflag;

    MyDouble a;
    /* See size1 above*/
    int size2;
};

/*Set the parameters of the BH module*/
void set_blackhole_params(ParameterSet * ps)
{
    int ThisTask;
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    if(ThisTask == 0) {
        blackhole_params.BlackHoleAccretionFactor = param_get_double(ps, "BlackHoleAccretionFactor");
        blackhole_params.BlackHoleEddingtonFactor = param_get_double(ps, "BlackHoleEddingtonFactor");

        blackhole_params.BlackHoleFeedbackFactor = param_get_double(ps, "BlackHoleFeedbackFactor");

        blackhole_params.BlackHoleFeedbackMethod = (enum BlackHoleFeedbackMethod) param_get_enum(ps, "BlackHoleFeedbackMethod");
        blackhole_params.BlackHoleRepositionEnabled = param_get_int(ps, "BlackHoleRepositionEnabled");
        
        blackhole_params.BlackHoleKineticOn = param_get_int(ps,"BlackHoleKineticOn");
        blackhole_params.BHKE_EddingtonThrFactor = param_get_double(ps, "BHKE_EddingtonThrFactor");
        blackhole_params.BHKE_EddingtonMFactor = param_get_double(ps, "BHKE_EddingtonMFactor");
        blackhole_params.BHKE_EddingtonMPivot = param_get_double(ps, "BHKE_EddingtonMPivot");
        blackhole_params.BHKE_EddingtonMIndex = param_get_double(ps, "BHKE_EddingtonMIndex");
        blackhole_params.BHKE_EffRhoFactor = param_get_double(ps, "BHKE_EffRhoFactor");
        blackhole_params.BHKE_EffCap = param_get_double(ps, "BHKE_EffCap");
        blackhole_params.BHKE_InjEnergyThr = param_get_double(ps, "BHKE_InjEnergyThr");
        blackhole_params.BHKE_SfrCritOverDensity = param_get_double(ps, "CritOverDensity");
        /***********************************************************************************/
        blackhole_params.BH_DynFrictionMethod = param_get_int(ps, "BH_DynFrictionMethod");
        blackhole_params.BH_DFBoostFactor = param_get_int(ps, "BH_DFBoostFactor");
        blackhole_params.BH_DFbmax = param_get_double(ps, "BH_DFbmax");
        blackhole_params.BH_DRAG = param_get_int(ps, "BH_DRAG");
        blackhole_params.MergeGravBound = param_get_int(ps, "MergeGravBound");
        blackhole_params.SeedBHDynMass = param_get_double(ps,"SeedBHDynMass");

        blackhole_params.SeedBlackHoleMass = param_get_double(ps, "SeedBlackHoleMass");
        blackhole_params.MaxSeedBlackHoleMass = param_get_double(ps,"MaxSeedBlackHoleMass");
        blackhole_params.SeedBlackHoleMassIndex = param_get_double(ps,"SeedBlackHoleMassIndex");
        /***********************************************************************************/
    }
    MPI_Bcast(&blackhole_params, sizeof(struct BlackholeParams), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/* accretion routines */
static void
blackhole_accretion_postprocess(int n, TreeWalk * tw);

static void
blackhole_accretion_reduce(int place, TreeWalkResultBHAccretion * remote, enum TreeWalkReduceMode mode, TreeWalk * tw);

static void
blackhole_accretion_copy(int place, TreeWalkQueryBHAccretion * I, TreeWalk * tw);

/* Initializes the minimum potentials*/
static void
blackhole_accretion_preprocess(int n, TreeWalk * tw);

static void
blackhole_accretion_ngbiter(TreeWalkQueryBHAccretion * I,
        TreeWalkResultBHAccretion * O,
        TreeWalkNgbIterBHAccretion * iter,
        LocalTreeWalk * lv);



/*************************************************************************************/
/* DF routines */
static void
blackhole_dynfric_postprocess(int n, TreeWalk * tw);

static int
blackhole_dynfric_haswork(int n, TreeWalk * tw);

static void
blackhole_dynfric_reduce(int place, TreeWalkResultBHDynfric * remote, enum TreeWalkReduceMode mode, TreeWalk * tw);

static void
blackhole_dynfric_copy(int place, TreeWalkQueryBHDynfric * I, TreeWalk * tw);

static void
blackhole_dynfric_ngbiter(TreeWalkQueryBHDynfric * I,
        TreeWalkResultBHDynfric * O,
        TreeWalkNgbIterBHDynfric * iter,
        LocalTreeWalk * lv);

/*************************************************************************************/



/* feedback routines */

static void
blackhole_feedback_postprocess(int n, TreeWalk * tw);

static int
blackhole_feedback_haswork(int n, TreeWalk * tw);

static void
blackhole_feedback_reduce(int place, TreeWalkResultBHFeedback * remote, enum TreeWalkReduceMode mode, TreeWalk * tw);

static void
blackhole_feedback_copy(int place, TreeWalkQueryBHFeedback * I, TreeWalk * tw);

static void
blackhole_feedback_ngbiter(TreeWalkQueryBHFeedback * I,
        TreeWalkResultBHFeedback * O,
        TreeWalkNgbIterBHFeedback * iter,
        LocalTreeWalk * lv);

#define BHPOTVALUEINIT 1.0e29

static double blackhole_soundspeed(double entropy, double rho, const double atime) {
    /* rho is comoving !*/
    double cs = sqrt(GAMMA * entropy * pow(rho, GAMMA_MINUS1));

    cs *= pow(atime, -1.5 * GAMMA_MINUS1);

    return cs;
}

/* Adds the injected black hole energy to an internal energy and caps it at a maximum temperature*/
static double
add_injected_BH_energy(double unew, double injected_BH_energy, double mass, double uu_in_cgs)
{
    unew += injected_BH_energy / mass;
    const double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * uu_in_cgs;
    /* Cap temperature*/
    if(unew > 5.0e8 / u_to_temp_fac)
        unew = 5.0e8 / u_to_temp_fac;

    return unew;
}

static int
get_random_dir(int i, double dir[3])
{
    double theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
    double phi = 2 * M_PI * get_random_number(P[i].ID + 4);

    dir[0] = sin(theta) * cos(phi);
    dir[1] = sin(theta) * sin(phi);
    dir[2] = cos(theta);
    return 0;
}

/* check if two BHs are gravitationally bounded, input dv, da, dx in code unit */
/* same as Bellovary2011, Tremmel2017 */
static int
check_grav_bound(double dx[3], double dv[3], double da[3], const double atime)
{
    int j;
    double KE = 0;
    double PE = 0;

    for(j = 0; j < 3; j++){
        KE += 0.5 * pow(dv[j], 2);
        PE += da[j] * dx[j];
    }

    KE /= (atime * atime); /* convert to proper velocity */
    PE /= atime; /* convert to proper unit */

    /* The gravitationally bound condition is PE + KE < 0.
     * Still merge if it is marginally bound so that we merge
     * particles at zero distance and velocity from each other.*/
    return (PE + KE <= 0);
}



static void
collect_BH_info(int * ActiveBlackHoles, int NumActiveBlackHoles, struct BHPriv *priv, FILE * FdBlackholeDetails)
{
    int i;

    struct BHinfo * infos = (struct BHinfo *) mymalloc("BHDetailCache", NumActiveBlackHoles * sizeof(struct BHinfo));
    report_memory_usage("BLACKHOLE");

    const int size = sizeof(struct BHinfo) - sizeof(infos[0].size1) - sizeof(infos[0].size2);

    #pragma omp parallel for
    for(i = 0; i < NumActiveBlackHoles; i++)
    {
        const int p_i = ActiveBlackHoles ? ActiveBlackHoles[i] : i;
        int PI = P[p_i].PI;

        struct BHinfo *info = &infos[i];
        info->size1 = size;
        info->size2 = size;
        info->ID = P[p_i].ID;
        info->Mass = BHP(p_i).Mass;
        info->Mdot = BHP(p_i).Mdot;
        info->Density = BHP(p_i).Density;
        info->minTimeBin = BHP(p_i).minTimeBin;
        info->encounter = BHP(p_i).encounter;

        if(priv->MinPot) {
            info->MinPot = priv->MinPot[PI];
        }
        info->BH_Entropy = priv->BH_Entropy[PI];
        int k;
        for(k=0; k < 3; k++) {
            info->MinPotPos[k] = BHP(p_i).MinPotPos[k] - PartManager->CurrentParticleOffset[k];
            info->BH_SurroundingGasVel[k] = priv->BH_SurroundingGasVel[PI][k];
            info->BH_accreted_momentum[k] = priv->BH_accreted_momentum[PI][k];
            info->BH_DragAccel[k] = BHP(p_i).DragAccel[k];
            info->BH_GravAccel[k] = P[p_i].GravAccel[k];
            info->Pos[k] = P[p_i].Pos[k] - PartManager->CurrentParticleOffset[k];
            info->Velocity[k] = P[p_i].Vel[k];
            info->BH_DFAccel[k] = BHP(p_i).DFAccel[k];
        }

        /****************************************************************************/
        /* Output some DF info for debugging */
        info->BH_SurroundingDensity = priv->BH_SurroundingDensity[PI];
        info->BH_SurroundingRmsVel = priv->BH_SurroundingRmsVel[PI];
        info->BH_SurroundingParticles = priv->BH_SurroundingParticles[PI];
        info->BH_SurroundingVel[0] = priv->BH_SurroundingVel[PI][0];
        info->BH_SurroundingVel[1] = priv->BH_SurroundingVel[PI][1];
        info->BH_SurroundingVel[2] = priv->BH_SurroundingVel[PI][2];

        /****************************************************************************/
        info->BH_accreted_BHMass = priv->BH_accreted_BHMass[PI];
        info->BH_accreted_Mass = priv->BH_accreted_Mass[PI];
        info->BH_FeedbackWeightSum = priv->BH_FeedbackWeightSum[PI];

        info->SPH_SwallowID = priv->SPH_SwallowID[PI];
        info->SwallowID =  BHP(p_i).SwallowID;
        info->CountProgs = BHP(p_i).CountProgs;
        info->Swallowed =  P[p_i].Swallowed;
        /************************************************************************************************/
        /* When SeedBHDynMass is larger than gas particle mass, we have three mass tracer of blackhole. */
        /* BHP(p_i).Mass : intrinsic mass of BH, accreted every (active) time step.                     */
        /* P[p_i].Mass :  Dynamic mass of BH, used for gravitational interaction.                       */
        /*                Starts to accrete gas particle when BHP(p_i).Mass > SeedBHDynMass             */
        /* BHP(p_i).Mtrack: Initialized as gas particle mass, and is capped at SeedBHDynMass,           */
        /*                 it traces BHP(p_i).Mass by swallowing gas when BHP(p_i).Mass < SeedBHDynMass */
        /************************************************************************************************/
        info->Mtrack = BHP(p_i).Mtrack;
        info->Mdyn = P[p_i].Mass;
        
        info->KineticFdbkEnergy = BHP(p_i).KineticFdbkEnergy;
        info->NumDM = priv->NumDM[PI];
        info->V1sumDM[0] = priv->V1sumDM[PI][0];
        info->V1sumDM[1] = priv->V1sumDM[PI][1];
        info->V1sumDM[2] = priv->V1sumDM[PI][2];
        info->V2sumDM = priv->V2sumDM[PI];
        info->MgasEnc = priv->MgasEnc[PI];
        info->KEflag = priv->KEflag[PI];
        
        info->a = priv->atime;
    }

    fwrite(infos,sizeof(struct BHinfo),NumActiveBlackHoles,FdBlackholeDetails);
    fflush(FdBlackholeDetails);
    myfree(infos);
    int64_t totalN;

    sumup_large_ints(1, &NumActiveBlackHoles, &totalN);
    message(0, "Written details of %ld blackholes in %lu bytes each.\n", totalN, sizeof(struct BHinfo));
}


void
blackhole(const ActiveParticles * act, double atime, Cosmology * CP, ForceTree * tree, const struct UnitSystem units, FILE * FdBlackHoles, FILE * FdBlackholeDetails)
{
    /* Do nothing if no black holes*/
    int64_t totbh;
    MPI_Allreduce(&SlotsManager->info[5].size, &totbh, 1, MPI_INT64, MPI_SUM, MPI_COMM_WORLD);
    if(totbh == 0)
        return;
    int i;

    walltime_measure("/Misc");
    struct BHPriv priv[1] = {0};
    priv->units = units;

    /*************************************************************************/
    TreeWalk tw_dynfric[1] = {{0}};
    tw_dynfric->ev_label = "BH_DYNFRIC";
    tw_dynfric->visit = (TreeWalkVisitFunction) treewalk_visit_ngbiter;
    tw_dynfric->ngbiter_type_elsize = sizeof(TreeWalkNgbIterBHDynfric);
    tw_dynfric->ngbiter = (TreeWalkNgbIterFunction) blackhole_dynfric_ngbiter;
    tw_dynfric->haswork = blackhole_dynfric_haswork;
    tw_dynfric->postprocess = (TreeWalkProcessFunction) blackhole_dynfric_postprocess;
    tw_dynfric->fill = (TreeWalkFillQueryFunction) blackhole_dynfric_copy;
    tw_dynfric->reduce = (TreeWalkReduceResultFunction) blackhole_dynfric_reduce;
    tw_dynfric->query_type_elsize = sizeof(TreeWalkQueryBHDynfric);
    tw_dynfric->result_type_elsize = sizeof(TreeWalkResultBHDynfric);
    tw_dynfric->tree = tree;
    tw_dynfric->priv = priv;

    /*************************************************************************/
    TreeWalk tw_accretion[1] = {{0}};

    tw_accretion->ev_label = "BH_ACCRETION";
    tw_accretion->visit = (TreeWalkVisitFunction) treewalk_visit_ngbiter;
    tw_accretion->ngbiter_type_elsize = sizeof(TreeWalkNgbIterBHAccretion);
    tw_accretion->ngbiter = (TreeWalkNgbIterFunction) blackhole_accretion_ngbiter;
    tw_accretion->haswork = blackhole_dynfric_haswork;
    tw_accretion->postprocess = (TreeWalkProcessFunction) blackhole_accretion_postprocess;
    tw_accretion->preprocess = (TreeWalkProcessFunction) blackhole_accretion_preprocess;
    tw_accretion->fill = (TreeWalkFillQueryFunction) blackhole_accretion_copy;
    tw_accretion->reduce = (TreeWalkReduceResultFunction) blackhole_accretion_reduce;
    tw_accretion->query_type_elsize = sizeof(TreeWalkQueryBHAccretion);
    tw_accretion->result_type_elsize = sizeof(TreeWalkResultBHAccretion);
    tw_accretion->tree = tree;
    tw_accretion->priv = priv;

    /*************************************************************************/

    TreeWalk tw_feedback[1] = {{0}};
    tw_feedback->ev_label = "BH_FEEDBACK";
    tw_feedback->visit = (TreeWalkVisitFunction) treewalk_visit_ngbiter;
    tw_feedback->ngbiter_type_elsize = sizeof(TreeWalkNgbIterBHFeedback);
    tw_feedback->ngbiter = (TreeWalkNgbIterFunction) blackhole_feedback_ngbiter;
    tw_feedback->haswork = blackhole_feedback_haswork;
    tw_feedback->fill = (TreeWalkFillQueryFunction) blackhole_feedback_copy;
    tw_feedback->postprocess = (TreeWalkProcessFunction) blackhole_feedback_postprocess;
    tw_feedback->reduce = (TreeWalkReduceResultFunction) blackhole_feedback_reduce;
    tw_feedback->query_type_elsize = sizeof(TreeWalkQueryBHFeedback);
    tw_feedback->result_type_elsize = sizeof(TreeWalkResultBHFeedback);
    tw_feedback->tree = tree;
    tw_feedback->priv = priv;
    tw_feedback->repeatdisallowed = 1;

    priv->atime = atime;
    priv->a3inv = 1./(atime * atime * atime);
    priv->hubble = hubble_function(CP, atime);
    priv->CP = CP;

    /* Build the queue once, since it is really 'all black holes' and similar for all treewalks*/
    treewalk_build_queue(tw_dynfric, act->ActiveParticle, act->NumActiveParticle, 0);
    /* Now we have a BH queue and we can re-use it*/
    int * ActiveBlackHoles = tw_dynfric->WorkSet;
    int64_t NumActiveBlackHoles = tw_dynfric->WorkSetSize;
    /* If this queue is empty, nothing to do.*/
    MPI_Allreduce(&NumActiveBlackHoles, &totbh, 1, MPI_INT64, MPI_SUM, MPI_COMM_WORLD);
    if(totbh == 0) {
        myfree(ActiveBlackHoles);
        return;
    }

    /* We can re-use the current queue for these treewalks*/
    tw_accretion->haswork = NULL;
    tw_dynfric->haswork = NULL;

    /*************************************************************************/
    /*  Dynamical Friction Treewalk */

    /* Environment variables for DF */
    priv->BH_SurroundingRmsVel = (MyFloat *) mymalloc("BH_SurroundingRmsVel", SlotsManager->info[5].size * sizeof(priv->BH_SurroundingRmsVel));
    priv->BH_SurroundingVel = (MyFloat (*) [3]) mymalloc("BH_SurroundingVel", 3* SlotsManager->info[5].size * sizeof(priv->BH_SurroundingVel[0]));
    priv->BH_SurroundingParticles = (int *)mymalloc("BH_SurroundingParticles", SlotsManager->info[5].size * sizeof(priv->BH_SurroundingParticles));
    priv->BH_SurroundingDensity = (MyFloat *) mymalloc("BH_SurroundingDensity", SlotsManager->info[5].size * sizeof(priv->BH_SurroundingDensity));
    /* guard treewalk */
    if (blackhole_params.BH_DynFrictionMethod > 0)
        treewalk_run(tw_dynfric, ActiveBlackHoles, NumActiveBlackHoles);

    /*************************************************************************/

    walltime_measure("/BH/DynFric");

    /* Let's determine which particles may be swallowed and calculate total feedback weights */
    priv->SPH_SwallowID = (MyIDType *) mymalloc("SPH_SwallowID", SlotsManager->info[0].size * sizeof(MyIDType));
    memset(priv->SPH_SwallowID, 0, SlotsManager->info[0].size * sizeof(MyIDType));

    /* Computed in accretion, used in feedback*/
    priv->BH_FeedbackWeightSum = (MyFloat *) mymalloc("BH_FeedbackWeightSum", SlotsManager->info[5].size * sizeof(MyFloat));

    /* These are initialized in preprocess and used to reposition the BH in postprocess*/
    priv->MinPot = (MyFloat *) mymalloc("BH_MinPot", SlotsManager->info[5].size * sizeof(MyFloat));

    /* Local to this treewalk*/
    priv->BH_Entropy = (MyFloat *) mymalloc("BH_Entropy", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->BH_SurroundingGasVel = (MyFloat (*) [3]) mymalloc("BH_SurroundVel", 3* SlotsManager->info[5].size * sizeof(priv->BH_SurroundingGasVel[0]));
    
    /* For AGN kinetic feedback */
    priv->NumDM = mymalloc("NumDM", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->V2sumDM = mymalloc("V2sumDM", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->V1sumDM = (MyFloat (*) [3]) mymalloc("V1sumDM", 3* SlotsManager->info[5].size * sizeof(priv->V1sumDM[0]));
    priv->MgasEnc = mymalloc("MgasEnc", SlotsManager->info[5].size * sizeof(MyFloat));
    /* mark the state of AGN kinetic feedback */
    priv->KEflag = mymalloc("KEflag", SlotsManager->info[5].size * sizeof(int));

    /* This allocates memory*/
    treewalk_run(tw_accretion, ActiveBlackHoles, NumActiveBlackHoles);

    /*************************************************************************/

    walltime_measure("/BH/Accretion");
    MPIU_Barrier(MPI_COMM_WORLD);

    /* Now do the swallowing of particles and dump feedback energy */

    /* Ionization counters*/
    priv[0].N_sph_swallowed = ta_malloc("n_sph_swallowed", int64_t, omp_get_max_threads());
    priv[0].N_BH_swallowed = ta_malloc("n_BH_swallowed", int64_t, omp_get_max_threads());
    memset(priv[0].N_sph_swallowed, 0, sizeof(int64_t) * omp_get_max_threads());
    memset(priv[0].N_BH_swallowed, 0, sizeof(int64_t) * omp_get_max_threads());

    priv->BH_accreted_Mass = (MyFloat *) mymalloc("BH_accretedmass", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->BH_accreted_BHMass = (MyFloat *) mymalloc("BH_accreted_BHMass", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->BH_accreted_Mtrack = (MyFloat *) mymalloc("BH_accreted_Mtrack", SlotsManager->info[5].size * sizeof(MyFloat));
    priv->BH_accreted_momentum = (MyFloat (*) [3]) mymalloc("BH_accretemom", 3* SlotsManager->info[5].size * sizeof(priv->BH_accreted_momentum[0]));

    treewalk_run(tw_feedback, ActiveBlackHoles, NumActiveBlackHoles);

    /*************************************************************************/
    walltime_measure("/BH/Feedback");

    if(FdBlackholeDetails){
        collect_BH_info(ActiveBlackHoles, NumActiveBlackHoles, priv, FdBlackholeDetails);
    }

    myfree(priv->BH_accreted_momentum);
    myfree(priv->BH_accreted_Mtrack);
    myfree(priv->BH_accreted_BHMass);
    myfree(priv->BH_accreted_Mass);

    /*****************************************************************/
    myfree(priv->KEflag);
    myfree(priv->MgasEnc);
    myfree(priv->V1sumDM);
    myfree(priv->V2sumDM);
    myfree(priv->NumDM);
    
    myfree(priv->BH_SurroundingGasVel);
    myfree(priv->BH_Entropy);
    myfree(priv->MinPot);

    myfree(priv->BH_FeedbackWeightSum);
    myfree(priv->SPH_SwallowID);

    /*****************************************************************/
    myfree(priv->BH_SurroundingDensity);
    myfree(priv->BH_SurroundingParticles);
    myfree(priv->BH_SurroundingVel);
    myfree(priv->BH_SurroundingRmsVel);

    /*****************************************************************/
    myfree(ActiveBlackHoles);

    int64_t Ntot_gas_swallowed, Ntot_BH_swallowed;
    int64_t N_sph_swallowed = 0, N_BH_swallowed = 0;
    for(i = 0; i < omp_get_max_threads(); i++) {
        N_sph_swallowed += priv[0].N_sph_swallowed[i];
        N_BH_swallowed += priv[0].N_BH_swallowed[i];
    }
    ta_free(priv[0].N_BH_swallowed);
    ta_free(priv[0].N_sph_swallowed);

    MPI_Reduce(&N_sph_swallowed, &Ntot_gas_swallowed, 1, MPI_INT64, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT64, MPI_SUM, 0, MPI_COMM_WORLD);

    MPIU_Barrier(MPI_COMM_WORLD);
    message(0, "Accretion done: %d gas particles swallowed, %d BH particles swallowed\n",
                Ntot_gas_swallowed, Ntot_BH_swallowed);

    int total_bh;
    double total_mdoteddington;
    double total_mass_holes, total_mdot;

    double Local_BH_mass = 0;
    double Local_BH_Mdot = 0;
    double Local_BH_Medd = 0;
    int Local_BH_num = 0;
    /* Compute total mass of black holes
     * present by summing contents of black hole array*/
    #pragma omp parallel for reduction(+ : Local_BH_num) reduction(+: Local_BH_mass) reduction(+: Local_BH_Mdot) reduction(+: Local_BH_Medd)
    for(i = 0; i < SlotsManager->info[5].size; i ++)
    {
        if(BhP[i].SwallowID != (MyIDType) -1)
            continue;
        Local_BH_num++;
        Local_BH_mass += BhP[i].Mass;
        Local_BH_Mdot += BhP[i].Mdot;
        Local_BH_Medd += BhP[i].Mdot/BhP[i].Mass;
    }

    MPI_Reduce(&Local_BH_mass, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Local_BH_Mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Local_BH_Medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Local_BH_num, &total_bh, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(FdBlackHoles)
    {
        /* convert to solar masses per yr */
        double mdot_in_msun_per_year =
            total_mdot * (units.UnitMass_in_g / SOLAR_MASS) / (units.UnitTime_in_s / SEC_PER_YEAR);

        total_mdoteddington *= 1.0 / ((4 * M_PI * GRAVITY * LIGHTCGS * PROTONMASS /
                    (0.1 * LIGHTCGS * LIGHTCGS * THOMPSON)) * units.UnitTime_in_s);

        fprintf(FdBlackHoles, "%g %d %g %g %g %g\n",
                atime, total_bh, total_mass_holes, total_mdot, mdot_in_msun_per_year, total_mdoteddington);
        fflush(FdBlackHoles);
    }
    walltime_measure("/BH/Info");
}


/*************************************************************************************/
/* DF routines */
static void
blackhole_dynfric_postprocess(int n, TreeWalk * tw){

    int PI = P[n].PI;
    int j;

    /***********************************************************************************/
    /* This is Gizmo's implementation of dynamic friction                              */
    /* c.f. section 3.1 in http://www.tapir.caltech.edu/~phopkins/public/notes_blackholes.pdf */
    /* Compute dynamic friction accel when DF turned on                                */
    /* averaged value for colomb logarithm and integral over the distribution function */
    /* acc_friction = -4*pi*G^2 * Mbh * log(lambda) * rho * f_of_x * bhvel / |bhvel^3| */
    /*       f_of_x = [erf(x) - 2*x*exp(-x^2)/sqrt(pi)]                                */
    /*       lambda = b_max * v^2 / G / (M+m)                                          */
    /*        b_max = Size of system (e.g. Rvir)                                       */
    /*            v = Relative velocity of BH with respect to the environment          */
    /*            M = Mass of BH                                                       */
    /*            m = individual mass elements composing the large system (e.g. m<<M)  */
    /*            x = v/sqrt(2)/sigma                                                  */
    /*        sigma = width of the max. distr. of the host system                      */
    /*                (e.g. sigma = v_disp / 3                                         */

    if(BH_GET_PRIV(tw)->BH_SurroundingDensity[PI] > 0){
        double bhvel;
        double lambda, x, f_of_x;
        const double a_erf = 8 * (M_PI - 3) / (3 * M_PI * (4. - M_PI));

        /* normalize velocity/dispersion */
        BH_GET_PRIV(tw)->BH_SurroundingRmsVel[PI] /= BH_GET_PRIV(tw)->BH_SurroundingDensity[PI];
        BH_GET_PRIV(tw)->BH_SurroundingRmsVel[PI] = sqrt(BH_GET_PRIV(tw)->BH_SurroundingRmsVel[PI]);
        for(j = 0; j < 3; j++)
            BH_GET_PRIV(tw)->BH_SurroundingVel[PI][j] /= BH_GET_PRIV(tw)->BH_SurroundingDensity[PI];

        /* Calculate Coulumb Logarithm */
        bhvel = 0;
        for(j = 0; j < 3; j++)
        {
            bhvel += pow(P[n].Vel[j] - BH_GET_PRIV(tw)->BH_SurroundingVel[PI][j], 2);
        }
        bhvel = sqrt(bhvel);

        /* There is no parameter in physical unit, so I kept everything in code unit */

        x = bhvel / sqrt(2) / (BH_GET_PRIV(tw)->BH_SurroundingRmsVel[PI] / 3);
        /* First term is aproximation of the error function */
        f_of_x = x / fabs(x) * sqrt(1 - exp(-x * x * (4 / M_PI + a_erf * x * x)
            / (1 + a_erf * x * x))) - 2 * x / sqrt(M_PI) * exp(-x * x);
        /* Floor at zero */
        if (f_of_x < 0)
            f_of_x = 0;

        lambda = 1. + blackhole_params.BH_DFbmax * pow((bhvel/BH_GET_PRIV(tw)->atime),2) / BH_GET_PRIV(tw)->CP->GravInternal / P[n].Mass;

        for(j = 0; j < 3; j++)
        {
            /* prevent DFAccel from exploding */
            if(bhvel > 0){
                BHP(n).DFAccel[j] = - 4. * M_PI * BH_GET_PRIV(tw)->CP->GravInternal * BH_GET_PRIV(tw)->CP->GravInternal * P[n].Mass * BH_GET_PRIV(tw)->BH_SurroundingDensity[PI] *
                log(lambda) * f_of_x * (P[n].Vel[j] - BH_GET_PRIV(tw)->BH_SurroundingVel[PI][j]) / pow(bhvel, 3);
                BHP(n).DFAccel[j] *= BH_GET_PRIV(tw)->atime;  // convert to code unit of acceleration
                BHP(n).DFAccel[j] *= blackhole_params.BH_DFBoostFactor; // Add a boost factor
            }
            else{
                BHP(n).DFAccel[j] = 0;
            }
        }
#ifdef DEBUG
        message(2,"x=%e, log(lambda)=%e, fof_x=%e, Mbh=%e, ratio=%e \n",
           x,log(lambda),f_of_x,P[n].Mass,BHP(n).DFAccel[0]/P[n].GravAccel[0]);
#endif
    }
    else
    {
        message(2, "Dynamic Friction density is zero for BH %ld. Surroundingpart %d, mass %g, hsml %g, dens %g, pos %g %g %g.\n",
            P[n].ID, BH_GET_PRIV(tw)->BH_SurroundingParticles[PI], BHP(n).Mass, P[n].Hsml, BHP(n).Density, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]);
        for(j = 0; j < 3; j++)
        {
            BHP(n).DFAccel[j] = 0;
        }
    }
}

/*******************************************************************/
static int
blackhole_dynfric_haswork(int n, TreeWalk * tw){
    /*Black hole not being swallowed*/
    return (P[n].Type == 5) && (!P[n].Swallowed);
}

static void
blackhole_dynfric_reduce(int place, TreeWalkResultBHDynfric * remote, enum TreeWalkReduceMode mode, TreeWalk * tw){
    int PI = P[place].PI;

    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingDensity[PI], remote->SurroundingDensity);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingParticles[PI], remote->SurroundingParticles);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingVel[PI][0], remote->SurroundingVel[0]);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingVel[PI][1], remote->SurroundingVel[1]);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingVel[PI][2], remote->SurroundingVel[2]);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingRmsVel[PI], remote->SurroundingRmsVel);

}

static void
blackhole_dynfric_copy(int place, TreeWalkQueryBHDynfric * I, TreeWalk * tw){
    /* SPH kernel width should be the only thing needed */
    I->Hsml = P[place].Hsml;
}


static void
blackhole_dynfric_ngbiter(TreeWalkQueryBHDynfric * I,
        TreeWalkResultBHDynfric * O,
        TreeWalkNgbIterBHDynfric * iter,
        LocalTreeWalk * lv){

   if(iter->base.other == -1) {
        iter->base.mask = 1 + 2 + 4 + 8 + 16 + 32;
        iter->base.Hsml = I->Hsml;
        iter->base.symmetric = NGB_TREEFIND_ASYMMETRIC;
        density_kernel_init(&iter->dynfric_kernel, I->Hsml, GetDensityKernelType());
        return;
    }

    int other = iter->base.other;
    double r = iter->base.r;
    double r2 = iter->base.r2;

    /* Collect Star/+DM/+Gas density/velocity for DF computation */
    if(P[other].Type == 4 || (P[other].Type == 1 && blackhole_params.BH_DynFrictionMethod > 1) ||
        (P[other].Type == 0 && blackhole_params.BH_DynFrictionMethod == 3) ){
        if(r2 < iter->dynfric_kernel.HH) {
            double u = r * iter->dynfric_kernel.Hinv;
            double wk = density_kernel_wk(&iter->dynfric_kernel, u);
            float mass_j = P[other].Mass;
            int k;
            O->SurroundingParticles += 1;
            O->SurroundingDensity += (mass_j * wk);
            for (k = 0; k < 3; k++){
                O->SurroundingVel[k] += (mass_j * wk * P[other].Vel[k]);
                O->SurroundingRmsVel += (mass_j * wk * pow(P[other].Vel[k], 2));
            }
        }
    }
}

/*************************************************************************************/


static void
blackhole_accretion_postprocess(int i, TreeWalk * tw)
{
    int k;
    int PI = P[i].PI;
    if(BHP(i).Density > 0)
    {
        BH_GET_PRIV(tw)->BH_Entropy[PI] /= BHP(i).Density;
        for(k = 0; k < 3; k++)
            BH_GET_PRIV(tw)->BH_SurroundingGasVel[PI][k] /= BHP(i).Density;
    }

    double mdot = 0;		/* if no accretion model is enabled, we have mdot=0 */

    double rho = BHP(i).Density;
    double bhvel = 0;
    for(k = 0; k < 3; k++)
        bhvel += pow(P[i].Vel[k] - BH_GET_PRIV(tw)->BH_SurroundingGasVel[PI][k], 2);

    bhvel = sqrt(bhvel);
    bhvel /= BH_GET_PRIV(tw)->atime;
    double rho_proper = rho * BH_GET_PRIV(tw)->a3inv;

    double soundspeed = blackhole_soundspeed(BH_GET_PRIV(tw)->BH_Entropy[PI], rho, BH_GET_PRIV(tw)->atime);

    /* Note: we take here a radiative efficiency of 0.1 for Eddington accretion */
    double meddington = (4 * M_PI * GRAVITY * LIGHTCGS * PROTONMASS / (0.1 * LIGHTCGS * LIGHTCGS * THOMPSON)) * BHP(i).Mass
        * BH_GET_PRIV(tw)->units.UnitTime_in_s / BH_GET_PRIV(tw)->CP->HubbleParam;

    double norm = pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);

    if(norm > 0)
        mdot = 4. * M_PI * blackhole_params.BlackHoleAccretionFactor * BH_GET_PRIV(tw)->CP->GravInternal * BH_GET_PRIV(tw)->CP->GravInternal *
            BHP(i).Mass * BHP(i).Mass * rho_proper / norm;

    if(blackhole_params.BlackHoleEddingtonFactor > 0.0 &&
        mdot > blackhole_params.BlackHoleEddingtonFactor * meddington) {
        mdot = blackhole_params.BlackHoleEddingtonFactor * meddington;
    }
    BHP(i).Mdot = mdot;

    double dtime = get_dloga_for_bin(P[i].TimeBin, P[i].Ti_drift) / BH_GET_PRIV(tw)->hubble;

    BHP(i).Mass += BHP(i).Mdot * dtime;

    /*************************************************************************/

    if(blackhole_params.BH_DRAG > 0){
        /* a_BH = (v_gas - v_BH) Mdot/M_BH                                   */
        /* motivated by BH gaining momentum from the accreted gas            */
        /*c.f.section 3.2,in http://www.tapir.caltech.edu/~phopkins/public/notes_blackholes.pdf */
        double fac = 0;
        if (blackhole_params.BH_DRAG == 1) fac = BHP(i).Mdot/P[i].Mass;
        if (blackhole_params.BH_DRAG == 2) fac = blackhole_params.BlackHoleEddingtonFactor * meddington/BHP(i).Mass;
        fac *= BH_GET_PRIV(tw)->atime; /* dv = acc * kick_fac = acc * a^{-1}dt, therefore acc = a*dv/dt  */
        for(k = 0; k < 3; k++) {
            BHP(i).DragAccel[k] = -(P[i].Vel[k] - BH_GET_PRIV(tw)->BH_SurroundingGasVel[PI][k])*fac;
        }
    }
    else{
        for(k = 0; k < 3; k++){
            BHP(i).DragAccel[k] = 0;
        }
    }
    /*************************************************************************/
    
    if(blackhole_params.BlackHoleKineticOn == 1){
        /* accumulate kenetic feedback energy by dE = epsilon x mdot x c^2 */
        /* epsilon = Min(rho_BH/(BHKE_EffRhoFactor*rho_sfr),BHKE_EffCap)   */
        /* KE is released when exceeding injection energy threshold        */
        BH_GET_PRIV(tw)->KEflag[PI] = 0;
        double Edd_ratio = BHP(i).Mdot/meddington;
        double lam_thresh = blackhole_params.BHKE_EddingtonThrFactor;
        double x = blackhole_params.BHKE_EddingtonMFactor * pow(BHP(i).Mass/blackhole_params.BHKE_EddingtonMPivot, blackhole_params.BHKE_EddingtonMIndex);
        if (lam_thresh > x)
            lam_thresh = x;
        if (Edd_ratio < lam_thresh){
            /* mark this timestep is accumulating KE feedback energy */
            BH_GET_PRIV(tw)->KEflag[PI] = 1;
            const double rho_crit_baryon = BH_GET_PRIV(tw)->CP->OmegaBaryon * 3 * pow(BH_GET_PRIV(tw)->CP->Hubble, 2) / (8 * M_PI * BH_GET_PRIV(tw)->CP->GravInternal);
            const double rho_sfr = blackhole_params.BHKE_SfrCritOverDensity * rho_crit_baryon;
            double epsilon = (BHP(i).Density/rho_sfr)/blackhole_params.BHKE_EffRhoFactor;
            if (epsilon > blackhole_params.BHKE_EffCap){
                epsilon = blackhole_params.BHKE_EffCap;
            }
            
            BHP(i).KineticFdbkEnergy += epsilon * (BHP(i).Mdot * dtime * pow(LIGHTCGS / BH_GET_PRIV(tw)->units.UnitVelocity_in_cm_per_s, 2));
        }
        
        /* decide whether to release KineticFdbkEnergy*/
        double vdisp = 0;
        double numdm = BH_GET_PRIV(tw)->NumDM[PI];
        if (numdm>0){
            vdisp = BH_GET_PRIV(tw)->V2sumDM[PI]/numdm;
            int d;
            for(d = 0; d<3; d++){
                vdisp -= pow(BH_GET_PRIV(tw)->V1sumDM[PI][d]/numdm,2);
            }
            if(vdisp > 0)
                vdisp = sqrt(vdisp / 3);
        }
        
        double KE_thresh = 0.5 * vdisp * vdisp * BH_GET_PRIV(tw)->MgasEnc[PI];
        KE_thresh *= blackhole_params.BHKE_InjEnergyThr;
        
        if (BHP(i).KineticFdbkEnergy > KE_thresh){
            /* mark KineticFdbkEnergy is ready to be released in the feedback treewalk */
            BH_GET_PRIV(tw)->KEflag[PI] = 2;
        }
    }
}

static void
blackhole_accretion_preprocess(int n, TreeWalk * tw)
{
    int j;
    BH_GET_PRIV(tw)->MinPot[P[n].PI] = P[n].Potential;

    for(j = 0; j < 3; j++) {
        BHP(n).MinPotPos[j] = P[n].Pos[j];
    }

}

static void
blackhole_feedback_postprocess(int n, TreeWalk * tw)
{
    const int PI = P[n].PI;
    if(BH_GET_PRIV(tw)->BH_accreted_BHMass[PI] > 0){
       BHP(n).Mass += BH_GET_PRIV(tw)->BH_accreted_BHMass[PI];
    }
    if(BH_GET_PRIV(tw)->BH_accreted_Mass[PI] > 0)
    {
        /* velocity feedback due to accretion; momentum conservation.
         * This does nothing with repositioning on.*/
        const MyFloat accmass = BH_GET_PRIV(tw)->BH_accreted_Mass[PI];
        int k;
        /* Need to add the momentum from Mtrack as well*/
        for(k = 0; k < 3; k++)
            P[n].Vel[k] = (P[n].Vel[k] * P[n].Mass + BH_GET_PRIV(tw)->BH_accreted_momentum[PI][k]) /
                    (P[n].Mass + accmass + BH_GET_PRIV(tw)->BH_accreted_Mtrack[PI]);
        P[n].Mass += accmass;
    }

    if(blackhole_params.SeedBHDynMass>0){
        if(BH_GET_PRIV(tw)->BH_accreted_Mtrack[PI] > 0){
            BHP(n).Mtrack += BH_GET_PRIV(tw)->BH_accreted_Mtrack[PI];
        }
        if(BHP(n).Mtrack > blackhole_params.SeedBHDynMass){
            BHP(n).Mtrack = blackhole_params.SeedBHDynMass; /*cap Mtrack at SeedBHDynMass*/
        }
    }
    /* Reset KineticFdbkEnerg to 0 after released */
    if(BH_GET_PRIV(tw)->KEflag[PI] == 2){
        BHP(n).KineticFdbkEnergy = 0;
    }
}

static void
blackhole_accretion_ngbiter(TreeWalkQueryBHAccretion * I,
        TreeWalkResultBHAccretion * O,
        TreeWalkNgbIterBHAccretion * iter,
        LocalTreeWalk * lv)
{

    if(iter->base.other == -1) {
        O->BH_minTimeBin = TIMEBINS;
        O->encounter = 0;

        O->BH_MinPot = BHPOTVALUEINIT;

        int d;
        for(d = 0; d < 3; d++) {
            O->BH_MinPotPos[d] = I->base.Pos[d];
        }
        iter->base.mask = 1 + 2 + 4 + 8 + 16 + 32;
        iter->base.Hsml = I->Hsml;
        iter->base.symmetric = NGB_TREEFIND_ASYMMETRIC;

        density_kernel_init(&iter->accretion_kernel, I->Hsml, GetDensityKernelType());
        density_kernel_init(&iter->feedback_kernel, I->Hsml, GetDensityKernelType());
        return;
    }

    int other = iter->base.other;
    double r = iter->base.r;
    double r2 = iter->base.r2;

    if(P[other].Mass < 0) return;

    if(P[other].Type != 5) {
        if (O->BH_minTimeBin > P[other].TimeBin)
            O->BH_minTimeBin = P[other].TimeBin;
    }

     /* BH does not accrete wind */
    if(winds_is_particle_decoupled(other)) return;

    /* Find the black hole potential minimum. */
    if(r2 < iter->accretion_kernel.HH)
    {
        if(P[other].Potential < O->BH_MinPot)
        {
            int d;
            O->BH_MinPot = P[other].Potential;
            for(d = 0; d < 3; d++) {
                O->BH_MinPotPos[d] = P[other].Pos[d];
                O->BH_MinPotVel[d] = P[other].Vel[d];
            }
        }
    }

    /* Accretion / merger doesn't do self interaction */
    if(P[other].ID == I->ID) return;

    /* we have a black hole merger. Now we use 2 times GravitationalSoftening as merging criteria, previously we used the SPH smoothing length. */
    if(P[other].Type == 5 && r < (2*FORCE_SOFTENING(0,1)/2.8))
    {
        O->encounter = 1; // mark the event when two BHs encounter each other

        int flag = 0; // the flag for BH merge

        if(blackhole_params.BlackHoleRepositionEnabled == 1) // directly merge if reposition is enabled
            flag = 1;
        if(blackhole_params.MergeGravBound == 0)
            flag = 1;
        /* apply Grav Bound check only when Reposition is disabled, otherwise BHs would be repositioned to the same location but not merge */
        if(blackhole_params.MergeGravBound == 1 && blackhole_params.BlackHoleRepositionEnabled == 0){

            double dx[3];
            double dv[3];
            double da[3];
            int d;

            for(d = 0; d < 3; d++){
                dx[d] = NEAREST(I->base.Pos[d] - P[other].Pos[d], PartManager->BoxSize);
                dv[d] = I->Vel[d] - P[other].Vel[d];
                /* we include long range PM force, short range force and DF */
                da[d] = (I->Accel[d] - P[other].GravAccel[d] - P[other].GravPM[d] - BHP(other).DFAccel[d]);
            }
            flag = check_grav_bound(dx,dv,da, BH_GET_PRIV(lv->tw)->atime);
            /*if(flag == 0)
                message(0, "dx %g %g %g dv %g %g %g da %g %g %g\n",dx[0], dx[1], dx[2], dv[0], dv[1], dv[2], da[0], da[1], da[2]);*/
        }

        /* do the merge */
        if(flag == 1)
        {
            O->encounter = 0;
            MyIDType readid, newswallowid;

            #pragma omp atomic read
            readid = (BHP(other).SwallowID);

            /* Here we mark the black hole as "ready to be swallowed" using the SwallowID.
             * The actual swallowing is done in the feedback treewalk by setting Swallowed = 1
             * and merging the masses.*/
            do {
                /* Generate the new ID from the old*/
                if(readid != (MyIDType) -1 && readid < I->ID ) {
                    /* Already marked, prefer to be swallowed by a bigger ID */
                    newswallowid = I->ID;
                } else if(readid == (MyIDType) -1 && P[other].ID < I->ID) {
                    /* Unmarked, the BH with bigger ID swallows */
                    newswallowid = I->ID;
                }
                else
                    break;
            /* Swap in the new id only if the old one hasn't changed:
             * in principle an extension, but supported on at least clang >= 9, gcc >= 5 and icc >= 18.*/
            } while(!__atomic_compare_exchange_n(&(BHP(other).SwallowID), &readid, newswallowid, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
        }
    }


    if(P[other].Type == 0) {
        if(r2 < iter->accretion_kernel.HH) {
            double u = r * iter->accretion_kernel.Hinv;
            double wk = density_kernel_wk(&iter->accretion_kernel, u);
            float mass_j = P[other].Mass;

            O->SmoothedEntropy += (mass_j * wk * SPHP(other).Entropy);
            O->GasVel[0] += (mass_j * wk * P[other].Vel[0]);
            O->GasVel[1] += (mass_j * wk * P[other].Vel[1]);
            O->GasVel[2] += (mass_j * wk * P[other].Vel[2]);

            /* here we have a gas particle; check for swallowing */

            /* compute accretion probability */
            double p = 0;

            MyFloat BHPartMass = I->Mass;
            /* If SeedBHDynMass is larger than gas paricle mass, we use Mtrack to do the gas accretion
             * when BHP.Mass < SeedBHDynMass. Mtrack is initialized as gas particle mass and is capped
             * at SeedBHDynMass. Mtrack traces the BHP.Mass by stochastically swallowing gas and
             * therefore ensures mass conservation.*/
            if(blackhole_params.SeedBHDynMass > 0 && I->Mtrack < blackhole_params.SeedBHDynMass)
                BHPartMass = I->Mtrack;

            /* This is an averaged Mdot, because Mdot increases BH_Mass but not Mass.
             * So if the total accretion is significantly above the dynamical mass,
             * a particle is swallowed. */
            if((I->BH_Mass - BHPartMass) > 0 && I->Density > 0)
                p = (I->BH_Mass - BHPartMass) * wk / I->Density;

            /* compute random number, uniform in [0,1] */
            const double w = get_random_number(P[other].ID);
            if(w < p)
            {
                MyIDType * SPH_SwallowID = BH_GET_PRIV(lv->tw)->SPH_SwallowID;
                MyIDType readid, newswallowid;
                #pragma omp atomic read
                readid = SPH_SwallowID[P[other].PI];
                do {
                    /* Already marked, prefer to be swallowed by a bigger ID.
                     * Not marked, the SwallowID is 0 */
                    if(readid < I->ID + 1) {
                        newswallowid = I->ID + 1;
                    }
                    else
                        break;
                    /* Swap in the new id only if the old one hasn't changed*/
                } while(!__atomic_compare_exchange_n(&SPH_SwallowID[P[other].PI], &readid, newswallowid, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
            }
        }

        if(r2 < iter->feedback_kernel.HH) {
            /* update the feedback weighting */
            double mass_j;
            if(HAS(blackhole_params.BlackHoleFeedbackMethod, BH_FEEDBACK_OPTTHIN)) {
                double redshift = 1./BH_GET_PRIV(lv->tw)->atime - 1;
                double nh0 = get_neutral_fraction_sfreff(redshift, BH_GET_PRIV(lv->tw)->hubble, &P[other], &SPHP(other));
                if(r2 > 0)
                    O->FeedbackWeightSum += (P[other].Mass * nh0) / r2;
            } else {
                if(HAS(blackhole_params.BlackHoleFeedbackMethod, BH_FEEDBACK_MASS)) {
                    mass_j = P[other].Mass;
                } else {
                    mass_j = P[other].Hsml * P[other].Hsml * P[other].Hsml;
                }
                if(HAS(blackhole_params.BlackHoleFeedbackMethod, BH_FEEDBACK_SPLINE)) {
                    double u = r * iter->feedback_kernel.Hinv;
                    O->FeedbackWeightSum += (mass_j *
                          density_kernel_wk(&iter->feedback_kernel, u)
                           );
                } else {
                    O->FeedbackWeightSum += (mass_j);
                }
            }
        }
    }
    
    /* collect info for sigmaDM and Menc for kinetic feedback */
    if(blackhole_params.BlackHoleKineticOn == 1){
        if(P[other].Type == 1 && r2 < iter->feedback_kernel.HH){
            O->NumDM += 1;
            int d;
            for(d = 0; d < 3; d++){
                double vel = P[other].Vel[d] - I->Vel[d];
                O->V1sumDM[d] += vel;
                O->V2sumDM += vel * vel;
            }
        }
        if(P[other].Type == 0 && r2 < iter->feedback_kernel.HH){
            O->MgasEnc += P[other].Mass;
        }
    }
}


/**
 * perform blackhole swallow / merger;
 */
static void
blackhole_feedback_ngbiter(TreeWalkQueryBHFeedback * I,
        TreeWalkResultBHFeedback * O,
        TreeWalkNgbIterBHFeedback * iter,
        LocalTreeWalk * lv)
{

    if(iter->base.other == -1) {

        iter->base.mask = 1 + 32;
        iter->base.Hsml = I->Hsml;
        /* Swallow is symmetric, but feedback dumping is asymetric;
         * we apply a cut in r to break the symmetry. */
        iter->base.symmetric = NGB_TREEFIND_SYMMETRIC;
        density_kernel_init(&iter->feedback_kernel, I->Hsml, GetDensityKernelType());
        return;
    }

    int other = iter->base.other;
    double r2 = iter->base.r2;
    double r = iter->base.r;
    /* Exclude self interaction */

    if(P[other].ID == I->ID) return;

     /* BH does not accrete wind */
    if(winds_is_particle_decoupled(other))
        return;


     /* we have a black hole merger! */
    if(P[other].Type == 5 && BHP(other).SwallowID != (MyIDType) -1)
    {
        if(BHP(other).SwallowID != I->ID) return;

        /* Swallow the particle*/
        /* A note on Swallowed vs SwallowID: black hole particles which have been completely swallowed
         * (ie, their mass has been added to another particle) have Swallowed = 1.
         * These particles are ignored in future tree walks. We set Swallowed here so that this process is atomic:
         * the total mass in the tree is always conserved.
         *
         * We also set SwallowID != -1 in the accretion treewalk. This marks the black hole as ready to be swallowed
         * by something. It is actually swallowed only by the nearby black hole with the largest ID. In rare cases
         * it may happen that the swallower is itself swallowed before swallowing the marked black hole. However,
         * in practice the new swallower should also take the marked black hole next timestep.
         */

        BHP(other).SwallowTime = BH_GET_PRIV(lv->tw)->atime;
        P[other].Swallowed = 1;
        O->BH_CountProgs += BHP(other).CountProgs;
        O->BH_Mass += (BHP(other).Mass);

        if (blackhole_params.SeedBHDynMass>0 && I->Mtrack>0){
        /* Make sure that the final dynamic mass (I->Mass + O->Mass) = MAX(SeedDynMass, total_gas_accreted),
           I->Mtrack only need to be updated when I->Mtrack < SeedBHDynMass, */
            if(I->Mtrack < blackhole_params.SeedBHDynMass && BHP(other).Mtrack < blackhole_params.SeedBHDynMass){
            /* I->Mass = SeedBHDynMass, total_gas_accreted = I->Mtrack + BHP(other).Mtrack */
                O->acMtrack += BHP(other).Mtrack;
                double delta_m = I->Mtrack + BHP(other).Mtrack - blackhole_params.SeedBHDynMass;
                O->Mass += DMAX(0,delta_m);
            }
            if(I->Mtrack >= blackhole_params.SeedBHDynMass && BHP(other).Mtrack < blackhole_params.SeedBHDynMass){
            /* I->Mass = gas_accreted, total_gas_accreted = I->Mass + BHP(other).Mtrack */
                O->Mass += BHP(other).Mtrack;
            }
            if(I->Mtrack < blackhole_params.SeedBHDynMass && BHP(other).Mtrack >= blackhole_params.SeedBHDynMass){
            /* I->Mass = SeedBHDynMass, P[other].Mass = gas_accreted,
               total_gas_accreted = I->track + P[other].Mass */
                O->acMtrack += BHP(other).Mtrack;
                O->Mass += (P[other].Mass + I->Mtrack - blackhole_params.SeedBHDynMass);
            }
            if(I->Mtrack >= blackhole_params.SeedBHDynMass && BHP(other).Mtrack >= blackhole_params.SeedBHDynMass){
            /* trivial case, total_gas_accreted = I->Mass + P[other].Mass */
                O->Mass += P[other].Mass;
            }
        }
        else{
            O->Mass += P[other].Mass;
        }

        /* Conserve momentum during accretion*/
        int d;
        for(d = 0; d < 3; d++)
            O->AccretedMomentum[d] += (P[other].Mass * P[other].Vel[d]);

        if(BHP(other).SwallowTime < BH_GET_PRIV(lv->tw)->atime)
            endrun(2, "Encountered BH %i swallowed at earlier time %g\n", other, BHP(other).SwallowTime);

        int tid = omp_get_thread_num();
        BH_GET_PRIV(lv->tw)->N_BH_swallowed[tid]++;

    }

    MyIDType * SPH_SwallowID = BH_GET_PRIV(lv->tw)->SPH_SwallowID;

    /* perform thermal or kinetic feedback into non-swallowed particles. */
    if(P[other].Type == 0 && SPH_SwallowID[P[other].PI] == 0 &&
        (r2 < iter->feedback_kernel.HH && P[other].Mass > 0) &&
            (I->FeedbackWeightSum > 0))
    {
        double u = r * iter->feedback_kernel.Hinv;
        double wk = 1.0;
        double mass_j;

        if(HAS(blackhole_params.BlackHoleFeedbackMethod, BH_FEEDBACK_MASS)) {
            mass_j = P[other].Mass;
        } else {
            mass_j = P[other].Hsml * P[other].Hsml * P[other].Hsml;
        }
        if(HAS(blackhole_params.BlackHoleFeedbackMethod, BH_FEEDBACK_SPLINE))
            wk = density_kernel_wk(&iter->feedback_kernel, u);
        
        /* thermal feedback */
        if(I->FeedbackEnergy > 0 && I->FdbkChannel == 0){
            const double injected_BH = I->FeedbackEnergy * mass_j * wk / I->FeedbackWeightSum;
            /* Set a flag for star-forming particles:
                * we want these to cool to the EEQOS via
                * tcool rather than trelax.*/
            if(sfreff_on_eeqos(&SPHP(other), BH_GET_PRIV(lv->tw)->a3inv)) {
                /* We cannot atomically set a bitfield.
                 * This flag is never read in this thread loop, and we are careful not to
                 * do this with a swallowed particle (as this can race with IsGarbage being set).
                 * So lack of atomicity is (I think) not a problem.*/
                //#pragma omp atomic write
                P[other].BHHeated = 1;
            }
            const double enttou = pow(SPHP(other).Density * BH_GET_PRIV(lv->tw)->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
            const double uu_in_cgs = BH_GET_PRIV(lv->tw)->units.UnitEnergy_in_cgs / BH_GET_PRIV(lv->tw)->units.UnitMass_in_g;

            double entold, entnew;
            double * entptr = &(SPHP(other).Entropy);
            #pragma omp atomic read
            entold = *entptr;
            do {
                entnew = add_injected_BH_energy(entold * enttou, injected_BH, P[other].Mass, uu_in_cgs) / enttou;
                /* Swap in the new gas entropy only if the old one hasn't changed.*/
            } while(!__atomic_compare_exchange(entptr, &entold, &entnew, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
        }
        
        /* kinetic feedback */
        if(I->KEFeedbackEnergy > 0 && I->FdbkChannel == 1){
            /* Kick the gas particle*/
            double dvel = sqrt(2 * I->KEFeedbackEnergy * wk / I->Density);
            double dir[3];
            get_random_dir(other, dir);
            int j;
            for(j = 0; j < 3; j++){
                #pragma omp atomic update
                P[other].Vel[j] += (dvel*dir[j]);
            }
        }
    }

    /* Swallowing a gas */
    /* This will only be true on one thread so we do not need a lock here*/
    /* Note that it will rarely happen that gas is swallowed by a BH which is itself swallowed.
     * In that case we do not swallow this particle: all swallowing changes before this are temporary*/
    if(P[other].Type == 0 && SPH_SwallowID[P[other].PI] == I->ID+1)
    {
        /* We do not know how to notify the tree of mass changes. so
         * enforce a mass conservation. */
        if(blackhole_params.SeedBHDynMass > 0 && I->Mtrack < blackhole_params.SeedBHDynMass) {
            /* we just add gas mass to Mtrack instead of dynMass */
            O->acMtrack += P[other].Mass;
        } else
            O->Mass += P[other].Mass;
        P[other].Mass = 0;
        /* Conserve momentum during accretion*/
        int d;
        for(d = 0; d < 3; d++)
            O->AccretedMomentum[d] += (P[other].Mass * P[other].Vel[d]);

        slots_mark_garbage(other, PartManager, SlotsManager);

        int tid = omp_get_thread_num();
        BH_GET_PRIV(lv->tw)->N_sph_swallowed[tid]++;
    }
}

static void
blackhole_accretion_reduce(int place, TreeWalkResultBHAccretion * remote, enum TreeWalkReduceMode mode, TreeWalk * tw)
{
    int k;
    MyFloat * MinPot = BH_GET_PRIV(tw)->MinPot;

    int PI = P[place].PI;
    if(MinPot[PI] > remote->BH_MinPot)
    {
        BHP(place).JumpToMinPot = blackhole_params.BlackHoleRepositionEnabled;
        MinPot[PI] = remote->BH_MinPot;
        for(k = 0; k < 3; k++) {
            /* Movement occurs in drift.c */
            BHP(place).MinPotPos[k] = remote->BH_MinPotPos[k];
            BHP(place).MinPotVel[k] = remote->BH_MinPotVel[k];
        }
    }

    BHP(place).encounter = remote->encounter;

    if (mode == 0 || BHP(place).minTimeBin > remote->BH_minTimeBin) {
        BHP(place).minTimeBin = remote->BH_minTimeBin;
    }

    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_FeedbackWeightSum[PI], remote->FeedbackWeightSum);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_Entropy[PI], remote->SmoothedEntropy);
    for (k = 0; k < 3; k++){
        TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_SurroundingGasVel[PI][k], remote->GasVel[k]);
        TREEWALK_REDUCE(BH_GET_PRIV(tw)->V1sumDM[PI][k], remote->V1sumDM[k]);
    }
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->NumDM[PI], remote->NumDM);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->V2sumDM[PI], remote->V2sumDM);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->MgasEnc[PI], remote->MgasEnc);
}

static void
blackhole_accretion_copy(int place, TreeWalkQueryBHAccretion * I, TreeWalk * tw)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        I->Vel[k] = P[place].Vel[k];
        I->Accel[k] = P[place].GravAccel[k] + P[place].GravPM[k] + BHP(place).DFAccel[k];
    }
    I->Hsml = P[place].Hsml;
    I->Mass = P[place].Mass;
    I->BH_Mass = BHP(place).Mass;
    I->Density = BHP(place).Density;
    I->ID = P[place].ID;
    I->Mtrack = BHP(place).Mtrack;
}

static int
blackhole_feedback_haswork(int n, TreeWalk * tw)
{
    /*Black hole not being swallowed*/
    return (P[n].Type == 5) && (!P[n].Swallowed) && (BHP(n).SwallowID == (MyIDType) -1);
}

static void
blackhole_feedback_copy(int i, TreeWalkQueryBHFeedback * I, TreeWalk * tw)
{
    I->Hsml = P[i].Hsml;
    I->BH_Mass = BHP(i).Mass;
    I->ID = P[i].ID;
    I->Mtrack = BHP(i).Mtrack;
    I->Density = BHP(i).Density;
    int PI = P[i].PI;

    I->FeedbackWeightSum = BH_GET_PRIV(tw)->BH_FeedbackWeightSum[PI];
    I->FdbkChannel = 0; /* thermal feedback mode */

    double dtime = get_dloga_for_bin(P[i].TimeBin, P[i].Ti_drift) / BH_GET_PRIV(tw)->hubble;

    I->FeedbackEnergy = blackhole_params.BlackHoleFeedbackFactor * 0.1 * BHP(i).Mdot * dtime *
                pow(LIGHTCGS / BH_GET_PRIV(tw)->units.UnitVelocity_in_cm_per_s, 2);
    I->KEFeedbackEnergy = 0;
    if (blackhole_params.BlackHoleKineticOn == 1 && BH_GET_PRIV(tw)->KEflag[PI] > 0){
        I->FdbkChannel = 1; /* kinetic feedback mode, (no thermal feedback for this timestep) */
        /* KEflag = 1: KEFeedbackEnergy is accumulating; KEflag = 2: KEFeedbackEnergy is released. */
        if (BH_GET_PRIV(tw)->KEflag[PI] == 2){
            I->KEFeedbackEnergy = BHP(i).KineticFdbkEnergy;
        }
    }
}

static void
blackhole_feedback_reduce(int place, TreeWalkResultBHFeedback * remote, enum TreeWalkReduceMode mode, TreeWalk * tw)
{
    int k;
    int PI = P[place].PI;

    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_accreted_Mass[PI], remote->Mass);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_accreted_BHMass[PI], remote->BH_Mass);
    TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_accreted_Mtrack[PI], remote->acMtrack);
    for(k = 0; k < 3; k++) {
        TREEWALK_REDUCE(BH_GET_PRIV(tw)->BH_accreted_momentum[PI][k], remote->AccretedMomentum[k]);
    }

    TREEWALK_REDUCE(BHP(place).CountProgs, remote->BH_CountProgs);
}

/* Sample from a power law to get the initial BH mass*/
static double
bh_powerlaw_seed_mass(MyIDType ID)
{
    /* compute random number, uniform in [0,1] */
    const double w = get_random_number(ID+23);
    /* Normalisation for this power law index*/
    double norm = pow(blackhole_params.MaxSeedBlackHoleMass, 1+blackhole_params.SeedBlackHoleMassIndex)
                - pow(blackhole_params.SeedBlackHoleMass, 1+blackhole_params.SeedBlackHoleMassIndex);
    /* Sample from the CDF:
     * w  = [M^(1+x) - M_0^(1+x)]/[M_1^(1+x) - M_0^(1+x)]
     * w [M_1^(1+x) - M_0^(1+x)] + M_0^(1+x) = M^(1+x)
     * M = pow((w [M_1^(1+x) - M_0^(1+x)] + M_0^(1+x)), 1/(1+x))*/
    double mass = pow(w * norm + pow(blackhole_params.SeedBlackHoleMass, 1+blackhole_params.SeedBlackHoleMassIndex),
                      1./(1+blackhole_params.SeedBlackHoleMassIndex));
    return mass;
}

void blackhole_make_one(int index, const double atime) {
    if(P[index].Type != 0)
        endrun(7772, "Only Gas turns into blackholes, what's wrong?");

    int child = index;

    /* Make the new particle a black hole: use all the P[i].Mass
     * so we don't have lots of low mass tracers.
     * If the BH seed mass is small this may lead to a mismatch
     * between the gas and BH mass. */
    child = slots_convert(child, 5, -1, PartManager, SlotsManager);

    BHP(child).base.ID = P[child].ID;
    /* The accretion mass should always be the seed black hole mass,
     * irrespective of the gravitational mass of the particle.*/
    if(blackhole_params.MaxSeedBlackHoleMass > 0)
        BHP(child).Mass = bh_powerlaw_seed_mass(P[child].ID);
    else
        BHP(child).Mass = blackhole_params.SeedBlackHoleMass;

    BHP(child).Mseed = BHP(child).Mass;
    BHP(child).Mdot = 0;
    BHP(child).FormationTime = atime;
    BHP(child).SwallowID = (MyIDType) -1;
    BHP(child).Density = 0;

    /* It is important to initialize MinPotPos to the current position of
     * a BH to avoid drifting to unknown locations (0,0,0) immediately
     * after the BH is created. */
    int j;
    for(j = 0; j < 3; j++) {
        BHP(child).MinPotPos[j] = P[child].Pos[j];
        BHP(child).DFAccel[j] = 0;
        BHP(child).DragAccel[j] = 0;
    }
    BHP(child).JumpToMinPot = 0;
    BHP(child).CountProgs = 1;

    if (blackhole_params.SeedBHDynMass>0){
        BHP(child).Mtrack = P[child].Mass;
        P[child].Mass = blackhole_params.SeedBHDynMass;
    }
    else{
        BHP(child).Mtrack = -1; /* This column is not used then. */
    }
    /* Initialize KineticFdbkEnergy, keep zero if BlackHoleKineticOn is not turned on */
    BHP(child).KineticFdbkEnergy = 0;
}
