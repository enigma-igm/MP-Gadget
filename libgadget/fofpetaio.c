#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include <bigfile-mpi.h>

#include "utils.h"
#include "utils/mpsort.h"

#include "partmanager.h"
#include "slotsmanager.h"
#include "petaio.h"
#include "exchange.h"
#include "fof.h"
#include "walltime.h"

static void fof_register_io_blocks(int MetalReturnOn, int BlackHoleOn, struct IOTable * IOTable);
static void fof_write_header(BigFile * bf, int64_t TotNgroups, const double atime, const double * MassTable, Cosmology * CP, MPI_Comm Comm);
static void build_buffer_fof(FOFGroups * fof, BigArray * array, IOTableEntry * ent, struct conversions * conv);

static int fof_distribute_particles(struct part_manager_type * halo_pman, struct slots_manager_type * halo_sman, MPI_Comm Comm);

static void fof_radix_Group_GrNr(const void * a, void * radix, void * arg) {
    uint64_t * u = (uint64_t *) radix;
    struct BaseGroup * f = (struct BaseGroup*) a;
    u[0] = f->GrNr;
}

void fof_save_particles(FOFGroups * fof, const char * OutputDir, const char * FOFFileBase, int num, int SaveParticles, Cosmology * CP, double atime, const double * MassTable, int MetalReturnOn, int BlackholeOn, MPI_Comm Comm) {
    int i;
    struct IOTable FOFIOTable = {0};
    char * fname = fastpm_strdup_printf("%s/%s_%03d", OutputDir, FOFFileBase, num);
    message(0, "Saving particle groups into %s\n", fname);

    fof_register_io_blocks(MetalReturnOn, BlackholeOn, &FOFIOTable);
    /* sort the groups according to group-number */
    mpsort_mpi(fof->Group, fof->Ngroups, sizeof(struct Group),
            fof_radix_Group_GrNr, 8, NULL, Comm);

    BigFile bf = {0};
    if(0 != big_file_mpi_create(&bf, fname, Comm)) {
        endrun(0, "Failed to open file at %s\n", fname);
    }
    myfree(fname);
    struct conversions conv = {0};
    conv.atime = atime;
    conv.hubble = hubble_function(CP, atime);

    MPIU_Barrier(Comm);
    fof_write_header(&bf, fof->TotNgroups, atime, MassTable, CP, Comm);

    for(i = 0; i < FOFIOTable.used; i ++) {
        /* only process the particle blocks */
        char blockname[128];
        int ptype = FOFIOTable.ent[i].ptype;
        BigArray array = {0};
        if(ptype == PTYPE_FOF_GROUP) {
            sprintf(blockname, "FOFGroups/%s", FOFIOTable.ent[i].name);
            build_buffer_fof(fof, &array, &FOFIOTable.ent[i], &conv);
            message(0, "Writing Block %s\n", blockname);

            petaio_save_block(&bf, blockname, &array, 1);
            petaio_destroy_buffer(&array);
        }
    }
    destroy_io_blocks(&FOFIOTable);
    walltime_measure("/FOF/IO/WriteFOF");

    if(SaveParticles) {
        struct IOTable IOTable = {0};
        register_io_blocks(&IOTable, 1, MetalReturnOn);
        struct part_manager_type halo_pman = {0};
        struct slots_manager_type halo_sman = {0};
        if(fof_distribute_particles(&halo_pman, &halo_sman, Comm)) {
            myfree(halo_sman.Base);
            myfree(halo_pman.Base);
            destroy_io_blocks(&IOTable);
            return;
        }

        int * selection = (int *) mymalloc("Selection", sizeof(int) * halo_pman.NumPart);

        int ptype_offset[6]={0};
        int ptype_count[6]={0};
        petaio_build_selection(selection, ptype_offset, ptype_count, halo_pman.Base, halo_pman.NumPart, NULL);

        walltime_measure("/FOF/IO/argind");

        for(i = 0; i < IOTable.used; i ++) {
            /* only process the particle blocks */
            char blockname[128];
            int ptype = IOTable.ent[i].ptype;
            BigArray array = {0};
            if(ptype < 6 && ptype >= 0) {
                sprintf(blockname, "%d/%s", ptype, IOTable.ent[i].name);
                petaio_build_buffer(&array, &IOTable.ent[i], selection + ptype_offset[ptype], ptype_count[ptype], halo_pman.Base, &halo_sman, &conv);

                message(0, "Writing Block %s\n", blockname);

                petaio_save_block(&bf, blockname, &array, 1);
                petaio_destroy_buffer(&array);
            }
        }
        myfree(selection);
        myfree(halo_sman.Base);
        myfree(halo_pman.Base);
        walltime_measure("/FOF/IO/WriteParticles");
        destroy_io_blocks(&IOTable);
    }

    big_file_mpi_close(&bf, Comm);

    message(0, "Group catalogues saved.\n");
}

struct PartIndex {
    uint64_t origin;
    union {
        int64_t sortKey;
        int targetTask;
    };
};

static int fof_sorted_layout(int i, const void * userdata) {
    struct part_manager_type * halo_pman = (struct part_manager_type *) userdata;
    return halo_pman->Base[i].TargetTask;
}

static void fof_radix_sortkey(const void * c1, void * out, void * arg) {
    uint64_t * u = (uint64_t *) out;
    const struct PartIndex * pi = (const struct PartIndex *) c1;
    *u = pi->sortKey;
}
static void fof_radix_origin(const void * c1, void * out, void * arg) {
    uint64_t * u = (uint64_t *) out;
    const struct PartIndex * pi = (const struct PartIndex *) c1;
    *u = pi->origin;
}
#if 0
/*Unused functions*/
static int fof_select_particle(int i) {
    return P[i].GrNr > 0;
}
static int fof_cmp_sortkey(const void * c1, const void * c2) {
    const struct PartIndex * p1 = c1;
    const struct PartIndex * p2 = c2;
    return (p1->sortKey > p2->sortKey) - (p1->sortKey < p2->sortKey);
}
static int fof_cmp_origin(const void * c1, const void * c2) {
    const struct PartIndex * p1 = c1;
    const struct PartIndex * p2 = c2;
    return (p1->origin > p2->origin) - (p1->origin < p2->origin);
}
#endif

static int
order_by_type_and_grnr(const void *a, const void *b)
{
    const struct particle_data * pa  = (const struct particle_data *) a;
    const struct particle_data * pb  = (const struct particle_data *) b;

    if(pa->Type < pb->Type)
        return -1;
    if(pa->Type > pb->Type)
        return +1;
    if(pa->GrNr < pb->GrNr)
        return -1;
    if(pa->GrNr > pb->GrNr)
        return +1;

    return 0;
}

static int
fof_try_particle_exchange(struct part_manager_type * halo_pman, struct slots_manager_type * halo_sman, MPI_Comm Comm)
{
    struct PartIndex * pi = (struct PartIndex *) mymalloc("PartIndex", sizeof(struct PartIndex) * halo_pman->NumPart);
    int ThisTask;
    MPI_Comm_rank(Comm, &ThisTask);

    int64_t i = 0;
    /* Build the index: needs to be done each time we loop as may have changed*/
    /* Yu: found it! this shall be int64 */
    const uint64_t task_origin_offset = PartManager->MaxPart + 1Lu;
    #pragma omp parallel for
    for(i = 0; i < halo_pman->NumPart; i ++) {
        pi[i].origin = task_origin_offset * ((uint64_t) ThisTask) + i;
        pi[i].sortKey = halo_pman->Base[i].GrNr;
    }
    /* sort pi to decide targetTask */
    mpsort_mpi(pi, halo_pman->NumPart, sizeof(struct PartIndex),
            fof_radix_sortkey, 8, NULL, Comm);

    //int64_t Npig = count_sum(NpigLocal);
    //int64_t offsetLocal = MPIU_cumsum(NpigLocal, Comm);

    //size_t chunksize = (Npig / NTask) + (Npig % NTask != 0);

    #pragma omp parallel for
    for(i = 0; i < halo_pman->NumPart; i ++) {
/* YU: A typo error here, should be IMIN, DMIN is for double but this should have tainted TargetTask,
   offset and chunksize are int  */
        //ptrdiff_t offset = offsetLocal + i;
        //pi[i].targetTask = IMIN(offset / chunksize, NTask - 1);
    /* YU: let's see if we keep the FOF particle load on the processes, IO would be faster
           (as at high z many ranks has no FOF), communication becomes sparse. */
        pi[i].targetTask = ThisTask;
    }
    /* return pi to the original processors */
    mpsort_mpi(pi, halo_pman->NumPart, sizeof(struct PartIndex), fof_radix_origin, 8, NULL, Comm);
    /* Target task is copied into the particle table, unioned with Dthsml.
     * This is a bit of a hack: probably the elegant thing to do is to unify slot
     * and main structure, then mpsort the combination directly. */
#ifdef DEBUG
    int NTask;
    MPI_Comm_size(Comm, &NTask);
    #pragma omp parallel for
    for(i = 0; i < halo_pman->NumPart; i ++) {
        halo_pman->Base[i].TargetTask = -1;
        if(pi[i].targetTask >= NTask || pi[i].targetTask < 0)
            endrun(23, "pi %d is impossible %d of %d tasks\n",i,pi[i].targetTask, NTask);
    }
#endif
    #pragma omp parallel for
    for(i = 0; i < halo_pman->NumPart; i ++) {
        size_t index = pi[i].origin % task_origin_offset;
        if(index >= (size_t) halo_pman->NumPart)
            endrun(23, "entry %d has index %lu (npiglocal %d)\n", i, index, halo_pman->NumPart);
        halo_pman->Base[index].TargetTask = pi[i].targetTask;
    }
    myfree(pi);
#ifdef DEBUG
    #pragma omp parallel for
    for(i = 0; i < halo_pman->NumPart; i ++) {
        if(halo_pman->Base[i].TargetTask < 0)
            endrun(4, "TargetTask %d not changed %d! neighbours: %d %d\n", i, halo_pman->Base[i].TargetTask, halo_pman->Base[i-1].TargetTask, halo_pman->Base[i+1].TargetTask);
    }
#endif

    walltime_measure("/FOF/IO/Distribute");

    return domain_exchange(fof_sorted_layout, halo_pman, 1, NULL, halo_pman, halo_sman, 10000, Comm);
}

static int
fof_distribute_particles(struct part_manager_type * halo_pman, struct slots_manager_type * halo_sman, MPI_Comm Comm)
{
    int64_t i, NpigLocal = 0;
    int64_t GrNrMax = -1;   /* will mark particles that are not in any group */
    int64_t GrNrMaxGlobal = 0;

    /* SlotsManager Needs initializing!*/
    memcpy(halo_sman, SlotsManager, sizeof(struct slots_manager_type));

    halo_sman->Base = NULL;
    for(i = 0; i < 6; i ++) {
        halo_sman->info[i].size = 0;
        halo_sman->info[i].maxsize = 0;
    }

    int64_t atleast[6]={0};
    /* Count how many particles we have: array reductions are an openmp 4.5 feature.*/
#if (defined _OPENMP) && (_OPENMP >= 201511)
    #pragma omp parallel for reduction(+: NpigLocal, atleast)
#endif
    for(i = 0; i < PartManager->NumPart; i ++) {
        if(P[i].GrNr >= 0) {
            NpigLocal++;
            int type = P[i].Type;
            /* How many of slot type?*/
            if(type < 6 && type >= 0 && halo_sman->info[type].enabled)
                atleast[type]++;
        }
    }
    double FOFPartAllocFactor = (double) PartManager->MaxPart / PartManager->NumPart;
    halo_pman->MaxPart = NpigLocal * FOFPartAllocFactor;
    struct particle_data * halopart = (struct particle_data *) mymalloc("HaloParticle", sizeof(struct particle_data) * halo_pman->MaxPart);
    halo_pman->Base = halopart;
    halo_pman->NumPart = NpigLocal;
    halo_pman->BoxSize = PartManager->BoxSize;
    memcpy(halo_pman->CurrentParticleOffset, PartManager->CurrentParticleOffset, 3 * sizeof(PartManager->CurrentParticleOffset[0]));

    /* We leave extra space in the hope that we can avoid compacting slots in the fof exchange*/
    for(i = 0; i < 6; i ++)
        atleast[i]*= FOFPartAllocFactor;

    slots_reserve(0, atleast, halo_sman);

    NpigLocal = 0;

    for(i = 0; i < PartManager->NumPart; i ++) {
        if(P[i].GrNr < 0)
            continue;
        if(P[i].GrNr > GrNrMax)
            GrNrMax = P[i].GrNr;
        memcpy(&halo_pman->Base[NpigLocal], &P[i], sizeof(P[i]));
        struct slot_info * info = &(halo_sman->info[P[i].Type]);
        char * oldslotptr = SlotsManager->info[P[i].Type].ptr;
        if(info->enabled) {
            memcpy(info->ptr + info->size * info->elsize, oldslotptr+P[i].PI * info->elsize, info->elsize);
            halo_pman->Base[NpigLocal].PI = info->size;
            info->size++;
        }
        NpigLocal ++;
    }
    if(NpigLocal != halo_pman->NumPart)
        endrun(3, "Error in NpigLocal %ld != %ld!\n", NpigLocal, halo_pman->NumPart);
    MPI_Allreduce(&GrNrMax, &GrNrMaxGlobal, 1, MPI_INT, MPI_MAX, Comm);
    message(0, "GrNrMax before exchange is %d\n", GrNrMaxGlobal);

    if(fof_try_particle_exchange(halo_pman, halo_sman, Comm)) {
        message(1930, "Failed to exchange and write particles for the FOF. This is non-fatal, continuing\n");
        return 1;
    }

    /* Sort locally by group number*/
    qsort_openmp(halopart, halo_pman->NumPart, sizeof(struct particle_data), order_by_type_and_grnr);
    GrNrMax = -1;
    #pragma omp parallel for reduction(max: GrNrMax)
    for(i = 0; i < halo_pman->NumPart; i ++) {
        if(halopart[i].GrNr > GrNrMax)
            GrNrMax = halopart[i].GrNr;
    }

    MPI_Allreduce(&GrNrMax, &GrNrMaxGlobal, 1, MPI_INT, MPI_MAX, Comm);
    message(0, "GrNrMax after exchange is %d\n", GrNrMaxGlobal);
    return 0;
}

static void build_buffer_fof(FOFGroups * fof, BigArray * array, IOTableEntry * ent, struct conversions * conv) {

    int64_t npartLocal = fof->Ngroups;

    petaio_alloc_buffer(array, ent, npartLocal);
    /* fill the buffer */
    char * p = (char *) array->data;
    int i;
    for(i = 0; i < fof->Ngroups; i ++) {
        ent->getter(i, p, fof->Group, NULL, conv);
        p += array->strides[0];
    }
}

static void fof_write_header(BigFile * bf, int64_t TotNgroups, const double atime, const double * MassTable, Cosmology * CP, MPI_Comm Comm) {
    BigBlock bh;
    if(0 != big_file_mpi_create_block(bf, &bh, "Header", NULL, 0, 0, 0, Comm)) {
        endrun(0, "Failed to create header\n");
    }
    int i;
    int k;
    int64_t npartLocal[6];
    int64_t npartTotal[6];

    for (k = 0; k < 6; k ++) {
        npartLocal[k] = 0;
    }
#if (defined _OPENMP) && (_OPENMP >= 201511)
    #pragma omp parallel for reduction(+: npartLocal)
#endif
    for (i = 0; i < PartManager->NumPart; i ++) {
        if(P[i].GrNr < 0) continue; /* skip those not in groups */
        npartLocal[P[i].Type] ++;
    }

    MPI_Allreduce(npartLocal, npartTotal, 6, MPI_INT64, MPI_SUM, Comm);

    /* conversion from peculiar velocity to RSD */
    const double hubble = hubble_function(CP, atime);
    double RSD = 1.0 / (atime * hubble);

    int pecvel = GetUsePeculiarVelocity();
    if(!pecvel) {
        RSD /= atime; /* Conversion from internal velocity to RSD */
    }
    big_block_set_attr(&bh, "NumPartInGroupTotal", npartTotal, "u8", 6);
    big_block_set_attr(&bh, "NumFOFGroupsTotal", &TotNgroups, "u8", 1);
    big_block_set_attr(&bh, "RSDFactor", &RSD, "f8", 1);
    big_block_set_attr(&bh, "MassTable", MassTable, "f8", 6);
    big_block_set_attr(&bh, "Time", &atime, "f8", 1);
    big_block_set_attr(&bh, "BoxSize", &PartManager->BoxSize, "f8", 1);
    big_block_set_attr(&bh, "OmegaLambda", &CP->OmegaLambda, "f8", 1);
    big_block_set_attr(&bh, "Omega0", &CP->Omega0, "f8", 1);
    big_block_set_attr(&bh, "HubbleParam", &CP->HubbleParam, "f8", 1);
    big_block_set_attr(&bh, "CMBTemperature", &CP->CMBTemperature, "f8", 1);
    big_block_set_attr(&bh, "OmegaBaryon", &CP->OmegaBaryon, "f8", 1);
    big_block_set_attr(&bh, "UsePeculiarVelocity", &pecvel, "i4", 1);
    big_block_mpi_close(&bh, Comm);
}


#define SIMPLE_PROPERTY_FOF(name, field, type, items) \
    SIMPLE_GETTER(GT ## name , field, type, items, struct Group ) \
    SIMPLE_SETTER(ST ## name , field, type, items, struct Group) \

SIMPLE_PROPERTY_FOF(GroupID, base.GrNr, uint32_t, 1)
SIMPLE_PROPERTY_FOF(MinID, base.MinID, uint64_t, 1)
SIMPLE_PROPERTY_FOF(Imom, Imom[0][0], float, 9)
/* FIXME: set Jmom to use peculiar velocity */
SIMPLE_PROPERTY_FOF(Jmom, Jmom[0], float, 3)

static void GTFirstPos(int i, float * out, void * baseptr, void * smanptr, const struct conversions * params) {
    /* Remove the particle offset before saving*/
    struct Group * grp = (struct Group *) baseptr;
    int d;
    for(d = 0; d < 3; d ++) {
        out[d] = grp[i].base.FirstPos[d] - PartManager->CurrentParticleOffset[d];
        while(out[d] > PartManager->BoxSize) out[d] -= PartManager->BoxSize;
        while(out[d] <= 0) out[d] += PartManager->BoxSize;
    }
}

static void STFirstPos(int i, float * out, void * baseptr, void * smanptr, const struct conversions * params) {
    int d;
    struct Group * grp = (struct Group *) baseptr;
    for(d = 0; d < 3; d ++) {
        grp->base.FirstPos[i] = out[d];
    }
}

static void GTMassCenterPosition(int i, double * out, void * baseptr, void * smanptr, const struct conversions * params) {
    /* Remove the particle offset before saving*/
    struct Group * grp = (struct Group *) baseptr;
    int d;
    for(d = 0; d < 3; d ++) {
        out[d] = grp[i].CM[d] - PartManager->CurrentParticleOffset[d];
        while(out[d] > PartManager->BoxSize) out[d] -= PartManager->BoxSize;
        while(out[d] <= 0) out[d] += PartManager->BoxSize;
    }
}

static void STMassCenterPosition(int i, double * out, void * baseptr, void * smanptr, const struct conversions * params) {
    int d;
    struct Group * grp = (struct Group *) baseptr;
    for(d = 0; d < 3; d ++) {
        grp->CM[d] = out[d];
    }
}

static void GTMassCenterVelocity(int i, float * out, void * baseptr, void * slotptr, const struct conversions * params) {
    double fac;
    struct Group * Group = (struct Group *) baseptr;
    if (GetUsePeculiarVelocity()) {
        fac = 1.0 / params->atime;
    } else {
        fac = 1.0;
    }

    int d;
    for(d = 0; d < 3; d ++) {
        out[d] = fac * Group[i].Vel[d];
    }
}
SIMPLE_PROPERTY_FOF(Mass, Mass, float, 1)
SIMPLE_PROPERTY_FOF(MassByType, MassType[0], float, 6)
SIMPLE_PROPERTY_FOF(LengthByType, LenType[0], uint32_t , 6)
SIMPLE_PROPERTY_FOF(StarFormationRate, Sfr, float, 1)
SIMPLE_PROPERTY_FOF(GasMetalMass, GasMetalMass, float, 1)
SIMPLE_PROPERTY_FOF(StellarMetalMass, StellarMetalMass, float, 1)
SIMPLE_PROPERTY_FOF(GasMetalElemMass, GasMetalElemMass[0], float, NMETALS)
SIMPLE_PROPERTY_FOF(StellarMetalElemMass, StellarMetalElemMass[0], float, NMETALS)
SIMPLE_PROPERTY_FOF(BlackholeMass, BH_Mass, float, 1)
SIMPLE_PROPERTY_FOF(BlackholeAccretionRate, BH_Mdot, float, 1)
SIMPLE_PROPERTY_FOF(MassHeIonized, MassHeIonized, float, 1)

static void fof_register_io_blocks(int MetalReturnOn, int BlackholeOn, struct IOTable * IOTable) {
    IOTable->used = 0;
    IOTable->allocated = 100;
    /* Allocate high so we can do a domain exchange,
     * potentially increasing the slots, around this*/
    IOTable->ent = (struct IOTableEntry *) mymalloc2("IOTable", IOTable->allocated* sizeof(IOTableEntry));

    IO_REG(GroupID, "u4", 1, PTYPE_FOF_GROUP, IOTable);
    IO_REG(Mass, "f4", 1, PTYPE_FOF_GROUP, IOTable);
    IO_REG(MassCenterPosition, "f8", 3, PTYPE_FOF_GROUP, IOTable);
    IO_REG(FirstPos, "f4", 3, PTYPE_FOF_GROUP, IOTable);
    IO_REG(MinID, "u8", 1, PTYPE_FOF_GROUP, IOTable);
    IO_REG(Imom, "f4", 9, PTYPE_FOF_GROUP, IOTable);
    IO_REG(Jmom, "f4", 3, PTYPE_FOF_GROUP, IOTable);
    IO_REG_WRONLY(MassCenterVelocity, "f4", 3, PTYPE_FOF_GROUP, IOTable);
    IO_REG(LengthByType, "u4", 6, PTYPE_FOF_GROUP, IOTable);
    IO_REG(MassByType, "f4", 6, PTYPE_FOF_GROUP, IOTable);
    IO_REG(MassHeIonized, "f4", 1, PTYPE_FOF_GROUP, IOTable);
    /* Zero if star formation is not on*/
    IO_REG(StarFormationRate, "f4", 1, PTYPE_FOF_GROUP, IOTable);
    if(MetalReturnOn) {
        IO_REG(GasMetalMass, "f4", 1, PTYPE_FOF_GROUP, IOTable);
        IO_REG(StellarMetalMass, "f4", 1, PTYPE_FOF_GROUP, IOTable);
        /* Zero if metal return is not on*/
        IO_REG(GasMetalElemMass, "f4", NMETALS, PTYPE_FOF_GROUP, IOTable);
        IO_REG(StellarMetalElemMass, "f4", NMETALS, PTYPE_FOF_GROUP, IOTable);
    }
    /* Zero if black hole is not on*/
    if(BlackholeOn) {
        IO_REG(BlackholeMass, "f4", 1, PTYPE_FOF_GROUP, IOTable);
        IO_REG(BlackholeAccretionRate, "f4", 1, PTYPE_FOF_GROUP, IOTable);
    }
}
