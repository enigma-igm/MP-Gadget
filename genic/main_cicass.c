#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#include <bigfile-mpi.h>
#include <libgenic/allvars.h>
#include <libgenic/proto.h>
#include <libgenic/thermal.h>
#include <libgadget/walltime.h>
#include <libgadget/petapm.h>
#include <libgadget/utils.h>

struct io_header_1
{
    unsigned int npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
    double mass[6];          /*!< mass[1] gives the particle mass */
    double time;             /*!< time (=cosmological scale factor) of snapshot */
    double redshift;         /*!< redshift of snapshot */
    int flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
    int flag_feedback;  /*!< flags whether feedback from star formation is included */
    unsigned int npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                 the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
    int flag_cooling;   /*!< flags whether radiative cooling is included */
    int num_files;      /*!< determines the number of files that are used for a snapshot */
    double BoxSize;          /*!< Simulation box size (in code units) */
    double Omega0;           /*!< matter density */
    double OmegaLambda;      /*!< vacuum energy density */
    double HubbleParam;      /*!< little 'h' */
    int flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
    int flag_metals;         /*!< flags whether metal enrichment is included */
    int hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
    char fill[84];        /*!< fills to 256 Bytes */
}
header1;

int main(int argc, char **argv)
{
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    init_endrun();
    
    if(argc < 2)
        endrun(0,"Please pass a parameter file.\n");
    
    tamalloc_init();
    
    read_parameterfile(argv[1]);
    
    // FBD HACK:
    // READ IN GADGET FORMAT CICASS OUTPUT HERE.
    // CODE COPIED FROM VOLKER'S io_input.c: load_snapshot_gas
    FILE *fd;
    char *fname = All2.PowerP.FileWithInputSpectrum;//"/mnt/quasar/fdavies/cdm_1Mpc_128_v30_IC.dat"; // CICASS IC filename
    int i, j, k, dummy, dummy2;
    int n;
    
    if(!(fd = fopen(fname, "r")))
    {
        printf("can't open file `%s`\n", fname);
        exit(1);
    }
    
    printf("reading `%s' ...\n", fname);
    fflush(stdout);
    
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&header1, sizeof(header1), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
    
    int NumPart = header1.npart[0];
    int NumPartSph = NumPart;
    printf("%i %i\n",NumPart,NumPartSph);
    
    float **dm_pos = (float**)malloc((size_t)NumPart*sizeof(float*));
    float **dm_vel = (float**)malloc((size_t)NumPart*sizeof(float*));
    float **gas_pos = (float**)malloc((size_t)NumPartSph*sizeof(float*));
    float **gas_vel = (float**)malloc((size_t)NumPartSph*sizeof(float*));
    for (i=0; i<NumPart; i++) {
        dm_pos[i] = (float*)malloc((size_t)3*sizeof(float));
        dm_vel[i] = (float*)malloc((size_t)3*sizeof(float));
        gas_pos[i] = (float*)malloc((size_t)3*sizeof(float));
        gas_vel[i] = (float*)malloc((size_t)3*sizeof(float));
    }
    
    /* Pos */
    fread(&dummy, sizeof(dummy), 1, fd);
    //printf("dummy position read");
    for(k = 0; k < 6; k++)
    {
        if(k == 0)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(gas_pos[n], sizeof(float), 3, fd);
            }
            printf("Read gadget Pos1 %f %f %f\n",gas_pos[1][0],gas_pos[1][1],gas_pos[1][2]);
        }
        else if (k == 1)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(dm_pos[n], sizeof(float), 3, fd);
            }
            printf("Read gadget Pos2 %f %f %f\n",dm_pos[1][0],dm_pos[1][1],dm_pos[1][2]);
        }
        else {
            fseek(fd, header1.npart[k] * sizeof(float) * 3, SEEK_CUR); // moves on in the file without reading
        }
    }
    fread(&dummy2, sizeof(dummy2), 1, fd);
    if(dummy2 != dummy)
    {
        printf("Pos: I/O Panic! dummy=%d dumy2=%d\n",dummy,dummy2);
        exit(1);
    }
    printf("Read gadget Pos done\n");
    
    /* Vel */
    fread(&dummy, sizeof(dummy), 1, fd);
    for(k = 0; k < 6; k++)
    {
        if(k == 0)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(gas_vel[n], sizeof(float), 3, fd);
            }
            printf("Read gadget Vel1\n");
        }
        else if (k == 1)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(dm_vel[n], sizeof(float), 3, fd);
            }
            printf("Read gadget Vel2\n");
        }
        else
            fseek(fd, header1.npart[k] * sizeof(float) * 3, SEEK_CUR);
    }
    fread(&dummy2, sizeof(dummy2), 1, fd);
    if(dummy2 != dummy)
    {
        printf("Vel: I/O Panic! dummy=%d dumy2=%d\n",dummy,dummy2);
        exit(1);
    }
    printf("Read gadget Vel done\n");
    
    fclose(fd);
    
    
    mymalloc_init(All.MaxMemSizePerNode);
    
    walltime_init(&All.CT);
    
    int64_t TotNumPart = (int64_t) All2.Ngrid*All2.Ngrid*All2.Ngrid;
    int64_t TotNumPartGas = (int64_t) All2.ProduceGas*All2.NgridGas*All2.NgridGas*All2.NgridGas;
    
    init_cosmology(&All.CP, All.TimeIC);
    
    //init_powerspectrum(ThisTask, All.TimeIC, All.UnitLength_in_cm, &All.CP, &All2.PowerP);
    All.NumThreads = omp_get_max_threads();
    
    petapm_init(All.BoxSize, All.Nmesh, All.NumThreads);
    /*Initialise particle spacings*/
    const double meanspacing = All.BoxSize / DMAX(All2.Ngrid, All2.NgridGas);
    const double shift_gas = -All2.ProduceGas * 0.5 * (All.CP.Omega0 - All.CP.OmegaBaryon) / All.CP.Omega0 * meanspacing;
    double shift_dm = All2.ProduceGas * 0.5 * All.CP.OmegaBaryon / All.CP.Omega0 * meanspacing;
    double shift_nu = 0;
    if(!All2.ProduceGas && All2.NGridNu > 0) {
        double OmegaNu = get_omega_nu(&All.CP.ONu, 1);
        shift_nu = -0.5 * (All.CP.Omega0 - OmegaNu) / All.CP.Omega0 * meanspacing;
        shift_dm = 0.5 * OmegaNu / All.CP.Omega0 * meanspacing;
    }
    
    /*Write the header*/
    char buf[4096];
    snprintf(buf, 4096, "%s/%s", All.OutputDir, All.InitCondFile);
    BigFile bf;
    if(0 != big_file_mpi_create(&bf, buf, MPI_COMM_WORLD)) {
        endrun(0, "%s\n", big_file_get_error_message());
    }
    /*Massive neutrinos*/
    
    const int64_t TotNu = (int64_t) All2.NGridNu*All2.NGridNu*All2.NGridNu;
    double total_nufrac = 0;
    struct thermalvel nu_therm;
    saveheader(&bf, TotNumPart, TotNumPartGas, TotNu, total_nufrac);
    
    /*Save the transfer functions*/
    //save_all_transfer_tables(&bf, ThisTask);
    
    /*Use 'total' (CDM + baryon) transfer function
     * unless DifferentTransferFunctions are on.
     */
    enum TransferType DMType = DELTA_CB, GasType = DELTA_CB, NuType = DELTA_NU;
    if(All2.ProduceGas && All2.PowerP.DifferentTransferFunctions) {
        DMType = DELTA_CDM;
        GasType = DELTA_BAR;
    }
    
    /*First compute and write CDM*/
    double mass[6] = {0};
    /*Can neglect neutrinos since this only matters for the glass force.*/
    compute_mass(mass, TotNumPart, TotNumPartGas, 0, 0);
    /*Not used*/
    int size[3], offset[3];
    int NumPartCDM = get_size_offset(size, offset, All2.Ngrid);
    int NumPartGas = get_size_offset(size, offset, All2.NgridGas);
    /*Space for both CDM and baryons*/
    struct ic_part_data * ICP = (struct ic_part_data *) mymalloc("PartTable", (NumPartCDM + All2.ProduceGas * NumPartGas)*sizeof(struct ic_part_data));
    
    for (i = 0; i < NumPartCDM; i++) {
        for (j = 0; j < 3; j++) {
            ICP[i].Pos[j] = dm_pos[i][j];
            ICP[i].Vel[j] = dm_vel[i][j]*sqrt(All.TimeIC);
        }
    }
    
    if(All2.WDM_therm_mass > 0){
        //int i;
        double v_th = WDM_V0(All.TimeIC, All2.WDM_therm_mass, All.CP.Omega0 - All.CP.OmegaBaryon - get_omega_nu(&All.CP.ONu, 1), All.CP.HubbleParam, All.UnitVelocity_in_cm_per_s);
        printf("Adding residual WDM thermal velocity of %E km/s\n",v_th);
        if(!All.IO.UsePeculiarVelocity)
            v_th /= sqrt(All.TimeIC);
        struct thermalvel WDM;
        init_thermalvel(&WDM, v_th, 10000/v_th, 0);
        unsigned int * seedtable = init_rng(All2.Seed+1,All2.Ngrid);
        gsl_rng * g_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
        /*Seed the random number table with the Id.*/
        gsl_rng_set(g_rng, seedtable[0]);
        
        for(i = 0; i < NumPartCDM; i++) {
            /*Find the slab, and reseed if it has zero z rank*/
            if(i % All2.Ngrid == 0) {
                uint64_t id = id_offset_from_index(i, All2.Ngrid);
                /*Seed the random number table with x,y index.*/
                gsl_rng_set(g_rng, seedtable[id / All2.Ngrid]);
            }
            add_thermal_speeds(&WDM, g_rng, ICP[i].Vel);
        }
        gsl_rng_free(g_rng);
        myfree(seedtable);
    }
    
    write_particle_data(1, &bf, 0, All2.Ngrid, ICP, NumPartCDM);

    for (i = NumPartCDM; i < NumPartCDM+NumPartGas; i++) {
        for (j = 0; j < 3; j++) {
            ICP[i].Pos[j] = gas_pos[i-NumPartCDM][j];
            ICP[i].Vel[j] = gas_vel[i-NumPartCDM][j]*sqrt(All.TimeIC);
        }
    }
    
    write_particle_data(0, &bf, TotNumPart, All2.NgridGas, ICP+NumPartCDM, NumPartGas);

    myfree(ICP);
    
    big_file_mpi_close(&bf, MPI_COMM_WORLD);
    
    walltime_summary(0, MPI_COMM_WORLD);
    walltime_report(stdout, 0, MPI_COMM_WORLD);
    
    message(0, "IC's copied from CICASS.\n");
    message(0, "Initial scale factor = %g\n", All.TimeIC);
    
    MPI_Finalize();        /* clean up & finalize MPI */
    return 0;
}
