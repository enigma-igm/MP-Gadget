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
#include <libgadget/physconst.h>
#include <libgadget/petapm.h>
#include <libgadget/utils.h>
#include <libgadget/partmanager.h>

#define GLASS_SEED_HASH(seed) ((seed) * 9999721L)

static void print_spec(int ThisTask, const int Ngrid, struct genic_config All2, Cosmology * CP);

// Structure of CICASS IC binary output
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
  int thread_provided, ThisTask, NumTask;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_provided);
  if(thread_provided != MPI_THREAD_FUNNELED)
    message(1, "MPI_Init_thread returned %d != MPI_THREAD_FUNNELED\n", thread_provided);

  if(argc < 2)
    endrun(0,"Please pass a parameter file.\n");

  tamalloc_init();

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NumTask);
  printf("%i %i\n",ThisTask,NumTask);

  /* Genic Specific configuration structure*/
  struct genic_config All2 = {0};

  Cosmology CP;
  int ShowBacktrace;
  double MaxMemSizePerNode;
  read_parameterfile(argv[1], &All2, &ShowBacktrace, &MaxMemSizePerNode, &CP);

  mymalloc_init(MaxMemSizePerNode);

  init_endrun(ShowBacktrace);
    
  // END OF CICASS SECTION

  struct ClockTable Clocks;
  walltime_init(&Clocks);

  int64_t TotNumPart = (int64_t) All2.Ngrid*All2.Ngrid*All2.Ngrid;
  int64_t TotNumPartGas = (int64_t) All2.ProduceGas*All2.NgridGas*All2.NgridGas*All2.NgridGas;

  init_cosmology(&CP, All2.TimeIC);

  //MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  //init_powerspectrum(ThisTask, All2.TimeIC, All2.UnitLength_in_cm, &CP, &All2.PowerP);

  petapm_module_init(omp_get_max_threads());

  /*Initialise particle spacings*/
  const double meanspacing = All2.BoxSize / DMAX(All2.Ngrid, All2.NgridGas);
  double shift_gas = -All2.ProduceGas * 0.5 * (CP.Omega0 - CP.OmegaBaryon) / CP.Omega0 * meanspacing;
  double shift_dm = All2.ProduceGas * 0.5 * CP.OmegaBaryon / CP.Omega0 * meanspacing;
  
  double shift_nu = 0;
  if(!All2.ProduceGas && All2.NGridNu > 0) {
      double OmegaNu = get_omega_nu(&CP.ONu, 1);
      shift_nu = -0.5 * (CP.Omega0 - OmegaNu) / CP.Omega0 * meanspacing;
      shift_dm = 0.5 * OmegaNu / CP.Omega0 * meanspacing;
  }
    
  if(All2.PrePosGridCenter){
      shift_dm += 0.5 * meanspacing;
      shift_gas += 0.5 * meanspacing;
      shift_nu += 0.5 * meanspacing;
  }

  /*Write the header*/
  char buf[4096];
  snprintf(buf, 4096, "%s/%s", All2.OutputDir, All2.InitCondFile);
  BigFile bf;
  if(0 != big_file_mpi_create(&bf, buf, MPI_COMM_WORLD)) {
      endrun(0, "%s\n", big_file_get_error_message());
  }
  /*Massive neutrinos*/

  const int64_t TotNu = (int64_t) All2.NGridNu*All2.NGridNu*All2.NGridNu;
  double total_nufrac = 0;
  struct thermalvel nu_therm;
  if(TotNu > 0) {
    const double kBMNu = 3*CP.ONu.kBtnu / (CP.MNu[0]+CP.MNu[1]+CP.MNu[2]);
    double v_th = NU_V0(All2.TimeIC, kBMNu, All2.UnitVelocity_in_cm_per_s);
    if(!All2.UsePeculiarVelocity)
        v_th /= sqrt(All2.TimeIC);
    total_nufrac = init_thermalvel(&nu_therm, v_th, All2.Max_nuvel/v_th, 0);
    message(0,"F-D velocity scale: %g. Max particle vel: %g. Fraction of mass in particles: %g\n",v_th*sqrt(All2.TimeIC), All2.Max_nuvel*sqrt(All2.TimeIC), total_nufrac);
  }
  saveheader(&bf, TotNumPart, TotNumPartGas, TotNu, total_nufrac, All2.BoxSize, &CP, All2);

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
  PetaPM pm[1];

  double UnitTime_in_s = All2.UnitLength_in_cm / All2.UnitVelocity_in_cm_per_s;
  double Grav = GRAVITY / pow(All2.UnitLength_in_cm, 3) * All2.UnitMass_in_g * pow(UnitTime_in_s, 2);

  petapm_init(pm, All2.BoxSize, 0, All2.Nmesh, Grav, MPI_COMM_WORLD);

  /*First compute and write CDM*/
  double mass[6] = {0};
  /*Can neglect neutrinos since this only matters for the glass force.*/
  compute_mass(mass, TotNumPart, TotNumPartGas, 0, 0, All2.BoxSize, &CP, All2);
  /*Not used*/
  IDGenerator idgen_cdm[1];
  IDGenerator idgen_gas[1];

  idgen_init(idgen_cdm, pm, All2.Ngrid, All2.BoxSize);
  idgen_init(idgen_gas, pm, All2.NgridGas, All2.BoxSize);

  int NumPartCDM = idgen_cdm->NumPart;
  int NumPartGas = idgen_gas->NumPart;
    
    printf("from igden: %i %i\n",NumPartCDM,NumPartGas);

  /*Space for both CDM and baryons*/
  struct ic_part_data * ICP = (struct ic_part_data *) mymalloc("PartTable", (NumPartCDM + All2.ProduceGas * NumPartGas)*sizeof(struct ic_part_data));

  /* If we have incoherent glass files, we need to store both the particle tables
   * to ensure that there are no close particle pairs*/
  /*Make the table for the CDM*/
//  if(!All2.MakeGlassCDM) {
//      setup_grid(idgen_cdm, shift_dm, mass[1], ICP);
//  } else {
//      setup_glass(idgen_cdm, pm, 0, GLASS_SEED_HASH(All2.Seed), mass[1], ICP, All2.UnitLength_in_cm, All2.OutputDir);
//  }
//
//  /*Make the table for the baryons if we need, using the second half of the memory.*/
//  if(All2.ProduceGas) {
//    if(!All2.MakeGlassGas) {
//        setup_grid(idgen_gas, shift_gas, mass[0], ICP+NumPartCDM);
//    } else {
//        setup_glass(idgen_gas, pm, 0, GLASS_SEED_HASH(All2.Seed + 1), mass[0], ICP+NumPartCDM, All2.UnitLength_in_cm, All2.OutputDir);
//    }
//    /*Do coherent glass evolution to avoid close pairs*/
//    if(All2.MakeGlassGas || All2.MakeGlassCDM)
//        glass_evolve(pm, 14, "powerspectrum-glass-tot", ICP, NumPartCDM+NumPartGas, All2.UnitLength_in_cm, All2.OutputDir);
//  }
    
    
    // FBD HACK:
    // READ IN GADGET FORMAT CICASS OUTPUT HERE.
    // CODE COPIED FROM VOLKER'S io_input.c: load_snapshot_gas
    FILE *fd;
    char *fname = All2.PowerP.FileWithInputSpectrum;//"/mnt/quasar/fdavies/cdm_1Mpc_128_v30_IC.dat"; // CICASS IC filename
    int i, j, k, dummy, dummy2;
    int n;
    // only read in the file on the primary task
    
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
    int NumPartTask = NumPart/NumTask;
    float skip[3];
    float invec[3];
    printf("%i %i\n",NumPart,NumPartSph);
    
    /* Pos */
    fread(&dummy, sizeof(dummy), 1, fd);
    //printf("dummy position read");
    for(k = 0; k < 6; k++)
    {
        if(k == 0)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(skip, sizeof(float), 3, fd);
            }
            //for (n=0; n<8; n++) {
            //    printf("Read gadget Pos1 %f %f %f\n",gas_pos[n][0],gas_pos[n][1],gas_pos[n][2]);
            //}
        }
        else if (k == 1)
        {
            for (n = 0; n < ThisTask*NumPartTask; n++) {
                fread(skip, sizeof(float), 3, fd);
            }
            for(n = ThisTask*NumPartTask; n < (ThisTask+1)*NumPartTask; n++)
            {
                fread(invec, sizeof(float), 3, fd);
                ICP[n-ThisTask*NumPartTask].Pos[0] = invec[0]; ICP[n-ThisTask*NumPartTask].Pos[1] = invec[1]; ICP[n-ThisTask*NumPartTask].Pos[2] = invec[2];
            }
            for (n = (ThisTask+1)*NumPartTask; n < NumPart; n++) {
                fread(skip, sizeof(float), 3, fd);
            }
            printf("Read gadget DM Pos %f %f %f\n",ICP[1].Pos[0],ICP[1].Pos[1],ICP[1].Pos[2]);
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
    printf("Read gadget DM Pos done\n");
    
    /* Vel */
    fread(&dummy, sizeof(dummy), 1, fd);
    for(k = 0; k < 6; k++)
    {
        if(k == 0)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(skip, sizeof(float), 3, fd);
            }
            //printf("Read gadget Vel1 %f %f %f\n",gas_vel[1][0],gas_vel[1][1],gas_vel[1][1]);
        }
        else if (k == 1)
        {
            for (n = 0; n < ThisTask*NumPartTask; n++) {
                fread(skip, sizeof(float), 3, fd);
            }
            for(n = ThisTask*NumPartTask; n < (ThisTask+1)*NumPartTask; n++)
            {
                fread(invec, sizeof(float), 3, fd);
                ICP[n-ThisTask*NumPartTask].Vel[0] = invec[0]*sqrt(All2.TimeIC);
                ICP[n-ThisTask*NumPartTask].Vel[1] = invec[1]*sqrt(All2.TimeIC);
                ICP[n-ThisTask*NumPartTask].Vel[2] = invec[2]*sqrt(All2.TimeIC);
            }
            for (n = (ThisTask+1)*NumPartTask; n < NumPart; n++) {
                fread(skip, sizeof(float), 3, fd);
            }
            printf("Read gadget DM Vel %f %f %f\n",ICP[1].Vel[0],ICP[1].Vel[1],ICP[1].Vel[2]);
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
    printf("Read gadget DM Vel done\n");
    
    fclose(fd);

  /*Write initial positions into ICP struct (for CDM and gas)*/
  //int j,k;
  for(j=0; j<NumPartCDM+NumPartGas; j++)
      for(k=0; k<3; k++)
          ICP[j].PrePos[k] = ICP[j].Pos[k];

  if(NumPartCDM > 0) {
    //displacement_fields(pm, DMType, ICP, NumPartCDM, &CP, All2);

    /*Add a thermal velocity to WDM particles*/
    if(All2.WDM_therm_mass > 0){
        //int i;
        double v_th = WDM_V0(All2.TimeIC, All2.WDM_therm_mass, CP.Omega0 - CP.OmegaBaryon - get_omega_nu(&CP.ONu, 1), CP.HubbleParam, All2.UnitVelocity_in_cm_per_s);
        if(!All2.UsePeculiarVelocity)
           v_th /= sqrt(All2.TimeIC);
        struct thermalvel WDM;
        init_thermalvel(&WDM, v_th, 10000/v_th, 0);
        unsigned int * seedtable = init_rng(All2.Seed+1,All2.Ngrid);
        gsl_rng * g_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
        /*Seed the random number table with the Id.*/
        gsl_rng_set(g_rng, seedtable[0]);

        for(i = 0; i < NumPartCDM; i++) {
             /*Find the slab, and reseed if it has zero z rank*/
             if(i % All2.Ngrid == 0) {
                  uint64_t id = idgen_create_id_from_index(idgen_cdm, i);
                  /*Seed the random number table with x,y index.*/
                  gsl_rng_set(g_rng, seedtable[id / All2.Ngrid]);
             }
             add_thermal_speeds(&WDM, g_rng, ICP[i].Vel);
        }
        gsl_rng_free(g_rng);
        myfree(seedtable);
    }
      
    MPI_Barrier(MPI_COMM_WORLD);
    write_particle_data(idgen_cdm, 1, &bf, 0, All2.SavePrePos, All2.NumFiles, All2.NumWriters, ICP);
  }

  /*Now make the gas if required*/
  if(All2.ProduceGas) {
      
      // Re-open the file
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
      
      
          /* Pos */
          fread(&dummy, sizeof(dummy), 1, fd);
          //printf("dummy position read");
          for(k = 0; k < 6; k++)
          {
              if(k == 0)
              {
                  for (n = 0; n < ThisTask*NumPartTask; n++) {
                      fread(skip, sizeof(float), 3, fd);
                  }
                  for(n = ThisTask*NumPartTask; n < (ThisTask+1)*NumPartTask; n++)
                  {
                      fread(invec, sizeof(float), 3, fd);
                      ICP[n-ThisTask*NumPartTask].Pos[0] = invec[0]; ICP[n-ThisTask*NumPartTask].Pos[1] = invec[1]; ICP[n-ThisTask*NumPartTask].Pos[2] = invec[2];
                  }
                  for (n = (ThisTask+1)*NumPartTask; n < NumPart; n++) {
                      fread(skip, sizeof(float), 3, fd);
                  }
                  printf("Read gadget gas Pos %f %f %f\n",ICP[1].Pos[0],ICP[1].Pos[1],ICP[1].Pos[2]);
              }
              else if (k == 1)
              {
                  for(n = 0; n < header1.npart[k]; n++)
                  {
                      fread(skip, sizeof(float), 3, fd);
                  }
                  
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
          printf("Read gadget gas Pos done\n");
          
          /* Vel */
          fread(&dummy, sizeof(dummy), 1, fd);
          for(k = 0; k < 6; k++)
          {
              if(k == 0)
              {
                  for (n = 0; n < ThisTask*NumPartTask; n++) {
                      fread(skip, sizeof(float), 3, fd);
                  }
                  for(n = ThisTask*NumPartTask; n < (ThisTask+1)*NumPartTask; n++)
                  {
                      fread(invec, sizeof(float), 3, fd);
                      ICP[n-ThisTask*NumPartTask].Vel[0] = invec[0]*sqrt(All2.TimeIC);
                      ICP[n-ThisTask*NumPartTask].Vel[1] = invec[1]*sqrt(All2.TimeIC);
                      ICP[n-ThisTask*NumPartTask].Vel[2] = invec[2]*sqrt(All2.TimeIC);
                  }
                  for (n = (ThisTask+1)*NumPartTask; n < NumPart; n++) {
                      fread(skip, sizeof(float), 3, fd);
                  }
                  printf("Read gadget gas Vel %f %f %f\n",ICP[1].Vel[0],ICP[1].Vel[1],ICP[1].Vel[2]);
                  //printf("Read gadget Vel1 %f %f %f\n",gas_vel[1][0],gas_vel[1][1],gas_vel[1][1]);
              }
              else if (k == 1)
              {
                  for(n = 0; n < header1.npart[k]; n++)
                  {
                      fread(skip, sizeof(float), 3, fd);
                  }
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
          printf("Read gadget gas Vel done\n");
          
          fclose(fd);
      
    //displacement_fields(pm, GasType, ICP+NumPartCDM, NumPartGas, &CP, All2);

      MPI_Barrier(MPI_COMM_WORLD);
      
    write_particle_data(idgen_gas, 0, &bf, TotNumPart, All2.SavePrePos, All2.NumFiles, All2.NumWriters, ICP+NumPartCDM);
  }
  myfree(ICP);

  /*Now add random velocity neutrino particles*/
  if(All2.NGridNu > 0) {
      int i;
      IDGenerator idgen_nu[1];
      idgen_init(idgen_nu, pm, All2.NGridNu, All2.BoxSize);

      int NumPartNu = idgen_nu->NumPart;
      ICP = (struct ic_part_data *) mymalloc("PartTable", NumPartNu*sizeof(struct ic_part_data));

      NumPartNu = setup_grid(idgen_nu, shift_nu, mass[2], ICP);

	  /*Write initial positions into ICP struct (for neutrinos)*/
	  for(j=0; j<NumPartNu; j++)
		  for(k=0; k<3; k++)
		      ICP[j].PrePos[k] = ICP[j].Pos[k];

      displacement_fields(pm, NuType, ICP, NumPartNu, &CP, All2);
      unsigned int * seedtable = init_rng(All2.Seed+2,All2.NGridNu);
      gsl_rng * g_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
      /*Just in case*/
      gsl_rng_set(g_rng, seedtable[0]);
      for(i = 0; i < NumPartNu; i++) {
           /*Find the slab, and reseed if it has zero z rank*/
           if(i % All2.NGridNu == 0) {
                uint64_t id = idgen_create_id_from_index(idgen_nu, i);
                /*Seed the random number table with x,y index.*/
                gsl_rng_set(g_rng, seedtable[id / All2.NGridNu]);
           }
           add_thermal_speeds(&nu_therm, g_rng, ICP[i].Vel);
      }
      gsl_rng_free(g_rng);
      myfree(seedtable);

      write_particle_data(idgen_nu, 2, &bf, TotNumPart+TotNumPartGas, All2.SavePrePos, All2.NumFiles, All2.NumWriters, ICP);
      myfree(ICP);
  }

  petapm_destroy(pm);
  big_file_mpi_close(&bf, MPI_COMM_WORLD);

  walltime_summary(0, MPI_COMM_WORLD);
  walltime_report(stdout, 0, MPI_COMM_WORLD);

  message(0, "IC's copied from CICASS.\n");
  message(0, "Initial scale factor = %g\n", All2.TimeIC);
    
  //print_spec(ThisTask, All2.Ngrid, All2, &CP);

  MPI_Finalize();		/* clean up & finalize MPI */
  return 0;
}

void print_spec(int ThisTask, const int Ngrid, struct genic_config All2, Cosmology * CP)
{
  if(ThisTask == 0)
    {
      double k, kstart, kend, DDD;
      char buf[1000];
      FILE *fd;

      sprintf(buf, "%s/inputspec_%s.txt", All2.OutputDir, All2.InitCondFile);

      fd = fopen(buf, "w");
      if (fd == NULL) {
        message(1, "Failed to create powerspec file at:%s\n", buf);
        return;
      }
      DDD = GrowthFactor(CP, All2.TimeIC, 1.0);

      fprintf(fd, "# %12g %12g\n", 1/All2.TimeIC-1, DDD);
      /* print actual starting redshift and linear growth factor for this cosmology */
      kstart = 2 * M_PI / (2*All2.BoxSize * (CM_PER_MPC / All2.UnitLength_in_cm));	/* 2x box size Mpc/h */
      kend = 2 * M_PI / (All2.BoxSize/(8*Ngrid) * (CM_PER_MPC / All2.UnitLength_in_cm));	/* 1/8 mean spacing Mpc/h */

      message(1,"kstart=%lg kend=%lg\n",kstart,kend);

      for(k = kstart; k < kend; k *= 1.025)
	  {
	    double po = pow(DeltaSpec(k, DELTA_TOT),2);
	    fprintf(fd, "%12g %12g\n", k, po);
	  }
      fclose(fd);
    }
}

