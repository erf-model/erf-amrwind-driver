#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF.H>

#include <MultiBlockContainer.H>

std::string inputs_name;

using namespace amrex;

/**
 * Function to set the refine_grid_layout flags to (1,1,0) by default
 * since the ERF default is different from the amrex default (1,1,1)
 * Also set max_grid_size to very large since the only reason for
 * chopping grids is if Nprocs > Ngrids
*/
void add_par () {
   ParmParse pp("amr");

   // Set the refine_grid_layout flags to (1,1,0) by default
   pp.add("refine_grid_layout_x",1);
   pp.add("refine_grid_layout_y",1);
   pp.add("refine_grid_layout_z",0);

   // n_proper is the minimum number of coarse cells between coarse-fine boundaries
   // between levels (ell and ell+1) and levels (ell-1 and ell).   We want this to be
   // greater than or equal to the stencil width (a function of spatial order) divided by
   // ref_ratio (which can be 2,3 or 4).  This ensures that fillpatch at level (ell)
   // does not need to reach beyond level (ell-1). Here to be conservative we set this to 2
   // (rather than the amrex default of 1).
   pp.add("n_proper",2);

   int max_grid_size = 2048;
   pp.queryAdd("max_grid_size",max_grid_size);

   // This will set the default value of blocking_factor to be 1, but will allow
   //     the user to override it in the inputs file or on command line
   int blocking_factor = 1;
   pp.queryAdd("blocking_factor",blocking_factor);

   int n_error_buf = 0;
   pp.queryAdd("n_error_buf",n_error_buf);
}

/**
 * Main driver -- creates the ERF object, calls ERF.InitData() and ERF.Evolve()
 * Also includes the multiblock interface in the case where there is more than one ERF object
*/
int main (int argc, char* argv[])
{

#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    if (argc < 2) {
        // Print usage and exit with error code if no input file was provided.
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }

    // Look for "-h" or "--help" flag and print usage
    for (auto i = 1; i < argc; i++) {
        const std::string param(argv[i]);
        if ((param == "--help") || (param == "-h") || (param == "--usage")) {
            ERF::print_banner(MPI_COMM_WORLD, std::cout);
            ERF::print_usage(MPI_COMM_WORLD, std::cout);
            return 0;
        }
    }

    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        // Print usage and exit with error code if we cannot find the input file
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD, "Input file does not exist = " +
                                std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

  //  print_banner(MPI_COMM_WORLD, std::cout);
    // Check to see if the command line contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                ERF::writeBuildInfo(std::cout);
                return 0;
            }
        }
    }
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);

    // Save the inputs file name for later.
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = amrex::second();

    {
        // Vector of constructor parameters for MultiBlock
        std::vector<amrex::RealBox> rb_v;
        std::vector<int> max_level_v;
        std::vector<int> coord_v;
        std::vector<amrex::Vector<int>> n_cell_v;
        std::vector<amrex::Array<int,AMREX_SPACEDIM>> is_per_v;
        std::vector<amrex::Vector<amrex::IntVect>> ref_rat_v;
        std::vector<std::string> prefix_v;
        int max_step{1};

        // Local constructor parameters for vector
        amrex::RealBox rb;
        int max_level{0};
        int coord{0};
        amrex::Vector<int> n_cell = {1,1,1};
        amrex::Array<int,AMREX_SPACEDIM> is_per = {1,1,1};
        amrex::Vector<amrex::IntVect> ref_rat = {amrex::IntVect(1,1,1)};

        // Parse max steps for the block
        {
            ParmParse pp;
            pp.query("max_step", max_step);
        }

        // Parse data for erf constructor
        {
            ParmParse pp("erf");
            amrex::Vector<amrex::Real> lo  = {0.,0.,0.};
            amrex::Vector<amrex::Real> hi  = {0.,0.,0.};
            amrex::Vector<int> periodicity = {1,1,1};
            pp.queryarr("prob_lo",lo);
            pp.queryarr("prob_hi",hi);
            rb.setLo(lo);
            rb.setHi(hi);
            pp.query("max_level",max_level);
            pp.query("coord",coord);
            pp.queryarr("n_cell",n_cell);
            pp.queryarr("is_periodic",periodicity);
            {
                for( int i(0); i<AMREX_SPACEDIM; i++ ) is_per[i] = periodicity[i];
            }
            pp.queryarr("ref_ratio",ref_rat);

            rb_v.push_back(rb);
            max_level_v.push_back(max_level);
            coord_v.push_back(coord);
            n_cell_v.push_back(n_cell);
            is_per_v.push_back(is_per);
            ref_rat_v.push_back(ref_rat);
            prefix_v.push_back("erf");
        }

        // Construct a MultiBlockContainer
        MultiBlockContainer mbc(rb_v, max_level_v, n_cell_v,
                                coord_v, ref_rat_v, is_per_v,
                                prefix_v, max_step);

        // Initialize data
        mbc.InitializeBlocks();

        // Advance blocks a timestep
        mbc.AdvanceBlocks();

        // wallclock time
        Real end_total = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "\nTotal Time: " << end_total << '\n';
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
