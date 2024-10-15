#include "ERF.H"

using namespace amrex;

/** Perform time-integration with multiblock coupled with ERF
 */
// advance solution over specified block steps
void
ERF::Evolve_MB (int MBstep, int max_block_step)
{
    Real cur_time = t_new[0];

    int step;

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    // for finer levels (with or without subcycling)
    for (int Bstep(0); Bstep < max_block_step && cur_time < stop_time; ++Bstep)
    {
        step = Bstep + MBstep - 1;

        Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt(step);

        // Make sure we have read enough of the boundary plane data to make it through this timestep
        if (input_bndry_planes)
        {
            m_r2d->read_input_files(cur_time,dt[0],m_bc_extdir_vals);
        }

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time  += dt[0];

        Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        post_timestep(step, cur_time, dt[0]);

        if (writeNow(cur_time, dt[0], step+1, m_plot_int_1, m_plot_per_1)) {
            last_plot_file_step_1 = step+1;
            WritePlotFile(1,plot_var_names_1);
        }

        if (writeNow(cur_time, dt[0], step+1, m_plot_int_2, m_plot_per_2)) {
            last_plot_file_step_2 = step+1;
            WritePlotFile(2,plot_var_names_2);
        }

        if (writeNow(cur_time, dt[0], step+1, m_check_int, m_check_per)) {
            last_check_file_step = step+1;
#ifdef ERF_USE_NETCDF
            if (check_type == "netcdf") {
               WriteNCCheckpointFile();
            }
#endif
            if (check_type == "native") {
               WriteCheckpointFile();
            }
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }
}
