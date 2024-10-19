#include "amr-wind/incflo.H"

using namespace amrex;

/** Perform time-integration with multiblock coupled with ERF
 */
void incflo::Evolve_MultiBlock(int MBstep, int max_block_step)
{
    BL_PROFILE("amr-wind::incflo::Evolve()");

    int step = 0;

    amrex::Print() << "\n======================================================"
      "========================\n";
    amrex::Print() << "AMR_WIND STEP" << std::endl;
    while (step < max_block_step) {
        m_time.new_timestep();
        step+=1;
        amrex::Real time0 = amrex::ParallelDescriptor::second();

        regrid_and_update();

        if (m_prescribe_vel) {
            pre_advance_stage2();
            ComputePrescribeDt();
        } else {
            pre_advance_stage1();
            pre_advance_stage2();
        }

        amrex::Real time1 = amrex::ParallelDescriptor::second();
        // Advance to time t + dt
        for (int fixed_point_iteration = 0;
             fixed_point_iteration < m_fixed_point_iterations;
             ++fixed_point_iteration)
            do_advance(fixed_point_iteration);

        amrex::Print() << std::endl;
        amrex::Real time2 = amrex::ParallelDescriptor::second();
        post_advance_work();
        amrex::Real time3 = amrex::ParallelDescriptor::second();

        amrex::Print() << "WallClockTime: " << m_time.time_index()
                       << " Pre: " << std::setprecision(3) << (time1 - time0)
                       << " Solve: " << std::setprecision(4) << (time2 - time1)
                       << " Post: " << std::setprecision(3) << (time3 - time2)
                       << " Total: " << std::setprecision(4) << (time3 - time0)
                       << std::endl;

        amrex::Print() << "Solve time per cell: " << std::setprecision(4)
                       << amrex::ParallelDescriptor::NProcs() *
                              (time2 - time1) /
                              static_cast<amrex::Real>(m_cell_count)
                       << std::endl;
    }
    amrex::Print() << "\n======================================================"
                      "========================\n"
                   << std::endl;

}
