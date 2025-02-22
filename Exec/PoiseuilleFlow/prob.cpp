#include "prob.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);

  pp.query("prob_type", parms.prob_type);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const amrex::Box& /*bx*/,
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& /*state_pert*/,
    amrex::Array4<amrex::Real      > const& x_vel,
    amrex::Array4<amrex::Real      > const& y_vel,
    amrex::Array4<amrex::Real      > const& z_vel,
    amrex::Array4<amrex::Real      > const& /*r_hse*/,
    amrex::Array4<amrex::Real      > const& /*p_hse*/,
    amrex::Array4<amrex::Real const> const& /*z_nd*/,
    amrex::Array4<amrex::Real const> const& /*z_cc*/,
    amrex::GeometryData const& geomdata,
    amrex::Array4<amrex::Real const> const& /*mf_m*/,
    amrex::Array4<amrex::Real const> const& /*mf_u*/,
    amrex::Array4<amrex::Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/)
{
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx      = geomdata.CellSize();
      const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

      // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
      if (parms.prob_type == 10)
          //x_vel(i, j, k) = 1.0 - z_h * z_h;
          x_vel(i, j, k) = 0.95 * z_h * (2 - z_h);
      else if (parms.prob_type == 11)
	        x_vel(i, j, k) = z_h;
      else
          x_vel(i, j, k) = 0.0;
  });

  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    y_vel(i, j, k) = 0.0;
  });

  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    z_vel(i, j, k) = 0.0;
  });
}
