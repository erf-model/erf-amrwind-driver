#include "prob.H"
#include "AMReX_Random.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("A_0", parms.A_0);
  pp.query("KE_0", parms.KE_0);
  pp.query("KE_decay_height", parms.KE_decay_height);
  pp.query("KE_decay_order", parms.KE_decay_order);

  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
  pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
  pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
  pp.query("T_0_Pert_Mag", parms.T_0_Pert_Mag);
  pp.query("pert_rhotheta", parms.pert_rhotheta);

  pp.query("pert_deltaU", parms.pert_deltaU);
  pp.query("pert_deltaV", parms.pert_deltaV);
  pp.query("pert_periods_U", parms.pert_periods_U);
  pp.query("pert_periods_V", parms.pert_periods_V);
  pp.query("pert_ref_height", parms.pert_ref_height);
  parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
  parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
  parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
  parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

  // ABL-with-bubble
  pp.query("use_bubble", parms.use_bubble);
  if (parms.use_bubble) {
    amrex::Print() << "Initializing with bubble" << std::endl;
    pp.query("bubble_loc_x", parms.bubble_loc_x);
    pp.query("bubble_loc_y", parms.bubble_loc_y);
    pp.query("bubble_loc_z", parms.bubble_loc_z);
    pp.query("bubble_radius", parms.bubble_radius);
    pp.query("bubble_temp_ratio", parms.bubble_temp_ratio);
  } else {
    amrex::Print() << "Not initializing with bubble" << std::endl;
  }

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const amrex::Box&  bx,
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& state_pert,
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
    const SolverChoice& sc)
{

  const bool use_moisture = (sc.moisture_type != MoistureType::None);

  ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    //const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Add temperature perturbations
    if ((z <= parms.pert_ref_height) && (parms.T_0_Pert_Mag != 0.0)) {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*parms.T_0_Pert_Mag;
    }

    // For bubble
    amrex::Real ratio = 1.0;
    const amrex::Real bcx = parms.bubble_loc_x;
    const amrex::Real bcy = parms.bubble_loc_y;
    const amrex::Real bcz = parms.bubble_loc_z;
    if (parms.use_bubble) {
      amrex::Real radius = std::sqrt((x-bcx)*(x-bcx) + (y-bcy)*(y-bcy) + (z-bcz)*(z-bcz));
      ratio = 1.0 + (parms.bubble_temp_ratio - 1.0) * exp(-0.5 * radius* radius / (parms.bubble_radius * parms.bubble_radius));
      state_pert(i, j, k, RhoScalar_comp) = 300.0 * ratio * parms.rho_0;
    }
    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    //  state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r*r);

    // Set an initial value for QKE
    state_pert(i, j, k, RhoKE_comp) = parms.KE_0;

    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = 0.0;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }
  });

  // Set the x-velocity
  ParallelForRNG(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the x-velocity as a constant gradient in z
    //x_vel(i, j, k) = parms.U_0;
    x_vel(i, j, k) = parms.U_0 * z / prob_hi[2];

    if ((z <= parms.pert_ref_height) && (parms.U_0_Pert_Mag != 0.0))
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real x_vel_prime = (rand_double*2.0 - 1.0)*parms.U_0_Pert_Mag;
        x_vel(i, j, k) += x_vel_prime;
    }
    if (parms.pert_deltaU != 0.0)
    {
        const amrex::Real yl = y - prob_lo[1];
        const amrex::Real zl = z / parms.pert_ref_height;
        const amrex::Real damp = std::exp(-0.5 * zl * zl);
        x_vel(i, j, k) += parms.ufac * damp * z * std::cos(parms.aval * yl);
    }
  });

  // Set the y-velocity
  ParallelForRNG(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel(i, j, k) = parms.V_0;
    if ((z <= parms.pert_ref_height) && (parms.V_0_Pert_Mag != 0.0))
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real y_vel_prime = (rand_double*2.0 - 1.0)*parms.V_0_Pert_Mag;
        y_vel(i, j, k) += y_vel_prime;
    }
    if (parms.pert_deltaV != 0.0)
    {
        const amrex::Real xl = x - prob_lo[0];
        const amrex::Real zl = z / parms.pert_ref_height;
        const amrex::Real damp = std::exp(-0.5 * zl * zl);
        y_vel(i, j, k) += parms.vfac * damp * z * std::cos(parms.bval * xl);
    }
  });

  // Set the z-velocity
  ParallelForRNG(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const int dom_lo_z = geomdata.Domain().smallEnd()[2];
    const int dom_hi_z = geomdata.Domain().bigEnd()[2];

    // Set the z-velocity
    if (k == dom_lo_z || k == dom_hi_z+1)
    {
        z_vel(i, j, k) = 0.0;
    }
    else if (parms.W_0_Pert_Mag != 0.0)
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real z_vel_prime = (rand_double*2.0 - 1.0)*parms.W_0_Pert_Mag;
        z_vel(i, j, k) = parms.W_0 + z_vel_prime;
    }
  });
}
