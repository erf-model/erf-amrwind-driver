#include <MultiBlockContainer.H>
#include <AMReX_NonLocalBC.H>
#include <ERF.H>
#include <ABLReadERF.H>

using namespace amrex;

// Vector input constructor
MultiBlockContainer::MultiBlockContainer (
  const std::vector<amrex::RealBox>& rb_v,
  std::vector<int> max_level_in_v,
  const std::vector<amrex::Vector<int>>& n_cell_in_v,
  std::vector<int> coord_v,
  const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
  const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
  std::vector<std::string> prefix_v,
  int max_step)
: erf1(rb_v[0],max_level_in_v[0],n_cell_in_v[0],coord_v[0],ref_ratios_v[0],is_per_v[0],prefix_v[0]),
  m_max_step(max_step)
{
  // Store ptr to container to call member functions
  amrwind.SetMultiBlockPointer(this);
  amrwind.set_read_erf(read_erf);
  erf1.SetMultiBlockPointer(this);

  // Set the permutation/sign of dtos
  dtos_efroma.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
  dtos_efroma.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};
  dtos_afrome.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
  dtos_afrome.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};

  // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
  {
    const amrex::Geometry geom = amrwind.repo().mesh().Geom(0);
    amrex::Real dx = geom.CellSize(0);
    amrex::Real dy = geom.CellSize(1);
    amrex::Real dz = geom.CellSize(2);
    int offx = static_cast<int>(amrex::Math::floor((geom.ProbLo(0) - rb_v[0].lo(0)) / dx));
    int offy = static_cast<int>(amrex::Math::floor((geom.ProbLo(1) - rb_v[0].lo(1)) / dy));
    int offz = static_cast<int>(amrex::Math::floor((geom.ProbLo(2) - rb_v[0].lo(2)) / dz));
    dtos_efroma.offset = amrex::IntVect{AMREX_D_DECL(offx, offy, offz)};
    dtos_afrome.offset = amrex::IntVect{AMREX_D_DECL(-offx, -offy, -offz)};
  }
}

// Destructor
MultiBlockContainer::~MultiBlockContainer ()
{
}

// Initialize block data
void
MultiBlockContainer::InitializeBlocks ()
{
  amrex::Print() << "======== MultiBlock Intitialization I ========"  << "\n";
  amrex::ParmParse pp("mbc");
  pp.query("erf_to_amrwind_dl_ratio", erf_to_aw_dl_ratio);
  pp.query("do_two_way_coupling", do_two_way_coupling);
  pp.query("two_way_coupling_frequency", two_way_coupling_frequency);
  SetBoxLists();

  amrex::Print() << "======== ERF1 INITIALIZATION ========"  << "\n";
  erf1.InitData();

  amrex::Print() << "======== AMRWIND INITIALIZATION ========"  << "\n";
  amrwind.InitData();

  // Allocate boundary registers for AMR-Wind
  // Only two fields supported for now -- velocity and temperature
  bndry1.resize(num_fields);
  bndry2.resize(num_fields);

  amrex::BoxArray ba(amrwind.boxArray(0));
  amrex::DistributionMapping dm{ba};
  const int in_rad     = 1;
  const int out_rad    = 1;
  const int extent_rad = 0;

  // for (i = 0; i < nfields; i++)
  // TODO: COULD THIS BE LOOP OVER NFIELDS?
  {
    // Velocity
    bndry1[0] = new BndryRegister(ba, dm, in_rad, out_rad, extent_rad, 3);
    bndry2[0] = new BndryRegister(ba, dm, in_rad, out_rad, extent_rad, 3);

    // Temperature
    bndry1[1] = new BndryRegister(ba, dm, in_rad, out_rad, extent_rad, 1);
    bndry2[1] = new BndryRegister(ba, dm, in_rad, out_rad, extent_rad, 1);
  }

  amrex::Print() << "======== MultiBlock Intitialization II ========"  << "\n";
  SetBlockCommMetaData();

  if (do_two_way_coupling) {
    amrex::Print() << "======== FILLPATCH A->E ========"  << "\n";
    FillPatchBlocksAE ();
  }
}

// Set up BoxList vector for use with Communication Meta Data
void
MultiBlockContainer::SetBoxLists ()
{
  // when copying from amr-wind to erf, data from the entire amr-wind domain is used
  amrex::Box awbox = amrwind.repo().mesh().Geom(0).Domain();
  abox_efroma = amrex::Box(awbox.smallEnd(), awbox.bigEnd());
  amrex::Dim3 se = dtos_efroma.Inverse(amrex::lbound(abox_efroma));
  amrex::Dim3 be = dtos_efroma.Inverse(amrex::ubound(abox_efroma));
  ebox_efroma = amrex::Box({se.x, se.y, se.z}, {be.x, be.y, be.z});

  // when copying from erf to amr-wind, we only copy data at the boundaries of the amr-wind domain
  bool ok_to_continue = true;
  // refine ERF domain for error checking later
  const amrex::Box erf_refined_domain{amrex::refine(erf1.domain_p[0], erf_to_aw_dl_ratio)};
  //amrex::Print() << "ERF refined domain: " << erf_refined_domain << std::endl;
  for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
    auto ori = oit();
    amrex::IntVect sev = awbox.smallEnd();
    amrex::IntVect bev = awbox.bigEnd();
    if (ori.faceDir() == 1) {
      sev[ori.coordDir()] = bev[ori.coordDir()];
      bev[ori.coordDir()] += 1;
    }
    else {
      bev[ori.coordDir()] = sev[ori.coordDir()];
      sev[ori.coordDir()] -= 1;
    }
    amrex::Box abx(sev,bev);
    aboxvec_afrome.push_back(abx);
    amrex::Dim3 ese = dtos_afrome(amrex::lbound(abx));
    amrex::Dim3 ebe = dtos_afrome(amrex::ubound(abx));
    amrex::Box ebx({ese.x, ese.y, ese.z}, {ebe.x, ebe.y, ebe.z});
    eboxvec_afrome.push_back(ebx);

    //amrex::Print() << "Ori: " << ori << " A-W BNDRY in A-W Coords: " << abx << std::endl;
    //amrex::Print() << "Ori: " << ori << " A-W BNDRY in ERF Coords: " << ebx << std::endl;

    // Do some error checking to ensure that the AMR-Wind is fully within the ERF domain (not touching any boundaries)
    // Note: in theory this error check could trip for a non-boundary plane mass_inflow in AMR-Wind
    auto bctype = amrwind.repo().get_field("velocity").bc_type()[ori];
    bool need_bndry = (bctype == BC::mass_inflow) || (bctype == BC::mass_inflow_outflow);
    if ( !(erf_refined_domain.contains(ebx)) and need_bndry ) {
      amrex::Print() << "ERF domain must fully contain the AMR-Wind boundary planes, does not in direction "
                     << ori.coordDir() << " on face " << ori.faceDir() << std::endl;
      ok_to_continue = false;
    }
  }
  if (!ok_to_continue) amrex::Abort("Invalid domains for ERF/AMR-Wind coupling");
}

// Set up MB Communication Meta Data
void
MultiBlockContainer::SetBlockCommMetaData ()
{
  amrex::IntVect nghost(0);
  amrex::NonLocalBC::MultiBlockCommMetaData *cmd_efroma_tmp =
    new amrex::NonLocalBC::MultiBlockCommMetaData(
      erf1.vars_new[0][Vars::cons], ebox_efroma,
      amrwind.repo().get_field("temperature")(0), nghost, dtos_efroma);
  cmd_efroma.push_back(cmd_efroma_tmp);
}

// Advance blocks
void
MultiBlockContainer::AdvanceBlocks()
{
  amrex::Print() << "======== STARTING MAIN DRIVER FOR: "
                 << m_max_step << " STEPS" << "\n";

  // Calculate ratio of AMR-Wind to ERF deltaT
  // and verify that it is an integer
  int aw_to_erf_dt_ratio;
  {
    amrex::Real erf_dt = erf1.get_dt(0);
    amrex::Real amrwind_dt = amrwind.time().delta_t();
    amrex::Real aw_to_erf_dt = amrwind_dt / erf_dt;
    aw_to_erf_dt_ratio = static_cast<int>(std::round(aw_to_erf_dt));
    const amrex::Real eps = 1e-8;
    AMREX_ALWAYS_ASSERT(std::abs(aw_to_erf_dt - aw_to_erf_dt_ratio) <= eps);
    amrex::Print() << "ERF will make " << aw_to_erf_dt_ratio
                   << " steps for each AMR-Wind step." << std::endl;
  }

  // NOTE: step here corresponds to the number of AMR-Wind steps taken
  for (int step(0); step < m_max_step; ++step) {

    if (step == 0) {
      fill_amrwind_bndry(bndry1, this);//, true);
    }

    amrex::Print() << "======== ERF EVOLVE ========"  << "\n";
    erf1.Evolve_MB(aw_to_erf_dt_ratio*step+1, aw_to_erf_dt_ratio);

    fill_amrwind_bndry(bndry2, this);

    amrex::Print() << "======== AMRWIND EVOLVE ========"  << "\n";
    amrwind.Evolve_MultiBlock(step+1,1);

    if (do_two_way_coupling && (((step+1) % two_way_coupling_frequency) == 0)) {
      amrex::Print() << "======== FILLPATCH A->E ========"  << "\n";
      FillPatchBlocksAE ();
    }
    amrex::Print() << "======== DRIVER STEP COMPLETE ========"  << "\n";

    // swap old and new boundary fields
    // TODO: COULD BE LOOP OVER NFIELDS?
    std::swap(bndry2[0],bndry1[0]); // Velocity
    std::swap(bndry2[1],bndry1[1]); // Temperature
  }
}

// Fill AMR-Wind Boundary Regsiter from ERF1
void
MultiBlockContainer::CopyERFtoAMRWindBoundaryReg (
  amrex::BndryRegister& receive_br,
  amrex::Orientation ori,
  amrex::Real time,
  const std::string &field) {

  // Need a ghost cell in case AMR-Wind boundary to be filled coincides with ERF boundary
  // FOR NOW - don't support this, ERF must be interior of AMR-Wind
  // nghost[ori.coordDir()] = 1;
  amrex::IntVect nghost(0);

  // ERF level where we are getting the data
  const int erf_source_level = 0;

  bool on_old_time{time == erf1.get_t_old(0)};
  bool on_new_time{time == erf1.get_t_new(0)};
  AMREX_ALWAYS_ASSERT(on_new_time || on_old_time);
  /*
  amrex::Print() << "IN COPY ERF TO AWBR " << std::endl;
  amrex::Print() << "    TIME IS " <<  time << std::endl;
  amrex::Print() << "OLD TIME IS " <<  erf1.get_t_old(0) << std::endl;
  amrex::Print() << "NEW TIME IS " <<  erf1.get_t_new(0) << std::endl;
  */
  if (on_old_time) {
    amrex::Print() << std::endl << "FILLPATCHING _ " << field
                   << " _ FROM ERF TO AMR WIND ON _ old _ TIME, ORIENTATION " << ori << std::endl;
  }
  else {
    amrex::Print() << std::endl << "FILLPATCHING _ " << field
                   << " _ FROM ERF TO AMR WIND ON _ new _ TIME, ORIENTATION " << ori << std::endl;
  }
  amrex::Vector<amrex::MultiFab>& erf_data = on_old_time ? erf1.vars_old[erf_source_level] : erf1.vars_new[erf_source_level];

  // For selecting only subdomain of erf domain for data processing
  std::map<int,int> mfmap;
  amrex::Vector<int> new_dl;
  amrex::BoxList new_bl;
  const amrex::BoxArray& old_ba = erf_data[Vars::cons].boxArray();
  const amrex::DistributionMapping& old_dm = erf_data[Vars::cons].DistributionMap();
  for (int i = 0; i < old_ba.size(); i++) {
    // refine ERF boxes as newmf needs to be at the resolution of AMR-Wind
    amrex::Box isect = refine(old_ba[i], erf_to_aw_dl_ratio) & eboxvec_afrome[ori];
    if (isect.ok()) {
      new_bl.push_back(isect);
      new_dl.push_back(old_dm[i]);
      mfmap.insert({new_dl.size() - 1, i});
    }
  }
  amrex::BoxArray new_ba(new_bl);
  amrex::DistributionMapping new_dm(new_dl);

  if (field == "temperature") {
    amrex::MultiFab newmf(new_ba, new_dm, 1, 0); //nghost);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      for (amrex::MFIter mfi(newmf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Array4<const amrex::Real> erf_arr = erf_data[Vars::cons].const_array(mfmap[mfi.index()]);
        const amrex::Array4<      amrex::Real> erf_arr_copy = newmf.array(mfi);
        int r{erf_to_aw_dl_ratio};
        amrex::ParallelFor(mfi.growntilebox(),[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          erf_arr_copy(i, j, k) = erf_arr(i/r, j/r, k/r, RhoScalar_comp)
                                / erf_arr(i/r, j/r, k/r, Rho_comp);
        });
      }
    }

    // Copy data
    amrex::NonLocalBC::MultiBlockCommMetaData cmd =
      amrex::NonLocalBC::MultiBlockCommMetaData(
        receive_br[ori].multiFab(), aboxvec_afrome[ori], newmf, nghost, dtos_afrome);

    amrex::NonLocalBC::ParallelCopy(
      receive_br[ori].multiFab(), newmf, cmd, 0, 0, 1, dtos_afrome );

  } else if (field == "velocity") {
    amrex::MultiFab newmf(new_ba, new_dm, 3, 0);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      for (amrex::MFIter mfi(newmf,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Array4<const amrex::Real> erf_vel_arr[3] = {
          erf_data[Vars::xvel].const_array(mfmap[mfi.index()]),
          erf_data[Vars::yvel].const_array(mfmap[mfi.index()]),
          erf_data[Vars::zvel].const_array(mfmap[mfi.index()])  };

        const amrex::Array4<amrex::Real> erf_arr_copy = newmf.array(mfi);
        const int ndir  = ori.coordDir(); // Normal Dir
        //const int nface_idx = ori.isLow() ? eboxvec_afrome[ori].smallEnd(ndir) : eboxvec_afrome[ori].bigEnd(ndir);
        const int tdir1 = (ndir + 1) % 3; // First tangenential dir
        const int tdir2 = (ndir + 2) % 3; // second tangential dir

        int r{erf_to_aw_dl_ratio};
        amrex::ParallelFor(mfi.growntilebox(),[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          amrex::IntVect idx {i, j, k};

          // idx_c are the coarse indices mapped from the refined idx
          amrex::IntVect idx_c {idx};
          idx_c[tdir1] = idx[tdir1] / r;
          idx_c[tdir2] = idx[tdir2] / r;
          idx_c[ndir]  = idx[ndir]  / r;
          // interpolate tangential velocity faces to center
          amrex::IntVect idx_tp1 {idx_c};
          idx_tp1[tdir1] += 1;
          erf_arr_copy(idx, tdir1) = 0.5 * (erf_vel_arr[tdir1](idx_c) + erf_vel_arr[tdir1](idx_tp1));

          idx_tp1 = idx_c;
          idx_tp1[tdir2] += 1;
          erf_arr_copy(idx, tdir2) = 0.5 * (erf_vel_arr[tdir2](idx_c) + erf_vel_arr[tdir2](idx_tp1));

          // take normal velocity from face
          idx_c[ndir] = (idx[ndir] + 1) / r;
          // idx_n[ndir] = nface_idx; // earlier solution
          erf_arr_copy(idx, ndir) = erf_vel_arr[ndir](idx_c);
        });
      }
    }

    // Copy data
    amrex::NonLocalBC::MultiBlockCommMetaData cmd =
      amrex::NonLocalBC::MultiBlockCommMetaData(
        receive_br[ori].multiFab(), aboxvec_afrome[ori], newmf, nghost, dtos_afrome);

    amrex::NonLocalBC::ParallelCopy(
      receive_br[ori].multiFab(), newmf, cmd, 0, 0, 3, dtos_afrome );
  } else {
    amrex::Abort("ERF to AMR-Wind copying only supported for fields: temperature, velocity");
  }
}

void
MultiBlockContainer::PopulateErfTimesteps (amrex::Real* tsteps) {
  tsteps[0] = erf1.get_t_old(0);
  tsteps[1] = erf1.get_t_new(0);
}

// Wrapper for ParallelCopy between ERF and AMRWIND
void
MultiBlockContainer::FillPatchBlocksAE()
{
  // Create temporary multifabs to store data from AMR-Wind
  // Note: AMR-Wind data is cell centered for all variables, and we only care about the valid data
  //       so no ghost cells for Temp or Dens. But we will be interpolating velocity to faces so we
  //       need one ghost cell for velocity
  // TODO: make these only correspond to the needed region (box_efroma) rather than the full ERF domain
  const amrex::BoxArray& ba            = erf1.vars_new[0][Vars::cons].boxArray();
  const amrex::DistributionMapping& dm = erf1.vars_new[0][Vars::cons].DistributionMap();
  amrex::MultiFab Temp_AW{ba, dm, 1, 0};
  amrex::MultiFab Dens_AW{ba, dm, 1, 0};
  amrex::MultiFab Vel_AW {ba, dm, 3, 1};

  // Bring AMR-Wind data to temporary multifabs
  amrex::NonLocalBC::ParallelCopy(
    Temp_AW, amrwind.repo().get_field("temperature")(0), *(cmd_efroma[0]), 0, 0, 1, dtos_efroma);
  amrex::NonLocalBC::ParallelCopy(
    Dens_AW, amrwind.repo().get_field("density")(0),     *(cmd_efroma[0]), 0, 0, 1, dtos_efroma);
  amrex::NonLocalBC::ParallelCopy(
    Vel_AW,  amrwind.repo().get_field("velocity")(0),    *(cmd_efroma[0]), 0, 0, 3, dtos_efroma);

  // Compute ERF variables from AMR-Wind variables and store in ERF data structures

  // Cell centered quantities
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    // compute cell-centered scalar and velocities
    auto cons_arr    = erf1.vars_new[0][Vars::cons][mfi].array();
    auto vel_aw_arr  = Vel_AW[mfi].array();
    auto dens_aw_arr = Dens_AW[mfi].array();
    auto temp_aw_arr = Temp_AW[mfi].array();
    amrex::Box ibox = box & ebox_efroma; // intersection of boxes
    amrex::ParallelFor(ibox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
      // Save AMR-Wind Scalar into ERF data
      cons_arr(i, j, k, RhoScalar_comp) = cons_arr(i, j, k, Rho_comp) * temp_aw_arr(i, j, k);
      // Cell-centered velocity using energy preserving correction from Sprague & Satkauskas 2015
      amrex::Real dens_correction = std::sqrt(dens_aw_arr(i, j, k) / cons_arr(i, j, k, Rho_comp));
      vel_aw_arr(i, j, k, 0) *= dens_correction;
      vel_aw_arr(i, j, k, 1) *= dens_correction;
      vel_aw_arr(i, j, k, 2) *= dens_correction;
    });
  }

  // We need to fill the interior boundary cells of velocity so we can interpolate to faces
  // Notes: Because the AMR-Wind domain is assumed to not exceed the extent of the ERF domain
  //         and we will only be filling velocity at face centers within (rather than on the
  //         boundary of) the AMR-Wind domain, we can get away with just doing a FillBoundary
  //        This fills all boundary cells but in theory we only need to fill boundary cells
  //         surrounding ebox_efroma
  Vel_AW.FillBoundary(erf1.Geom(0).periodicity());

  // Move cell centered velocity to faces and store for ERF
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    auto vel_aw_arr = Vel_AW[mfi].array();
    // interpolate corrected velocities to face centers (only on interior)
    auto velx_arr = erf1.vars_new[0][Vars::xvel][mfi].array();
    amrex::Box xbox = mfi.nodaltilebox(0) & surroundingNodes(ebox_efroma, 0).grow(0, -1);
    amrex::ParallelFor(xbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
      velx_arr(i, j, k) = 0.5 * (vel_aw_arr(i, j, k, 0) + vel_aw_arr(i-1, j, k, 0));
    });

    auto vely_arr = erf1.vars_new[0][Vars::yvel][mfi].array();
    amrex::Box ybox =  mfi.nodaltilebox(1) & surroundingNodes(ebox_efroma, 1).grow(1, -1);
    amrex::ParallelFor(ybox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
      vely_arr(i, j, k) = 0.5 * (vel_aw_arr(i, j, k, 1) + vel_aw_arr(i, j-1, k, 1));
    });

    auto velz_arr = erf1.vars_new[0][Vars::zvel][mfi].array();
    amrex::Box zbox =  mfi.nodaltilebox(2) & surroundingNodes(ebox_efroma, 2).grow(2, -1);
    amrex::ParallelFor(zbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
      velz_arr(i, j, k) = 0.5 * (vel_aw_arr(i, j, k, 2) + vel_aw_arr(i, j, k-1, 2));
    });
  }
}
