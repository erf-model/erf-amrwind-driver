#ifndef MultiBlockContainer_H_
#define MultiBlockContainer_H_

#include <AMReX_NonLocalBC.H>
#include <ERF.H>
#include <amr-wind/incflo.H>

class MultiBlockContainer
{
public:

    // Constructor with vector of parameters
    MultiBlockContainer (const std::vector<amrex::RealBox>& rb_v,
                         std::vector<int> max_level_in_v,
                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                         std::vector<int> coord_v,
                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
                         std::vector<std::string> prefix_v,
                         int max_step);

    // Destructor
    ~MultiBlockContainer ();

    // Initialize block data
    void InitializeBlocks ();

    // Set up BoxList vector for use with Communication Meta Data
    void SetBoxLists ();

    // Set up MB Communication Meta Data
    void SetBlockCommMetaData ();

    // Advance blocks
    void AdvanceBlocks ();

    // Fill MFs from other block(s)
    void FillPatchBlocksAE ();

    // Fill AMR-Wind Boundary Regsiter from ERF1
  void CopyERFtoAMRWindBoundaryReg (amrex::BndryRegister& receive_br,  amrex::Orientation ori, amrex::Real time, const std::string& field);

  //void PopulateErfTimesteps (amrex::Vector<amrex::Real> tsteps) const {erf1.CoarseTimeState(/*tsteps*/);}
  void PopulateErfTimesteps(amrex::Real* tsteps);

    // Constructor with vector of parameters
    MultiBlockContainer (const std::vector<amrex::RealBox>& rb_v,
                         std::vector<int> max_level_in_v,
                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                         std::vector<int> coord_v,
                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
                         std::vector<std::string> prefix_v,
                         int max_step, int is_erf);
    // Destructor
//    ~MultiBlockContainer ();

    // Initialize block data
    void ERF_InitializeBlocks ();

    // Set up BoxList vector for use with Communication Meta Data
    void ERF_SetBoxLists ();

    // Set up MB Communication Meta Data
    void ERF_SetBlockCommMetaData ();

    // Advance blocks
    void ERF_AdvanceBlocks ();

    // Fill MFs from other block(s)
    void ERF_FillPatchBlocks (int src_ind, int dst_ind);

private:
  std::vector<amrex::Box> eboxvec_afrome;
  std::vector<amrex::Box> aboxvec_afrome;
    int m_max_step;
  // TODO: why does amrwind have to come before erf?
    incflo amrwind;
    ERF erf1;

    bool do_two_way_coupling;
  int two_way_coupling_frequency;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos_efroma;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos_afrome;
  std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*> cmd_efroma;
    amrex::Box ebox_efroma;
    amrex::Box abox_efroma;
  //std::vector<std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>> cmd_afrome;

    ERF erf2;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos;
    std::vector<std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>> cmd;
    std::vector<amrex::BoxList> blv;
};

#endif
