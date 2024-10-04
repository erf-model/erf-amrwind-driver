#include "ABLReadERF.H"

void read_erf(const amr_wind::SimTime& m_time,
              amrex::Vector<amrex::Real>& m_in_times,
              amr_wind::InletData& m_in_data,
              const amrex::Vector<amr_wind::Field*>& m_fields,
              MultiBlockContainer* mbc)
{
    amrex::Real time = m_time.new_time();
    AMREX_ALWAYS_ASSERT(m_in_times[0] <= time); // Can't go back in time for ERF data
    // return early if current erf data can still be interpolated in time
    if ((m_in_data.tn() <= time) && (time < m_in_data.tnp1())) {
        m_in_data.interpolate(time);
        return;
    }

    // Get current ERF time values
    mbc->PopulateErfTimesteps(m_in_times.data());
    AMREX_ALWAYS_ASSERT((m_in_times[0]<= time) && (time <= m_in_times[1]));
    const int index = 0;
    const int lev = 0;
    // FIXME FIXME TODO DELETE
    mbc->SetBoxLists();

    for (auto* fld : m_fields) {

      auto& field = *fld;
      const auto& geom = field.repo().mesh().Geom();

      amrex::Box domain = geom[lev].Domain();
      amrex::BoxArray ba(domain);
      amrex::DistributionMapping dm{ba};

      const int in_rad = 1;
      const int out_rad = 1;
      const int extent_rad = 0;
      amrex::BndryRegister bndry1(ba, dm, in_rad, out_rad, extent_rad, field.num_comp());
      amrex::BndryRegister bndry2(ba, dm, in_rad, out_rad, extent_rad, field.num_comp());

      if (field.name() == "velocity") {
        bndry1.setVal(1.0e13); // 1.0e13
        bndry2.setVal(1.0e13); // 1.0e13
      } else if (field.name() == "temperature") {
        bndry1.setVal(1.0e13); // 1.0e13
        bndry2.setVal(1.0e13); // 1.0e13
      }

      for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((!m_in_data.is_populated(ori)) ||
            (field.bc_type()[ori] != BC::mass_inflow)) {
          //continue;
        }
        if (field.bc_type()[ori] == BC::mass_inflow and time >= 0.0) {
          if (field.name() == "temperature") {
            amrex::Print() << "COPY ERF TO AMRWIND ... " << std::endl;
            mbc->CopyERFtoAMRWindBoundaryReg(bndry1, ori, m_in_times[0], field.name());
            mbc->CopyERFtoAMRWindBoundaryReg(bndry2, ori, m_in_times[1], field.name());
            /*
            if ( field.name() == "temperature") {
                mbc->CopyToBoundaryRegister(bndry1, ori);
            } else if ( field.name() == "velocity") {
                // mbc->CopyToBoundaryRegister(bndry1, bndry2, ori);
                } */
          } else if (field.name() == "velocity") {
            mbc->CopyERFtoAMRWindBoundaryReg(bndry1, ori, m_in_times[0], field.name());
            mbc->CopyERFtoAMRWindBoundaryReg(bndry2, ori, m_in_times[1], field.name());
          }
        } else {
          if (field.name() == "temperature") {
            bndry1[ori].setVal(300.0);
            bndry2[ori].setVal(300.0);
          } else if (field.name() == "velocity") {
            bndry1[ori].setVal(10.0, 0, 1);
            bndry2[ori].setVal(10.0, 0, 1);
            bndry1[ori].setVal( 0.0, 1, 1);
            bndry2[ori].setVal( 0.0, 1, 1);
            bndry1[ori].setVal( 0.0, 2, 1);
            bndry2[ori].setVal( 0.0, 2, 1);
          }
        }
        /*
        amrex::IntVect nghost(0);
        amrex::NonLocalBC::MultiBlockCommMetaData *cmd_full_tmp =
          new amrex::NonLocalBC::MultiBlockCommMetaData(bndry1[ori].multiFab(), domain,
                                                        bndry2[ori].multiFab(), nghost, mbc->dtos_etoa);
        */

        //std::cout << ori << " after break " << std::endl;

        //std::cout << "BNDRY REG " << field.name() << " " << ori << " " << std::endl;
        /*
          std::string facename1 =
          amrex::Concatenate(filename1 + '_', ori, 1);
          std::string facename2 =
          amrex::Concatenate(filename2 + '_', ori, 1);
          bndry1[ori].read(facename1);
          bndry2[ori].read(facename2);
        */
       //amrex::Print() << "*****" << m_in_times[0] << " " << m_in_times[1] << " " << time << std::endl;
        m_in_data.read_data_native(oit, bndry1, bndry2, lev, fld, time, m_in_times, true);

      }
    }
    m_in_data.interpolate(time);
}
