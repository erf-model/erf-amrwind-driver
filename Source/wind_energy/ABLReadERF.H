#ifndef ABLREADERF_H
#define ABLREADERF_H

#include <amr-wind/wind_energy/ABLBoundaryPlane.H>
#include "MultiBlockContainer.H"

void fill_old_bndry(amrex::Vector<amrex::BndryRegister*>& bndry1,
                    MultiBlockContainer* mbc);
void fill_new_bndry(amrex::Vector<amrex::BndryRegister*>& bndry2,
                    MultiBlockContainer* mbc);

void read_erf(const amrex::Real time,
              amrex::Vector<amrex::Real>& m_in_times,
              amr_wind::InletData& m_in_data,
              const amrex::Vector<amr_wind::Field*>& m_fields,
              MultiBlockContainer* mbc);

#endif /* ABLREADERF_H */
