#ifndef ABLREADERF_H
#define ABLREADERF_H

#include <amr-wind/wind_energy/ABLBoundaryPlane.H>
#include "MultiBlockContainer.H"

void read_erf(const amr_wind::SimTime& m_time,
              amrex::Vector<amrex::Real>& m_in_times,
              amr_wind::InletData& m_in_data,
              const amrex::Vector<amr_wind::Field*>& m_fields,
              MultiBlockContainer* mbc);

#endif /* ABLREADERF_H */
