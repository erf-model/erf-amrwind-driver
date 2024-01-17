function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_erf_lib erf_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_lib_name})
  set(ERF_SRC_DIR  ${ERF_HOME}/Source)
  set(AMRW_SRC_DIR  ${AMRWIND_HOME}/amr-wind)

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_lib_name})

  target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MOISTURE)

  if(ERF_ENABLE_MULTIBLOCK)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/MultiBlock/MultiBlockContainer.cpp
                   ${SRC_DIR}/wind_energy/ABLBoundaryPlane.cpp
                   ${SRC_DIR}/wind_energy/ABLFieldInit.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MULTIBLOCK)
    target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/MultiBlock)
  endif()

  # Generate AMR-Wind version header
  configure_file("${AMRWIND_HOME}/cmake/AMRWindVersion.H.in"
  "${CMAKE_BINARY_DIR}/amr-wind/AMRWindVersion.H" @ONLY)

  target_sources(${erf_lib_name} PRIVATE
                ${AMRW_SRC_DIR}/incflo_advance.cpp
                ${AMRW_SRC_DIR}/incflo.cpp
                ${AMRW_SRC_DIR}/incflo_compute_dt.cpp
                ${AMRW_SRC_DIR}/incflo_regrid.cpp
                ${AMRW_SRC_DIR}/helics.cpp
                ${AMRW_SRC_DIR}/CFDSim.cpp
                ${AMRW_SRC_DIR}/core/SimTime.cpp
                ${AMRW_SRC_DIR}/core/Field.cpp
                ${AMRW_SRC_DIR}/core/IntField.cpp
                ${AMRW_SRC_DIR}/core/FieldRepo.cpp
                ${AMRW_SRC_DIR}/core/ScratchField.cpp
                ${AMRW_SRC_DIR}/core/IntScratchField.cpp
                ${AMRW_SRC_DIR}/core/ViewField.cpp
                ${AMRW_SRC_DIR}/core/MLMGOptions.cpp
                ${AMRW_SRC_DIR}/core/MeshMap.cpp
                ${AMRW_SRC_DIR}/boundary_conditions/BCInterface.cpp
                ${AMRW_SRC_DIR}/boundary_conditions/FixedGradientBC.cpp
                ${AMRW_SRC_DIR}/boundary_conditions/wall_models/WallFunction.cpp
                ${AMRW_SRC_DIR}/convection/incflo_godunov_weno.cpp
                ${AMRW_SRC_DIR}/convection/incflo_godunov_ppm.cpp
                #${AMRW_SRC_DIR}/convection/incflo_godunov_plm.cpp
                ${AMRW_SRC_DIR}/convection/incflo_godunov_predict.cpp
                ${AMRW_SRC_DIR}/convection/incflo_godunov_advection.cpp
                ${AMRW_SRC_DIR}/convection/incflo_mol_fluxes.cpp
                ${AMRW_SRC_DIR}/derive/field_algebra.cpp
                ${AMRW_SRC_DIR}/diffusion/incflo_diffusion.cpp
                ${AMRW_SRC_DIR}/projection/incflo_apply_nodal_projection.cpp
                ${AMRW_SRC_DIR}/setup/init.cpp
                ${AMRW_SRC_DIR}/utilities/diagnostics.cpp
                ${AMRW_SRC_DIR}/utilities/io.cpp
                ${AMRW_SRC_DIR}/utilities/bc_ops.cpp
                ${AMRW_SRC_DIR}/utilities/console_io.cpp
                ${AMRW_SRC_DIR}/utilities/IOManager.cpp
                ${AMRW_SRC_DIR}/utilities/FieldPlaneAveraging.cpp
                ${AMRW_SRC_DIR}/utilities/FieldPlaneAveragingFine.cpp
                ${AMRW_SRC_DIR}/utilities/SecondMomentAveraging.cpp
                ${AMRW_SRC_DIR}/utilities/ThirdMomentAveraging.cpp
                ${AMRW_SRC_DIR}/utilities/PostProcessing.cpp
                ${AMRW_SRC_DIR}/utilities/DerivedQuantity.cpp
                ${AMRW_SRC_DIR}/utilities/DerivedQtyDefs.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/RefinementCriteria.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/CartBoxRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/FieldRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/GradientMagRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/CurvatureRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/QCriterionRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/VorticityMagRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/OversetRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/GeometryRefinement.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/BoxRefiner.cpp
                ${AMRW_SRC_DIR}/utilities/tagging/CylinderRefiner.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/Sampling.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/SamplingContainer.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/LineSampler.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/LidarSampler.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/DTUSpinnerSampler.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/PlaneSampler.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/ProbeSampler.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/FieldNorms.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/KineticEnergy.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/Enstrophy.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/FreeSurface.cpp
                ${AMRW_SRC_DIR}/utilities/sampling/WaveEnergy.cpp
                ${AMRW_SRC_DIR}/utilities/averaging/TimeAveraging.cpp
                ${AMRW_SRC_DIR}/utilities/averaging/ReAveraging.cpp
                ${AMRW_SRC_DIR}/utilities/averaging/ReynoldsStress.cpp
                #${AMRW_SRC_DIR}/utilities/ncutils/nc_interface.cpp
                #${AMRW_SRC_DIR}/utilities/ascent/ascent.H
                #${AMRW_SRC_DIR}/utilities/ascent/ascent.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABL.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLStats.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLFieldInitFile.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLMesoscaleForcing.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLWallFunction.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLFillInflow.cpp
                #${AMRW_SRC_DIR}/wind_energy/ABLBoundaryPlane.cpp
                ${AMRW_SRC_DIR}/wind_energy/MOData.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLFillMPL.cpp
                ${AMRW_SRC_DIR}/wind_energy/ABLModulatedPowerLaw.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/actuator_utils.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/Actuator.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/ActuatorContainer.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/FLLC.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/aero/AirfoilTable.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/wing/wing_ops.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/wing/FlatPlate.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/wing/FixedWing.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/turbine/turbine_utils.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/turbine/fast/FastIface.cpp
                #${AMRW_SRC_DIR}/wind_energy/actuator/turbine/fast/TurbineFast.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/disk/ActuatorDisk.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/disk/disk_ops.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/disk/Joukowsky_ops.cpp
                ${AMRW_SRC_DIR}/wind_energy/actuator/disk/uniform_ct_ops.cpp
                ${AMRW_SRC_DIR}/equation_systems/PDEBase.cpp
                ${AMRW_SRC_DIR}/equation_systems/DiffusionOps.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/icns_advection.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/icns.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ABLForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ABLMesoForcingMom.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/GeostrophicForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/HurricaneForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/BoussinesqBuoyancy.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/DensityBuoyancy.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/GravityForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/CoriolisForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/BodyForce.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ABLMeanBoussinesq.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ActuatorForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/icns/source_terms/SynthTurbForcing.cpp
                ${AMRW_SRC_DIR}/equation_systems/temperature/temperature.cpp
                ${AMRW_SRC_DIR}/equation_systems/temperature/source_terms/ABLMesoForcingTemp.cpp
                ${AMRW_SRC_DIR}/equation_systems/density/density.cpp
                ${AMRW_SRC_DIR}/equation_systems/tke/TKE.cpp
                ${AMRW_SRC_DIR}/equation_systems/tke/source_terms/KsgsM84Src.cpp
                ${AMRW_SRC_DIR}/equation_systems/tke/source_terms/KwSSTSrc.cpp
                ${AMRW_SRC_DIR}/equation_systems/sdr/SDR.cpp
                ${AMRW_SRC_DIR}/equation_systems/sdr/source_terms/SDRSrc.cpp
                ${AMRW_SRC_DIR}/equation_systems/levelset/levelset.cpp
                ${AMRW_SRC_DIR}/equation_systems/vof/SplitAdvection.cpp
                ${AMRW_SRC_DIR}/equation_systems/vof/vof.cpp
                ${AMRW_SRC_DIR}/turbulence/LaminarModel.cpp
                ${AMRW_SRC_DIR}/turbulence/turb_utils.cpp
                ${AMRW_SRC_DIR}/turbulence/LES/Smagorinsky.cpp
                ${AMRW_SRC_DIR}/turbulence/LES/OneEqKsgs.cpp
                ${AMRW_SRC_DIR}/turbulence/RANS/KOmegaSST.cpp
                ${AMRW_SRC_DIR}/turbulence/RANS/KOmegaSSTIDDES.cpp
                ${AMRW_SRC_DIR}/physics/BoussinesqBubble.cpp
                ${AMRW_SRC_DIR}/physics/BoussinesqBubbleFieldInit.cpp
                ${AMRW_SRC_DIR}/physics/ChannelFlow.cpp
                ${AMRW_SRC_DIR}/physics/RayleighTaylor.cpp
                ${AMRW_SRC_DIR}/physics/RayleighTaylorFieldInit.cpp
                ${AMRW_SRC_DIR}/physics/TaylorGreenVortex.cpp
                ${AMRW_SRC_DIR}/physics/FreeStream.cpp
                ${AMRW_SRC_DIR}/physics/ConvectingTaylorVortex.cpp
                ${AMRW_SRC_DIR}/physics/EkmanSpiral.cpp
                ${AMRW_SRC_DIR}/physics/SyntheticTurbulence.cpp
                ${AMRW_SRC_DIR}/physics/HybridRANSLESABL.cpp
                ${AMRW_SRC_DIR}/physics/VortexRing.cpp
                ${AMRW_SRC_DIR}/physics/BurggrafFlow.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/MultiPhase.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/VortexPatch.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/VortexPatchScalarVel.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/ZalesakDisk.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/ZalesakDiskScalarVel.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/DamBreak.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/SloshingTank.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/RainDrop.cpp
                ${AMRW_SRC_DIR}/physics/multiphase/BreakingWaves.cpp
                ${AMRW_SRC_DIR}/physics/udfs/UDF.cpp
                ${AMRW_SRC_DIR}/physics/udfs/LinearProfile.cpp
                ${AMRW_SRC_DIR}/physics/udfs/PowerLawProfile.cpp
                ${AMRW_SRC_DIR}/physics/udfs/BurggrafLid.cpp
                #${AMRW_SRC_DIR}/physics/mms/MMS.cpp
                #${AMRW_SRC_DIR}/physics/mms/MMSForcing.cpp
                ${AMRW_SRC_DIR}/overset/TiogaInterface.cpp
                ${AMRW_SRC_DIR}/immersed_boundary/IB.cpp
                ${AMRW_SRC_DIR}/immersed_boundary/bluff_body/bluff_body_ops.cpp
                ${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Box.cpp
                ${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Cylinder.cpp
                ${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Sphere.cpp
                ${AMRW_SRC_DIR}/mesh_mapping_models/ChannelFlowMap.cpp
                ${AMRW_SRC_DIR}/mesh_mapping_models/ConstantMap.cpp
                ${AMRW_SRC_DIR}/ocean_waves/OceanWaves.cpp
                ${AMRW_SRC_DIR}/ocean_waves/relaxation_zones/relaxation_zones_ops.cpp
                ${AMRW_SRC_DIR}/ocean_waves/relaxation_zones/LinearWaves.cpp
                ${AMRW_SRC_DIR}/ocean_waves/relaxation_zones/StokesWaves.cpp
                ${AMRW_SRC_DIR}/ocean_waves/relaxation_zones/HOSWaves.cpp
                )

  target_include_directories(${erf_lib_name} SYSTEM PUBLIC
                            ${AMRW_SRC_DIR}/../
                            )

  if(ERF_ENABLE_WARM_NO_PRECIP)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_WARM_NO_PRECIP)
  endif()

  if(ERF_ENABLE_POISSON_SOLVE)
    target_sources(${erf_lib_name} PRIVATE
                   ${ERF_SRC_DIR}/Utils/ERF_PoissonSolve.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_POISSON_SOLVE)
  endif()

  if(ERF_ENABLE_PARTICLES)
    target_sources(${erf_lib_name} PRIVATE
                   ${ERF_SRC_DIR}/Particles/TracerPC.cpp
                   ${ERF_SRC_DIR}/Particles/HydroPC.cpp)
    target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Particles)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_PARTICLES)
  endif()

  if(ERF_ENABLE_NETCDF)
    target_sources(${erf_lib_name} PRIVATE
                   ${ERF_SRC_DIR}/IO/NCInterface.cpp
                   ${ERF_SRC_DIR}/IO/NCPlotFile.cpp
                   ${ERF_SRC_DIR}/IO/NCCheckpoint.cpp
                   ${ERF_SRC_DIR}/IO/NCMultiFabFile.cpp
                   ${ERF_SRC_DIR}/IO/ReadFromMetgrid.cpp
                   ${ERF_SRC_DIR}/IO/ReadFromWRFBdy.cpp
                   ${ERF_SRC_DIR}/IO/ReadFromWRFInput.cpp
                   ${ERF_SRC_DIR}/IO/NCColumnFile.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NETCDF)
  endif()

  if(ERF_ENABLE_RRTMGP)
    target_sources(${erf_lib_name} PRIVATE
                   ${ERF_SRC_DIR}/Radiation/Init_rrtmgp.cpp
                   ${ERF_SRC_DIR}/Radiation/Finalize_rrtmgp.cpp
                   ${ERF_SRC_DIR}/Radiation/Run_longwave_rrtmgp.cpp
                   ${ERF_SRC_DIR}/Radiation/Run_shortwave_rrtmgp.cpp
                   ${ERF_SRC_DIR}/Radiation/Cloud_rad_props.cpp
                   ${ERF_SRC_DIR}/Radiation/Aero_rad_props.cpp
                   ${ERF_SRC_DIR}/Radiation/Optics.cpp
                   ${ERF_SRC_DIR}/Radiation/Radiation.cpp
                   ${ERF_SRC_DIR}/Radiation/Albedo.cpp
                   ${ERF_HOME}/Submodules/RRTMGP/cpp/examples/mo_load_coefficients.cpp
                   ${ERF_HOME}/Submodules/RRTMGP/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
                  )

    # The interface code needs to know about the RRTMGP includes
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_RRTMGP)

    target_include_directories(${erf_lib_name} SYSTEM PUBLIC
                               ${ERF_HOME}/Submodules/RRTMGP/cpp/extensions/fluxes_byband
                               ${ERF_HOME}/Submodules/RRTMGP/cpp/extensions/cloud_optics
                               ${ERF_HOME}/Submodules/RRTMGP/cpp/examples
                              )
  endif()

  if(ERF_ENABLE_HDF5)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_HDF5)
  endif()

  target_sources(${erf_lib_name}
     PRIVATE
       ${ERF_SRC_DIR}/Derive.cpp
       ${ERF_SRC_DIR}/ERF.cpp
       ${ERF_SRC_DIR}/ERF_make_new_level.cpp
       ${ERF_SRC_DIR}/ERF_Tagging.cpp
       ${ERF_SRC_DIR}/Advection/AdvectionSrcForMom.cpp
       ${ERF_SRC_DIR}/Advection/AdvectionSrcForState.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/ABLMost.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/MOSTAverage.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_cons.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_xvel.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_yvel.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_zvel.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_bndryreg.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_wrfbdy.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/BoundaryConditions_metgrid.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/ERF_FillPatch.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/ERF_FillPatcher.cpp
       ${ERF_SRC_DIR}/BoundaryConditions/ERF_PhysBCFunct.cpp
       ${ERF_SRC_DIR}/Diffusion/DiffusionSrcForMom_N.cpp
       ${ERF_SRC_DIR}/Diffusion/DiffusionSrcForMom_T.cpp
       ${ERF_SRC_DIR}/Diffusion/DiffusionSrcForState_N.cpp
       ${ERF_SRC_DIR}/Diffusion/DiffusionSrcForState_T.cpp
       ${ERF_SRC_DIR}/Diffusion/ComputeStress_N.cpp
       ${ERF_SRC_DIR}/Diffusion/ComputeStress_T.cpp
       ${ERF_SRC_DIR}/Diffusion/ComputeStrain_N.cpp
       ${ERF_SRC_DIR}/Diffusion/ComputeStrain_T.cpp
       ${ERF_SRC_DIR}/Diffusion/ComputeTurbulentViscosity.cpp
       ${ERF_SRC_DIR}/Diffusion/NumericalDiffusion.cpp
       ${ERF_SRC_DIR}/Diffusion/PBLModels.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_custom.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_from_hse.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_from_input_sounding.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_from_wrfinput.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_from_metgrid.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init_uniform.cpp
       ${ERF_SRC_DIR}/Initialization/ERF_init1d.cpp
       ${ERF_SRC_DIR}/IO/Checkpoint.cpp
       ${ERF_SRC_DIR}/IO/ERF_ReadBndryPlanes.cpp
       ${ERF_SRC_DIR}/IO/ERF_WriteBndryPlanes.cpp
       ${ERF_SRC_DIR}/IO/ERF_Write1DProfiles.cpp
       ${ERF_SRC_DIR}/IO/ERF_WriteScalarProfiles.cpp
       ${ERF_SRC_DIR}/IO/Plotfile.cpp
       ${ERF_SRC_DIR}/IO/writeJobInfo.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_ComputeTimestep.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_Advance.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_TimeStep.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_advance_dycore.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_advance_microphysics.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_advance_radiation.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_make_buoyancy.cpp
       #${ERF_SRC_DIR}/TimeIntegration/ERF_make_condensation_source.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_make_fast_coeffs.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_slow_rhs_pre.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_ApplySpongeZoneBCs.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_slow_rhs_post.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_fast_rhs_N.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_fast_rhs_T.cpp
       ${ERF_SRC_DIR}/TimeIntegration/ERF_fast_rhs_MT.cpp
       ${ERF_SRC_DIR}/Utils/MomentumToVelocity.cpp
       ${ERF_SRC_DIR}/Utils/TerrainMetrics.cpp
       ${ERF_SRC_DIR}/Utils/VelocityToMomentum.cpp
       ${ERF_SRC_DIR}/Utils/InteriorGhostCells.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/Init_SAM.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/Cloud_SAM.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/IceFall.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/Precip.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/PrecipFall.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/Diagnose_SAM.cpp
       ${ERF_SRC_DIR}/Microphysics/SAM/Update_SAM.cpp
       ${ERF_SRC_DIR}/Microphysics/Kessler/Init_Kessler.cpp
       ${ERF_SRC_DIR}/Microphysics/Kessler/Kessler.cpp
       ${ERF_SRC_DIR}/Microphysics/Kessler/Diagnose_Kessler.cpp
       ${ERF_SRC_DIR}/Microphysics/Kessler/Update_Kessler.cpp
       ${ERF_SRC_DIR}/Microphysics/FastEddy/Init_FE.cpp
       ${ERF_SRC_DIR}/Microphysics/FastEddy/FastEddy.cpp
       ${ERF_SRC_DIR}/Microphysics/FastEddy/Diagnose_FE.cpp
       ${ERF_SRC_DIR}/Microphysics/FastEddy/Update_FE.cpp
  )

  if(NOT "${erf_exe_name}" STREQUAL "erf_unit_tests")
    target_sources(${erf_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${erf_lib_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(ERF_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${erf_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${erf_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(ERF_ENABLE_RRTMGP)
    target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Radiation)
  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ERF include directories
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Advection)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/BoundaryConditions)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/DataStructs)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Diffusion)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Initialization)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/IO)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/TimeIntegration)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Utils)
  target_include_directories(${erf_lib_name} PUBLIC ${CMAKE_BINARY_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Microphysics)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Microphysics/Null)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Microphysics/SAM)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Microphysics/Kessler)
  target_include_directories(${erf_lib_name} PUBLIC ${ERF_SRC_DIR}/Microphysics/FastEddy)

  if(ERF_ENABLE_RRTMGP)
     target_link_libraries(${erf_lib_name} PUBLIC yakl)
     target_link_libraries(${erf_lib_name} PUBLIC rrtmgp)
  endif()

  #Link to amrex library
  target_link_libraries_system(${erf_lib_name} PUBLIC amrex)
  target_link_libraries_system(${erf_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api)
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  #Define what we want to be installed during a make install
  install(TARGETS ${erf_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_erf_lib)

function(build_erf_exe erf_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(ERF_SRC_DIR  ${ERF_HOME}/Source)

  target_link_libraries(${erf_exe_name} PRIVATE ${erf_lib_name} AMReX-Hydro::amrex_hydro_api)
  target_link_libraries(${erf_exe_name}  PUBLIC ${erf_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_exe_name})

  target_sources(${erf_exe_name} PRIVATE
                ${ERF_SRC_DIR}/Initialization/ERF_init_bcs.cpp
                )

  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${erf_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
