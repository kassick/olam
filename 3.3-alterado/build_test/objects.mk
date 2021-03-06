#Makefile objects.mk

# Define main source.

MAIN = $(OMODEL)/olammain.f90

MPIF = $(OMODEL)/olam_mpi.F90 \
       $(OMODEL)/olam_mpi_sea.F90 \
       $(OMODEL)/olam_mpi_land.F90

# Define objects.

OBJ = $(ARC)($(MODEL_MODS)/max_dims.o) \
      $(ARC)($(MODEL_MODS)/misc_coms.o) \
      $(ARC)($(MODEL_MODS)/consts_coms.o) \
      $(ARC)($(MODEL_MODS)/mem_para.o) \
      $(ARC)($(MODEL_MODS)/plotcolors.o) \
      $(ARC)($(OUTILS)/hdf5_utils.o) \
      $(ARC)($(MODEL_MODS)/micro_coms.o) \
      $(ARC)($(MODEL_MODS)/oplot_coms.o) \
      $(ARC)($(OISAN)/isan_coms.o) \
      $(ARC)($(MODEL_MODS)/mem_zonavg.o) \
      $(ARC)($(MODEL_MODS)/mem_rayf.o) \
      $(ARC)($(MODEL_MODS)/mem_mclat.o) \
      $(ARC)($(MODEL_MODS)/mem_grid.o) \
      $(ARC)($(MODEL_MODS)/ke_coms.o) \
      $(ARC)($(SEA)/sea_coms.o) \
      $(ARC)($(LEAF)/leaf_coms.o) \
      $(ARC)($(RADIATE)/mem_harr.o) \
      $(ARC)($(MODEL_MODS)/var_tables.o) \
      $(ARC)($(MODEL_MODS)/mem_addsc.o) \
      $(ARC)($(MODEL_MODS)/mem_basic.o)  \
      $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
      $(ARC)($(MODEL_MODS)/mem_micro.o) \
      $(ARC)($(MODEL_MODS)/mem_nudge.o) \
      $(ARC)($(MODEL_MODS)/mem_turb.o) \
      $(ARC)($(CONVECT)/mem_cuparm.o) \
      $(ARC)($(CONVECT)/kuo_coms.o) \
      $(ARC)($(CONVECT)/grell_coms.o) \
      $(ARC)($(CONVECT)/grell_deep_coms.o) \
      $(ARC)($(CONVECT)/grell_shallow_coms.o) \
      $(ARC)($(ED)/mem_ed.o) \
      $(ARC)($(ED)/disturbance_coms.o) \
      $(ARC)($(ED)/fusion_fission_coms.o) \
      $(ARC)($(ED)/c34constants.o) \
      $(ARC)($(ED)/offline_coms.o) \
      $(ARC)($(ED)/ed_structure_defs.o) \
      $(ARC)($(MODEL_MODS)/oname_coms.o) \
      $(ARC)($(ED)/pft_coms.o) \
      $(ARC)($(ED)/canopy_radiation_coms.o) \
      $(ARC)($(ED)/decomposition_coms.o) \
      $(ARC)($(ED)/ed_options.o) \
      $(ARC)($(ED)/phenology_coms.o) \
      $(ARC)($(ED)/lphys_interface.o) \
      $(ARC)($(MODEL_MODS)/mem_tend.o) \
      $(ARC)($(RADIATE)/mem_radiate.o) \
      $(ARC)($(RADIATE)/fu_liou/aero_coms.o) \
      $(ARC)($(LEAF)/mem_leaf.o) \
      $(ARC)($(LEAF)/mem_mksfc.o) \
      $(ARC)($(LEAF)/leaf3_interface.o) \
      $(ARC)($(SEA)/mem_sea.o) \
      $(ARC)($(MODEL_MODS)/mem_sflux.o) \
      $(ARC)($(LEAF)/leaf_database_read.o) \
      \
      $(ARC)($(OMODEL)/alloc.o)   \
      $(ARC)($(OMODEL)/para_init.o)   \
      $(ARC)($(LEAF)/para_init_land.o)   \
      $(ARC)($(SEA)/para_init_sea.o)   \
      \
      $(ARC)($(OMODEL)/fill_itabs.o)   \
      $(ARC)($(OMODEL)/triangle_geometry.o)   \
      $(ARC)($(OMODEL)/icosahedron.o)   \
      $(ARC)($(OMODEL)/spring_dynamics.o)   \
      $(ARC)($(OMODEL)/cartesian.o)   \
      $(ARC)($(OMODEL)/spawn_nest.o)   \
      $(ARC)($(OMODEL)/triangle_utils.o)   \
      $(ARC)($(OMODEL)/omic_coll.o) \
      $(ARC)($(OMODEL)/omic_driv.o) \
      $(ARC)($(OMODEL)/omic_gamma.o) \
      $(ARC)($(OMODEL)/omic_init.o) \
      $(ARC)($(OMODEL)/omic_misc.o) \
      $(ARC)($(OMODEL)/omic_nuc.o)  \
      $(ARC)($(OMODEL)/omic_tabs.o) \
      $(ARC)($(OMODEL)/omic_vap.o) \
      $(ARC)($(OMODEL)/obnd.o) \
      $(ARC)($(OMODEL)/ohhi.o)  \
      $(ARC)($(OMODEL)/olhi.o)  \
      $(ARC)($(OMODEL)/olam_grid.o) \
      $(ARC)($(OMODEL)/olam_run.o) \
      $(ARC)($(OMODEL)/modsched.o) \
      $(ARC)($(OMODEL)/oname.o) \
      $(ARC)($(OMODEL)/oname_check.o) \
      $(ARC)($(OMODEL)/othrm.o) \
      $(ARC)($(OMODEL)/timestep.o) \
      $(ARC)($(OMODEL)/ocio.o) \
      $(ARC)($(OMODEL)/history_write.o) \
      $(ARC)($(OMODEL)/history_start.o) \
      $(ARC)($(OMODEL)/scalar_transport.o) \
      $(ARC)($(OMODEL)/surface_fluxes.o) \
      $(ARC)($(OMODEL)/turb_k.o) \
      $(ARC)($(OMODEL)/thiltend_long.o) \
      $(ARC)($(OMODEL)/veltend_long.o) \
      $(ARC)($(OMODEL)/prog_wrtu.o) \
      $(ARC)($(OMODEL)/massflux.o) \
      $(ARC)($(OMODEL)/oplot_interface.o) \
      $(ARC)($(OMODEL)/oplot_lib.o) \
      $(ARC)($(OMODEL)/tileslab.o) \
      $(ARC)($(OMODEL)/contslab.o) \
      $(ARC)($(OMODEL)/vectslab.o) \
      $(ARC)($(OMODEL)/olamplot.o) \
      $(ARC)($(OMODEL)/o_ncar.o) \
      $(ARC)($(OMODEL)/topo_database_read.o) \
      \
      $(ARC)($(RADIATE)/rad_driv.o) \
      $(ARC)($(RADIATE)/ccmp_raddriv.o) \
      $(ARC)($(RADIATE)/cc_rad.o) \
      $(ARC)($(RADIATE)/mp_rad.o) \
      $(ARC)($(RADIATE)/harr_radinit.o) \
      $(ARC)($(RADIATE)/harr_raddriv.o) \
      $(ARC)($(RADIATE)/harr_rad.o) \
      $(ARC)($(RADIATE)/rad_mclat.o) \
      $(ARC)($(RADIATE)/rad_stable.o) \
      $(ARC)($(RADIATE)/fu_liou/extras.o) \
      $(ARC)($(RADIATE)/fu_liou/taucorr.o) \
      $(ARC)($(RADIATE)/fu_liou/fuinput.o) \
      $(ARC)($(RADIATE)/fu_liou/fuoutput.o) \
      $(ARC)($(RADIATE)/fu_liou/calipso_output.o) \
      $(ARC)($(RADIATE)/fu_liou/fuprint.o) \
      $(ARC)($(RADIATE)/fu_liou/ma_tip.o) \
      $(ARC)($(RADIATE)/fu_liou/gflq.o) \
      $(ARC)($(RADIATE)/fu_liou/uvcor_all.o) \
      $(ARC)($(RADIATE)/fu_liou/vla.o) \
      $(ARC)($(RADIATE)/fu_liou/icedirsfc.o) \
      $(ARC)($(RADIATE)/fu_liou/rad_multi_0403.o) \
      $(ARC)($(RADIATE)/fu_liou/sktbl_ht02a.o) \
      $(ARC)($(RADIATE)/fu_liou/seiji_k2b.o) \
      $(ARC)($(RADIATE)/fu_liou/seiji_solver_0403.o) \
      $(ARC)($(RADIATE)/fu_liou/aerosols_0403.o) \
      $(ARC)($(RADIATE)/fu_liou/aqua_wnflt_0404.o) \
      $(ARC)($(RADIATE)/fu_liou/misc_0403.o) \
      $(ARC)($(RADIATE)/fu_liou/seiji_twostreamsolv_sw_v21.09022004.o) \
      $(ARC)($(RADIATE)/fu_liou/fuliou_raddriv.o) \
      $(ARC)($(RADIATE)/fu_liou/fuliou_init.o) \
      $(ARC)($(RADIATE)/fu_liou/aerosols.o) \
      $(ARC)($(CONVECT)/cuparm_driver.o) \
      $(ARC)($(CONVECT)/kuo.o) \
      $(ARC)($(CONVECT)/grell_downdraft.o) \
      $(ARC)($(CONVECT)/grell_env_capmax_ens.o) \
      $(ARC)($(CONVECT)/grell_deep_capmax_ens.o) \
      $(ARC)($(CONVECT)/grell_shallow_capmax_ens.o) \
      $(ARC)($(CONVECT)/grell_updraft.o) \
      $(ARC)($(ED)/allometry.o) \
      $(ARC)($(ED)/decomposition.o) \
      $(ARC)($(ED)/ed_canopy_update.o) \
      $(ARC)($(ED)/ed_database.o) \
      $(ARC)($(ED)/edio.o) \
      $(ARC)($(ED)/ed_params.o) \
      $(ARC)($(ED)/ed_history.o) \
      $(ARC)($(ED)/ed_startup.o) \
      $(ARC)($(ED)/ed_init_state_vars.o) \
      $(ARC)($(ED)/vegetation_dynamics.o) \
      $(ARC)($(ED)/phenology.o) \
      $(ARC)($(ED)/fire.o) \
      $(ARC)($(ED)/forestry.o) \
      $(ARC)($(ED)/disturbance.o) \
      $(ARC)($(ED)/landuse_init.o) \
      $(ARC)($(ED)/growth_balive.o) \
      $(ARC)($(ED)/structural_growth.o) \
      $(ARC)($(ED)/reproduction.o) \
      $(ARC)($(ED)/ed_offline_met.o) \
      $(ARC)($(ED)/fis_fus_utils.o) \
      $(ARC)($(ED)/mortality.o) \
      $(ARC)($(LEAF)/leaf3_startup.o) \
      $(ARC)($(LEAF)/leaf3_init_atm.o) \
      $(ARC)($(LEAF)/leaf3.o) \
      $(ARC)($(LEAF)/leaf_plot.o) \
      $(ARC)($(ED)/lphys_full.o) \
      $(ARC)($(LEAF)/makesfc.o) \
      $(ARC)($(LEAF)/ndvi_database_read.o) \
      $(ARC)($(ED)/radiate_offline.o) \
      $(ARC)($(ED)/sfluxes_offline.o) \
      $(ARC)($(ED)/sfcrad_ed.o) \
      $(ARC)($(SEA)/sea_init_atm.o) \
      $(ARC)($(SEA)/sea.o) \
      $(ARC)($(SEA)/sea_startup.o) \
      $(ARC)($(SEA)/sst_database_read.o) \
      $(ARC)($(SEA)/seaice_database_read.o) \
      \
      $(ARC)($(OISAN)/isan_driver.o) \
      $(ARC)($(OISAN)/astp.o) \
      $(ARC)($(OISAN)/asti.o) \
      $(ARC)($(OISAN)/fldsisan.o) \
      $(ARC)($(OISAN)/file_inv.o) \
      \
      $(ARC)($(OUTILS)/charutils.o) \
      $(ARC)($(OUTILS)/dateutils.o) \
      $(ARC)($(OUTILS)/filelist.o) \
      $(ARC)($(OUTILS)/hdf5_f2c.o) \
      $(ARC)($(OUTILS)/read_cdc.o) \
      $(ARC)($(OUTILS)/interp_lib.o) \
      $(ARC)($(OUTILS)/map_proj.o) \
      $(ARC)($(OUTILS)/therm_lib.o) \
      $(ARC)($(OUTILS)/utils_f.o) \
      $(ARC)($(OUTILS)/quicksort.o)
