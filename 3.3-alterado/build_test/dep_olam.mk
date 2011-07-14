#Makefile include dep_olam.mk

$(ARC)($(MODEL_MODS)/mem_rayf.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(MODEL_MODS)/misc_coms.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(MODEL_MODS)/oname_coms.o): \
                     $(ARC)($(ED)/mem_ed.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(MODEL_MODS)/mem_grid.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(MODEL_MODS)/mem_ijtabs.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(MODEL_MODS)/mem_zonavg.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(MODEL_MODS)/mem_basic.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(CONVECT)/mem_cuparm.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(MODEL_MODS)/oplot_coms.o): \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(MODEL_MODS)/var_tables.o): \
                     $(ARC)($(MODEL_MODS)/max_dims.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(ED)/ed_structure_defs.o): \
                     $(ARC)($(ED)/mem_ed.o) \
                     $(ARC)($(ED)/offline_coms.o) \
                     $(ARC)($(ED)/fusion_fission_coms.o) \
                     $(ARC)($(ED)/c34constants.o) \
                     $(ARC)($(ED)/disturbance_coms.o)

$(ARC)($(ED)/canopy_radiation_coms.o): \
                     $(ARC)($(ED)/mem_ed.o)

$(ARC)($(ED)/lphys_interface.o): \
                     $(ARC)($(ED)/c34constants.o)  \
                     $(ARC)($(ED)/pft_coms.o) \
                     $(ARC)($(ED)/ed_options.o) \
                     $(ARC)($(ED)/ed_structure_defs.o)

$(ARC)($(LEAF)/leaf_coms.o): \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(LEAF)/mem_leaf.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(LEAF)/leaf3_interface.o): \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(SEA)/seaice_database_read.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(SEA)/sea_coms.o): \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(SEA)/mem_sea.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(SEA)/var_tables.o)

$(ARC)($(MODEL_MODS)/mem_micro.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(MODEL_MODS)/mem_radiate.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(MODEL_MODS)/mem_addsc.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(MODEL_MODS)/mem_sflux.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(MODEL_MODS)/mem_turb.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(MODEL_MODS)/mem_tend.o): \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/mem_addsc.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)

$(ARC)($(MODEL_MODS)/mem_nudge.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(OISAN)/isan_coms.o): \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(OMODEL)/alloc.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/max_dims.o)  \
                     $(ARC)($(RADIATE)/aero_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(CONVECT)/mem_cuparm.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/mem_radiate.o)  \
                     $(ARC)($(MODEL_MODS)/mem_addsc.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/fill_itabs.o): \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/triangle_geometry.o): \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/icosahedron.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)

$(ARC)($(OMODEL)/spring_dynamics.o): \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/cartesian.o): \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)

$(ARC)($(OMODEL)/spawn_nest.o): \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/triangle_utils.o): \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/fill_jtabs.o): \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/omic_coll.o): \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/omic_gamma.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/omic_driv.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/omic_init.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/omic_misc.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/omic_nuc.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/omic_tabs.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/omic_vap.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/oconv.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(CONVECT)/mem_cuparm.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/conv_coms.o)

$(ARC)($(OMODEL)/rshcupar.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)

$(ARC)($(OMODEL)/cup_dn.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/grell_deep_coms.o)

$(ARC)($(OMODEL)/cup_up.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/cup_env_capmax_ens.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/cup_grell_deep_capmax_ens.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/mem_turb.o) \
                     $(ARC)($(MODEL_MODS)/grell_shallow_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/grell_coms.o)

$(ARC)($(OMODEL)/cup_grell_shallow_capmax_ens.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/mem_basic.o) \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/grell_coms.o) \
                     $(ARC)($(MODEL_MODS)/grell_deep_coms.o)

$(ARC)($(OMODEL)/grell_master.o): \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/grell_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(CONVECT)/mem_cuparm.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)

$(ARC)($(OMODEL)/obnd.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/ohhi.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/olhi.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/mem_zonavg.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_mclat.o)  \
                     $(ARC)($(MODEL_MODS)/mem_radiate.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/olam_grid.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oname_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/olam_run.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_rayf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o) \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(OISAN)/isan_coms.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/modsched.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)

$(ARC)($(OMODEL)/oname.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(MODEL_MODS)/isan_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)  \
                     $(ARC)($(ED)/mem_ed.o)  \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_rayf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oname_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/oname_check.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/max_dims.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oname_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/othrm.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/timestep.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/para_wrtu.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(OMODEL)/para_init.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_para.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OUTILS)/hdf5_utils.o): \
		    $(ARC)($(MODEL_MODS)/rastro_evts.o)\
		    $(ARC)($(OUTILS)/rastro.o) \
		    $(ARC)($(OUTILS)/rastro_f.o)

$(ARC)($(OMODEL)/ocio.o): \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(OMODEL)/analysis_write.o): \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_para.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(OMODEL)/history_write.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(OMODEL)/history_start.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(RADIATE)/rad_driv.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_radiate.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(RADIATE)/mem_harr.o)  \
                     $(ARC)($(MODEL_MODS)/mem_mclat.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(RADIATE)/harr_radinit.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(RADIATE)/rad_mclat.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(RADIATE)/mem_mclat.o)  \
                     $(ARC)($(RADIATE)/mem_radiate.o)

$(ARC)($(OMODEL)/scalar_transport.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(OMODEL)/surface_fluxes.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(CONVECT)/mem_cuparm.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(OMODEL)/turb_k.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)

$(ARC)($(OMODEL)/thiltend_long.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/veltend_long.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/prog_wrtu.o): \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_rayf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/massflux.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)

$(ARC)($(OMODEL)/oplot_interface.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oname_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(MODEL_MODS)/plotcolors.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/oplot_lib.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(CONVECT)/mem_cuparm.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/mem_radiate.o)  \
                     $(ARC)($(MODEL_MODS)/mem_addsc.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/mem_turb.o)  \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(OMODEL)/tileslab.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/oname_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(MODEL_MODS)/plotcolors.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(OMODEL)/contslab.o): \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)

$(ARC)($(OMODEL)/vectslab.o): \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(OMODEL)/olam_mpi_land.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(OMODEL)/olam_mpi_sea.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(OMODEL)/olam_mpi.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_para.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)

$(ARC)($(OMODEL)/olamplot.o): \
                     $(ARC)($(MODEL_MODS)/plotcolors.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)

$(ARC)($(OMODEL)/o_ncar.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_para.o)  \
                     $(ARC)($(MODEL_MODS)/plotcolors.o)

$(ARC)($(OMODEL)/mksfc_topo.o): \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_para.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(OMODEL)/topo_database.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(ED)/pft_coms.o): \
                     $(ARC)($(ED)/mem_ed.o)  \
                     $(ARC)($(ED)/disturbance_coms.o)

$(ARC)($(ED)/allometry.o): \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(ED)/decomposition_coms.o): \
                     $(ARC)($(ED)/mem_ed.o)

$(ARC)($(ED)/decomposition.o): \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/leaf_coms.o) \
                     $(ARC)($(ED)/pft_coms.o) \
                     $(ARC)($(ED)/decomposition_coms.o) \
                     $(ARC)($(ED)/ed_options.o)

$(ARC)($(ED)/ed_canopy_update.o): \
                     $(ARC)($(ED)/canopy_radiation_coms.o)  \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)

$(ARC)($(ED)/ed_database.o): \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)

$(ARC)($(ED)/edio.o): \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/var_tables.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/decomposition_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)

$(ARC)($(ED)/ed_params.o): \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/canopy_radiation_coms.o) \
                     $(ARC)($(ED)/decomposition_coms.o)

$(ARC)($(ED)/ed_startup.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(ED)/mem_ed.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(ED)/ed_init_state_vars.o): \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(ED)/mem_ed.o) \
                     $(ARC)($(MODEL_MODS)/consts_coms.o) \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(ED)/phenology.o): \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/phenology_coms.o)  \
                     $(ARC)($(ED)/decomposition_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(ED)/mortality.o): \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(ED)/pft_coms.o)

$(ARC)($(ED)/fire.o): \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/leaf_coms.o) \
                     $(ARC)($(ED)/pft_coms.o)

$(ARC)($(ED)/forestry.o): \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/leaf_coms.o) \
                     $(ARC)($(ED)/mem_ed.o)

$(ARC)($(ED)/disturbance.o): \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(ED)/ed_options.o) \
                     $(ARC)($(ED)/pft_coms.o) \
                     $(ARC)($(LEAF)/leaf_coms.o) \
                     $(ARC)($(ED)/mem_ed.o) \
                     $(ARC)($(ED)/decomposition_coms.o) \
                     $(ARC)($(LEAF)/mem_leaf.o)

$(ARC)($(ED)/landuse_init.o): \
                     $(ARC)($(ED)/disturbance_coms.o)  \
                     $(ARC)($(ED)/ed_options.o) \
                     $(ARC)($(MODEL_MODS)/consts_coms.o) \
                     $(ARC)($(LEAF)/mem_leaf.o) \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(ED)/vegetation_dynamics.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(ED)/ed_options.o) \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(ED)/growth_balive.o): \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/decomposition_coms.o)  \
                     $(ARC)($(ED)/pft_coms.o)

$(ARC)($(ED)/structural_growth.o): \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(ED)/decomposition_coms.o) \
                     $(ARC)($(ED)/pft_coms.o)

$(ARC)($(ED)/reproduction.o): \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(ED)/pft_coms.o) \
                     $(ARC)($(ED)/decomposition_coms.o) \
                     $(ARC)($(LEAF)/leaf_coms.o) \
                     $(ARC)($(ED)/fusion_fission_coms.o)

$(ARC)($(ED)/fis_fus_utils.o): \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(ED)/pft_coms.o)  \
                     $(ARC)($(ED)/mem_ed.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(ED)/decomposition_coms.o)  \
                     $(ARC)($(ED)/fusion_fission_coms.o)

$(ARC)($(LEAF)/leaf3_startup.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(LEAF)/misc_coms.o)  \
                     $(ARC)($(LEAF)/consts_coms.o)  \
                     $(ARC)($(LEAF)/var_tables.o)  \
                     $(ARC)($(LEAF)/max_dims.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(LEAF)/leaf3_init_atm.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_micro.o) \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/ed_options.o)

$(ARC)($(LEAF)/leaf3.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(LEAF)/leaf3_interface.o) \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(LEAF)/leaf_plot.o): \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/oplot_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(ED)/lphys_full.o): \
                     $(ARC)($(ED)/c34constants.o)  \
                     $(ARC)($(ED)/pft_coms.o) \
                     $(ARC)($(ED)/ed_options.o) \
                     $(ARC)($(ED)/ed_structure_defs.o) \
                     $(ARC)($(ED)/lphys_interface.o)

$(ARC)($(LEAF)/makesfc.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(LEAF)/mem_mksfc.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(LEAF)/leaf_database.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(ED)/ed_offline_met.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(ED)/ed_options.o)  \
                     $(ARC)($(ED)/offline_coms.o)  \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)

$(ARC)($(LEAF)/ndvi_sfcfile_read.o): \
                     $(ARC)($(LEAF)/leaf_coms.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o) \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(ED)/sfcrad_ed.o): \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(ED)/canopy_radiation_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(LEAF)/leaf_coms.o)

$(ARC)($(ED)/radiate_offline.o): \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(ED)/canopy_radiation_coms.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(ED)/sfluxes_offline.o): \
                     $(ARC)($(ED)/ed_structure_defs.o)  \
                     $(ARC)($(LEAF)/mem_leaf.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(SEA)/sea_init_atm.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_sflux.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(SEA)/sea.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)

$(ARC)($(SEA)/sea_startup.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)

$(ARC)($(SEA)/sst_sfcfile_read.o): \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(SEA)/mem_sea.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(SEA)/sst_database_read.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)

$(ARC)($(SEA)/sst_database.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(OUTILS)/hdf5_utils.o)  \
                     $(ARC)($(LEAF)/mem_mksfc.o)  \
                     $(ARC)($(SEA)/sea_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(OISAN)/isan_driver.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_zonavg.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/isan_coms.o)

$(ARC)($(OISAN)/astp.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/isan_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(OISAN)/asti.o): \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/isan_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o) \
                     $(ARC)($(MODEL_MODS)/mem_zonavg.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(MODEL_MODS)/max_dims.o)

$(ARC)($(OISAN)/fldsisan.o): \
		     $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(MODEL_MODS)/mem_nudge.o)  \
                     $(ARC)($(MODEL_MODS)/mem_basic.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/consts_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_ijtabs.o)  \
                     $(ARC)($(MODEL_MODS)/mem_micro.o)  \
                     $(ARC)($(MODEL_MODS)/micro_coms.o)  \
                     $(ARC)($(MODEL_MODS)/mem_tend.o)  \
                     $(ARC)($(MODEL_MODS)/isan_coms.o)

$(ARC)($(OISAN)/file_inv.o): \
                     $(ARC)($(MODEL_MODS)/isan_coms.o)  \
                     $(ARC)($(MODEL_MODS)/misc_coms.o)  \
                     $(ARC)($(MODEL_MODS)/max_dims.o)


$(ARC)($(RADIATE)/calipso_output.o):  \
                     $(ARC)($(RADIATE)/fuinput.o) \
                     $(ARC)($(RADIATE)/fuoutput.o)

$(ARC)($(RADIATE)/fuinput.o):  \
                     $(ARC)($(RADIATE)/extras.o) \
                     $(ARC)($(RADIATE)/taucorr.o)

$(ARC)($(RADIATE)/fuoutput.o):  \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/fuprint.o):  \
                     $(ARC)($(RADIATE)/extras.o) \
                     $(ARC)($(RADIATE)/fuinput.o) \
                     $(ARC)($(RADIATE)/fuoutput.o)

$(ARC)($(RADIATE)/gflq.o):  \
                     $(ARC)($(RADIATE)/extras.o) \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/ma_tip.o):  \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/rad_multi_0403.o):  \
                     $(ARC)($(RADIATE)/uvcor_all.o) \
                     $(ARC)($(RADIATE)/vla.o) \
                     $(ARC)($(RADIATE)/icedirsfc.o) \
                     $(ARC)($(RADIATE)/fuinput.o) \
                     $(ARC)($(RADIATE)/fuoutput.o)

$(ARC)($(RADIATE)/seiji_k2b.o):  \
                     $(ARC)($(RADIATE)/sktbl_ht02a.o) \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/seiji_solver_0403.o):  \
                     $(ARC)($(RADIATE)/fuoutput.o) \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/uvcor_all.o):  \
                     $(ARC)($(RADIATE)/fuoutput.o) \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/vla.o):  \
                     $(ARC)($(RADIATE)/fuinput.o)

$(ARC)($(RADIATE)/fuliou_init.o):  \
		    $(ARC)($(MODEL_MODS)/rastro_evts.o)\
                     $(ARC)($(RADIATE)/aero_coms.o)

$(ARC)($(RADIATE)/aerosols.o):  \
                     $(ARC)($(RADIATE)/aero_coms.o)

$(ARC)($(RADIATE)/fuliou_raddriv.o):  \
                     $(ARC)($(RADIATE)/aero_coms.o) \
                     $(ARC)($(MODEL_MODS)/misc_coms.o) \
                     $(ARC)($(MODEL_MODS)/micro_coms.o) \
                     $(ARC)($(RADIATE)/mem_radiate.o) \
                     $(ARC)($(MODEL_MODS)/mem_micro.o) \
                     $(ARC)($(MODEL_MODS)/consts_coms.o) \
                     $(ARC)($(MODEL_MODS)/mem_basic.o) \
                     $(ARC)($(MODEL_MODS)/mem_grid.o) \
                     $(ARC)($(RADIATE)/extras.o) \
                     $(ARC)($(RADIATE)/rad_multi_0403.o) \
                     $(ARC)($(RADIATE)/fuinput.o) \
                     $(ARC)($(RADIATE)/fuoutput.o) \
                     $(ARC)($(RADIATE)/vla.o) \
                     $(ARC)($(RADIATE)/gflq.o) \
                     $(ARC)($(RADIATE)/calipso_output.o)


