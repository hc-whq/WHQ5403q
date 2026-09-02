cmake_minimum_required (VERSION 2.8)

add_library(hydro_SedCNPmodel STATIC
#=====||__WHQ5403q__||=====!
#
	module_SedCNP_config.f90
#
#=====||__WHQ5403q__||=====!
	module_SedCNPvariables.f90
	module_AlloDeallocate.f90
	module_overlandSed.f90
	module_channelSed.f90
	module_WriteCNPini.f90
	module_CNPfunctions.f90
	module_CNPparams.f90
	module_SedCNP_in.f90
        module_SedCNP_out.f90
	module_Carbon.f90
	module_Nitrogen.f90
        module_Phosphorus.f90
	module_CNPupdates.f90
	module_ReadCNPini.f90
	module_CNPmain.f90
	module_SedCNPmodel_driver.f90
)

