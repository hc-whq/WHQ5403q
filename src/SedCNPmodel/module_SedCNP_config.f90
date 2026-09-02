!=====||__WHQ5403q__||=====!
!
module SedCNP_config
   ! Single choke point through which every SedCNPmodel source file reaches
   ! config_base's SedCNPmodel namelist state. Every file that instead did its
   ! own independent "use config_base, only: SedCNPmodel" produced its own
   ! gfortran-derived snapshot of that derived type; combining two or more such
   ! independently-derived snapshots in one scope (e.g. inside
   ! module_SedCNPmodel_driver.f90, which pulls in several of these files at
   ! once) trips a gfortran diamond-import "Mismatch in components of derived
   ! type" error. Routing everyone through one shared .mod avoids that.
   use config_base, only: SedCNPmodel
   implicit none
end module SedCNP_config
!
!=====||__WHQ5403q__||=====!
