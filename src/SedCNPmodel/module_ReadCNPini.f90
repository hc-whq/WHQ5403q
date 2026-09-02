!> @file module_ReadCNPini.f90
!! @brief This module reads the initial conditions for Carbon (C), Nitrogen (N),
!! and Phosphorus (P) from specified files for all model components
!! (Soil, Aquifer, Groundwater, Channel).
module ReadCNPini

   use module_SedCNPvariables
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,         only: SedCNPmodel
!
!=====||__WHQ5403q__||=====!
   use CNPparams

   contains

   !> @brief Reads initial C, N, and P storages for the soil component.
   !! @details This subroutine reads initial conditions for all soil C, N, and P
   !! pools from specified input files. If a file is not specified, cannot be
   !! opened, or is invalid, it performs a "cold start" by initializing all
   !! relevant pools to zero. Units are converted from kg/ha to kg per grid cell.
   subroutine ReadCNPini_So()

      implicit none

      integer :: i,j, isl, idum
      real    :: dum
      real(8) :: areaxyha      ! Grid cell area in hectares [ha]
      logical :: read_successful
      logical :: message_printed
      integer :: io_status
      character(len=256)  :: filename

      areaxyha = domain%areaxy / 10000.0
      message_printed = .FALSE.

      ! ====================================================================
      !  Read Soil Carbon (C) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      IF (len(trim(SedCNPmodel%CiniSo_file)) > 0) THEN
         filename = trim(SedCNPmodel%CiniSo_file)
         open(591010, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (591010, *, iostat=io_status) ! Read header
            
            IF (io_status == 0) THEN
read_loop_c_so: DO i = 1, domain%ix
                  DO j = 1, domain%jx
                     DO isl = 1, domain%nsl
                        IF (isl .eq. 1) THEN
                           read (591010, *, iostat=io_status) idum, idum, &
                              idum, So_LPOCLIT0(i,j), So_LPOCEXC0(i,j), &
                              So_LPOCMAN0(i,j), So_LPOCRES0(i,isl,j), &
                              So_RPOC0(i,isl,j), So_LDOC0(i,isl,j), &
                              So_RDOC0(i,isl,j), So_MBMC0(i,isl,j)
                           
                           ! Convert from kg-C/ha to kg-C per cell
                           So_LPOCLIT0(i,j) = So_LPOCLIT0(i,j) * areaxyha
                           So_LPOCEXC0(i,j) = So_LPOCEXC0(i,j) * areaxyha
                           So_LPOCMAN0(i,j) = So_LPOCMAN0(i,j) * areaxyha
                           So_LPOCRES0(i,isl,j) = So_LPOCRES0(i,isl,j) * areaxyha
                           So_RPOC0(i,isl,j) = So_RPOC0(i,isl,j) * areaxyha
                           So_LDOC0(i,isl,j) = So_LDOC0(i,isl,j) * areaxyha
                           So_RDOC0(i,isl,j) = So_RDOC0(i,isl,j) * areaxyha
                           So_MBMC0(i,isl,j) = So_MBMC0(i,isl,j) * areaxyha
                        ELSE
                           read (591010, *, iostat=io_status) idum, idum, &
                              idum, dum, dum, dum, &
                              So_LPOCRES0(i,isl,j), So_RPOC0(i,isl,j), &
                              So_LDOC0(i,isl,j), So_RDOC0(i,isl,j), &
                              So_MBMC0(i,isl,j)
                           
                           ! Convert from kg-C/ha to kg-C per cell
                           So_LPOCRES0(i,isl,j) = So_LPOCRES0(i,isl,j) * areaxyha
                           So_RPOC0(i,isl,j) = So_RPOC0(i,isl,j) * areaxyha
                           So_LDOC0(i,isl,j) = So_LDOC0(i,isl,j) * areaxyha
                           So_RDOC0(i,isl,j) = So_RDOC0(i,isl,j) * areaxyha
                           So_MBMC0(i,isl,j) = So_MBMC0(i,isl,j) * areaxyha
                        END IF
                     END DO
                     IF (io_status /= 0) EXIT read_loop_c_so
                  END DO
               END DO read_loop_c_so
               
               IF (io_status == 0) read_successful = .TRUE.
               close (591010)
            END IF
         END IF
      END IF
      
      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Soil C initial conditions not read.', &
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         So_LPOCLIT0 = 0.0_8
         So_LPOCEXC0 = 0.0_8
         So_LPOCMAN0 = 0.0_8
         So_LPOCRES0 = 0.0_8
         So_RPOC0 = 0.0_8
         So_LDOC0 = 0.0_8
         So_RDOC0 = 0.0_8
         So_MBMC0 = 0.0_8
      END IF

      ! ====================================================================
      !  Read Soil Nitrogen (N) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%NiniSo_file)) > 0) THEN
         filename = trim(SedCNPmodel%NiniSo_file)
         open(592010, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (592010, *, iostat=io_status) ! Read header
            
            IF (io_status == 0) THEN
read_loop_n_so: DO i = 1, domain%ix
                  DO j = 1, domain%jx
                     DO isl = 1, domain%nsl
                        IF (isl .eq. 1) THEN
                           read (592010, *, iostat=io_status) idum, idum, &
                              idum, So_LPONLIT0(i,j), So_LPONEXC0(i,j), &
                              So_LPONMAN0(i,j), So_LPONRES0(i,isl,j), &
                              So_RPON0(i,isl,j), So_LDON0(i,isl,j), &
                              So_RDON0(i,isl,j), So_MBMN0(i,isl,j), & 
                              So_NH40(i,isl,j), So_NO30(i,isl,j)
                           
                           ! Convert from kg-N/ha to kg-N per cell
                           So_LPONLIT0(i,j) = So_LPONLIT0(i,j) * areaxyha
                           So_LPONEXC0(i,j) = So_LPONEXC0(i,j) * areaxyha
                           So_LPONMAN0(i,j) = So_LPONMAN0(i,j) * areaxyha
                           So_LPONRES0(i,isl,j) = So_LPONRES0(i,isl,j) * areaxyha
                           So_RPON0(i,isl,j) = So_RPON0(i,isl,j) * areaxyha
                           So_LDON0(i,isl,j) = So_LDON0(i,isl,j) * areaxyha
                           So_RDON0(i,isl,j) = So_RDON0(i,isl,j) * areaxyha
                           So_MBMN0(i,isl,j) = So_MBMN0(i,isl,j) * areaxyha
                           So_NH40(i,isl,j) = So_NH40(i,isl,j) * areaxyha
                           So_NO30(i,isl,j) = So_NO30(i,isl,j) * areaxyha
                        ELSE
                           read (592010, *, iostat=io_status) idum, idum, &
                              idum, dum, dum, dum, So_LPONRES0(i,isl,j), &
                              So_RPON0(i,isl,j), So_LDON0(i,isl,j), &
                              So_RDON0(i,isl,j), So_MBMN0(i,isl,j), &
                              So_NH40(i,isl,j), So_NO30(i,isl,j)
                           
                           ! Convert from kg-N/ha to kg-N per cell
                           So_LPONRES0(i,isl,j) = So_LPONRES0(i,isl,j) * areaxyha
                           So_RPON0(i,isl,j) = So_RPON0(i,isl,j) * areaxyha
                           So_LDON0(i,isl,j) = So_LDON0(i,isl,j) * areaxyha
                           So_RDON0(i,isl,j) = So_RDON0(i,isl,j) * areaxyha
                           So_MBMN0(i,isl,j) = So_MBMN0(i,isl,j) * areaxyha
                           So_NH40(i,isl,j) = So_NH40(i,isl,j) * areaxyha
                           So_NO30(i,isl,j) = So_NO30(i,isl,j) * areaxyha
                        END IF
                     END DO
                     IF (io_status /= 0) EXIT read_loop_n_so
                  END DO
               END DO read_loop_n_so
               
               IF (io_status == 0) read_successful = .TRUE.
               close (592010)
            END IF
         END IF
      END IF

      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Soil N initial conditions not read.', &
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         So_LPONLIT0 = 0.0_8
         So_LPONEXC0 = 0.0_8
         So_LPONMAN0 = 0.0_8
         So_LPONRES0 = 0.0_8
         So_RPON0 = 0.0_8
         So_LDON0 = 0.0_8
         So_RDON0 = 0.0_8
         So_MBMN0 = 0.0_8
         So_NH40 = 0.0_8
         So_NO30 = 0.0_8
      END IF

      ! ====================================================================
      !  Read Soil Phosphorus (P) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%PiniSo_file)) > 0) THEN
         filename = trim(SedCNPmodel%PiniSo_file)
         open(593010, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (593010, *, iostat=io_status) ! Read header
            
            IF (io_status == 0) THEN
read_loop_p_so: DO i = 1, domain%ix
                  DO j = 1, domain%jx
                     DO isl = 1, domain%nsl
                        IF (isl .eq. 1) THEN
                           read (593010, *, iostat=io_status) idum, idum, &
                              idum, So_LPOPLIT0(i,j), So_LPOPEXC0(i,j), &
                              So_LPOPMAN0(i,j), So_LPOPRES0(i,isl,j), &
                              So_RPOP0(i,isl,j), So_LDOP0(i,isl,j), &
                              So_RDOP0(i,isl,j), So_MBMP0(i,isl,j), &
                              So_PO40(i,isl,j), So_PIPA0(i,isl,j), &
                              So_PIPS0(i,isl,j)
                           
                           ! Convert from kg-P/ha to kg-P per cell
                           So_LPOPLIT0(i,j) = So_LPOPLIT0(i,j) * areaxyha
                           So_LPOPEXC0(i,j) = So_LPOPEXC0(i,j) * areaxyha
                           So_LPOPMAN0(i,j) = So_LPOPMAN0(i,j) * areaxyha
                           So_LPOPRES0(i,isl,j) = So_LPOPRES0(i,isl,j) * areaxyha
                           So_RPOP0(i,isl,j) = So_RPOP0(i,isl,j) * areaxyha
                           So_LDOP0(i,isl,j) = So_LDOP0(i,isl,j) * areaxyha
                           So_RDOP0(i,isl,j) = So_RDOP0(i,isl,j) * areaxyha
                           So_MBMP0(i,isl,j) = So_MBMP0(i,isl,j) * areaxyha
                           So_PO40(i,isl,j) = So_PO40(i,isl,j) * areaxyha
                           So_PIPA0(i,isl,j) = So_PIPA0(i,isl,j) * areaxyha
                           So_PIPS0(i,isl,j) = So_PIPS0(i,isl,j) * areaxyha
                        ELSE
                           read (593010, *, iostat=io_status) idum, idum, &
                              idum, dum, dum, dum, &
                              So_LPOPRES0(i,isl,j), So_RPOP0(i,isl,j), &
                              So_LDOP0(i,isl,j), So_RDOP0(i,isl,j), &
                              So_MBMP0(i,isl,j), So_PO40(i,isl,j), &
                              So_PIPA0(i,isl,j), So_PIPS0(i,isl,j)

                           ! Convert from kg-P/ha to kg-P per cell
                           So_LPOPRES0(i,isl,j) = So_LPOPRES0(i,isl,j) * areaxyha
                           So_RPOP0(i,isl,j) = So_RPOP0(i,isl,j) * areaxyha
                           So_LDOP0(i,isl,j) = So_LDOP0(i,isl,j) * areaxyha
                           So_RDOP0(i,isl,j) = So_RDOP0(i,isl,j) * areaxyha
                           So_MBMP0(i,isl,j) = So_MBMP0(i,isl,j) * areaxyha
                           So_PO40(i,isl,j) = So_PO40(i,isl,j) * areaxyha
                           So_PIPA0(i,isl,j) = So_PIPA0(i,isl,j) * areaxyha
                           So_PIPS0(i,isl,j) = So_PIPS0(i,isl,j) * areaxyha
                        END IF
                     END DO
                     IF (io_status /= 0) EXIT read_loop_p_so
                  END DO
               END DO read_loop_p_so
               
               IF (io_status == 0) read_successful = .TRUE.
               close (593010)
            END IF
         END IF
      END IF

      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Soil P initial conditions not read.', &
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         So_LPOPLIT0 = 0.0_8
         So_LPOPEXC0 = 0.0_8
         So_LPOPMAN0 = 0.0_8
         So_LPOPRES0 = 0.0_8
         So_RPOP0 = 0.0_8
         So_LDOP0 = 0.0_8
         So_RDOP0 = 0.0_8
         So_MBMP0 = 0.0_8
         So_PO40 = 0.0_8
         So_PIPA0 = 0.0_8
         So_PIPS0 = 0.0_8
      END IF

   end subroutine ReadCNPini_So


!    !> @brief Reads initial C, N, and P storages for the aquifer component.
!    !! @details This subroutine reads initial conditions for all aquifer C, N,
!    !! and P pools from specified files. If a file is not specified or is
!    !! invalid, it performs a "cold start" by initializing all pools to zero.
!    subroutine ReadCNPini_Aq ()

!       implicit none

!       integer :: i,j, idum
!       logical :: read_successful
!       logical :: message_printed
!       integer :: io_status
!       character(len=256)  :: filename

!       message_printed = .FALSE.

!       ! ====================================================================
!       !  Read Aquifer Carbon (C) Initial Conditions
!       ! ====================================================================
!       read_successful = .FALSE.
!       IF (len(trim(SedCNPmodel%CiniAq_file)) > 0) THEN
!          filename = trim(SedCNPmodel%CiniAq_file)
!          open(591020, file=filename, status='old', action='read', &
!               iostat=io_status)

!          IF (io_status == 0) THEN
!             read (591020, *, iostat=io_status) ! Read header
            
!             IF (io_status == 0) THEN
! read_loop_c_aq: DO i = 1, domain%ix
!                   DO j = 1, domain%jx
!                      read (591020, *, iostat=io_status) idum, idum, &
!                         Aq_LDOC0(i,j), Aq_RDOC0(i,j) ! kg-C
!                      IF (io_status /= 0) EXIT read_loop_c_aq
!                   END DO
!                END DO read_loop_c_aq
               
!                IF (io_status == 0) read_successful = .TRUE.
!                close (591020)
!             END IF
!          END IF
!       END IF

!       ! Cold start if file read failed or not specified
!       IF (.NOT. read_successful) THEN
!          IF (.NOT. message_printed) THEN
!             write(*,*) 'INFO: Aquifer C initial conditions not read.', &
!                       ' Performing cold start.'
!             message_printed = .TRUE.
!          END IF
!          Aq_LDOC0 = 0.0_8
!          Aq_RDOC0 = 0.0_8
!       END IF

!       ! ====================================================================
!       !  Read Aquifer Nitrogen (N) Initial Conditions
!       ! ====================================================================
!       read_successful = .FALSE.
!       message_printed = .FALSE.
!       IF (len(trim(SedCNPmodel%NiniAq_file)) > 0) THEN
!          filename = trim(SedCNPmodel%NiniAq_file)
!          open(592020, file=filename, status='old', action='read', &
!               iostat=io_status)

!          IF (io_status == 0) THEN
!             read (592020, *, iostat=io_status) ! Read header
            
!             IF (io_status == 0) THEN
! read_loop_n_aq: DO i = 1, domain%ix
!                   DO j = 1, domain%jx
!                      read (592020, *, iostat=io_status) idum, idum, &
!                         Aq_LDON0(i,j), Aq_RDON0(i,j), Aq_NH40(i,j), &
!                         Aq_NO30(i,j) ! kg-N 
!                      IF (io_status /= 0) EXIT read_loop_n_aq
!                   END DO
!                END DO read_loop_n_aq
               
!                IF (io_status == 0) read_successful = .TRUE.
!                close (592020)
!             END IF
!          END IF
!       END IF

!       ! Cold start if file read failed or not specified
!       IF (.NOT. read_successful) THEN
!          IF (.NOT. message_printed) THEN
!             write(*,*) 'INFO: Aquifer N initial conditions not read.', &
!                       ' Performing cold start.'
!             message_printed = .TRUE.
!          END IF
!          Aq_LDON0 = 0.0_8
!          Aq_RDON0 = 0.0_8
!          Aq_NH40 = 0.0_8
!          Aq_NO30 = 0.0_8
!       END IF

!       ! ====================================================================
!       !  Read Aquifer Phosphorus (P) Initial Conditions
!       ! ====================================================================
!       read_successful = .FALSE.
!       message_printed = .FALSE.
!       IF (len(trim(SedCNPmodel%PiniAq_file)) > 0) THEN
!          filename = trim(SedCNPmodel%PiniAq_file)
!          open(593020, file=filename, status='old', action='read', &
!               iostat=io_status)

!          IF (io_status == 0) THEN
!             read (593020, *, iostat=io_status) ! Read header
            
!             IF (io_status == 0) THEN
! read_loop_p_aq: DO i = 1, domain%ix
!                   DO j = 1, domain%jx
!                      read (593020, *, iostat=io_status) idum, idum, &
!                         Aq_LDOP0(i,j), Aq_RDOP0(i,j), Aq_PO40(i,j), &
!                         Aq_PIPA0(i,j), Aq_PIPS0(i,j) ! kg-P
!                      IF (io_status /= 0) EXIT read_loop_p_aq
!                   END DO
!                END DO read_loop_p_aq
               
!                IF (io_status == 0) read_successful = .TRUE.
!                close (593020)
!             END IF
!          END IF
!       END IF

!       ! Cold start if file read failed or not specified
!       IF (.NOT. read_successful) THEN
!          IF (.NOT. message_printed) THEN
!             write(*,*) 'INFO: Aquifer P initial conditions not read.', &
!                       ' Performing cold start.'
!             message_printed = .TRUE.
!          END IF
!          Aq_LDOP0 = 0.0_8
!          Aq_RDOP0 = 0.0_8
!          Aq_PO40 = 0.0_8
!          Aq_PIPA0 = 0.0_8
!          Aq_PIPS0 = 0.0_8
!       END IF

!    end subroutine ReadCNPini_Aq


   !> @brief Reads initial C, N, and P storages for the groundwater component.
   !! @details This subroutine reads initial conditions for all groundwater C,
   !! N, and P pools from specified files. It reads one record per groundwater
   !! basin. If a file is not specified or is invalid, it performs a "cold start"
   !! by initializing all pools to zero.
   subroutine ReadCNPini_Gw ()

      implicit none

      integer :: gwid, idum
      logical :: read_successful
      logical :: message_printed
      integer :: io_status
      character(len=256)  :: filename

      message_printed = .FALSE.

      ! ====================================================================
      !  Read Groundwater Carbon (C) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      IF (len(trim(SedCNPmodel%CiniGw_file)) > 0) THEN
         filename = trim(SedCNPmodel%CiniGw_file)
         open(591030, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (591030, *, iostat=io_status) ! Read header

            IF (io_status == 0) THEN
read_loop_c_gw: DO gwid = 1, domain%nbasin
                  read (591030, *, iostat=io_status) idum, Gw_LDOC0(gwid), &
                     Gw_RDOC0(gwid) ! kg-C
                  IF (io_status /= 0) EXIT read_loop_c_gw
               END DO read_loop_c_gw
               
               IF (io_status == 0) read_successful = .TRUE.
               close (591030)
            END IF
         END IF
      END IF

      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Groundwater C initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Gw_LDOC0 = 0.0_8
         Gw_RDOC0 = 0.0_8
      END IF

      ! ====================================================================
      !  Read Groundwater Nitrogen (N) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%NiniGw_file)) > 0) THEN
         filename = trim(SedCNPmodel%NiniGw_file)
         open(592030, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (592030, *, iostat=io_status) ! Read header

            IF (io_status == 0) THEN
read_loop_n_gw: DO gwid = 1, domain%nbasin
                  read (592030, *, iostat=io_status) idum, Gw_LDON0(gwid), &
                     Gw_RDON0(gwid), Gw_NH40(gwid), Gw_NO30(gwid) ! kg-N
                  IF (io_status /= 0) EXIT read_loop_n_gw
               END DO read_loop_n_gw
               
               IF (io_status == 0) read_successful = .TRUE.
               close (592030)
            END IF
         END IF
      END IF

      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Groundwater N initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Gw_LDON0 = 0.0_8
         Gw_RDON0 = 0.0_8
         Gw_NH40 = 0.0_8
         Gw_NO30 = 0.0_8
      END IF

      ! ====================================================================
      !  Read Groundwater Phosphorus (P) Initial Conditions
      ! ====================================================================
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%PiniGw_file)) > 0) THEN
         filename = trim(SedCNPmodel%PiniGw_file)
         open(593030, file=filename, status='old', action='read', &
              iostat=io_status)

         IF (io_status == 0) THEN
            read (593030, *, iostat=io_status) ! Read header

            IF (io_status == 0) THEN
read_loop_p_gw: DO gwid = 1, domain%nbasin
                  read (593030, *, iostat=io_status) idum, Gw_LDOP0(gwid), &
                     Gw_RDOP0(gwid), Gw_PO40(gwid)  ! kg-P
                  IF (io_status /= 0) EXIT read_loop_p_gw
               END DO read_loop_p_gw
               
               IF (io_status == 0) read_successful = .TRUE.
               close (593030)
            END IF
         END IF
      END IF

      ! Cold start if file read failed or not specified
      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Groundwater P initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Gw_LDOP0 = 0.0_8
         Gw_RDOP0 = 0.0_8
         Gw_PO40 = 0.0_8
      END IF

   end subroutine ReadCNPini_Gw


   subroutine ReadCNPini_Ch ()   !Y.Kwon

      !to read in initial CNP storages of the channels

      implicit none

      !integer :: k  !Y.Kwon
      integer              :: ich  !, chid   !BK20240512 !BK20241102
      integer              :: idum
      logical :: read_successful
      logical :: message_printed
      integer              :: io_status

      character(len=256)   :: filename   !BK20231002

      message_printed = .FALSE.

      ! --- C
      read_successful = .FALSE.
      IF (len(trim(SedCNPmodel%CiniCh_file)) > 0) THEN
         filename = trim(SedCNPmodel%CiniCh_file)
         open (591040, file=filename, status='old', action='read', iostat=io_status)

         IF (io_status == 0) THEN
            read (591040, *, iostat=io_status) ! Read header

            IF (io_status == 0) THEN
               read_loop_c_ch: DO ich = 1, domain%nch    !BK20240613
                  !chid  = SedCNP_hydro%linkID(ich)  !BK20240512 !BK20241102
                  read (591040, *, iostat=io_status) idum, &
                        Ch_ALGC0(ich,1), Ch_ALGC0(ich,2), Ch_ALGC0(ich,3), &
                        Ch_ZOOC0(ich), Ch_MBMC0(ich), Ch_DIC0(ich), &
                        Ch_LPOC0(ich), Ch_RPOC0(ich), Ch_LDOC0(ich), &
                        Ch_RDOC0(ich), Ch_POCDEPOSIT0(ich)     ! kg-C
                  !unit conversion   !BK20231002 --- commented out BK20231018
                  !Ch_ALGC0(k,1) =  Ch_ALGC0(k,1) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_ALGC0(k,2) =  Ch_ALGC0(k,2) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_ALGC0(k,3) =  Ch_ALGC0(k,3) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_ZOOC0(k) =  Ch_ZOOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_MBMC0(k) =  Ch_MBMC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_LPOC0(k) =  Ch_LPOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_RPOC0(k) =  Ch_RPOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_LDOC0(k) =  Ch_LDOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_RDOC0(k) =  Ch_RDOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C
                  !Ch_POCDEPOSIT0(k) =  Ch_POCDEPOSIT0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-C

                  ! put a small SEED values for planktons (Redfield ratio applied)
                  if (Ch_ALGC0(ich,1) <= 0.00106) Ch_ALGC0(ich,1) = 0.00106
                  if (Ch_ALGC0(ich,2) <= 0.00106) Ch_ALGC0(ich,2) = 0.00106
                  if (Ch_ALGC0(ich,3) <= 0.00106) Ch_ALGC0(ich,3) = 0.00106
                  if (Ch_ZOOC0(ich) <= 0.00106) Ch_ZOOC0(ich) = 0.00106

                  IF (io_status /= 0) EXIT read_loop_c_ch
               END DO read_loop_c_ch

               IF (io_status == 0) read_successful = .TRUE.
               close (591040)
            END IF
         END IF
      END IF

      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Channel C initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Ch_ALGC0 = 0.00001_8
         Ch_ZOOC0 = 0.00001_8
         Ch_MBMC0 = 0.0_8
         Ch_DIC0 = 0.0_8
         Ch_LPOC0 = 0.0_8
         Ch_RPOC0 = 0.0_8
         Ch_LDOC0 = 0.0_8
         Ch_RDOC0 = 0.0_8
         Ch_POCDEPOSIT0 = 0.0_8
      END IF

      ! --- N
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%NiniCh_file)) > 0) THEN
         filename = trim(SedCNPmodel%NiniCh_file)
         open (592040, file=filename, status='old', action='read', iostat=io_status)

         IF (io_status == 0) THEN
            read (592040, *, iostat=io_status) ! Read header
            
            IF (io_status == 0) THEN
               read_loop_n_ch: DO ich = 1, domain%nch   !BK20240613
                  !chid  = SedCNP_hydro%linkID(ich)  !BK20240512 !BK20241102
                  read (592040, *, iostat=io_status) idum, &
                        Ch_ALGN0(ich,1), Ch_ALGN0(ich,2), Ch_ALGN0(ich,3), &
                        Ch_ZOON0(ich), Ch_MBMN0(ich), Ch_LPON0(ich), &
                        Ch_RPON0(ich), Ch_LDON0(ich), Ch_RDON0(ich), &
                        Ch_NH40(ich), Ch_NO30(ich), Ch_PONDEPOSIT0(ich)     ! kg-N

                  !apply CNR  !BK20231115 --- commented out BK20231122
                  !Ch_ALGN0(k,1) =  Ch_ALGC0(k,1) / CNRALG
                  !Ch_ALGN0(k,2) =  Ch_ALGC0(k,2) / CNRALG
                  !Ch_ALGN0(k,3) =  Ch_ALGC0(k,3) / CNRALG
                  !Ch_ZOON0(k) =  Ch_ZOOC0(k) / CNRZOO
                  !Ch_MBMN0(k) =  Ch_MBMC0(k) / CNRMBM

                  !unit conversion   !BK20231002 --- commented out BK20231018
                  !Ch_ALGN0(k,1) =  Ch_ALGN0(k,1) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_ALGN0(k,2) =  Ch_ALGN0(k,2) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_ALGN0(k,3) =  Ch_ALGN0(k,3) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_ZOON0(k) =  Ch_ZOON0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_MBMN0(k) =  Ch_MBMN0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_LPON0(k) =  Ch_LPON0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_RPON0(k) =  Ch_RPON0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_LDON0(k) =  Ch_LDON0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_RDON0(k) =  Ch_RDON0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_NH40(k) =  Ch_NH40(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_NO30(k) =  Ch_NO30(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N
                  !Ch_PONDEPOSIT0(k) =  Ch_PONDEPOSIT0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-N

                  ! put a small SEED values for planktons (Redfield ratio applied)
                  if (Ch_ALGN0(ich,1) <= 0.00016) Ch_ALGN0(ich,1) = 0.00016
                  if (Ch_ALGN0(ich,2) <= 0.00016) Ch_ALGN0(ich,2) = 0.00016
                  if (Ch_ALGN0(ich,3) <= 0.00016) Ch_ALGN0(ich,3) = 0.00016
                  if (Ch_ZOON0(ich) <= 0.00016) Ch_ZOON0(ich) = 0.00016

                  IF (io_status /= 0) EXIT read_loop_n_ch
               END DO read_loop_n_ch

               IF (io_status == 0) read_successful = .TRUE.
               close (592040)
            END IF
         END IF
      END IF

      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Channel N initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Ch_ALGN0 = 0.00001_8
         Ch_ZOON0 = 0.00001_8
         Ch_MBMN0 = 0.0_8
         Ch_LPON0 = 0.0_8
         Ch_RPON0 = 0.0_8
         Ch_LDON0 = 0.0_8
         Ch_RDON0 = 0.0_8
         Ch_NH40 = 0.0_8
         Ch_NO30 = 0.0_8
         Ch_PONDEPOSIT0 = 0.0_8
      END IF

      ! --- p
      read_successful = .FALSE.
      message_printed = .FALSE.
      IF (len(trim(SedCNPmodel%PiniCh_file)) > 0) THEN
         filename = trim(SedCNPmodel%PiniCh_file)
         open (593040, file=filename, status='old', action='read', iostat=io_status)

         IF (io_status == 0) THEN
            read (593040, *, iostat=io_status) ! Read header

            IF (io_status == 0) THEN
               read_loop_p_ch: DO ich = 1, domain%nch  !BK20240613
                  !chid  = SedCNP_hydro%linkID(ich)  !BK20240512 !BK20241102
                  read (593040, *, iostat=io_status) idum, &
                        Ch_ALGP0(ich,1), Ch_ALGP0(ich,2), Ch_ALGP0(ich,3), &
                        Ch_ZOOP0(ich), Ch_MBMP0(ich), Ch_LPOP0(ich), &
                        Ch_RPOP0(ich), Ch_LDOP0(ich), Ch_RDOP0(ich), &
                        Ch_PO40(ich), Ch_PIPA0(ich), Ch_PIPS0(ich), &
                        Ch_POPDEPOSIT0(ich)  ! kg-P
                  
                  !apply CPR  !BK20231115 --- commented out BK20231122
                  !Ch_ALGP0(k,1) =  Ch_ALGC0(k,1) / CPRALG
                  !Ch_ALGP0(k,2) =  Ch_ALGC0(k,2) / CPRALG
                  !Ch_ALGP0(k,3) =  Ch_ALGC0(k,3) / CPRALG
                  !Ch_ZOOP0(k) =  Ch_ZOOC0(k) / CPRZOO
                  !Ch_MBMP0(k) =  Ch_MBMC0(k) / CPRMBM
                  
                  !unit conversion   !BK20231002 --- commented out BK20231018
                  !Ch_ALGP0(k,1) =  Ch_ALGC0(k,1) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_ALGP0(k,2) =  Ch_ALGC0(k,2) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_ALGP0(k,3) =  Ch_ALGC0(k,3) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_ZOOP0(k) =  Ch_ZOOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_MBMP0(k) =  Ch_MBMC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_LPOP0(k) =  Ch_LPOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_RPOP0(k) =  Ch_RPOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_LDOP0(k) =  Ch_LDOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_RDOP0(k) =  Ch_RDOC0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_PO40(k) =  Ch_PO40(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_PIPA0(k) =  Ch_PIPA0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_PIPS0(k) =  Ch_PIPS0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P
                  !Ch_POPDEPOSIT0(k) =  Ch_POCDEPOSIT0(k) * SedCNP_hydro%Wstg_st(k)   ! kg-P

                  ! put a small SEED values for planktons (Redfield ratio applied)
                  if (Ch_ALGP0(ich,1) <= 0.00001) Ch_ALGP0(ich,1) = 0.00001
                  if (Ch_ALGP0(ich,2) <= 0.00001) Ch_ALGP0(ich,2) = 0.00001
                  if (Ch_ALGP0(ich,3) <= 0.00001) Ch_ALGP0(ich,3) = 0.00001
                  if (Ch_ZOOP0(ich) <= 0.00001) Ch_ZOOP0(ich) = 0.00001

                  IF (io_status /= 0) EXIT read_loop_p_ch
               END DO read_loop_p_ch

               IF (io_status == 0) read_successful = .TRUE.
               close (593040)
            END IF
         END IF
      END IF

      IF (.NOT. read_successful) THEN
         IF (.NOT. message_printed) THEN
            write(*,*) 'INFO: Channel P initial conditions not read.',&
                      ' Performing cold start.'
            message_printed = .TRUE.
         END IF
         Ch_ALGP0 = 0.00001_8
         Ch_ZOOP0 = 0.00001_8
         Ch_MBMP0 = 0.0_8
         Ch_LPOP0 = 0.0_8
         Ch_RPOP0 = 0.0_8
         Ch_LDOP0 = 0.0_8
         Ch_RDOP0 = 0.0_8
         Ch_PO40 = 0.0_8
         Ch_PIPA0 = 0.0_8
         Ch_PIPS0 = 0.0_8
         Ch_POPDEPOSIT0 = 0.0_8
      END IF

   end subroutine ReadCNPini_Ch

end module ReadCNPini
