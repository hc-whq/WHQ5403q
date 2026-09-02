!  Program Name:
!  Author(s)/Contact(s): Y.Kwon
!  Abstract:
!  History Log:
!
!  Usage:
!  Parameters: <Specify typical arguments passed>
!  Input Files:
!        <list file names and briefly describe the data they include>
!  Output Files:
!        <list file names and briefly describe the information they include>
!
!  Condition codes:
!        <list exit condition or error codes returned >
!        If appropriate, descriptive troubleshooting instructions or
!        likely causes for failures could be mentioned here with the
!        appropriate error code
!
!  User controllable options: <if applicable>

module WriteCNPini

   !use CNPvariables  !Y.Kwon
   use module_SedCNPvariables    !BK20230316
!=====||__WHQ5403q__||=====!
!
   use SedCNP_config,         only: SedCNPmodel   !BK20231002
!
!=====||__WHQ5403q__||=====!

   contains

   subroutine WriteCNPini_So() !Y.Kwon

      !to write CNP storages of the soil at the end of the simulation

      implicit none

      integer :: i, j, isl
      real(8) :: zero = 0.0
      real(8) :: areaxyha      !CNP grid cell area [ha]    !BK20231002
      character(len=256)  :: filename   !BK20230316
      character(len=256)  :: outdir    !BK20230316

      areaxyha = domain%areaxy / 10000.0
      outdir = trim(SedCNPmodel%WHQOUT_dir)  !BK20240102

      !--- C
      !BK20230316
      filename = trim(outdir)//'/'//'CiniSo_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (691010, file=trim(filename), status='unknown')      !BK20230316
      write (691010, '(3A6, 9A14)') 'I', 'J', 'ISL', 'So_LPOCLIT', 'So_LPOCEXC', 'So_LPOCMAN', &
              'So_LPOCRES', 'So_RPOC', 'So_LDOC', 'So_RDOC', 'So_MBMC', '(kg-C ha-1)'    !BK20231002

      DO i = 1, domain%ix !m
         DO j = 1, domain%jx !n
            DO isl = 1, domain%nsl
               IF (isl .eq. 1) THEN
                  ! unit conversion
                  So_LPOCLIT(i,j) = So_LPOCLIT(i,j) / areaxyha   ! kg-C ha-1
                  So_LPOCEXC(i,j) = So_LPOCEXC(i,j) / areaxyha   ! kg-C ha-1
                  So_LPOCMAN(i,j) = So_LPOCMAN(i,j) / areaxyha   ! kg-C ha-1
                  So_LPOCRES(i,isl,j) = So_LPOCRES(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_RPOC(i,isl,j) = So_RPOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_LDOC(i,isl,j) = So_LDOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_RDOC(i,isl,j) = So_RDOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_MBMC(i,isl,j) = So_MBMC(i,isl,j) / areaxyha   ! kg-C ha-1

                  write (691010, '(3I6, 8E14.6)') i, j, isl, So_LPOCLIT(i,j), So_LPOCEXC(i,j), &  !BK20231115
                          So_LPOCMAN(i,j), So_LPOCRES(i,isl,j), So_RPOC(i,isl,j), &
                          So_LDOC(i,isl,j), So_RDOC(i,isl,j), So_MBMC(i,isl,j)    !BK20231002
               ELSE
                  ! unit conversion
                  So_LPOCRES(i,isl,j) = So_LPOCRES(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_RPOC(i,isl,j) = So_RPOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_LDOC(i,isl,j) = So_LDOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_RDOC(i,isl,j) = So_RDOC(i,isl,j) / areaxyha   ! kg-C ha-1
                  So_MBMC(i,isl,j) = So_MBMC(i,isl,j) / areaxyha   ! kg-C ha-1
                  
                  write (691010, '(3I6, 8E14.6)') i, j, isl, zero, zero, zero, So_LPOCRES(i,isl,j), & !BK20231115
                          So_RPOC(i,isl,j), So_LDOC(i,isl,j), So_RDOC(i,isl,j), So_MBMC(i,isl,j)    !BK20231002
               END IF
            END DO
         END DO
      END DO
      close (691010)

      !--- N
      !BK20230316
      filename = trim(outdir)//'/'//'NiniSo_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (692010, file=trim(filename), status='unknown')      !BK20230316
      write (692010, '(3A6, 11A14)') 'I', 'J', 'ISL', 'So_LPONLIT', 'So_LPONEXC', 'So_LPONMAN', &
              'So_LPONRES', 'So_RPON', 'So_LDON', 'So_RDON', 'So_MBMN', 'So_NH4','So_NO3', '(kg-N ha-1)'    !BK20231002
      DO i = 1, domain%ix !m
         DO j = 1, domain%jx !n
            DO isl = 1, domain%nsl
               IF (isl .eq. 1) THEN
                  ! unit conversion
                  So_LPONLIT(i,j) = So_LPONLIT(i,j) / areaxyha   ! kg-N ha-1
                  So_LPONEXC(i,j) = So_LPONEXC(i,j) / areaxyha   ! kg-N ha-1
                  So_LPONMAN(i,j) = So_LPONMAN(i,j) / areaxyha   ! kg-N ha-1
                  So_LPONRES(i,isl,j) = So_LPONRES(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_RPON(i,isl,j) = So_RPON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_LDON(i,isl,j) = So_LDON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_RDON(i,isl,j) = So_RDON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_MBMN(i,isl,j) = So_MBMN(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_NH4(i,isl,j) = So_NH4(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_NO3(i,isl,j) = So_NO3(i,isl,j) / areaxyha   ! kg-N ha-1

                  write (692010, '(3I6, 10E14.6)') i, j, isl, So_LPONLIT(i,j), So_LPONEXC(i,j), & !BK20231115
                          So_LPONMAN(i,j), So_LPONRES(i,isl,j), So_RPON(i,isl,j), &
                          So_LDON(i,isl,j), So_RDON(i,isl,j), So_MBMN(i,isl,j), &
                          So_NH4(i,isl,j), So_NO3(i,isl,j)    !BK20231002
               ELSE
                  ! unit conversion
                  So_LPONRES(i,isl,j) = So_LPONRES(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_RPON(i,isl,j) = So_RPON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_LDON(i,isl,j) = So_LDON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_RDON(i,isl,j) = So_RDON(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_MBMN(i,isl,j) = So_MBMN(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_NH4(i,isl,j) = So_NH4(i,isl,j) / areaxyha   ! kg-N ha-1
                  So_NO3(i,isl,j) = So_NO3(i,isl,j) / areaxyha   ! kg-N ha-1

                  write (692010, '(3I6, 10E14.6)') i, j, isl, zero, zero, zero, So_LPONRES(i,isl,j), &
                          So_RPON(i,isl,j), So_LDON(i,isl,j), So_RDON(i,isl,j), &
                          So_MBMN(i,isl,j), So_NH4(i,isl,j), So_NO3(i,isl,j)    !BK20231002 !BK20231115
               END IF
            END DO
         END DO
      END DO
      close (692010)

      !--- P
      !BK20230316
      filename = trim(outdir)//'/'//'PiniSo_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (693010, file=trim(filename), status='unknown')      !BK20230316
      write (693010, '(3A6, 12A14)') 'I', 'J', 'ISL', 'So_LPOPLIT', 'So_LPOPEXC', 'So_LPOPMAN', &
              'So_LPOPRES', 'So_RPOP', 'So_LDOP', 'So_RDOP', 'So_MBMP', 'So_PO4', 'So_PIPA', &
              'So_PIPS', '(kg-P ha-1)'    !BK20231002
      DO i = 1, domain%ix !m
         DO j = 1, domain%jx !n
            DO isl = 1, domain%nsl
               IF (isl .eq. 1) THEN
                  ! unit conversion
                  So_LPOPLIT(i,j) = So_LPOPLIT(i,j) / areaxyha   ! kg-P ha-1
                  So_LPOPEXC(i,j) = So_LPOPEXC(i,j) / areaxyha   ! kg-P ha-1
                  So_LPOPMAN(i,j) = So_LPOPMAN(i,j) / areaxyha   ! kg-P ha-1
                  So_LPOPRES(i,isl,j) = So_LPOPRES(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_RPOP(i,isl,j) = So_RPOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_LDOP(i,isl,j) = So_LDOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_RDOP(i,isl,j) = So_RDOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_MBMP(i,isl,j) = So_MBMP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PO4(i,isl,j) = So_PO4(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PIPA(i,isl,j) = So_PIPA(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PIPS(i,isl,j) = So_PIPS(i,isl,j) / areaxyha   ! kg-P ha-1

                  write (693010, '(3I6, 11E14.6)') i, j, isl, So_LPOPLIT(i,j), So_LPOPEXC(i,j), & !BK20231115
                          So_LPOPMAN(i,j), So_LPOPRES(i,isl,j), So_RPOP(i,isl,j), &
                          So_LDOP(i,isl,j), So_RDOP(i,isl,j), So_MBMP(i,isl,j), &
                          So_PO4(i,isl,j), So_PIPA(i,isl,j), So_PIPS(i,isl,j)    !BK20231002
               ELSE
                  ! unit conversion
                  So_LPOPRES(i,isl,j) = So_LPOPRES(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_RPOP(i,isl,j) = So_RPOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_LDOP(i,isl,j) = So_LDOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_RDOP(i,isl,j) = So_RDOP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_MBMP(i,isl,j) = So_MBMP(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PO4(i,isl,j) = So_PO4(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PIPA(i,isl,j) = So_PIPA(i,isl,j) / areaxyha   ! kg-P ha-1
                  So_PIPS(i,isl,j) = So_PIPS(i,isl,j) / areaxyha   ! kg-P ha-1

                  write (693010, '(3I6, 11E14.6)') i, j, isl, zero, zero, zero, So_LPOPRES(i,isl,j), & !BK20231115
                          So_RPOP(i,isl,j), So_LDOP(i,isl,j), So_RDOP(i,isl,j), &
                          So_MBMP(i,isl,j), So_PO4(i,isl,j), So_PIPA(i,isl,j), So_PIPS(i,isl,j)    !BK20231002
               END IF
            END DO
         END DO
      END DO
      close (693010)

   end subroutine WriteCNPini_So


   ! subroutine WriteCNPini_Aq()  !Y.Kwon

   !    !to write CNP storages of the aquifer at the end of the simulation

   !    use config_base,           only: SedCNPmodel

   !    implicit none

   !    integer :: i, j
   !    character(len=256)  :: filename   !BK20230316
   !    character(len=256)  :: outdir    !BK20230316
   !    real(8)             :: Aq_WStorage !water storage in aquifer [m3]  !Y.Kwon

   !    outdir = trim(SedCNPmodel%WHQOUT_dir)  !BK20240102

   !    !--- C
   !    !BK20230316
   !    filename = trim(outdir)//'/'//'CiniAq_TEST_'//trim(dateSedCNP%olddate(1:4))//&
   !              trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
   !              trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
   !              '.dat'
   !    open (691020, file=trim(filename), status='unknown')      !BK20230316
   !    write (691020, '(2A6, 3A14)') 'I', 'J', 'Aq_LDOC', 'Aq_RDOC', '(kg-C)'   !BK20231002 !BK20231018
   !    DO i = 1, domain%ix !m
   !       DO j = 1, domain%jx !n
   !          ! unit conversion --- commented out BK20231018
   !          !Aq_WStorage= (SedCNP_hydro%WA(i,j)/1000.0) * domain%areaxy    ! water storage in aquifer [m3] !Y.Kwon  !BK20231002
   !          !
   !          !if (Aq_WStorage > 0.0) then  !Y.Kwon
   !          !   Aq_LDOC(i,j) = Aq_LDOC(i,j) / Aq_WStorage   ! kg-C m-3  !BK20230808
   !          !   Aq_RDOC(i,j) = Aq_RDOC(i,j) / Aq_WStorage   ! kg-C m-3  !BK20230808
   !          !else
   !          !   Aq_LDOC(i,j) = 0.0
   !          !   Aq_RDOC(i,j) = 0.0
   !          !endif
   !          write (691020, '(2I6, 2E14.6)') i, j, Aq_LDOC(i,j), Aq_RDOC(i,j)    !BK20231002  !BK20231115
   !       END DO
   !    END DO
   !    close (691020)

   !    !--- N
   !    !BK20230316
   !    filename = trim(outdir)//'/'//'NiniAq_TEST_'//trim(dateSedCNP%olddate(1:4))//&
   !              trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
   !              trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
   !              '.dat'
   !    open (692020, file=trim(filename), status='unknown')      !BK20230316
   !    write (692020, '(2A6, 5A14)') 'I', 'J', 'Aq_LDON', 'Aq_RDON', 'Aq_NH4', 'Aq_NO3', '(kg-N)'    !BK20231002 !BK20231018
   !    DO i = 1, domain%ix !m
   !       DO j = 1, domain%jx !n

   !          !Aq_WStorage= (SedCNP_hydro%WA(i,j)/1000.0) * domain%areaxy    ! water storage in aquifer [m3] !Y.Kwon  !BK20231002

   !          ! unit conversion --- commented out BK20231018
   !          !if (Aq_WStorage > 0.0) then  !Y.Kwon
   !          !   Aq_LDON(i,j) = Aq_LDON(i,j) / Aq_WStorage   ! kg-N m-3  !BK20230808
   !          !   Aq_RDON(i,j) = Aq_RDON(i,j) / Aq_WStorage   ! kg-N m-3  !BK20230808
   !          !   Aq_NH4(i,j) = Aq_NH4(i,j) / Aq_WStorage   ! kg-N m-3  !BK20230808
   !          !   Aq_NO3(i,j) = Aq_NO3(i,j) / Aq_WStorage   ! kg-N m-3  !BK20230808
   !          !else
   !          !   Aq_LDON(i,j) = 0.0
   !          !   Aq_RDON(i,j) = 0.0
   !          !   Aq_NH4(i,j)  = 0.0
   !          !   Aq_NO3(i,j)  = 0.0
   !          !endif
   !          write (692020, '(2I6, 4E14.6)') i, j, Aq_LDON(i,j), Aq_RDON(i,j), Aq_NH4(i,j), Aq_NO3(i,j)    !BK20231002 !BK20231115
   !       END DO
   !    END DO
   !    close (692020)

   !    !--- P
   !    !BK20230316
   !    filename = trim(outdir)//'/'//'PiniAq_TEST_'//trim(dateSedCNP%olddate(1:4))//&
   !              trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
   !              trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
   !              '.dat'
   !    open (693020, file=trim(filename), status='unknown')      !BK20230316
   !    write (693020, '(2A6, 6A14)') 'I', 'J', 'Aq_LDOP', 'Aq_RDOP', 'Aq_PO4', 'Aq_PIPA', &
   !            'Aq_PIPS', '(kg-P)'    !BK20231002 !BK20231018
   !    DO i = 1, domain%ix !m
   !       DO j = 1, domain%jx !n

   !          !Aq_WStorage= (SedCNP_hydro%WA(i,j)/1000.0) * domain%areaxy    ! water storage in aquifer [m3] !Y.Kwon  !BK20231002

   !          ! unit conversion --- commented out BK20231018
   !          !if (Aq_WStorage > 0.0) then  !Y.Kwon
   !          !   Aq_LDOP(i,j) = Aq_LDOP(i,j) / Aq_WStorage   ! kg-P m-3  !BK20230808
   !          !   Aq_RDOP(i,j) = Aq_RDOP(i,j) / Aq_WStorage   ! kg-P m-3  !BK20230808
   !          !   Aq_PO4(i,j) = Aq_PO4(i,j) / Aq_WStorage   ! kg-P m-3  !BK20230808
   !          !   Aq_PIPA(i,j) = Aq_PIPA(i,j) / Aq_WStorage   ! kg-P m-3  !BK20230808
   !          !   Aq_PIPS(i,j) = Aq_PIPS(i,j) / Aq_WStorage   ! kg-P m-3  !BK20230808
   !          !else
   !          !   Aq_LDOP(i,j) = 0.0
   !          !   Aq_RDOP(i,j) = 0.0
   !          !   Aq_PO4(i,j) = 0.0
   !          !   Aq_PIPA(i,j) = 0.0
   !          !   Aq_PIPS(i,j) = 0.0
   !          !endif
   !          write (693020, '(2I6, 5E14.6)') i, j, Aq_LDOP(i,j), Aq_RDOP(i,j), Aq_PO4(i,j), & !BK20231115
   !                  Aq_PIPA(i,j), Aq_PIPS(i,j)    !BK20231002
   !       END DO
   !    END DO
   !    close (693020)


   ! end subroutine WriteCNPini_Aq


   subroutine WriteCNPini_Gw()   !Y.Kwon

      !to write CNP storages of the subwatershed groundwater at the end of the simulation

      implicit none

      !integer :: ibasin, sbid           !BK20240513 
      integer :: gwid                   !BK20240613 
      character(len=256)  :: filename   !BK20230316
      character(len=256)  :: outdir     !BK20230316

      outdir = trim(SedCNPmodel%WHQOUT_dir)  !BK20240102

      !--- C
      !BK20230316
      filename = trim(outdir)//'/'//'CiniGw_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (691030, file=trim(filename), status='unknown')      !BK20230316
      write (691030, '(A6, 3A14)') 'GWID', 'Gw_LDOC', 'Gw_RDOC', '(kg-C)'    !BK20231002 !BK20231024
      DO gwid = 1, domain%nbasin   !BK20240613 
         write (691030, '(I6, 2E14.6)') gwid, Gw_LDOC(gwid), Gw_RDOC(gwid)    !BK20240613
      END DO
      close (691030)

      !--- N
      !BK20230316
      filename = trim(outdir)//'/'//'NiniGw_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (692030, file=trim(filename), status='unknown')      !BK20230316
      write (692030, '(A6, 5A14)') 'GWID', 'Gw_LDON', 'Gw_RDON','Gw_NH4','Gw_NO3', '(kg-N)'    !BK20231002 !BK20231024
      DO gwid = 1, domain%nbasin   !BK20240613
         write (692030, '(I6, 4E14.6)') gwid, Gw_LDON(gwid), Gw_RDON(gwid), Gw_NH4(gwid), Gw_NO3(gwid)    !BK20240613
      END DO
      close (692030)

      !--- P
      !BK20230316
      filename = trim(outdir)//'/'//'PiniGw_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (693030, file=trim(filename), status='unknown')      !BK20230316
      write (693030, '(A6, 6A14)') 'GWID', 'Gw_LDOP', 'Gw_RDOP','Gw_PO4', '(kg-P)'    !BK20231002 !BK20231024 !BK20250809
      DO gwid = 1, domain%nbasin   !BK20240613
         write (693030, '(I6, 5E14.6)') gwid, Gw_LDOC(gwid), Gw_RDOC(gwid), Gw_PO4(gwid)  !BK20240613 !BK20250809
      END DO
      close (693030)

   end subroutine WriteCNPini_Gw


   subroutine WriteCNPini_Ch()   !Y.Kwon

      !to write CNP storages of the channels at the end of the simulation

      implicit none

      integer              :: ich, ia  !, chid   !BK20240513  !BK20241102
      character(len=256)   :: filename   !BK20230316
      character(len=256)   :: outdir     !BK20230316

      outdir = trim(SedCNPmodel%WHQOUT_dir)  !BK20240102

      !--- C
      !BK20230316
      filename = trim(outdir)//'/'//'CiniCh_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (691040, file=trim(filename), status='unknown')      !BK20230316
      write (691040, '(A6, 12A14)') 'Ch_Index', 'Ch_ALGC(1)', 'Ch_ALGC(2)', 'Ch_ALGC(3)', &
              'Ch_ZOOC', 'Ch_MBMC', 'Ch_DIC', 'Ch_LPOC', 'Ch_RPOC', 'Ch_LDOC', 'Ch_RDOC', &
              'Ch_POCDEPOSIT', '(kg-C)'    !BK20231002 !BK20241102 !BK20251120
      DO ich = 1, domain%nch  !BK20240613
         write (691040, '(I6, 11E14.6)') ich, Ch_ALGC(ich,1), Ch_ALGC(ich,2), Ch_ALGC(ich,3), Ch_ZOOC(ich), & 
                 Ch_MBMC(ich), Ch_DIC(ich), Ch_LPOC(ich), Ch_RPOC(ich), Ch_LDOC(ich), Ch_RDOC(ich), & 
                 Ch_POCDEPOSIT(ich)   !BK20231002 !BK20240513 !BK20241102 !BK20251120
      END DO
      close (691040)

      !--- N
      !BK20230316
      filename = trim(outdir)//'/'//'NiniCh_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (692040, file=trim(filename), status='unknown')      !BK20230316
      write (692040, '(A6, 13A14)') 'Ch_Index', 'Ch_ALGN(1)', 'Ch_ALGN(2)', 'Ch_ALGN(3)', &
              'Ch_ZOON', 'Ch_MBMN', 'Ch_LPON', 'Ch_RPON', 'Ch_LDON', 'Ch_RDON', &
              'Ch_NH4', 'Ch_NO3','Ch_PONDEPOSIT', '(kg-N)'    !BK20231002
      DO ich = 1, domain%nch  !BK20240613
         write (692040, '(I6, 12E14.6)') ich, Ch_ALGN(ich,1), Ch_ALGN(ich,2), Ch_ALGN(ich,3), Ch_ZOON(ich), & 
                 Ch_MBMN(ich), Ch_LPON(ich), Ch_RPON(ich), Ch_LDON(ich), Ch_RDON(ich), Ch_NH4(ich), & 
                 Ch_NO3(ich), Ch_PONDEPOSIT(ich)      !BK20231002 !BK20240513 !BK20241102
      END DO
      close (692040)

      !--- P
      !BK20230316
      filename = trim(outdir)//'/'//'PiniCh_TEST_'//trim(dateSedCNP%olddate(1:4))//&
                trim(dateSedCNP%olddate(6:7))//trim(dateSedCNP%olddate(9:10))//&
                trim(dateSedCNP%olddate(12:13))//trim(dateSedCNP%olddate(15:16))//&
                '.dat'
      open (693040, file=trim(filename), status='unknown')      !BK20230316
      write (693040, '(A6, 14A14)') 'Ch_Index', 'Ch_ALGP(1)', 'Ch_ALGP(2)', 'Ch_ALGP(3)', &
              'Ch_ZOOP', 'Ch_MBMP', 'Ch_LPOP', 'Ch_RPOP', 'Ch_LDOP', 'Ch_RDOP', &
              'Ch_PO4', 'Ch_PIPA', 'Ch_PIPS', 'Ch_POPDEPOSIT', '(kg-P)'    !BK20231002
      DO ich = 1, domain%nch   !BK20240613
         write (693040, '(I6, 13E14.6)') ich, Ch_ALGP(ich,1), Ch_ALGP(ich,2), Ch_ALGP(ich,3), Ch_ZOOP(ich), & 
                 Ch_MBMP(ich), Ch_LPOP(ich), Ch_RPOP(ich), Ch_LDOP(ich), Ch_RDOP(ich), Ch_PO4(ich), & 
                 Ch_PIPA(ich), Ch_PIPS(ich), Ch_POPDEPOSIT(ich)    !BK20231002 !BK20240513 !BK20241102
      END DO
      close (693040)

   end subroutine WriteCNPini_Ch

end module WriteCNPini

