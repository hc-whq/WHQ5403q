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

module CNPparams

   contains

   subroutine ReadCNPparams()

      !to read in CNP parameters (CNPparams.dat)

      !use CNPvariables
      use module_SedCNPvariables !Y.Kwon
!=====||__WHQ5403q__||=====!
!
      use SedCNP_config, only: SedCNPmodel !BK20251130
!
!=====||__WHQ5403q__||=====!
      use module_hydro_stop, only:HYDRO_stop !BK20251130

      implicit none
      integer :: ilc, isl, imp, ia
      
      !BK20251130 
      logical :: file_exists
      integer :: ierr

      if (len_trim(SedCNPmodel%CNPparams_file) > 0) then
         inquire(file=trim(SedCNPmodel%CNPparams_file), exist=file_exists)
         if (file_exists) then
            open (5901, file=trim(SedCNPmodel%CNPparams_file), status='unknown', iostat=ierr)
         else
             ierr = 1 ! Flag to try defaults
         endif
      else
         ierr = 1 ! Flag to try defaults
      endif

      if (ierr /= 0) then
          !BK20251130 Check for default files
          inquire(file='./WHQIN/CNPparams.dat', exist=file_exists)
          if (file_exists) then
             open (5901, file='./WHQIN/CNPparams.dat', status='unknown', iostat=ierr)
          else
             open (5901, file='./WHQIN/CNPparams_TEST.dat', status='unknown', iostat=ierr)
          endif
      endif
      
      if (ierr /= 0) then
         call HYDRO_stop("Failed to open CNPparams data file")
      endif

      !open (5901, file='./WHQIN/CNPparams_TEST.dat', status='unknown')  !BK20231104
      read (5901, *)
      read (5901, *) FBIO, FMET, FDEC, FREF      ! FBIO + FMET + FDEC = 1,  0 < FREF < 1
      !read (5901, *)
      !read (5901, *) FMBMC, FMBMN, FMBMP  !BK20240114
      !read (5901, *)
      !read (5901, *) FALGC, FALGN, FALGP
      !read (5901, *)
      !read (5901, *) FZOOC, FZOON, FZOOP
        !calculate CNR and CPR for organisms (-):
        !CNRMBM = FMBMC / FMBMN
        !CNRALG = FALGC / FALGN
        !CNRZOO = FZOOC / FZOON
        !CPRMBM = FMBMC / FMBMP
        !CPRALG = FALGC / FALGP
        !CPRZOO = FZOOC / FZOOP

      DO ilc = 1, nlc
         read (5901, *)
         read (5901, *)
         read (5901, *) (FROOT(ilc, isl), isl = 1, domain%nsl)  ! (-)
         read (5901, *)
         read (5901, *) KLIT(ilc), KRES(ilc), KEXC(ilc), KMAN(ilc), KLDOM(ilc), KRPOM(ilc), &  !BK typo fixed 20230308
           KRDOM(ilc), KMBM(ilc), KDEN(ilc), KNIT(ilc), KVOL(ilc)    ! read in (day-1), converted to (sec-1)     !BK typo fixed 20230308
           ! unit conversion
           KLIT(ilc) = KLIT(ilc) / 86400.0 ! (sec-1)
           KRES(ilc) = KRES(ilc) / 86400.0 ! (sec-1)
           KEXC(ilc) = KEXC(ilc) / 86400.0 ! (sec-1)
           KMAN(ilc) = KMAN(ilc) / 86400.0 ! (sec-1)
           KLDOM(ilc) = KLDOM(ilc) / 86400.0 ! (sec-1)
           KRPOM(ilc) = KRPOM(ilc) / 86400.0 ! (sec-1)
           KRDOM(ilc) = KRDOM(ilc) / 86400.0 ! (sec-1)
           KMBM(ilc) = KMBM(ilc) / 86400.0 ! (sec-1)
           KDEN(ilc) = KDEN(ilc) / 86400.0 ! (sec-1)
           KNIT(ilc) = KNIT(ilc) / 86400.0 ! (sec-1) 
           KVOL(ilc) = KVOL(ilc) / 86400.0 ! (sec-1)
      END DO
      read (5901, *)
      read (5901, *)
      read (5901, *) KLPOMW, KLDOMW, KRPOMW, KRDOMW, KMBMW, KDENW, KNITW, KVOLW   ! read in (day-1), converted to (sec-1)
        ! unit conversion
        KLPOMW = KLPOMW / 86400.0 ! (sec-1)
        KLDOMW = KLDOMW / 86400.0 ! (sec-1)
        KRPOMW = KRPOMW / 86400.0 ! (sec-1)
        KRDOMW = KRDOMW / 86400.0 ! (sec-1)
        KMBMW = KMBMW / 86400.0 ! (sec-1)
        KDENW = KDENW / 86400.0 ! (sec-1)
        KNITW = KNITW / 86400.0 ! (sec-1)
        KVOLW = KVOLW / 86400.0 ! (sec-1)
      read (5901, *)
      !read (5901, *)
      !read (5901, *) KADSPO4, KADSPIPA   ! read in (day-1), converted to (sec-1)
      !  ! unit conversion
      !  KADSPO4 = KADSPO4 / 86400.0 ! (sec-1)
      !  KADSPIPA = KADSPIPA / 86400.0 ! (sec-1)
      KADSPO4 = 1.0 / 86400.0     ! (sec-1)             !BK20240623 hardcoded (as in SWAT)
      KADSPIPA = 0.0006 / 86400.0 ! (sec-1)             !BK20240623 hardcoded (as in SWAT)
      !read (5901, *)
      !read (5901, *) SIGMAPO4, SIGMAPIPA  ! (-)
      SIGMAPO4 = 0.667  !(-)             !BK20240623 hardcoded (as in SWAT)
      SIGMAPIPA = 0.25  !(-)             !BK20240623 hardcoded (as in SWAT)
      read (5901, *)
      read (5901, *) ALPHADOC, ALPHACHL, ALPHATSS  ! (-)
      read (5901, *)
      read (5901, *) ALBEDO  ! (%)
      read (5901, *)
      read (5901, *) SOLRADMAX  ! (w m-2)
      read (5901, *)
      read (5901, *) MONOD_N, MONOD_P  ! read in (g m-3), converted to (kg m-3)
         ! unit conversion
         MONOD_N = MONOD_N / 1000.0  ! (kg m-3)
         MONOD_P = MONOD_P / 1000.0  ! (kg m-3)
      read (5901, *)
      DO ia = 1, nalg
         read (5901, *) KAGMAX(ia), KARMAX(ia), KAEMAX(ia), KAMMAX(ia)  ! read in (day-1), converted to (sec-1) 
            ! unit conversion
            KAGMAX(ia) = KAGMAX(ia) / 86400.0 ! (sec-1)
            KARMAX(ia) = KARMAX(ia) / 86400.0 ! (sec-1)
            KAEMAX(ia) = KAEMAX(ia) / 86400.0 ! (sec-1)
            KAMMAX(ia) = KAMMAX(ia) / 86400.0 ! (sec-1)
      END DO
      read (5901, *)
      read (5901, *) EZI ! 0.0 - 1.0 (-)
      read (5901, *)
      read (5901, *) ZLOW, ZHALF  ! read in (g m-3), converted to (kg m-3)
         ! unit conversion
         ZLOW = ZLOW / 1000.0  ! (kg m-3)
         ZHALF = ZHALF / 1000.0  ! (kg m-3)
      read (5901, *)
      read (5901, *) KZIMAX, KZRMAX, KZMMAX  ! read in (day-1), converted to (sec-1) 
         ! unit conversion
         KZIMAX = KZIMAX / 86400.0 ! (sec-1)
         KZRMAX = KZRMAX / 86400.0 ! (sec-1)
         KZMMAX = KZMMAX / 86400.0 ! (sec-1)
      DO ia = 1, nalg
         read (5901, *)
         read (5901, *) (TALG(ia,imp), imp = 1, ntm)   ! (deg_C)
         read (5901, *) (MALG(ia,imp), imp = 1, ntm)   ! (-)
      END DO
      read (5901, *)
      read (5901, *) (TMBM(imp), imp = 1, ntm)
      read (5901, *) (MMBM(imp), imp = 1, ntm)
      read (5901, *)
      read (5901, *) (TZOO(imp), imp = 1, ntm)
      read (5901, *) (MZOO(imp), imp = 1, ntm)
      read (5901, *)
      read (5901, *) OMEGAPOM, OMEGAALG, OMEGAZOO  ! read in (m day-1), converted to (m sec-1)
         ! unit conversion
         OMEGAPOM = OMEGAPOM / 86400.0 ! (m sec-1)
         OMEGAALG = OMEGAALG / 86400.0 ! (m sec-1)
         OMEGAZOO = OMEGAZOO / 86400.0 ! (m sec-1)
      read (5901, *)
      read (5901, *) TAQUIFER  ! (deg-C)
      read (5901, *)

      !BK20260302
      read (5901, *) etaPOMsurf, etaDOMsurf, etaNH4surf, etaNO3surf, etaPO4surf
      read (5901, *)
      read (5901, *)             etaDOMintf, etaNH4intf, etaNO3intf, etaPO4intf
      read (5901, *)
      read (5901, *)             etaDOMperc, etaNH4perc, etaNO3perc, etaPO4perc
      read (5901, *)
      read (5901, *)             etaDOMsogw, etaNH4sogw, etaNO3sogw, etaPO4sogw
      read (5901, *)
      read (5901, *)             etaDOMgwch, etaNH4gwch, etaNO3gwch, etaPO4gwch
      read (5901, *)
      read (5901, *) etaPOMdsch, etaDOMdsch, etaNH4dsch, etaNO3dsch, etaPO4dsch

      close (5901)

   end subroutine ReadCNPparams         

end module CNPparams
