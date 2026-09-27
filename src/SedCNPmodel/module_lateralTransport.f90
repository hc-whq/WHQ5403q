!=====||__WHQ5403__||=====!
!
!> @brief Lateral loads (surface runoff, interflow) from the land cells to the channel links.
!! @details
!! SedCNP_lateral_option = 1 (default): the loads follow the cell water paths of the hydro
!! model, using the cell water fluxes of the time-step written to RTOUT (module_UDMAP):
!!   1. interflow: the interflow load of a cell (So_*intf) moves with q_intf to the receiving
!!      neighbour (sub_dir) and enters its soil; the fraction sub_exfil/(interflow received)
!!      exfiltrates with the water to the surface of the receiving cell (as in SUBSFC_RTNG);
!!      loads moving to a domain-edge cell leave the domain.
!!   2. surface: the surface loads of a cell (So_*surf, overland sediment Ssurf) and the
!!      exfiltrated interflow loads are mixed with the surface water of the cell (sfc_in plus the
!!      inflow from upslope cells) and leave with it: to the neighbour (sfc_out, sfc_dir), to the
!!      channel (sfc_chan; link of the stream pixel, internal lake links -> lake outlet), to the
!!      lake (sfc_lake; lake outlet link), across the domain boundary (sfc_bdry), or remain on the
!!      cell (sfc_rem), which returns to the LSM in the hydro model: the remaining loads are
!!      returned to the top soil layer (sediment to the overland sediment storage Sol).
!!      The cells are processed from upslope to downslope; the rare closed loops of the main
!!      flow directions are solved exactly.
!!   Groundwater loads are distributed by gw_to_ch (hydro-model weights, Lateral_Map_Init).
!! SedCNP_lateral_option = 0: previous method; the (previous time-step) loads of all cells of the
!! gw basin of a link are summed (lat_sum*, module_SedCNPvariables).
!! Both options fill latS_ch / latI_ch [kg dt-1], used by the channel sediment/CNP routines.
!
module module_lateralTransport

   use module_SedCNPvariables
   use SedCNP_config,     only: SedCNPmodel
   use module_SedCNP_in,  only: get3d_lsm_real, get2d_int, is_north_up
   use module_hydro_stop, only: hydro_stop

   implicit none

   private
   public :: Lateral_Transport_Init, Lateral_Loads

   integer, parameter   :: ntr = ntr_s + ntr_i      ! transported tracers (surface + exfiltrated interflow)
   integer, allocatable :: cell_link(:,:)           ! channel index receiving sfc_chan of the cell (0: none)
   integer, allocatable :: cell_lake(:,:)           ! lake outlet channel index receiving sfc_lake (0: none)
   integer              :: ubud = 5120              ! unit of debug/lateral_budget.csv

   contains

   !> @brief Allocates the link loads and, for option 1, maps the stream pixels and lake cells to
   !! channel links. Called once, after Lateral_Map_Init.
   subroutine Lateral_Transport_Init()

      use hashtable,       only: hash_t
      use ISO_FORTRAN_ENV, only: int64

      implicit none

      type(hash_t)         :: hash_table
      integer(kind=int64)  :: val
      logical              :: found
      integer, allocatable :: lakegrid(:,:)
      integer              :: i, j, il, status

      allocate(latS_ch(ntr_s, domain%nch), latI_ch(ntr_i, domain%nch))
      latS_ch = 0.0d0
      latI_ch = 0.0d0
      if (SedCNPmodel%SedCNP_lateral_option /= 1) return

      if (nps /= 4) call hydro_stop("Lateral_Transport_Init: 4 sediment grain sizes are assumed")
      if (domain%ix /= domain%ixrt .or. domain%jx /= domain%jxrt) &
         call hydro_stop("Lateral_Transport_Init: SedCNP_lateral_option = 1 requires AGGFACTRT = 1")

      ! --- stream pixels -> channel index (internal lake links -> lake outlet)
      allocate(cell_link(domain%ix, domain%jx), cell_lake(domain%ix, domain%jx))
      cell_link = 0
      cell_lake = 0
      call hash_table%set_all_idx(int(SedCNP_hydro%linkID, int64), domain%nch)
      do j = 1, domain%jx
         do i = 1, domain%ix
            if (SedCNP_hydro%linkID_grid(i,j) >= 1) then
               call hash_table%get(int(SedCNP_hydro%linkID_grid(i,j), int64), val, found)
               if (found) cell_link(i,j) = lake_ich(int(val))
            endif
         enddo
      enddo
      call hash_table%clear()

      ! --- lake cells (Fulldom LAKEGRID = LAKEPARM lake_id) -> lake outlet
      if (lake%active) then
         allocate(lakegrid(domain%ix, domain%jx))
         lakegrid = -9999
         status = get2d_int("LAKEGRID", lakegrid, domain%ix, domain%jx, './DOMAIN/Fulldom_hires.nc')
         if (status == 0) then
            if (is_north_up('./DOMAIN/Fulldom_hires.nc')) lakegrid = lakegrid(:, domain%jx:1:-1)
            do j = 1, domain%jx
               do i = 1, domain%ix
                  if (lakegrid(i,j) <= 0) cycle
                  do il = 1, lake%nlake
                     if (lake%lake_id(il) == lakegrid(i,j)) cell_lake(i,j) = lake%outlet_ich(il)
                  enddo
               enddo
            enddo
         else
            call progress_clear()   !WHQ5403
            write(6,*) 'WARNING: Lateral_Transport_Init: LAKEGRID not found in Fulldom_hires.nc;', &
                       ' surface loads to lakes are counted as losses'
         endif
         deallocate(lakegrid)
      endif

      call system('mkdir -p ./debug')
      open(ubud, file='./debug/lateral_budget.csv', status='replace')
      write(ubud,'(A)') 'itime,group,source_surface,interflow_moved,interflow_to_soil,interflow_exfiltrated,'// &
                        'interflow_boundary,to_channel,to_lake,boundary,unlinked,remaining_to_soil,residual'
      call progress_clear()   !WHQ5403
      write(6,'(A,I8,A,I8)') ' INFO: SedCNP lateral transport: stream pixels mapped ', count(cell_link > 0), &
                             '   lake cells mapped ', count(cell_lake > 0)

   end subroutine Lateral_Transport_Init


   !> @brief Lateral loads of the time-step to the channel links (latS_ch, latI_ch).
   !! Called every time-step after the soil/overland and groundwater calculations and before the
   !! channel calculations.
   subroutine Lateral_Loads(itime)

      implicit none

      integer, intent(in) :: itime
      integer             :: ich, ips
      real(8)             :: dt

      dt = real(SedCNPmodel%SedCNP_timestep, 8)
      latS_ch = 0.0d0
      latI_ch = 0.0d0

      if (SedCNPmodel%SedCNP_lateral_option == 1) then
         call Lateral_Transport(itime)
         return
      endif

      ! --- option 0: sum over the gw basin(s) of each link (loads of the previous time-step)
      do ich = 1, domain%nch
         if (lake%active) then
            if (lake%typel(ich) == 2) cycle
         endif
         do ips = 1, nps
            latS_ch(ips,ich) = lat_sum_sed(overSed%Ssurf0, ips, ich) * dt
         enddo
         latS_ch(5,ich)  = lat_sum2(So_LPOCsurf0, ich)
         latS_ch(6,ich)  = lat_sum2(So_RPOCsurf0, ich)
         latS_ch(7,ich)  = lat_sum2(So_LDOCsurf0, ich)
         latS_ch(8,ich)  = lat_sum2(So_RDOCsurf0, ich)
         latS_ch(9,ich)  = lat_sum2(So_MBMCsurf0, ich)
         latS_ch(10,ich) = lat_sum2(So_LPONsurf0, ich)
         latS_ch(11,ich) = lat_sum2(So_RPONsurf0, ich)
         latS_ch(12,ich) = lat_sum2(So_LDONsurf0, ich)
         latS_ch(13,ich) = lat_sum2(So_RDONsurf0, ich)
         latS_ch(14,ich) = lat_sum2(So_MBMNsurf0, ich)
         latS_ch(15,ich) = lat_sum2(So_NH4surf0, ich)
         latS_ch(16,ich) = lat_sum2(So_NO3surf0, ich)
         latS_ch(17,ich) = lat_sum2(So_LPOPsurf0, ich)
         latS_ch(18,ich) = lat_sum2(So_RPOPsurf0, ich)
         latS_ch(19,ich) = lat_sum2(So_LDOPsurf0, ich)
         latS_ch(20,ich) = lat_sum2(So_RDOPsurf0, ich)
         latS_ch(21,ich) = lat_sum2(So_MBMPsurf0, ich)
         latS_ch(22,ich) = lat_sum2(So_PO4surf0, ich)
         latS_ch(23,ich) = lat_sum2(So_PIPAsurf0, ich)
         latS_ch(24,ich) = lat_sum2(So_PIPSsurf0, ich)
         latI_ch(1,ich)  = lat_sum3l(So_LDOCintf0, ich)
         latI_ch(2,ich)  = lat_sum3l(So_RDOCintf0, ich)
         latI_ch(3,ich)  = lat_sum3l(So_LDONintf0, ich)
         latI_ch(4,ich)  = lat_sum3l(So_RDONintf0, ich)
         latI_ch(5,ich)  = lat_sum3l(So_NH4intf0, ich)
         latI_ch(6,ich)  = lat_sum3l(So_NO3intf0, ich)
         latI_ch(7,ich)  = lat_sum3l(So_LDOPintf0, ich)
         latI_ch(8,ich)  = lat_sum3l(So_RDOPintf0, ich)
         latI_ch(9,ich)  = lat_sum3l(So_PO4intf0, ich)
      enddo

   end subroutine Lateral_Loads


   !> @brief Option 1: interflow and surface load transport along the hydro-model cell paths.
   subroutine Lateral_Transport(itime)

      implicit none

      integer, intent(in)  :: itime
      integer              :: ix, jx, ncell
      real(8), allocatable :: sfc_in(:,:), sfc_out(:,:), sfc_dir(:,:), sfc_chan(:,:), sfc_lake(:,:)
      real(8), allocatable :: sfc_bdry(:,:), sfc_rem(:,:), sub_dir(:,:), sub_exfil(:,:), q_intf(:,:)
      real(8), allocatable :: src(:,:)       ! (ntr,ncell) surface loads of the cells [kg dt-1]
      real(8), allocatable :: inflow(:,:)    ! (ntr,ncell) loads from upslope cells
      real(8), allocatable :: remain(:,:)    ! (ntr,ncell) loads remaining on the cells
      real(8), allocatable :: wout(:)        ! total surface water leaving/remaining in the cell [mm]
      real(8), allocatable :: gross(:,:), fex(:,:)
      integer, allocatable :: nxt(:), indeg(:), order(:), walk(:), cyc(:)
      logical, allocatable :: done(:)
      real(8)              :: bud(ntr,11)   ! budget terms, see the csv header
      real(8)              :: tot(ntr), bvec(ntr), y(ntr), dt, a, amul
      character(len=256)   :: filename_RT
      integer              :: i, j, k, m, n, ri, rj, nord, head, ncyc, kk, status, ips
      logical              :: degenerate

      ix = domain%ix; jx = domain%jx; ncell = ix*jx
      dt = real(SedCNPmodel%SedCNP_timestep, 8)
      bud = 0.0d0

      ! --- cell water fluxes of the time-step (RTOUT, model grid)
      filename_RT = trim(SedCNPmodel%LDASOUT_dir)//'/'//&
                    trim(dateSedCNP%olddate(1:4))//trim(dateSedCNP%olddate(6:7))//&
                    trim(dateSedCNP%olddate(9:10))//trim(dateSedCNP%olddate(12:13))//&
                    trim(dateSedCNP%olddate(15:16))//'.RTOUT_DOMAIN1'
      allocate(sfc_in(ix,jx), sfc_out(ix,jx), sfc_dir(ix,jx), sfc_chan(ix,jx), sfc_lake(ix,jx), &
               sfc_bdry(ix,jx), sfc_rem(ix,jx), sub_dir(ix,jx), sub_exfil(ix,jx), q_intf(ix,jx))
      status = 0
      status = status + abs(get3d_lsm_real("sfc_in",    sfc_in,    ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_out",   sfc_out,   ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_dir",   sfc_dir,   ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_chan",  sfc_chan,  ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_lake",  sfc_lake,  ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_bdry",  sfc_bdry,  ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sfc_rem",   sfc_rem,   ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sub_dir",   sub_dir,   ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("sub_exfil", sub_exfil, ix, jx, trim(filename_RT)))
      status = status + abs(get3d_lsm_real("q_intf",    q_intf,    ix, jx, trim(filename_RT)))
      if (status /= 0) then
         call progress_clear()   !WHQ5403
         write(6,*) 'ERROR: Lateral_Transport: cell water fluxes (sfc_*, sub_*) not found in ', trim(filename_RT)
         call hydro_stop("SedCNP_lateral_option = 1 requires RTOUT from the hydro model with the cell water "// &
                         "fluxes (rerun the hydro model), or set SedCNP_lateral_option = 0")
      endif
      where (sfc_in    < 0.0d0) sfc_in    = 0.0d0   ! fill values
      where (sfc_out   < 0.0d0) sfc_out   = 0.0d0
      where (sfc_chan  < 0.0d0) sfc_chan  = 0.0d0
      where (sfc_lake  < 0.0d0) sfc_lake  = 0.0d0
      where (sfc_bdry  < 0.0d0) sfc_bdry  = 0.0d0
      where (sfc_rem   < 0.0d0) sfc_rem   = 0.0d0
      where (sub_exfil < 0.0d0) sub_exfil = 0.0d0
      where (q_intf    < 0.0d0) q_intf    = 0.0d0

      allocate(src(ntr,ncell), inflow(ntr,ncell), remain(ntr,ncell), wout(ncell))
      src = 0.0d0; inflow = 0.0d0; remain = 0.0d0

      ! --- surface loads of the cells (this time-step)
      do j = 1, jx
         do i = 1, ix
            k = i + (j-1)*ix
            do ips = 1, nps
               src(ips,k) = max(0.0d0, overSed%Ssurf(ips,i,j)) * dt
            enddo
            src(5,k)  = So_LPOCsurf(i,j); src(6,k)  = So_RPOCsurf(i,j); src(7,k)  = So_LDOCsurf(i,j)
            src(8,k)  = So_RDOCsurf(i,j); src(9,k)  = So_MBMCsurf(i,j)
            src(10,k) = So_LPONsurf(i,j); src(11,k) = So_RPONsurf(i,j); src(12,k) = So_LDONsurf(i,j)
            src(13,k) = So_RDONsurf(i,j); src(14,k) = So_MBMNsurf(i,j); src(15,k) = So_NH4surf(i,j)
            src(16,k) = So_NO3surf(i,j)
            src(17,k) = So_LPOPsurf(i,j); src(18,k) = So_RPOPsurf(i,j); src(19,k) = So_LDOPsurf(i,j)
            src(20,k) = So_RDOPsurf(i,j); src(21,k) = So_MBMPsurf(i,j); src(22,k) = So_PO4surf(i,j)
            src(23,k) = So_PIPAsurf(i,j); src(24,k) = So_PIPSsurf(i,j)
         enddo
      enddo
      where (src < 0.0d0) src = 0.0d0
      bud(:,1) = sum(src, dim=2)

      ! --- 1. interflow: to the soil of the receiving cell, exfiltrated part to its surface
      allocate(gross(ix,jx), fex(ix,jx))
      gross = 0.0d0
      do j = 1, jx
         do i = 1, ix
            if (q_intf(i,j) <= 0.0d0) cycle
            if (.not. receiver(sub_dir(i,j), i, j, ri, rj)) cycle
            if (ri == 1 .or. ri == ix .or. rj == 1 .or. rj == jx) cycle   ! leaves the domain
            gross(ri,rj) = gross(ri,rj) + q_intf(i,j)
         enddo
      enddo
      fex = 0.0d0
      where (gross > 0.0d0) fex = min(1.0d0, sub_exfil / gross)
      call move_intf(So_LDOCintf, So_LDOC, ntr_s+1)
      call move_intf(So_RDOCintf, So_RDOC, ntr_s+2)
      call move_intf(So_LDONintf, So_LDON, ntr_s+3)
      call move_intf(So_RDONintf, So_RDON, ntr_s+4)
      call move_intf(So_NH4intf,  So_NH4,  ntr_s+5)
      call move_intf(So_NO3intf,  So_NO3,  ntr_s+6)
      call move_intf(So_LDOPintf, So_LDOP, ntr_s+7)
      call move_intf(So_RDOPintf, So_RDOP, ntr_s+8)
      call move_intf(So_PO4intf,  So_PO4,  ntr_s+9)

      ! --- 2. surface: flow graph of the main overland directions
      allocate(nxt(ncell), indeg(ncell), order(ncell), done(ncell), walk(ncell), cyc(ncell))
      nxt = 0; indeg = 0
      do j = 1, jx
         do i = 1, ix
            k = i + (j-1)*ix
            wout(k) = sfc_out(i,j) + sfc_chan(i,j) + sfc_lake(i,j) + sfc_bdry(i,j) + sfc_rem(i,j)
            if (sfc_out(i,j) > 0.0d0) then
               if (receiver(sfc_dir(i,j), i, j, ri, rj)) then
                  nxt(k) = ri + (rj-1)*ix
                  indeg(nxt(k)) = indeg(nxt(k)) + 1
               endif
            endif
         enddo
      enddo
      ! upslope-to-downslope order (Kahn); cells left over form closed loops
      nord = 0
      do k = 1, ncell
         if (indeg(k) == 0) then
            nord = nord + 1
            order(nord) = k
         endif
      enddo
      head = 1
      do while (head <= nord)
         k = order(head); head = head + 1
         if (nxt(k) > 0) then
            indeg(nxt(k)) = indeg(nxt(k)) - 1
            if (indeg(nxt(k)) == 0) then
               nord = nord + 1
               order(nord) = nxt(k)
            endif
         endif
      enddo
      done = .false.
      do n = 1, nord
         call cell_out(order(n), src(:,order(n)) + inflow(:,order(n)), .true.)
         done(order(n)) = .true.
      enddo
      ! closed loops c1 -> c2 -> ... -> cm -> c1: inflow x of c1 from cm solved exactly
      walk = 0
      do kk = 1, ncell
         if (done(kk)) cycle
         ! find the loop reached from kk
         k = kk
         do while (walk(k) /= kk)
            walk(k) = kk
            if (nxt(k) == 0) exit   ! not expected: cells outside the Kahn order lie on loops
            k = nxt(k)
         enddo
         if (nxt(k) == 0) then
            call cell_out(k, src(:,k) + inflow(:,k), .true.)
            done(k) = .true.
            cycle
         endif
         ncyc = 0
         m = k
         do
            ncyc = ncyc + 1
            cyc(ncyc) = m
            m = nxt(m)
            if (m == k) exit
         enddo
         ! out_i = a_i (T_i + y_i), y_1 = x, y_{i+1} = out_i, x = out_m = amul x + bvec
         amul = 1.0d0; bvec = 0.0d0
         do n = 1, ncyc
            m = cyc(n)
            a = 0.0d0
            if (wout(m) > 0.0d0) a = sfc_out(ci(m), cj(m)) / wout(m)
            amul = a * amul
            bvec = a * (src(:,m) + inflow(:,m) + bvec)
         enddo
         degenerate = (1.0d0 - amul < 1.0d-12)
         if (degenerate) then
            y = 0.0d0          ! no water leaves the loop: the load of the last cell remains there
         else
            y = bvec / (1.0d0 - amul)
         endif
         do n = 1, ncyc
            m = cyc(n)
            tot = src(:,m) + inflow(:,m) + y
            a = 0.0d0
            if (wout(m) > 0.0d0) a = sfc_out(ci(m), cj(m)) / wout(m)
            y = a * tot
            if (degenerate .and. n == ncyc) then
               call cell_out(m, tot, .false.)
               remain(:,m) = remain(:,m) + y
            else
               call cell_out(m, tot, .false.)
            endif
            done(m) = .true.
         enddo
      enddo

      ! --- 3. remaining loads back to the soil (top layer) / overland sediment storage
      do j = 1, jx
         do i = 1, ix
            k = i + (j-1)*ix
            if (all(remain(:,k) == 0.0d0)) cycle
            do ips = 1, nps
               overSed%Sol(ips,i,j) = overSed%Sol(ips,i,j) + remain(ips,k)
            enddo
            So_LPOCRES(i,1,j) = So_LPOCRES(i,1,j) + remain(5,k)
            So_RPOC(i,1,j)    = So_RPOC(i,1,j)    + remain(6,k)
            So_LDOC(i,1,j)    = So_LDOC(i,1,j)    + remain(7,k)  + remain(ntr_s+1,k)
            So_RDOC(i,1,j)    = So_RDOC(i,1,j)    + remain(8,k)  + remain(ntr_s+2,k)
            So_MBMC(i,1,j)    = So_MBMC(i,1,j)    + remain(9,k)
            So_LPONRES(i,1,j) = So_LPONRES(i,1,j) + remain(10,k)
            So_RPON(i,1,j)    = So_RPON(i,1,j)    + remain(11,k)
            So_LDON(i,1,j)    = So_LDON(i,1,j)    + remain(12,k) + remain(ntr_s+3,k)
            So_RDON(i,1,j)    = So_RDON(i,1,j)    + remain(13,k) + remain(ntr_s+4,k)
            So_MBMN(i,1,j)    = So_MBMN(i,1,j)    + remain(14,k)
            So_NH4(i,1,j)     = So_NH4(i,1,j)     + remain(15,k) + remain(ntr_s+5,k)
            So_NO3(i,1,j)     = So_NO3(i,1,j)     + remain(16,k) + remain(ntr_s+6,k)
            So_LPOPRES(i,1,j) = So_LPOPRES(i,1,j) + remain(17,k)
            So_RPOP(i,1,j)    = So_RPOP(i,1,j)    + remain(18,k)
            So_LDOP(i,1,j)    = So_LDOP(i,1,j)    + remain(19,k) + remain(ntr_s+7,k)
            So_RDOP(i,1,j)    = So_RDOP(i,1,j)    + remain(20,k) + remain(ntr_s+8,k)
            So_MBMP(i,1,j)    = So_MBMP(i,1,j)    + remain(21,k)
            So_PO4(i,1,j)     = So_PO4(i,1,j)     + remain(22,k) + remain(ntr_s+9,k)
            So_PIPA(i,1,j)    = So_PIPA(i,1,j)    + remain(23,k)
            So_PIPS(i,1,j)    = So_PIPS(i,1,j)    + remain(24,k)
         enddo
      enddo
      bud(:,10) = sum(remain, dim=2)

      ! --- 4. budget (debug/lateral_budget.csv): sediment, TC, TN, TP
      ! surface: source + exfiltrated = channel + lake + boundary + unlinked + remaining
      ! interflow: moved = to soil + exfiltrated + boundary
      bud(:,11) = bud(:,1) + bud(:,2) - bud(:,3) - bud(:,4) &
                - (bud(:,6) + bud(:,7) + bud(:,8) + bud(:,9) + bud(:,10))
      call write_budget(itime, 'SED', [(n, n = 1, 4)])
      call write_budget(itime, 'TC',  [5, 6, 7, 8, 9, ntr_s+1, ntr_s+2])
      call write_budget(itime, 'TN',  [10, 11, 12, 13, 14, 15, 16, ntr_s+3, ntr_s+4, ntr_s+5, ntr_s+6])
      call write_budget(itime, 'TP',  [17, 18, 19, 20, 21, 22, 23, 24, ntr_s+7, ntr_s+8, ntr_s+9])
      if (maxval(abs(bud(:,11))) > 1.0d-8 * max(1.0d-30, maxval(bud(:,1) + bud(:,2)))) &
         write(6,*) 'WARNING: Lateral_Transport: load balance residual ', maxval(abs(bud(:,11))), ' at itime ', itime

      deallocate(sfc_in, sfc_out, sfc_dir, sfc_chan, sfc_lake, sfc_bdry, sfc_rem, sub_dir, sub_exfil, q_intf)
      deallocate(src, inflow, remain, wout, gross, fex, nxt, indeg, order, done, walk, cyc)

   contains

      integer function ci(kk)
         integer, intent(in) :: kk
         ci = mod(kk-1, ix) + 1
      end function ci

      integer function cj(kk)
         integer, intent(in) :: kk
         cj = (kk-1) / ix + 1
      end function cj

      !> distributes the load tot of cell kk over its surface water outputs
      subroutine cell_out(kk, tot, to_neighbour)
         integer, intent(in) :: kk
         real(8), intent(in) :: tot(ntr)
         logical, intent(in) :: to_neighbour   ! .false. for loop cells (neighbour outflow handled by the loop solution)
         integer             :: ii, jj, lk
         real(8)             :: cc(ntr)

         ii = ci(kk); jj = cj(kk)
         if (wout(kk) <= 0.0d0) then
            remain(:,kk) = remain(:,kk) + tot
            return
         endif
         cc = tot / wout(kk)
         if (sfc_out(ii,jj) > 0.0d0 .and. to_neighbour) then
            if (nxt(kk) > 0) then
               inflow(:,nxt(kk)) = inflow(:,nxt(kk)) + cc * sfc_out(ii,jj)
            else
               bud(:,8) = bud(:,8) + cc * sfc_out(ii,jj)
            endif
         endif
         if (sfc_chan(ii,jj) > 0.0d0) then
            lk = cell_link(ii,jj)
            if (lk > 0) then
               latS_ch(:,lk) = latS_ch(:,lk) + cc(1:ntr_s) * sfc_chan(ii,jj)
               latI_ch(:,lk) = latI_ch(:,lk) + cc(ntr_s+1:ntr) * sfc_chan(ii,jj)
               bud(:,6) = bud(:,6) + cc * sfc_chan(ii,jj)
            else
               bud(:,9) = bud(:,9) + cc * sfc_chan(ii,jj)
            endif
         endif
         if (sfc_lake(ii,jj) > 0.0d0) then
            lk = cell_lake(ii,jj)
            if (lk > 0) then
               latS_ch(:,lk) = latS_ch(:,lk) + cc(1:ntr_s) * sfc_lake(ii,jj)
               latI_ch(:,lk) = latI_ch(:,lk) + cc(ntr_s+1:ntr) * sfc_lake(ii,jj)
               bud(:,7) = bud(:,7) + cc * sfc_lake(ii,jj)
            else
               bud(:,9) = bud(:,9) + cc * sfc_lake(ii,jj)
            endif
         endif
         bud(:,8) = bud(:,8) + cc * sfc_bdry(ii,jj)
         remain(:,kk) = remain(:,kk) + cc * sfc_rem(ii,jj)
      end subroutine cell_out

      !> moves the interflow load intf(i,isl,j) of each cell to its receiving cell: soil pool
      !! pool(ri,isl,rj) and the exfiltrated fraction to the surface tracer itr of the receiving cell
      subroutine move_intf(intf, pool, itr)
         real(8), intent(in)    :: intf(:,:,:)
         real(8), intent(inout) :: pool(:,:,:)
         integer, intent(in)    :: itr
         integer                :: ii, jj, ll, rri, rrj
         real(8)                :: load

         do jj = 1, jx
            do ii = 1, ix
               if (q_intf(ii,jj) <= 0.0d0) cycle
               do ll = 1, domain%nsl
                  load = max(0.0d0, intf(ii,ll,jj))
                  if (load <= 0.0d0) cycle
                  bud(itr,2) = bud(itr,2) + load
                  if (.not. receiver(sub_dir(ii,jj), ii, jj, rri, rrj)) then
                     pool(ii,ll,jj) = pool(ii,ll,jj) + load         ! no receiving cell: stays in the cell soil
                     bud(itr,3) = bud(itr,3) + load
                  elseif (rri == 1 .or. rri == ix .or. rrj == 1 .or. rrj == jx) then
                     bud(itr,4) = bud(itr,4) + load                 ! leaves the domain
                  else
                     pool(rri,ll,rrj) = pool(rri,ll,rrj) + (1.0d0 - fex(rri,rrj)) * load
                     src(itr, rri + (rrj-1)*ix) = src(itr, rri + (rrj-1)*ix) + fex(rri,rrj) * load
                     bud(itr,3) = bud(itr,3) + (1.0d0 - fex(rri,rrj)) * load
                     bud(itr,5) = bud(itr,5) + fex(rri,rrj) * load
                  endif
               enddo
            enddo
         enddo
      end subroutine move_intf

      subroutine write_budget(it, grp, idx)
         integer,          intent(in) :: it
         character(len=*), intent(in) :: grp
         integer,          intent(in) :: idx(:)
         real(8)                      :: v(11)
         integer                      :: q
         do q = 1, 11
            v(q) = sum(bud(idx,q))
         enddo
         ! columns: surface source, interflow moved, to soil, exfiltrated, interflow boundary loss,
         ! to channel, to lake, boundary, unlinked, remaining to soil, residual
         write(ubud,'(I8,",",A,11(",",ES14.6))') it, grp, v(1), v(2), v(3), v(5), v(4), v(6), v(7), &
                                                   v(8), v(9), v(10), v(11)
      end subroutine write_budget

   end subroutine Lateral_Transport


   !> @brief Receiving neighbour (ri,rj) of cell (i,j) for the direction code
   !! 3*(dj+1) + (di+1) + 1 (1..9, 5 = none); .false. if none or outside the grid.
   logical function receiver(dircode, i, j, ri, rj)
      implicit none
      real(8), intent(in)  :: dircode
      integer, intent(in)  :: i, j
      integer, intent(out) :: ri, rj
      integer              :: code

      receiver = .false.
      ri = i; rj = j
      code = nint(dircode)
      if (code < 1 .or. code > 9 .or. code == 5) return
      ri = i + mod(code-1, 3) - 1
      rj = j + (code-1) / 3 - 1
      if (ri < 1 .or. ri > domain%ix .or. rj < 1 .or. rj > domain%jx) return
      receiver = .true.
   end function receiver

end module module_lateralTransport
!
!=====||__WHQ5403__||=====!
