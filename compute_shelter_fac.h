! compute_shelter_fac.h — One-time computation of shelter factor from
! land mask.  Included in init_arrays_wec_tile.
!
! Physical motivation:
!   WW3 runs at coarser resolution with a smoother coastline than ROMS.
!   Fine-scale embayments in ROMS would shelter waves from reaching the
!   interior, but WW3 cannot see them.  The breaking dissipation WW3
!   sends into those embayments is unphysical.  shelter_fac corrects
!   for this cross-resolution mismatch by scaling down brk(i,j) in
!   concave coastline regions.
!
! Algorithm:
!   1. shelter_fac = 1 - rmask  (land=1, ocean=0: "land influence")
!   2. N_shelter passes of isotropic Laplacian smoother (same stencil
!      as smooth_brk.h) bleed land signal into nearby ocean cells.
!   3. In embayments (3-sided land), land influence accumulates from
!      multiple directions → high value.  On straight coast → moderate.
!   4. Convert to scaling factor with threshold Lf0 using cosine taper:
!        shelter_fac = 1  where land_influence <= Lf0 (no reduction)
!        shelter_fac = 0.5*(1+cos(pi*(LI-Lf0)/(1-Lf0)))  for LI > Lf0
!      Cosine taper has zero derivative at both endpoints, eliminating
!      the gradient discontinuity that caused geometric artifacts with
!      the original linear ramp.
!   5. N_post passes of masked Laplacian smoother on shelter_fac to
!      remove staircase artifacts from grid-aligned coastline.
!
! Parameters:
!   N_shelter    = 8    detection radius  σ = Δx√(3N/8) ≈ 52m at Δx=30m
!   Lf0          = 0.5  straight-coast threshold (no reduction there)
!   N_post_shelt = 4    post-smoothing passes on shelter_fac (with masking)
!
! Required local variables (declared in calling subroutine):
!   iter_s, cff_sh, cff1_sh, Lf0_sh, cff_sum_sh, pi_sh
! Required module arrays:
!   shelter_fac, FX_brk, FE_brk, FE1_brk (allocated before this include)
!
! CPP guard: BREAKING && SHELTER_BRK

! ---- Parameters ----
      Lf0_sh = 0.5             ! straight-coast land_influence threshold
      cff_sh = 1.D0/12.D0     ! isotropy correction coefficient
      cff1_sh = 3.D0/16.D0    ! Laplacian weight

! ---- Step 1: Initialize land influence field ----
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          shelter_fac(i,j) = 1.0 - rmask(i,j)
        enddo
      enddo

! ---- Step 2: Smooth with N_shelter Laplacian passes ----
      do iter_s = 1, 8   ! N_shelter = 8

#  ifdef EXCHANGE
        call exchange_xxx(shelter_fac)
#  endif

! Fill ghost points at physical boundaries (Neumann BC)
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr-1,jend+1
            shelter_fac(istr-1,j)=shelter_fac(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr-1,jend+1
            shelter_fac(iend+1,j)=shelter_fac(iend,j)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr-1,iend+1
            shelter_fac(i,jstr-1)=shelter_fac(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr-1,iend+1
            shelter_fac(i,jend+1)=shelter_fac(i,jend)
          enddo
        endif
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
          shelter_fac(istr-1,jstr-1)=shelter_fac(istr,jstr)
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
          shelter_fac(istr-1,jend+1)=shelter_fac(istr,jend)
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          shelter_fac(iend+1,jstr-1)=shelter_fac(iend,jstr)
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
          shelter_fac(iend+1,jend+1)=shelter_fac(iend,jend)
        endif
#   endif
#  endif

! Isotropic Laplacian stencil (same as smooth_brk.h)
! Note: NO masking here — we want land signal to bleed into ocean

        do j=jstr-1,jend+1
          do i=istr,iend+1
            FX_brk(i,j)=(shelter_fac(i,j)-shelter_fac(i-1,j))
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr-1,iend+1
            FE1_brk(i,j)=(shelter_fac(i,j)-shelter_fac(i,j-1))
          enddo
          do i=istr,iend
            FE_brk(i,j)=FE1_brk(i,j)
     &             + cff_sh*( FX_brk(i+1,j)+FX_brk(i  ,j-1)
     &                       -FX_brk(i  ,j)-FX_brk(i+1,j-1))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend+1
            FX_brk(i,j)=FX_brk(i,j)
     &             + cff_sh*( FE1_brk(i,j+1)+FE1_brk(i-1,j  )
     &                       -FE1_brk(i,j  )-FE1_brk(i-1,j+1))
          enddo
          do i=istr,iend
            shelter_fac(i,j)=shelter_fac(i,j)
     &             + cff1_sh*( FX_brk(i+1,j)-FX_brk(i,j)
     &                        +FE_brk(i,j+1)-FE_brk(i,j))
          enddo
        enddo

      enddo  ! iter_s

! ---- Step 3: Convert land_influence to scaling factor ----
! shelter_fac currently holds the smoothed land_influence.
! Cosine taper:  zero derivative at both endpoints → no gradient kink.
!   land_influence <= Lf0  →  shelter_fac = 1.0
!   Lf0 < LI < 1          →  shelter_fac = 0.5*(1+cos(pi*(LI-Lf0)/(1-Lf0)))
!   land_influence >= 1    →  shelter_fac = 0.0

      pi_sh = 4.D0*atan(1.D0)

#  ifdef EXCHANGE
      call exchange_xxx(shelter_fac)
#  endif

      do j=jstr,jend
        do i=istr,iend
          if (rmask(i,j) > 0.5) then
            cff_sum_sh = shelter_fac(i,j)  ! land_influence value
            if (cff_sum_sh .le. Lf0_sh) then
              shelter_fac(i,j) = 1.0
            elseif (cff_sum_sh .ge. 1.0) then
              shelter_fac(i,j) = 0.0
            else
              shelter_fac(i,j) = 0.5*(1.0
     &           + cos(pi_sh*(cff_sum_sh-Lf0_sh)/(1.0-Lf0_sh)))
            endif
          else
            shelter_fac(i,j) = 0.0   ! land cells: brk already 0
          endif
        enddo
      enddo

! ---- Step 4: Post-smooth shelter_fac to remove staircase artifacts ----
! Same Laplacian stencil but WITH masking: land stays 0, smoothing
! only acts on ocean cells.  Removes grid-aligned patterns from the
! discrete coastline mask.

      do iter_s = 1, 4   ! N_post_shelt = 4

#  ifdef EXCHANGE
        call exchange_xxx(shelter_fac)
#  endif

! Fill ghost points at physical boundaries (Neumann BC)
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr-1,jend+1
            shelter_fac(istr-1,j)=shelter_fac(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr-1,jend+1
            shelter_fac(iend+1,j)=shelter_fac(iend,j)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr-1,iend+1
            shelter_fac(i,jstr-1)=shelter_fac(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr-1,iend+1
            shelter_fac(i,jend+1)=shelter_fac(i,jend)
          enddo
        endif
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
          shelter_fac(istr-1,jstr-1)=shelter_fac(istr,jstr)
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
          shelter_fac(istr-1,jend+1)=shelter_fac(istr,jend)
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          shelter_fac(iend+1,jstr-1)=shelter_fac(iend,jstr)
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
          shelter_fac(iend+1,jend+1)=shelter_fac(iend,jend)
        endif
#   endif
#  endif

! Masked Laplacian: multiply differences by mask at flux points
! so that land cells act as Neumann boundaries (no flux into land).

        do j=jstr-1,jend+1
          do i=istr,iend+1
            FX_brk(i,j)=(shelter_fac(i,j)-shelter_fac(i-1,j))
#  ifdef MASKING
     &                  * umask(i,j)
#  endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr-1,iend+1
            FE1_brk(i,j)=(shelter_fac(i,j)-shelter_fac(i,j-1))
#  ifdef MASKING
     &                   * vmask(i,j)
#  endif
          enddo
          do i=istr,iend
            FE_brk(i,j)=FE1_brk(i,j)
     &             + cff_sh*( FX_brk(i+1,j)+FX_brk(i  ,j-1)
     &                       -FX_brk(i  ,j)-FX_brk(i+1,j-1))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend+1
            FX_brk(i,j)=FX_brk(i,j)
     &             + cff_sh*( FE1_brk(i,j+1)+FE1_brk(i-1,j  )
     &                       -FE1_brk(i,j  )-FE1_brk(i-1,j+1))
          enddo
          do i=istr,iend
            shelter_fac(i,j)=shelter_fac(i,j)
     &             + cff1_sh*( FX_brk(i+1,j)-FX_brk(i,j)
     &                        +FE_brk(i,j+1)-FE_brk(i,j))
          enddo
        enddo

      enddo  ! iter_s (post-smooth)

#  ifdef EXCHANGE
      call exchange_xxx(shelter_fac)
#  endif

      if (mynode==0) then
        write(*,'(7x,A,I2,A,F4.2,A,I2,A)')
     &   'SHELTER_BRK: computed shelter_fac (N_shelter=',8,
     &   ', Lf0=',Lf0_sh,', N_post=',4,')'
      endif
