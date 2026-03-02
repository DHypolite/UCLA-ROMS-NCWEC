! Iterative isotropic Laplacian smoothing of brk(i,j) at rho-points.
! Adapted from lmd_kpp_smooth_hbl.h.  N_brk_smooth passes with MPI
! exchange between each pass.  Mask-aware (Neumann BC at coastline).
! Conservative (Laplacian with Neumann BC preserves area integral).
!
! Stencil weights (isotropic Laplacian):
!
!       1/32  1/8  1/32
!       1/8   3/8  1/8
!       1/32  1/8  1/32
!
! Uses module-level scratch arrays FX_brk, FE_brk, FE1_brk.

      cff_s=1.D0/12.D0 ; cff1_s=3.D0/16.D0

      do iter_brk = 1, N_brk_smooth

#  ifdef EXCHANGE
        call exchange_xxx(brk)
#  endif

! Fill ghost points at physical boundaries (Neumann BC).
! For periodic directions the MPI exchange handles it.

#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr-1,jend+1
            brk(istr-1,j)=brk(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr-1,jend+1
            brk(iend+1,j)=brk(iend,j)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr-1,iend+1
            brk(i,jstr-1)=brk(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr-1,iend+1
            brk(i,jend+1)=brk(i,jend)
          enddo
        endif
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
          brk(istr-1,jstr-1)=brk(istr,jstr)
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
          brk(istr-1,jend+1)=brk(istr,jend)
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          brk(iend+1,jstr-1)=brk(iend,jstr)
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
          brk(iend+1,jend+1)=brk(iend,jend)
        endif
#   endif
#  endif

! Compute masked elementary differences and apply isotropic correction.

        do j=jstr-1,jend+1
          do i=istr,iend+1
            FX_brk(i,j)=(brk(i,j)-brk(i-1,j))
#  ifdef MASKING
     &                      *umask(i,j)
#  endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr-1,iend+1
            FE1_brk(i,j)=(brk(i,j)-brk(i,j-1))
#  ifdef MASKING
     &                      *vmask(i,j)
#  endif
          enddo
          do i=istr,iend
            FE_brk(i,j)=FE1_brk(i,j)
     &             + cff_s*( FX_brk(i+1,j)+FX_brk(i  ,j-1)
     &                      -FX_brk(i  ,j)-FX_brk(i+1,j-1))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend+1
            FX_brk(i,j)=FX_brk(i,j)
     &             + cff_s*( FE1_brk(i,j+1)+FE1_brk(i-1,j  )
     &                      -FE1_brk(i,j  )-FE1_brk(i-1,j+1))
          enddo
          do i=istr,iend
            brk(i,j)=brk(i,j)
     &             + cff1_s*( FX_brk(i+1,j)-FX_brk(i,j)
     &                       +FE_brk(i,j+1)-FE_brk(i,j))
#  ifdef MASKING
            brk(i,j)=brk(i,j)*rmask(i,j)
#  endif
          enddo
        enddo

      enddo  ! iter_brk
