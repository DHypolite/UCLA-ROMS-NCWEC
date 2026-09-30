      module wec_frc

      ! ===================================
      ! Wave Effect on Current (WEC) module
      ! ===================================
      ! read in surface gravity wave forcing 
      !  compute stokes drift and other wave-averaged term
      !                 For details see: 
      !             McWilliams, Restrepo, & Lane 2004 (JFM)
      !             Uchiyama et al. 2010 (Ocean Modeling)
      !             Romero et al. 2021 (Ocean Modeling)
      ! Translation from mrl_wci.F: Devin Dollery & Jeroen Molemaker (2020)
      ! re-coded by    Delphine Hypolite, Daniel Dauhajre (2022) 
      ! Broad Spectra terms (e.g., ust,vst) and code is courtesy of Leonel Romero.
      !                 --------------------------
      !         LATEST MERGE TO STABLE REPOSITORY: 
      !         Daniel Dauhajre and
      !         Delphine Hypolite (2024)
      !                ----------------------------


#include "cppdefs.opt"
! CPP flag WEC covering entire module so that the module is empty if
! WEC is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined WEC /* WEC for whole module! */

      use param
      use nc_read_write
      use roms_read_write
      use netcdf
      use grid
      use dimensions
      use ocean_vars   
      use scalars
      use mpi_exchanges
      !                 add brk2dx,brk2de to these
      use surf_flux, only: sustr, svstr
      implicit none

      ! ****************************
      ! MODULE SPECIFICATION SECTION
      !
      ! This part of the module (before the 'contains' section)
      ! is used to declare all variables that might be needed globally
      ! or if there are memory benefits to declaring only once, rather
      ! than for each call of a subroutine.

      ! ****************************

      ! Make all variable private to this module unless public specified
      private
      save     ! save -> Make all public variables global (might not be needed depending on compiler)

#include "wec_frc.opt"


!--------------------------------------------------------------------
!  WIND INDUCED WAVES: everything is defined at rho-point
!--------------------------------------------------------------------
! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
! uorb | BBL     | xi-component  of wave-induced bed orbital velocity [m/s]
! vorb | BBL     | eta-component of wave-induced bed orbital velocity [m/s]
! wdrx | MRL     | cosine of wave direction [non dimension]
! wdre | MRL     | sine of   wave direction [non dimension]
! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
! rdsp | ROLLER  | roller dissipation rate (\epsilon_r term) [m3/s3]
!--------------------------------------------------------------------

      real,allocatable,dimension(:,:) :: wfrq
      real,allocatable,dimension(:,:) :: whrm
# ifdef UORBOFF
      real,allocatable,dimension(:,:) :: uorb
      real,allocatable,dimension(:,:) :: vorb
# endif
# ifdef SURFACE_ROLLER
      real,allocatable,dimension(:,:) :: rdsp
# endif


!-------------------------------------------------------------------
!  WAVE AVEREAGED QUANTITIES AND TERMS
!--------------------------------------------------------------------
!  2D  |  brk2dx   |   xi-direciton 2D breaking dissipation (rho)
!  2D  |  brk2de   |  eta-direction 2D breaking dissipation (rho)
!  2D  |  frc2dx   |   xi-direciton 2D frictional dissipation (rho)
!  2D  |  frc2de   |  eta-direction 2D frictional dissipation (rho)
!  2D  |  ust2d    |   xi-direciton Stokes transport (u-point)
!  2D  |  vst2d    |  eta-direciton Stokes transport (v-point)
!  2D  |  ustr2d   |   xi-direciton Stokes transport (r-point)
!  2D  |  vstr2d   |  eta-direciton Stokes transport (r-point)
!  2D  |  ustr0    |   xi-direciton Surface Stokes drift (r-point)
!  2D  |  vstr0    |  eta-direciton Surface Stokes drift (r-point)
!  2D  |  sup      |  quasi-static wave set-up (rho-point)
!  2D  |  calP     |  pressure correction term (rho-point)
!  2D  |  Kapsrf   |  Bernoulli head terrm at the surface (rho-point)
!  2D  |  lm       |  mean wave length (rho-point)
!--------------------------------------------------------------------
!  3D  |  brk3dx   |   xi-direciton 3D breaking dissipation (rho)
!  3D  |  brk3de   |  eta-direction 3D breaking dissipation (rho)
!  3D  |  ust      |   xi-direciton 3D Stokes drift velocity (u-point)
!  3D  |  vst      |  eta-direciton 3D Stokes drift velocity (v-point)
!  3D  |  wst      |       vertical 3D Stokes drift velocity (rho-point)
!  3D  |  kvf      |  vertical vortex force term (K term, 3D, rho-point)
!  3D  |  Akb      |  breaking-wave-induced additional diffusivity (w-point)
!  3D  |  Akw      |  wave-induced additional diffusivity (rho-point)
!  3D  |  E_pre    |  previous time-step value for Akw estimation (rho)
!  3D  |  frc3dx   |   xi-direction 3D frictional dissipation (rho)
!  3D  |  frc3de   |  eta-direction 3D frictional dissipation (rho)
!--------------------------------------------------------------------

      ! Anything new that is public is used for prsgrd.F terms
      ! in compute_prsgrd_wec_terms.h

      real,allocatable,dimension(:,:),public :: brk2dx
      real,allocatable,dimension(:,:),public :: brk2de
      real,allocatable,dimension(:,:),public :: ust2d
      real,allocatable,dimension(:,:),public :: vst2d
      real,allocatable,dimension(:,:)        :: ustr2d ! At rho point
      real,allocatable,dimension(:,:)        :: vstr2d ! At rho point
      real,allocatable,dimension(:,:)        :: ustr0
      real,allocatable,dimension(:,:)        :: vstr0
      real,allocatable,dimension(:,:)        :: llmm
      real,allocatable,dimension(:,:),public :: frc2dx ! used in prestep, step2d, step3d
      real,allocatable,dimension(:,:),public :: frc2de ! used in prestep, step2d, step3d
      real,allocatable,dimension(:,:),public :: sup    ! used in step2d
      real, allocatable,dimension(:,:)       :: sup_zeta !sup + zeta (diagnostic)
# ifdef SOLVE3D
      real,allocatable,dimension(:,:),public :: calP   ! used in prsgrd.F
      real,allocatable,dimension(:,:),public :: Kapsrf ! used in prsgrd.F
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      real,allocatable, dimension(:,:,:) :: brk_r ! used in momentum eqn in wec_frc.F
      real,allocatable, dimension(:,:,:) :: brk3dx ! used in momentum eqn in wec_frc.F
      real,allocatable, dimension(:,:,:)  :: brk3de ! used in momentum eqn in wec_frc.F
#  endif
#  ifdef BODY_FRICTION
      real,allocatable,dimension(:,:,:)  :: frc_r  ! streaming profile workspace
      real,allocatable,dimension(:,:,:)  :: frc3dx ! used in momentum eqn in wec_frc.F
      real,allocatable, dimension(:,:,:) :: frc3de ! used in momentum eqn in wec_frc.F
#  endif
      real,allocatable,dimension(:,:,:),public :: ust
      real,allocatable,dimension(:,:,:),public :: vst
      real,allocatable,dimension(:,:,:),public :: wst    ! used in set_avg.F
      real,allocatable,dimension(:,:,:),public :: kvf    ! used in prsgrd32AC1.F
      real,allocatable,dimension(:,:,:),public :: Akb    ! used in lmd_kpp.F
      real,allocatable,dimension(:,:,:),public :: Akw    ! used in lmd_kpp.F
      real,allocatable,dimension(:,:,:)        :: E_pre
# endif  /* SOLVE3D */

      !--------------------------------------------------------------------------
      !         Pubic variables for all of ROMS code
      !         These are needed for inclusion of vertical vortex force,
      !         Bernoulli head, and surface pressure correction in
      !         pressure gradient (prsgrd.F)
      !     
      !         Some of these terms also used in calculation of stokes
      !         drift or other wave-averaged terms in wec_forces
      !----------------------------------------------------------------------
      !real, allocatable,dimension(:,:),public :: inv_ex
      real, allocatable,dimension(:,:),public :: kD
      real, allocatable,dimension(:,:),public :: kw
      real, allocatable,dimension(:,:),public :: Dstp
      real, allocatable,dimension(:,:),public :: actf
      real, allocatable,dimension(:,:),public :: act
      real, allocatable,dimension(:,:),public :: wdrx
      real, allocatable,dimension(:,:),public :: wdre
      real, allocatable,dimension(:,:,:),public :: ust_r
      real, allocatable,dimension(:,:,:),public :: vst_r
      real, allocatable,dimension(:,:,:),public :: ust_rm
      real, allocatable,dimension(:,:,:),public :: inv_dz
# ifdef BBL_S95
      ! Orbital velocity magnitude for r_D calculation in
      ! compute_rd_bott_drag.h
      real, allocatable,dimension(:,:),public :: Uorb_mag
      real, allocatable,dimension(:,:),public :: tauw_orb 
# endif




# if !defined MONOWAVE
      !This is used in cff4VF for calP calculation in prsgrd.F
      real, allocatable,dimension(:,:),public :: keff
# endif

      real, allocatable, dimension(:) :: kb0
      real, allocatable, dimension(:,:) :: wh
      real, allocatable, dimension(:,:) :: fr
      real, allocatable, dimension(:,:) :: brk
      real, allocatable, dimension(:,:) :: inv_d
      real, allocatable, dimension(:,:) :: inv_f
      real, allocatable, dimension(:,:) :: frc
      real, allocatable, dimension(:,:) :: stk
      real, allocatable, dimension(:,:) :: ebrk
      real, allocatable, dimension(:,:) :: erol
      real, allocatable, dimension(:,:), public :: hvf_int_xi  ! depth-integrated horiz vortex force xi [m^2/s^2] (normalized by dm*dn)
      real, allocatable, dimension(:,:), public :: hvf_int_eta ! depth-integrated horiz vortex force eta [m^2/s^2] (normalized by dm*dn)
      real, allocatable, dimension(:,:) :: hvf_int_xi_avg
      real, allocatable, dimension(:,:) :: hvf_int_eta_avg
# if defined WAVE_OFFLINE && !defined MONOWAVE
      real, allocatable, dimension(:,:) :: fac2d
      real, allocatable, dimension(:,:) :: inv_exef
# endif



      ! =======================================================
      ! REALISTIC WEC VARIABLES (FROM READ IN NETCDF FILE)
      ! =====================================================

# if !defined ANA_WEC_FRC
!--------------------------------------------------------------------
!wwv_time|                    | time of wind-induced waves
!--------------------------------------------------------------------
!  wwag  |                    | wave amplitude [m]
!  wwdg  |  (READ IN FROM     | wave direction [radians]
!  wwpg  |   INPUT FILE)      | wave period [s]
!  wwuob |  point data        | xi-orbital velocity [m/s]
!  wwvob |  for wind induced  | eta-orbital velocity [m/s]
!  wwdrx |  waves             | cosine wave direction [ND]
!  wwdre |                    | sine wave direction [ND]
!  wweb  |                    | breaking dissipation [m3/s3]
!  wwed  |                    | frictional dissipation [m3/s3]
!  wwUSg |                    | Depth averaged Stokes drift [m/s]
!  wwVSg |                    | Depth averaged Stokes drift [m/s]
!  wwus0g|                    | Surface Stokes drift [m/s]
!  wwvs0g|                    | Surface Stokes drift [m/s]
!  wwsup |                    | Set down [m]
!  wwlm  |                    | Mean wavelength [m] - for diffusivity approx.
!--------------------------------------------------------------------
!  wveb  |  interpolated      | breaking dissipation [m3/s3]
!  wved  |  for present time  | frictional dissipation [m3/s3]
!  wvqb  |  step              | fraction of breaking waves [ND]
!--------------------------------------------------------------------

      real,allocatable,dimension(:,:)        :: wwag
      real,allocatable,dimension(:,:)        :: wwdg
      real,allocatable,dimension(:,:)        :: wwpg
#   if defined WAVE_OFFLINE && defined UORBOFF
      real,allocatable,dimension(:,:)        :: wwuob
      real,allocatable,dimension(:,:)        :: wwvob
#   endif
      real,allocatable,dimension(:,:)        :: wwdrx
      real,allocatable,dimension(:,:)        :: wwdre
      real,allocatable,dimension(:,:)        :: wweb
      real,allocatable,dimension(:,:),public   :: wveb
      real,allocatable,dimension(:,:),public   :: wved
      real,allocatable,dimension(:,:)        :: wwed
      real,allocatable,dimension(:,:)        :: wwUSg
      real,allocatable,dimension(:,:)        :: wwVSg
      real,allocatable,dimension(:,:)        :: wwus0g
      real,allocatable,dimension(:,:)        :: wwvs0g
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      real,allocatable,dimension(:,:),public   :: wvqb
      real,allocatable,dimension(:,:)        :: wwqb
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      real,allocatable,dimension(:,:)        :: wwsup
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      real,allocatable,dimension(:,:),public :: wwlm
#   endif

      ! Grid rotation angles
      !         MAY NOT BE NEEDED ANYMORE
      real,allocatable,dimension(:,:)          :: angler_cos
      real,allocatable,dimension(:,:)          :: angler_sin

# elif defined ANA_WEC_FRC
!Daniel added all these below
      real,allocatable,dimension(:,:,:)        :: wwag
      real,allocatable,dimension(:,:,:)        :: wwdg
      real,allocatable,dimension(:,:,:)        :: wwdrx
      real,allocatable,dimension(:,:,:)        :: wwdre
      real,allocatable,dimension(:,:,:)        :: wweb
      real,allocatable,dimension(:,:),public   :: wveb
      real,allocatable,dimension(:,:),public   :: wved
      real,allocatable,dimension(:,:,:)        :: wwed
      real,allocatable,dimension(:,:,:)        :: wwUSg
      real,allocatable,dimension(:,:,:)        :: wwVSg
      real,allocatable,dimension(:,:,:)        :: wwus0g
      real,allocatable,dimension(:,:,:)        :: wwvs0g
      real,allocatable,dimension(:,:,:)        :: wwpg
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      real,allocatable,dimension(:,:),public   :: wvqb
      real,allocatable,dimension(:,:,:)        :: wwqb
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      real,allocatable,dimension(:,:,:)        :: wwsup
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      real,allocatable,dimension(:,:,:),public :: wwlm
#   endif

# endif /* !ANA_WEC_FRC */

 
#ifdef SOLVE3D
      ! For wstokes
      real, allocatable, dimension(:) :: Wrk
#endif
     ! For compute_wec_rhs_uv_terms
      real, allocatable, dimension(:,:) :: utmp
      real, allocatable, dimension(:,:) :: vtmp



      ! AVERAGED VARIABLES
      ! ==================

      ! from averages.h in MRL_WCI old roms code
      real,allocatable,dimension(:,:)   :: sup_avg
      real,allocatable,dimension(:,:)   :: ust2d_avg
      real,allocatable,dimension(:,:)   :: vst2d_avg
      real,allocatable,dimension(:,:,:) :: ust_avg
      real,allocatable,dimension(:,:,:) :: vst_avg
      real,allocatable,dimension(:,:,:) :: wst_avg
      real,allocatable,dimension(:,:,:) :: akb_avg
      real,allocatable,dimension(:,:,:) :: akw_avg
      real,allocatable,dimension(:,:,:) :: kvf_avg
      real,allocatable,dimension(:,:)   :: calp_avg
      real,allocatable,dimension(:,:)   :: kapsrf_avg
      real,allocatable,dimension(:,:)   :: sup_zeta_avg
      real,allocatable,dimension(:,:)   :: frc2dx_avg
      real,allocatable,dimension(:,:)   :: frc2de_avg
      real,allocatable,dimension(:,:)   :: brk2dx_avg
      real,allocatable,dimension(:,:)   :: brk2de_avg
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      real,allocatable,dimension(:,:,:) :: brk3dx_avg
      real,allocatable,dimension(:,:,:) :: brk3de_avg
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
      real,allocatable,dimension(:,:,:) :: frc3dx_avg
      real,allocatable,dimension(:,:,:) :: frc3de_avg
#  endif
      real,allocatable,dimension(:,:)   :: wveb_avg
      real,allocatable,dimension(:,:)   :: wved_avg











      ! ========================
      ! NETCDF RELATED VARIABLES
      ! ========================
      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_wec=0
      real    :: t_avg_dia_wec=0
      integer,save :: navg_wec = 0               ! number of samples in average
      !integer,save :: navg_dia_wec = 0           ! number of samples in average
      integer :: record_avg = nrpf_avg           ! Triggers making of initial file
      integer :: record_his = nrpf_his
      !integer :: record_dia_avg = nrpf_avg_dia   ! Triggers making of initial file
     ! integer :: record_dia_his = nrpf_his_dia
      real,save :: output_time_his = 0
      real,save :: output_time_avg = 0
      
      ! ATTRIBUTES FOR WEC OUTPUT (input netcdf declared as ncforce in
      ! wec_frc.opt)

      ! Vertical average of u_stokes
      character(len=5) :: ust2d_name      = 'ust2d'
      character(len=38):: ust2d_long_name =
     &                   'u depth-averaged Stokes drift velocity'
      character(len=14):: ust2d_units     = 'meter second-1'

      ! Vertical average of v_stokes
      character(len=5) :: vst2d_name      = 'vst2d'
      character(len=38):: vst2d_long_name =
     &                   'v depth-averaged Stokes drift velocity'
      character(len=14):: vst2d_units     = 'meter second-1'

# ifdef SOLVE3D

      ! u_stokes 3d
      character(len=3) :: ust_name        = 'ust'
      character(len=26):: ust_long_name   = 'u 3D Stokes drift velocity'
      character(len=14):: ust_units       = 'meter second-1'

      ! v_stokes 3d
      character(len=3) :: vst_name        = 'vst'
      character(len=26):: vst_long_name   = 'v 3D Stokes drift velocity'
      character(len=14):: vst_units       = 'meter second-1'

      ! Read in data surface u & v stokes drfit
!      character(len=4) :: ust0_name  = 'ust0'
!      character(len=4) :: vst0_name  = 'vst0'

      ! w_stokes 3d
      character(len=3) :: wst_name        = 'wst'
      character(len=35):: wst_long_name   =
     &                   'w vertical 3D Stokes drift velocity'
      character(len=14):: wst_units       = 'meter second-1'

      ! Akb
      character(len=3) :: Akb_name        = 'Akb'
      character(len=35):: Akb_long_name   =
     &                   'eddy viscosity due to wave breaking'
      character(len=15):: Akb_units       = 'meter2 second-1'

      ! Akw
      character(len=3) :: Akw_name        = 'Akw'
      character(len=37):: Akw_long_name   =
     &                   'eddy diffusivity due to primary waves'
      character(len=15):: Akw_units       = 'meter2 second-1'

      ! kvf
      character(len=3) :: kvf_name        = 'kvf'
      character(len=34):: kvf_long_name   =
     &                   'vertical vortex force (u^st du/dz)'
      character(len=14):: kvf_units       = 'meter second-2'

      ! calP
      character(len=4) :: calP_name        = 'calP'
      character(len=37):: calP_long_name   =
     &                   'surface pressure correction in prsgrd'
      character(len=15):: calP_units       = 'meter2 second-2'

      ! Kapsrf
      character(len=6) :: Kapsrf_name        = 'Kapsrf'
      character(len=37):: Kapsrf_long_name   =
     &                   'surface Bernoulli head in prsgrd term'
      character(len=15):: Kapsrf_units       = 'meter2 second-2'

# endif

      character(len=2) :: wveb_name = 'eb'     
      character(len=40) :: wveb_long_name = 
     &            "breaking dissipation"
      character(len=15) :: wveb_units = "meter3 second-3"

      character(len=2) :: wved_name = 'ed'     
      character(len=40) :: wved_long_name = 
     &            "frictional dissipation"
      character(len=15) :: wved_units = "meter3 second-3"

      character(len=7) :: frc2dx_name = 'frc2dx'      
      character(len=40) :: frc2dx_long_name = 
     &            "xi-dir frictional dissipation (stress)"
      character(len=15) :: frc2dx_units = "meter2 second-2"

      character(len=7) :: frc2de_name = 'frc2de'      
      character(len=40) :: frc2de_long_name = 
     &            "eta-dir frictional dissipation (stress)"
      character(len=15) :: frc2de_units = "meter2 second-2"

      character(len=7) :: brk2dx_name = 'brk2dx'      
      character(len=40) :: brk2dx_long_name = 
     &            "xi-dir breaking dissipation (stress)"
      character(len=15) :: brk2dx_units = "meter2 second-2"

      character(len=7) :: brk2de_name = 'brk2de'
      character(len=40) :: brk2de_long_name =
     &            "eta-dir breaking dissipation (stress)"
      character(len=15) :: brk2de_units = "meter2 second-2"

      character(len=7) :: brk3dx_name = 'brk3dx'
      character(len=40) :: brk3dx_long_name =
     &            "xi-dir 3D breaking body force"
      character(len=15) :: brk3dx_units = "meter second-2"

      character(len=7) :: brk3de_name = 'brk3de'
      character(len=40) :: brk3de_long_name =
     &            "eta-dir 3D breaking body force"
      character(len=15) :: brk3de_units = "meter second-2"

      character(len=7) :: frc3dx_name = 'frc3dx'
      character(len=40) :: frc3dx_long_name =
     &            "xi-dir 3D streaming body force"
      character(len=15) :: frc3dx_units = "meter second-2"

      character(len=7) :: frc3de_name = 'frc3de'
      character(len=40) :: frc3de_long_name =
     &            "eta-dir 3D streaming body force"
      character(len=15) :: frc3de_units = "meter second-2"

      character(len=10) :: hvf_int_xi_name = 'hvf_int_xi'
      character(len=40) :: hvf_int_xi_long_name =
     &          "Depth-int horiz vortex force, xi-comp"
      character(len=15) :: hvf_int_xi_units = "meter2 second-2"

      character(len=11) :: hvf_int_eta_name = 'hvf_int_eta'
      character(len=40) :: hvf_int_eta_long_name =
     &          "Depth-int horiz vortex force, eta-comp"
      character(len=15) :: hvf_int_eta_units = "meter2 second-2"

      character(len=3)  :: rD_name = 'r_D'
      character(len=40) :: rD_long_name =
     &          "Bottom drag coefficient"
      character(len=15) :: rD_units = "meter second-1"





      character(len=2) :: wec_edwv_name = 'ed'      !! wweb
      character(len=3) :: wec_lwv_name = 'lmw'      !! Mean Wave length
# ifdef UORBOFF
      character(len=4) :: uorb_name = 'uorb'
      character(len=4) :: vorb_name = 'vorb'
# endif      

      ! sup
      character(len=3) :: sup_name      = 'sup'
      character(len=40):: sup_long_name =
     &                   'quasi-static sea-level response (set-up)'
      character(len=5) :: sup_units     = 'meter'

      ! sup_zeta
      character(len=8) :: sup_zeta_name      = 'sup_zeta'
      character(len=40):: sup_zeta_long_name =
     &        'quasi-static sea-level response (set-up) + free surface'
      character(len=5) :: sup_zeta_units     = 'meter'





      !! Wave forcing time (Also its dimension name in netcdf file!!)
      character(len=8) :: wec_tim_name = 'wwv_time'


      ! Add netcdf WEC terms previously contained in ncvars
      ! ===================================================

      ! Array of bools for yes/no whether to output
      ! result variable to netcdf file
      ! ------------------------------

      ! Size of array hard coded, change as required
      !logical wrt_wec_vars(11) ! Previously part of wrthis in old roms
      ! same as wrt_wec_vars but for averaging of wec variables:
      !logical wrt_avg_vars(11) ! Previously part of wrtavg in old roms

      ! Index numbers are arbitrary, they are just place holders for
      ! the true or false bool for outputing WEC variables stored in
      ! wrt_wec_vars
!      integer, parameter :: indxSUP=1
!     &                    , indxUST2D=indxSUP+1, indxVST2D=indxSUP+2
!# ifdef SOLVE3D
!     &                    , indxUST=indxSUP+3, indxVST=indxSUP+4
!     &                    , indxAkb=indxSUP+6
!     &                    , indxAkw=indxSUP+7, indxKVF=indxSUP+8
!     &                    , indxCALP=indxSUP+9, indxKAPS=indxSUP+10
      ! indxWST used in set_avg.F
!      integer, parameter, public :: indxWST=indxSUP+5
!# endif

      ! ===============
      ! MISC. VARIABLES
      ! ===============

 
      ! coefficient variables for external calculations in prsgrd32AC1.F
      ! real, public :: cff_ext, cff_ext2 --> now just declare these in
      ! prsgrd.F in ifdef WEC

     ! Rename grid variables
      integer, parameter :: gnx = LLm
      integer, parameter :: nz  = N

      ! timestep variables for read in WEC data
      integer              :: it1 = 1  ! Placeholder for read in time
      integer              :: it2 = 2  ! Placeholder for read in time
      real, dimension(2)   :: wec_tim = [-99,-99] ! used to store both time steps (n <-> n+1)
      integer    :: ifile = 0 ! read in force file number
      ! current record number of data point in read in data - irec
      integer    :: irec  = 0


      logical :: wec_init_done = .false. ! T/F to confirm if wec has been initialized

      ! ================================
      ! WEC subroutines needed as public (accessed outside of module)
      ! ================================

!----------------------------------------------------------
      !         DEPRACATED FUNCTIONS 
      !public read_inp_wec         ! called from read_inp.F
      !public set_wec_frc_tile     ! called from set_forces.F
      !public wrt_his_wec          ! called from wrt_his.F
      !public def_his_def_var_wec  ! called from def_his.F --> OUT OF
! GET RID OF AVERAGES
!#ifdef AVERAGES
!      public read_inp_avg_wec     ! called from read_inp.F
!      public def_avg_def_var_wec  ! called from def_his.F
!      public set_avg_wec_tile     ! called from set_avg.F
!      public wrt_avg_wec          ! called from wrt_avg.F
!#endif
!----------------------------------------------------------



      public set_wec_frc          ! called from set_forces.F
      public init_arrays_wec_tile ! called from init_arrays.F
      public compute_wec_rhs_uv_terms ! called from step3d_uv1.F & pre_step3d4S.F
      !public wstokes_wec_tile     ! calle
      !public wec_forces           ! called from main.F (NOT ANYMORE,2024!)
      public wec_forces_tile      ! called from set_forces.F 
      public wrt_wec              ! called from main.F

      ! =======================
      ! END OF MODULE PRE-AMBLE
      ! =======================

      contains 
! ----------------------------------------------------------------------


      ! ######################
      ! ANALYTICAL SUBROUTINES
      ! ######################

#ifdef ANA_WEC_FRC
      subroutine ana_init_wec_frc_tile (istr,iend,jstr,jend)  ![

      use tracers, only: t

      implicit none

      integer, intent(in) :: istr, iend, jstr, jend

      integer :: i, j, k
        !Daniel changed to 90
      real    :: d_wave = 90. / 360.0 * 2 * pi
      real    :: cos_dw, sin_dw

      ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
      ! This had to appear after other variables declared else wont compile
#include "compute_auxiliary_bounds.h"

      cos_dw=cos(d_wave)
      sin_dw=cos(d_wave)
      if (mynode.eq.0) then
         print *, 'jstr-1=',jstr-1
         print *, 'jend+1=',jend+1
      endif
      ! removed -1 to try and debug heisenbug
      !do j=jstr,jend
      !  do i=istr,iend
      do j=jstr-1,jend+1
        do i=istr-1,iend+1

          !Chao's LES runs
          ! wave amplitude = 0.8 m
          ! wave length = 60 m
          ! turbulent langmuir number = 0.3
          whrm(i,j)= 0.8 * 2 ! 
          wfrq(i,j)= 1.01 !2*pi/6.2 sec (Tw = 6.2 sec from Yan et al 2021) 
                  !--> taken from deep water linear dispersion with D =120m, k = 1/60m 
          wdrx(i,j)=1.0       ! cosine wave direction (xi)
          wdre(i,j)= 0.0         ! sine wave direction (eta)
          wveb(i,j)= 0.0
          wved(i,j)= 0.0
        
          !ustr2d(i,j)=0.00297
          !vstr2d(i,j)=0.
          !ustr0(i,j) =0.068
          !vstr0(i,j) =0.
          !For monochromatic should not need these
          !ustr2d(i,j)=0.0
          !vstr2d(i,j)=0.
          !ustr0(i,j) =0.
          !vstr0(i,j) =0.
          sup(i,j)   =0.
          !llmm(i,j)  = 60. 
       
        enddo
      enddo

      wec_init_done = .true.
      if(mynode==0) print *, ' ana_init_wec done!'

      end subroutine ana_init_wec_frc_tile  !]
#endif /* ANA_WEC_FRC */


      !-------------------------------------------------
#ifdef ANA_WEC_FRC
      subroutine ana_wec_frc_tile (istr,iend,jstr,jend)  ![
      ! Analytical WEC forcing
      implicit none

      integer, intent(in) :: istr, iend, jstr, jend
      integer             :: i, j
      if (.not. wec_init_done) then
        call ana_init_wec_frc_tile(istr,iend,jstr,jend)
      endif

      end subroutine ana_wec_frc_tile  !]
#endif /* ANA_WEC_FRC */

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      ! ###################
      ! FORCING SUBROUTINES
      ! ###################

      subroutine init_arrays_wec_tile (istr,iend,jstr,jend)

      ! This is called from init_arrays_tile of init_arrays.F
      ! It is needed to "first touch" arrays to optimize
      ! NUMA shared memory allocation on linux.

      implicit none

      integer, intent(in) :: istr, iend, jstr, jend
      ! Declare looping indices
      integer :: i, j, k

      real, parameter :: init=0.0

      ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
      ! This had to declared last else won't compile
#include "compute_auxiliary_bounds.h"

      ! TESTING NCFORCE ALLOCATE
      allocate( nc_wwag%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwdg%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwpg%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wweb%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwed%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwqb%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwuob%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwvob%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwUsg%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwVsg%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwus0g%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwvs0g%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwsup%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_wwlmv%vdata( GLOBAL_2D_ARRAY,2) )










      !DPD adding, debug heisenbug
      !originally thought bug fix was to change ust_rm, inv_dz
      !allocation to N instead of 0:N, but then April 27 testing shows
      !that N or 0:N works fine. In true Heisenbug fashion, the error
      !appears and disappears with no explanation
      allocate( wfrq(GLOBAL_2D_ARRAY) )
      wfrq=0.
      allocate( ust_rm(GLOBAL_2D_ARRAY,N))
      ust_rm=0.
      allocate( ust_r(GLOBAL_2D_ARRAY,N))
      ust_r=0.
      allocate( vst_r(GLOBAL_2D_ARRAY,N))
      vst_r=0.
      allocate( inv_dz(GLOBAL_2D_ARRAY,N))
      inv_dz=0.
# ifdef BBL_S95
      allocate( Uorb_mag(GLOBAL_2D_ARRAY))
      Uorb_mag=0.
      allocate( tauw_orb(GLOBAL_2D_ARRAY))
      tauw_orb=0.
# endif
      !allocate( inv_ex(GLOBAL_2D_ARRAY))

      !Allocate what were A2d arrays for wec_force_tile 
      allocate( wh(GLOBAL_2D_ARRAY))
      wh=0.
      allocate( fr(GLOBAL_2D_ARRAY))
      fr=0.
      allocate( brk(GLOBAL_2D_ARRAY))
      brk=0.
      allocate( inv_d(GLOBAL_2D_ARRAY))
      inv_d=0.
      allocate( inv_f(GLOBAL_2D_ARRAY))
      inv_f=0.
      allocate( frc(GLOBAL_2D_ARRAY))
      frc=0.
      allocate( stk(GLOBAL_2D_ARRAY))
      stk=0.
      allocate( ebrk(GLOBAL_2D_ARRAY))
      ebrk=0.
      allocate( erol(GLOBAL_2D_ARRAY))
      erol=0.
      allocate( hvf_int_xi(GLOBAL_2D_ARRAY))
      hvf_int_xi=0.
      allocate( hvf_int_eta(GLOBAL_2D_ARRAY))
      hvf_int_eta=0.

# if defined WAVE_OFFLINE && !defined MONOWAVE
      allocate( fac2d(GLOBAL_2D_ARRAY))
      fac2d=0.
      allocate( inv_exef(GLOBAL_2D_ARRAY))
      inv_exef=0.
# endif

      allocate(kb0(0:N))
      kb0=0.
#ifdef UORBOFF
      allocate( uorb(GLOBAL_2D_ARRAY) )
      uorb=0.
      allocate( vorb(GLOBAL_2D_ARRAY) )
      vorb=0.
#endif
      allocate( whrm(GLOBAL_2D_ARRAY) )
      whrm=0.
      allocate( wdrx(GLOBAL_2D_ARRAY) )
      wdrx=0.
      allocate( wdre(GLOBAL_2D_ARRAY) )
      wdre=0.
# ifdef SURFACE_ROLLER
      allocate( rdsp(GLOBAL_2D_ARRAY) )
      rdsp=0.
# endif
      allocate( brk2dx(GLOBAL_2D_ARRAY) )
      brk2dx=0.
      allocate( brk2de(GLOBAL_2D_ARRAY) )
      brk2de=0.
      allocate( ust2d(GLOBAL_2D_ARRAY) )
      ust2d=0.
      allocate( vst2d(GLOBAL_2D_ARRAY) )
      vst2d=0.
      allocate( ustr2d(GLOBAL_2D_ARRAY) ) ! At rho point
      ustr2d=0.
      allocate( vstr2d(GLOBAL_2D_ARRAY) ) ! At rho point
      vstr2d=0.
      allocate( ustr0(GLOBAL_2D_ARRAY) )
      ustr0=0.
      allocate( vstr0(GLOBAL_2D_ARRAY) )
      vstr0=0.
      allocate( llmm(GLOBAL_2D_ARRAY) )
      llmm=0.
      allocate( frc2dx(GLOBAL_2D_ARRAY) ) ! used in prestep, step2d, step3d
      frc2dx=0.
      allocate( frc2de(GLOBAL_2D_ARRAY) ) ! used in prestep, step2d, step3d
      frc2de=0.
      allocate( sup(GLOBAL_2D_ARRAY) ) ! used in step2d
      sup=0.
      allocate( sup_zeta(GLOBAL_2D_ARRAY) ) ! used in step2d
      sup_zeta=0.


# ifdef SOLVE3D
      allocate( calP(GLOBAL_2D_ARRAY) ) ! used in prsgrd.F
      calP=0.
      allocate( Kapsrf(GLOBAL_2D_ARRAY) ) ! used in prsgrd.F
      Kapsrf=0.
      allocate( kD(GLOBAL_2D_ARRAY) )
      kD=0.
      allocate( kw(GLOBAL_2D_ARRAY) )
      kw=0.
      allocate( Dstp(GLOBAL_2D_ARRAY) )
      Dstp=0.
      allocate( actf(GLOBAL_2D_ARRAY) )
      actf=0.
      allocate( act(GLOBAL_2D_ARRAY) )
      act=0.
# if !defined MONOWAVE
      allocate( keff(GLOBAL_2D_ARRAY) )
      keff=0.0
# endif
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      allocate( brk_r(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
      brk_r=0.
      allocate( brk3dx(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
      brk3dx=0.
      allocate( brk3de(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
      brk3de=0.
#  endif
#  ifdef BODY_FRICTION
      allocate( frc_r(GLOBAL_2D_ARRAY,N) )
      frc_r=0.
      allocate( frc3dx(GLOBAL_2D_ARRAY,N) )
      frc3dx=0.
      allocate( frc3de(GLOBAL_2D_ARRAY,N) )
      frc3de=0.
#  endif
      allocate( ust(GLOBAL_2D_ARRAY,N) )
      ust=0.
      allocate( vst(GLOBAL_2D_ARRAY,N) )
      vst=0.
      allocate( wst(GLOBAL_2D_ARRAY,N) ) ! used in set_avg.F
      wst=0.
      allocate( kvf(GLOBAL_2D_ARRAY,N) ) ! used in prsgrd32AC1.F
      kvf=0.
      allocate( Akb(GLOBAL_2D_ARRAY,0:N) ) ! used in lmd_kpp.F
      Akb=0.
      allocate( Akw(GLOBAL_2D_ARRAY,0:N) ) ! used in lmd_kpp.F
      Akw=0.
      allocate( E_pre(GLOBAL_2D_ARRAY,0:N) )
      E_pre=0.
# endif  /* SOLVE3D */

      ! REALISTIC WEC VARIABLES (FROM READ IN NETCDF FILE)
      ! =======================

!#if defined BBL || defined WEC
# if !defined ANA_WEC_FRC && !defined WKB_WWAVE

      allocate( wwag(GLOBAL_2D_ARRAY) )
      wwag=0.
      allocate( wwdg(GLOBAL_2D_ARRAY) )
      wwdg=0.
      allocate( wwpg(GLOBAL_2D_ARRAY) )
      wwpg=0.
#  ifdef UORBOFF
      allocate( wwuob(GLOBAL_2D_ARRAY) )
      wwuob=0.
      allocate( wwvob(GLOBAL_2D_ARRAY) )
      wwvob=0.
#  endif
      allocate( wwdrx(GLOBAL_2D_ARRAY) )
      wwdrx=0.
      allocate( wwdre(GLOBAL_2D_ARRAY) )
      wwdre=0.
      allocate( wweb(GLOBAL_2D_ARRAY) )
      wweb=0.
      allocate( wveb(GLOBAL_2D_ARRAY) )
      wveb=0.
      allocate( wved(GLOBAL_2D_ARRAY) )
      wved=0.
      allocate( wwed(GLOBAL_2D_ARRAY) )
      wwed=0.
      allocate( wwUSg(GLOBAL_2D_ARRAY) )
      wwUSg=0.
      allocate( wwVSg(GLOBAL_2D_ARRAY) )
      wwVSg=0.
      allocate( wwus0g(GLOBAL_2D_ARRAY) )
      wwus0g=0.
      allocate( wwvs0g(GLOBAL_2D_ARRAY) )
      wwvs0g=0.
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      allocate( wvqb(GLOBAL_2D_ARRAY) )
      wvqb=0.
      allocate( wwqb(GLOBAL_2D_ARRAY) )
      wwqb=0.
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      allocate( wwsup(GLOBAL_2D_ARRAY) )
      wwsup=0.
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      allocate( wwlm(GLOBAL_2D_ARRAY) )
      wwlm=0.
#   endif

      ! Grid rotation angles
      allocate( angler_cos(GLOBAL_2D_ARRAY) )
      angler_cos=0.
      allocate( angler_sin(GLOBAL_2D_ARRAY) )
      angler_sin=0.

# elif defined ANA_WEC_FRC && defined WEC
     !Daniel add, we need more genral cpp since this is copied from
      ! above
      allocate( wwag(GLOBAL_2D_ARRAY,2) )
      wwag=0.
      allocate( wwdg(GLOBAL_2D_ARRAY,2) )
      wwdg=0.
      allocate( wwpg(GLOBAL_2D_ARRAY,2) )
      wwpg=0.
#  ifdef UORBOFF
      allocate( wwuob(GLOBAL_2D_ARRAY,2) )
      wwuob=0.
      allocate( wwvob(GLOBAL_2D_ARRAY,2) )
      wwvob=0.
#  endif
      allocate( wwdrx(GLOBAL_2D_ARRAY,2) )
      wwdrx=0.
      allocate( wwdre(GLOBAL_2D_ARRAY,2) )
      wwdre=0.
      allocate( wweb(GLOBAL_2D_ARRAY,2) )
      wweb=0.
      allocate( wveb(GLOBAL_2D_ARRAY) )
      wveb=0.
      allocate( wved(GLOBAL_2D_ARRAY) )
      wved=0.
      allocate( wwed(GLOBAL_2D_ARRAY,2) )
      wwed=0.
      allocate( wwUSg(GLOBAL_2D_ARRAY,2) )
      wwUSg=0.
      allocate( wwVSg(GLOBAL_2D_ARRAY,2) )
      wwVSg=0.
      allocate( wwus0g(GLOBAL_2D_ARRAY,2) )
      wwus0g=0.
      allocate( wwvs0g(GLOBAL_2D_ARRAY,2) )
      wwvs0g=0.
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      allocate( wvqb(GLOBAL_2D_ARRAY) )
      wvqb=0.
      allocate( wwqb(GLOBAL_2D_ARRAY,2) )
      wwqb=0.
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      allocate( wwsup(GLOBAL_2D_ARRAY,2) )
      wwsup=0.
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      allocate( wwlm(GLOBAL_2D_ARRAY,2) )
      wwlm=0.
#   endif

# endif /* !ANA_WEC_FRC && !WKB_WWAVE */
!#endif /* BBL || WEC */
#ifdef SOLVE3D
      allocate( Wrk(0:N))
      Wrk=0.
#endif
      !For compute_wec_rhs_uv_terms
      allocate( utmp(GLOBAL_2D_ARRAY))
      utmp=0.
      allocate( vtmp(GLOBAL_2D_ARRAY))
      vtmp=0.

      !#### end variables added as per forces.h #####

      ! AVERAGED VARIABLES
      ! ==================

      !ADD IF/DEF AROUND THIS      
      allocate( sup_avg(GLOBAL_2D_ARRAY) )
      sup_avg=0.
      allocate( ust2d_avg(GLOBAL_2D_ARRAY) )
      ust2d_avg=0.
      allocate( vst2d_avg(GLOBAL_2D_ARRAY) )
      vst2d_avg=0.
      allocate( ust_avg(GLOBAL_2D_ARRAY,N) )
      ust_avg=0.
      allocate( vst_avg(GLOBAL_2D_ARRAY,N) )
      vst_avg=0.
      allocate( wst_avg(GLOBAL_2D_ARRAY,N) )
      wst_avg=0.
      allocate( akb_avg(GLOBAL_2D_ARRAY,0:N) )
      akb_avg=0.
      allocate( akw_avg(GLOBAL_2D_ARRAY,0:N) )
      akw_avg=0.
      allocate( kvf_avg(GLOBAL_2D_ARRAY,N) )
      kvf_avg=0.
      allocate( calp_avg(GLOBAL_2D_ARRAY) )
      calp_avg=0.
      allocate( kapsrf_avg(GLOBAL_2D_ARRAY) )
      kapsrf_avg=0.
      allocate( sup_zeta_avg(GLOBAL_2D_ARRAY) )
      sup_zeta_avg=0.
      allocate( frc2dx_avg(GLOBAL_2D_ARRAY) )
      frc2dx_avg=0.
      allocate( frc2de_avg(GLOBAL_2D_ARRAY) )
      frc2de_avg=0.
      allocate( brk2dx_avg(GLOBAL_2D_ARRAY) )
      brk2dx_avg=0.
      allocate( brk2de_avg(GLOBAL_2D_ARRAY) )
      brk2de_avg=0.
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      allocate( brk3dx_avg(GLOBAL_2D_ARRAY,N) )
      brk3dx_avg=0.
      allocate( brk3de_avg(GLOBAL_2D_ARRAY,N) )
      brk3de_avg=0.
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
      allocate( frc3dx_avg(GLOBAL_2D_ARRAY,N) )
      frc3dx_avg=0.
      allocate( frc3de_avg(GLOBAL_2D_ARRAY,N) )
      frc3de_avg=0.
#  endif
      allocate( wveb_avg(GLOBAL_2D_ARRAY) )
      wveb_avg=0.
      allocate( wved_avg(GLOBAL_2D_ARRAY) )
      wved_avg=0.
      allocate( hvf_int_xi_avg(GLOBAL_2D_ARRAY) )
      hvf_int_xi_avg=0.
      allocate( hvf_int_eta_avg(GLOBAL_2D_ARRAY) )
      hvf_int_eta_avg=0.


      end subroutine init_arrays_wec_tile
      !---------------------------------------------------------------
      !---------------------------------------------------------------
  
#ifdef SOLVE3D
      subroutine wstokes_wec (tile)

      ! ===========================================
      ! Vertical Stokes drift velocity at rho-point
      ! ===========================================

      use hidden_mpi_vars
      use private_scratch

      implicit none
! # include "param" ! already included at top of module
      integer tile

# include "compute_tile_bounds.h"
      call wstokes_wec_tile(istr,iend,jstr,jend)

      end subroutine wstokes_wec
!-----------------------------------------------------------------------

      subroutine wstokes_wec_tile (istr,iend,jstr,jend)

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      ! define Wrk at top of module
      !real Wrk(0:N) 
      do j=jstr,jend
        do i=istr,iend
          Wrk(0)=0.D0
          do k=1,N,+1        !--> recursive
            Wrk(k) = Wrk(k-1) -Hz(i,j,k)
     &       *( pm(i,j)*(ust(i+1,j,k)-ust(i,j,k))
     &         +pn(i,j)*(vst(i,j+1,k)-vst(i,j,k)) )
          enddo

          wst(i,j,N)=+0.375*Wrk(N) +0.75*Wrk(N-1) -0.125*Wrk(N-2)
          do k=N-1,2,-1
            wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
          enddo
          wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
        enddo
      enddo
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then                       ! Set lateral
        do k=1,N                                   ! boundary
          do j=jstr,jend                           ! conditions
            wst(istr-1,j,k)=wst(istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            wst(iend+1,j,k)=wst(iend,j,k)
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jstr-1,k)=wst(i,jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jend+1,k)=wst(i,jend,k)
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jend+1,k)=wst(istr,jend,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(iend+1,jend+1,k)=wst(iend,jend,k)
        enddo
      endif
#  endif
# endif

      end subroutine wstokes_wec_tile
#endif  /* SOLVE3D */
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      subroutine compute_wec_rhs_uv_terms (istr,iend,jstr,jend,ru,rv)
!     &                                           wrkone, wrktwo)

      ! ========================================
      ! Add some of the WEC variables to the RHS
      ! of the equation system
      ! ======================
      
      ! DPD removed VF_ROBUST code since this is not actually 'robust'
      ! D.H. tested this in ~2020 and found that it led to instabilities 
      implicit none
      
      !Input/Output
      integer, intent(in) :: istr, iend, jstr, jend
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N), intent(inout) :: ru,rv
      !Local
      integer :: i, j, k
      real :: cff, cff1
# include "compute_auxiliary_bounds.h"
      

!      print *, 'IN compute_wec_rhs_uv_terms'
! ------------------------------------------------------------------------
!      Add horizontal vortex-force and advection terms & breaking terms
! ------------------------------------------------------------------------


! ======================================================================
!      xi-component horiz. vortex force is computed as follows:
!      hvfx  = ust*du/dx +  vst * dv/dx
!            = hvfx_u    +  hvfx_v
!            = hvfxu^{i-1/2,j,k} + hvfxv^{i-1/2,j,k}  
!            = cff1^{i-1/2,j,k}
!           
!            where i-1/2 is a u-pt corresponding to index=i in the loop
!            below for u-pt variables

!      eta-component vortex force is computed as follows
!      hvfy =  ust*du/dy +   vst * dv/dy
!            = hvfy_u    +   hvfy_v
!            = hvfyu^{i,j-1/2,k} + hvfyv^{i,j-1/2,k}  
!            = cff1^{i,j-1/2,,k}
!           
!            where j-1/2 is a v-pt corresponding to index=j in the loop
!            below for v-pt variables
! ======================================================================

      ! Reset depth-integrated vortex force accumulators.
      ! This subroutine is called twice per timestep (predictor+corrector);
      ! we keep only the last call's values (matching diagnostics timing).
      do j=jstr,jend
        do i=istrU,iend
          hvf_int_xi(i,j) = 0.
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          hvf_int_eta(i,j) = 0.
        enddo
      enddo

      do k=1,N                       ! <-- start k-loop
        do j=jstr,jend
          do i=istrU,iend

            ! Volume normalization dy * Hz at (i-1/2,j,k)
            cff  =0.5*dn_u(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
#ifdef UV_ADV
            cff1 = 0.25*(
               ! ust * du/dx at (i,j)
     &         (ust(i+1,j,k)+ust(i,j,k))*(u(i+1,j,k,nrhs)-u(i,j,k,nrhs))
               ! ust * du/dx at (i-1,j)
     &        +(ust(i,j,k)+ust(i-1,j,k))*(u(i,j,k,nrhs)-u(i-1,j,k,nrhs))
               ! vst * dv/dx at (i-1/2,j-1/2)
     &        +(vst(i,j,k)+vst(i-1,j,k))*(v(i,j,k,nrhs)-v(i-1,j,k,nrhs))
               ! vst * dvdx at (i-1/2,j+1/2)
     &        +(vst(i,j+1,k)+vst(i-1,j+1,k))
     &                              *(v(i,j+1,k,nrhs)-v(i-1,j+1,k,nrhs))

     &                                                                 )
#  else
           cff1 =0.D0
#  endif
            ! Accumulate depth-integrated horizontal vortex force (xi)
            hvf_int_xi(i,j) = hvf_int_xi(i,j) + cff*cff1

            ru(i,j,k)=ru(i,j,k) + cff*( cff1
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
                                                !3D breaking body-force
     &                                          +dm_u(i,j)*brk3dx(i,j,k)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
                                                !3D friction body-force
     &                                          +dm_u(i,j)*frc3dx(i,j,k)
#  endif
     &                                                                 )
          enddo
        enddo

        do j=jstrV,jend
          do i=istr,iend
            ! Volume normalization dx * Hz at (i,j-1/2,k)
            cff  =0.5*dm_v(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
#ifdef UV_ADV
            cff1 = 0.25*(
               ! ust * du/dy at (i-1/2,j-1/2)
     &         (ust(i,j,k)+ust(i,j-1,k))*(u(i,j,k,nrhs)-u(i,j-1,k,nrhs))
               ! ust * du/dy at (i+1/2,j-1/2)
     &        +(ust(i+1,j,k)+ust(i+1,j-1,k))
     &                              *(u(i+1,j,k,nrhs)-u(i+1,j-1,k,nrhs))
               ! vst * dv/dy at (i,j)
     &        +(vst(i,j+1,k)+vst(i,j,k))*(v(i,j+1,k,nrhs)-v(i,j,k,nrhs))
               ! vst * dv/dy at (i,j-1)
     &        +(vst(i,j,k)+vst(i,j-1,k))*(v(i,j,k,nrhs)-v(i,j-1,k,nrhs))
     &                                                                 )
# else
             cff1=0.D0
#endif
            ! Accumulate depth-integrated horizontal vortex force (eta)
            hvf_int_eta(i,j) = hvf_int_eta(i,j) + cff*cff1

            rv(i,j,k)=rv(i,j,k) + cff*( cff1
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
                                                !3D breaking body-force
     &                                          +dn_v(i,j)*brk3de(i,j,k)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
                                                !3D friction body-force
     &                                          +dn_v(i,j)*frc3de(i,j,k)
#  endif
     &                                                                 )
          enddo
        enddo
      enddo                ! <-- end k-loop

      ! Normalize hvf_int from ru-units [m^4/s^2] to diagnostic units
      ! [m^2/s^2] by dividing by cell area dm*dn = 1/(pm*pn).
      do j=jstr,jend
        do i=istrU,iend
          hvf_int_xi(i,j) = hvf_int_xi(i,j)
     &                    / (dm_u(i,j)*dn_u(i,j))
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          hvf_int_eta(i,j) = hvf_int_eta(i,j)
     &                     / (dm_v(i,j)*dn_v(i,j))
        enddo
      enddo

      end subroutine compute_wec_rhs_uv_terms

!-----------------------------------------------------------------------

      subroutine set_wec_frc![
      ! READ IN WEC FORCING
      ! modeled after set_bgc_surf_frc in bgc.F
      ! in updated repo as of June 2024


      implicit none
   

      ! IF ANALYTICAL FORCING
      ! =====================
# ifdef ANA_WEC_FRC
      ! NO ISTR,JSTR,IEND,JEND --> RE-MAKE THIS
      !call ana_wec_frc_tile(istr,iend,jstr,jend)
      ierr=-99
# else
 
      ! ===============================
      !         READ IN WW3 FORCING
      ! ================================
      !print *, 'SET WEC FRC HELLO' 
      ! wave amplitude [m]
      call set_frc_data(nc_wwag,wwag,'r')     
      !print *, 'success wwag!' 
      ! wave direction [degrees from true north] 
      call set_frc_data(nc_wwdg, wwdg, 'r')
      
      ! wave period [sec] 
      call set_frc_data(nc_wwpg, wwpg, 'r')     

      ! breaking dissipation [m^3/s^3] 
      call set_frc_data(nc_wweb, wweb, 'r')

      ! frictional dissipation [m^3 / s^3]
      call set_frc_data(nc_wwed, wwed, 'r')

# if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      ! fraction of breaking waves [N.D.]      
      call set_frc_data(nc_wwqb, wwqb, 'r')
# endif

# if defined WAVE_OFFLINE && defined UORBOFF
      ! u-orbital velocity [m/s]
      ! --> U-east (not xi!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwuob, wwuob, 'r')

      ! v-orbital velocity [m/s]     
      ! --> V-north (not eta!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwvob, wwvob, 'r')
#endif


#ifndef MONOWAVE
      ! FOR BROADBAND FORCING & LEONEL'S METHOD 

      ! depth-averaged u stokes drift [m/s] 
      ! --> U-east (not xi!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwUSg, wwUSg, 'r')

      ! depth-averaged v stokes drift [m/s]
      ! --> V-north (not eta!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwVSg, wwVSg, 'r')

      ! surface u stokes drift [m/s]
      ! --> U-east (not xi!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwus0g, wwus0g, 'r')

      ! surface v stokes drift [m/s]
      ! --> V-north (not eta!!)  converted in conv_wec_frc
      call set_frc_data(nc_wwvs0g, wwvs0g, 'r')

#endif

# if defined WAVE_OFFLINE && defined SUP_OFF
      ! wave set-down [m]
      call set_frc_data(nc_wwsup, wwsup, 'r')
#endif

# if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM

      ! mean wavelenght [m]
      call set_frc_data(nc_wwlmv, wwlm, 'r')
# endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MAKE CONVERSIONS AFTER READING IN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call conv_wec_frc

      !print *, 'SUCCESS conv_wec_frc'
#endif /* NOT DEF ANA_WEC_FRC */
       
      end subroutine set_wec_frc !]
      !-------------------------------------------

      subroutine conv_wec_frc
      ! Convert wave forcing data read in in set_wec_frc
      ! to appropriate variables / units conventions for WEC in ROMS
      ! this function is based off of set_wec_inp2vars in old version of
      ! code (pre-2024, pre-github repo)
  
      implicit none
      !integer :: istr,iend,jstr,jend
      integer ::  i,j
      real :: wdir=0. !interim variable to transform wave direction

!      integer :: tile=0
!#include "compute_tile_bounds.h"
!#include "compute_auxiliary_bounds.h"
!#include "compute_extended_bounds.h"


      ! CHANGE TO nx,ny convention once we are sure of how that works
      do j=0,ny+1
         do i=0,nx+1
       !do j=jstrR,jendR
       !   do i=istrR,iendR
      !do j=jstr-1,jend+1
      !   do i=istr-1,iend+1
          ! ------- Wave amplitude -----------! 
          ! Wave amplitude is read in, so convert to hrms for later
          ! hrms = 2 x amplitude
          ! keep same variable for efficiency (sasha-style)
          wwag(i,j) = 2 * wwag(i,j)
#ifdef MASKING
     &                   * rmask(i,j)
#endif
          !What is used in wec_forces
          whrm(i,j) = wwag(i,j)
          ! ------- Wave direction -----------! 
          !Calculate wave direction at i,j just for interim
          ! wdir --> angle from N to travel to 0 deg
          wdir = -(deg2rad*wwdg(i,j) - 1.5*pi)
#ifdef CURVGRID
          ! CODE IN CHECK FOR angler??
     &                     - angler(i,j) !angler MUST BE IN RADIANS!
#endif
          ! Direction (wwdg) is read in as one variable so can't be
          ! reused for dir-x and dir-y, need new interim variables
          !wwdrx(i,j) = cos(wdir)
          wdrx(i,j) = cos(wdir)
# ifdef MASKING
     &                             *rmask(i,j)
# endif
          wdre(i,j) = sin(wdir)
          !wwdre(i,j) = sin(wdir)
# ifdef MASKING
     &                             *rmask(i,j)
# endif

          ! ------- Wave period -----------! 
          ! conver period to frequency
          wwpg(i,j) = 2.0*pi / max(wwpg(i,j),0.1D0) 
# ifdef MASKING
     &                             *rmask(i,j)
# endif
          wfrq(i,j) = wwpg(i,j)

                    ! Set grid rotation angles (if first timestep)
          ! ------------------------
          ! Need once as doesn't change.

          !     DANIEL NOTE --> IF ONLY USED HERE, JUST USE cos,
          !     sin...why make up 2 2D arrays for this?
          if(FIRST_TIME_STEP) then
            angler_cos(i,j) = cos(angler(i,j))
            angler_sin(i,j) = sin(angler(i,j))
          endif
          
#ifndef MONOWAVE        
      ! ------- wave model read in stokes drifts ---!  
      ! use grid rotation angles to convert wwUSg, wwVSg, wwus0g,
      ! wwvs0g to xi-, eta-dir (from east, north dirs)

      ! FOR BROADBAND FORCING & LEONEL'S METHOD 

      !depth-averaged xi-dir ust
          ustr2d(i,j) = wwUSg(i,j) * angler_cos(i,j) 
     &                + wwVSg(i,j) * angler_sin(i,j)
      !depth-averaged eta-dir vst
          vstr2d(i,j) = -wwUSg(i,j) * angler_sin(i,j) 
     &                + wwVSg(i,j) * angler_cos(i,j)
      !surface  xi-dir ust
          ustr0(i,j) = wwus0g(i,j) * angler_cos(i,j) 
     &                + wwvs0g(i,j) * angler_sin(i,j)
      !surface eta-dir vst
          vstr0(i,j) = -wwus0g(i,j) * angler_sin(i,j) 
     &                + wwvs0g(i,j) * angler_cos(i,j)
#endif

#ifdef UORBOFF
      ! ------- wave orbital velocities ---!  

        uorb(i,j) = wwuob(i,j) * angler_cos(i,j) 
     &             + wwvob(i,j) * angler_sin(i,j)

        vorb(i,j) = -wwuob(i,j) * angler_sin(i,j) 
     &             + wwvob(i,j) * angler_cos(i,j)

#endif




          !---- remaining variables just get masked ---! 
# ifdef MASKING /* remaining variables not converted just rmask'ed */

# ifndef MONOWAVE
          ustr2d(i,j) = ustr2d(i,j)*rmask(i,j)
          vstr2d(i,j) = vstr2d(i,j)*rmask(i,j)
          ustr0(i,j) = ustr0(i,j)*rmask(i,j)
          vstr0(i,j) = vstr0(i,j)*rmask(i,j)
# endif
          !going from ww_ vars to wv_ vars (dumb, but keep it for now)
          wveb(i,j) = wweb(i,j)*rmask(i,j)
!#   if defined WAVE_OFFLINE && defined WAVE_FRICTION
          wved(i,j) = wwed(i,j)*rmask(i,j)
!#   endif
# ifndef MONOWAVE
          sup(i,j) = wwsup(i,j)*rmask(i,j)
#  ifdef WAVE_DIFF_FROM_LM
          llmm(i,j) = wwlm(i,j) * rmask(i,j)
#  endif
# endif
# ifdef UORBOFF
          uorb(i,j) = uorb(i,j)*rmask(i,j)
          vorb(i,j) = vorb(i,j)*rmask(i,j)
# endif

# endif

         enddo
      enddo




      end subroutine conv_wec_frc



!-------------------------------------------------------------------------

       ! THIS FUNCTION MAY BE OUTDATED?
      subroutine wec_forces (tile)
     ! This is mrl_wci (tile) from mrl_wci.F in old code.
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer  tile
!# include "param" ! DevinD - already included in top of module
# include "compute_tile_bounds.h"
      call wec_forces_tile (istr,iend,jstr,jend)
      call wstokes_wec_tile(istr,iend,jstr,jend)
      !print *, 'SUCCESS SUCCESS SUCCESS wec_forces'
      end subroutine wec_forces
! ----------------------------------------------------------------------

      subroutine wec_forces_tile (istr,iend,jstr,jend)

      ! NOTE: input wave data is at rho-points, thus conversion from
      ! whrm at rho-points, is needed to get ust2d and u-points, etc.

      ! This is taken from mrl_wci.F of old code.

!   Evaluating wave-averaged terms and Stokes drift based on MRL04
!   (see: McWilliams, Restrepo & Lane, 2004, JFM, 511, pp.135-178)
!   inputs: R.M.S. wave height (m);
!           peak wave freqiuency (rad/s);
!           mean wave direction: wdrx & wdre (non dimensional);
!           breaking dissipation   (\epsilon_b / \rho, m3/s3);
!           roller dissipation     (\epsilon_r / \rho, m3/s3); and
!           frictional dissipation (\epsilon_d / \rho, m3/s3)
!
      use coupling

      implicit none
      !Inputs/outputs
      integer, intent(in) :: istr,iend,jstr,jend
      !Local
      integer :: i,j,imin,jmin
      real ::   cff,cff1,cff2,cff3,cff4,inv_g,khd,kh,umag,
!    &   cff6,cff7,
     &   eps,wave_ramp,wramp2,kbrk,fb,fb0,fb1,fb2,inv_fbs,
     &   bconst,ka_f00,hz0,z_tide,
     !&   , Uorb1,tauc,tauw --> Uorb1 = public Uorb_mag, tauc, tauw in
     !pre_step3d for bottom drag calculation
     &   dtinv,inv_zb,a_brk,a_kv,c1o3,c4o3,inv_k,dd,intfb,!ust_ker,intfb,

     &   fbsrf,fn1,fn2,khmax,bz1,bz2,beta,abot,delta,vst_ker,
     &   a_frc,kfrc
# if defined WAVE_OFFLINE && !defined MONOWAVE 
     &         ,usec,vsec,USEZ,VSEZ,ustm0w
     &         ,FAC,FACo,COUNTR, USTERR, ustm0
     &         ,ust_ker
# endif
     &    ,L,hs,zz,zzm,zzb,hsfs
     &         ,wrk4(0:N)
# if defined SUP_OFF
     &         ,supm, Fsup
# endif
!# if defined WAVE_DIFF_FROM_LM
     &         ,keffkw
!# endif
# ifdef SOLVE3D
      integer ::  k,kk
# endif
      parameter ( eps=1.e-10,
     &            bconst=0.1D0, !0.03D0,               ! breaking contribution to KPP
     &            a_kv=1.2D0, !1.2D0                 ! breaking scale for eddy visc.
     &            a_brk=0.2D0,                ! breaking scale for body force
     &            a_frc=3.D0,                 ! friction scale for body force
     &            khmax=20.D0,                ! deep-water limit for k x dep
     &            c1o3=0.3333333333333333,    ! 1/3
     &            c4o3=1.3333333333333333     ! 4/3
     &           )
# ifdef WAVE_OFFLINE
      real roller_fraction
      parameter ( roller_fraction=0.0)
# endif
# include "compute_auxiliary_bounds.h"


!====================================================
!               October 2022 Update - Daniel D.
! calP, Kapsrf, and kvf now all calculated in compute_prsgrd_wec_terms.h
! which is called from prsgrd.F  this is done so that these terms are
! added consistently to the predictor and corrector step using index
! nrhs for Eulerian u,v.

! Previously this was done here with index nstp, which does not equal
! nrhs during the corrector stage (see set_global_definitions.h)
!===================================================



!=========================================================================
! explicit wavenumber estimator

! initial ramping coefficients

      wave_ramp = 1.0
      ! get rid of wramp2 ??
      wramp2= wave_ramp**2
!
! Evaluate 2DH wave-current intreraction variable.
! ================================================
!
! 2DH Stokes velocities, breaking, roller and bottom-friction
! dissipation terms, defined at horizontal rho-points.
!
     ! if (mynode.eq.0) then
     !    print *, 'knew=', knew
     ! endif
      inv_g = 1.0/g
      do j=0,ny+1
        do i=0,nx+1
      !do j=jstr-1,jend+1
      !  do i=istr-1,iend+1
          Dstp(i,j)=max(h(i,j)+zeta(i,j,knew),eps)
          inv_d(i,j)=1.0/Dstp(i,j)

          fr(i,j) =wfrq(i,j)
          !if (mynode.eq.0) then
          !   if (i.eq.10) then
          !      if (j.eq.10) then
          !         print *, 'whrm(i,j)=', whrm(i,j)
          !      endif
          !   endif
          !endif
          wh(i,j) =wave_ramp*whrm(i,j)
          inv_f(i,j) =1.0/max(fr(i,j),eps)
          khd = Dstp(i,j)*(fr(i,j)**2)*inv_g

          kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &               +khd*(0.3555555555 + khd*(0.1608465608
     &               +khd*(0.0632098765 + khd*(0.0217540484
     &                            +khd*0.0065407983)))))) )

          kD(i,j) =kh
          kw(i,j) =kh*inv_d(i,j)
#ifdef BREAKING
          ebrk(i,j) =wramp2*wveb(i,j)            ! epsilon_b (m3/s3)
#else
          ebrk(i,j) =0.0
#endif
          !if (mynode==0 .and. j==1 .and. i==1) then
          !     print*,'ebrk(i,j)=',ebrk(i,j)
          !endif
          cff =ebrk(i,j)

#ifdef WAVE_FRICTION
          frc(i,j) =wramp2*wved(i,j)*kw(i,j)
     &                           *inv_f(i,j)     ! 2DH friction term
#else
          frc(i,j) =0.0
#endif
#ifdef BREAKING
          brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
#else
          brk(i,j) =0.0
#endif
          act(i,j) =inv_f(i,j)*(
     &                      0.125*g*(wh(i,j)**2) ) ! wave action density

          stk(i,j) =act(i,j)*kw(i,j)*inv_d(i,j)  ! depth-averaged Stokes
# ifdef MASKING
          kw(i,j) = kw(i,j)*rmask(i,j)
          kD(i,j) = kD(i,j)*rmask(i,j)
          act(i,j) =act(i,j)*rmask(i,j)
          stk(i,j) =stk(i,j)*rmask(i,j)
          brk(i,j) =brk(i,j)*rmask(i,j)
          frc(i,j) =frc(i,j)*rmask(i,j)
# endif   /* apply land mask */
        enddo
      enddo              ! <-- discard inv_d (ifndef BBL_F00)


# if defined BREAKING && defined BRK0
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          brk(i,j)=0.D0
        enddo
      enddo
# endif

# if defined WAVE_FRICTION && defined FRC0
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          frc(i,j)=0.D0
        enddo
      enddo
# endif

!
! sup : quasi-static sea-level response, set-up. no interaction.
! ==============================================================
!
# if !defined SUP_OFF
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstrR,jendR
      !  do i=istrR,iendR
          sup(i,j) =-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
#  ifdef MASKING
     &                                        *rmask(i,j)
#  endif
        sup_zeta(i,j) = sup(i,j) + zeta(i,j,knew)
        enddo
      enddo
# endif

# if defined SUP_OFF
      actf(:,:)=0.
! determine the effective filtered action density g a^2/w
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstrR,jendR
      !  do i=istrR,iendR
          supm=-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kw(i,j)*h(i,j),khmax)),eps)
! re calculate sup with \hat{\eta}
          if (kD(i,j).lt.6.28 .and.supm.ne.0..and. sup(i,j)*supm .gt. 0.0) then
          Fsup=min(1.,sqrt(sup(i,j)/supm))
          else
          Fsup=1.
          endif
          actf(i,j)=act(i,j)*Fsup**2
        enddo
      enddo
# endif
!
!
! 2DH depth-averaged Stokes drift and breaking acceleration
! =========================================================
!
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstrR,jendR
      !  do i=istr,iendR
# if defined MONOWAVE
          !Monochromatic
          ust2d(i,j)  =0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
# else 
          !Read in from WW3
          ust2d(i,j)=0.5*(ustr2d(i-1,j)+ustr2d(i,j))

# endif
          brk2dx(i,j) = 0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
          frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
#  ifdef MASKING
          brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
#  endif
#  if defined BREAKING && defined SURFACE_BREAK && !defined BODY_BREAK
          !ADD brk2dx to sustr (add breaking stress to global surface
          !(wind) stress, frc2d (bottom streaming) added explicitly in
          !step3d
          sustr(i,j) = sustr(i,j) + brk2dx(i,j)
#  endif
        enddo
      enddo
      
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstr,jendR
      !  do i=istrR,iendR
# if defined MONOWAVE 
          !Monochromatic
          vst2d(i,j)  =0.5*(stk(i,j-1)*wdre(i,j-1)+stk(i,j)*wdre(i,j))
# else
          !Read in from WW3
          vst2d(i,j)=0.5*(vstr2d(i,j-1)+vstr2d(i,j))

# endif
          brk2de(i,j) = 0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
          frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
#  ifdef MASKING
          brk2de(i,j) =brk2de(i,j)*vmask(i,j)
          frc2de(i,j) =frc2de(i,j)*vmask(i,j)
#  endif
#  if defined BREAKING && defined SURFACE_BREAK && !defined BODY_BREAK
          !ADD brk2de to svstr (add breaking stress to global surface
          !(wind) stress, frc2d (bottom streaming) added explicitly in
          !step3d
          svstr(i,j) = svstr(i,j) + brk2de(i,j)
#  endif
 
        enddo
      enddo    

!
! Combined wave-current bottom drag coefficient r_D at rho-points
! ===============================================================
!
# if defined BBL_S95
      !inv_zb=1.0/max(Zob,eps)
# endif
      do j=1,ny
         do i=1,nx
      !do j=jstr,jend
      !  do i=istr,iend
# if defined BBL_S95
!          cff1 = 0.5*(u(i,j,1,nrhs)+u(i+1,j,1,nrhs))
!          cff2 = 0.5*(v(i,j,1,nrhs)+v(i,j+1,1,nrhs))
!          umag = sqrt( cff1*cff1+cff2*cff2 )

#  ifndef UORBOFF
          Uorb_mag(i,j) = 0.5*fr(i,j)*wh(i,j)
     &              /max(sinh(min(kD(i,j),khmax)),eps)
#  else
          Uorb_mag(i,j) = sqrt(uorb(i,j)**2+vorb(i,j)**2)
#  endif
          !Define public tauw_orb for use in compute_rd_bott_drag.h
          tauw_orb(i,j)=0.695*(Uorb_mag(i,j)**1.48)*((Zob*fr(i,j))**0.52)

!#    ifdef SOLVE3D
!          cff1=umag*(vonKar/log(1.+0.5*Hz(i,j,1)*inv_zb))**2
!#    else
!          cff1=umag*(vonKar/log(1.+Dstp(i,j)*inv_zb))**2
!#    endif

!          tauc=cff1*umag
!          tauw=0.695*(Uorb1**1.48)*((Zob*fr(i,j))**0.52)
!          cff2=1.0 + 1.2*((tauw/max(tauw+tauc,eps))**3.2)
!          r_D(i,j)=cff1*cff2
!#    ifdef MASKING
!          r_D(i,j)=r_D(i,j)*rmask(i,j)
!#     endif


!       DPD COMMENTED THIS OUT, r_D shoudl just be regular old r_D
!# else /* BBL_S95 */
!           r_D(i,j)=0.0

# endif /* BBL_S95 */
        enddo
      enddo

      ! DPD also commented this out, again wec_frc should not do
      ! anything to r_D, this is done compute_rd_bott_drag.h
!#  ifndef EW_PERIODIC
!      if (WESTERN_EDGE) then
!        do j=jstr,jend
!          r_D(istr-1,j)=r_D(istr,j)
!        enddo
!      endif
!      if (EASTERN_EDGE) then
!        do j=jstr,jend
!          r_D(iend+1,j)=r_D(iend,j)
!        enddo
!      endif
!#  endif
!#  ifndef NS_PERIODIC
!      if (SOUTHERN_EDGE) then
!        do i=istr,iend
!          r_D(i,jstr-1)=r_D(i,jstr)
!        enddo
!      endif
!      if (NORTHERN_EDGE) then
!        do i=istr,iend
!          r_D(i,jend+1)=r_D(i,jend)
!        enddo
!      endif
!#   ifndef EW_PERIODIC
!      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
!        r_D(istr-1,jstr-1)=r_D(istr,jstr)
!      endif
!      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
!        r_D(istr-1,jend+1)=r_D(istr,jend)
!      endif
!      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
!        r_D(iend+1,jstr-1)=r_D(iend,jstr)
!      endif
!      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
!        r_D(iend+1,jend+1)=r_D(iend,jend)
!      endif
!#   endif
!#   endif












# ifdef SOLVE3D
!
! Conservative 3D wave-current interaction variables.
! ===================================================
!
! We estimate cell-averaged 3D Stokes drift velocity (ust_r) by
! presuming that roller Stokes drift has the same z-dependency as
! the primary waves.  Analytical volume-averaging operation is
! applied to 3D Stokes drift velocity such that:
!
!   u^st = 0.5*g*A^2/frq*k/dz/(1-exp(-4*k*D))*
!         ( exp(2*k*(z_up +h-D)) -exp(-2*k*(z_up +h+D))
!          -exp(2*k*(z_low+h-D)) +exp(-2*k*(z_low+h+D)) )
!
! in which sinh, cosh, and tanh functions in u^st are expanded to
! exponential function so as not to be singular when kD gets very
! large.  Here 0.5*g*A^2/frq is equivalent to wave action density.
! The similar rule is also adapted in computation of CalP, Kapsrf,
! and Akw terms (see the code below).
!
!
! 3D depth-dependent Stokes drift velocities
! ==========================================
!
      !print *, 'i0=', i0 
      !print *, 'i1=', i1
      !print *, 'istr-1=', istr-1
      !print *, 'iend+1=', iend+1
!Should this if be here? 
# if defined WAVE_OFFLINE
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstr-1,jend+1
      !  do i=istr-1,iend+1
           !inv_ex(i,j)=1.0/max(1.0-exp(-4.*kD(i,j)),eps)
           !ust_ker=act(i,j)*kw(i,j)*inv_ex(i,j)
          
           cff1 = exp( 2.*kw(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
           !cff6  = exp( 2.*kw(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
           do k=1,N,+1              ! <-- irreversible
              inv_dz(i,j,k)=1./Hz(i,j,k)
              !cff7= exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
              cff2 = exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &               -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
              !ust_rm(i,j,k) =ust_ker*inv_dz(i,j,k)*( cff7-cff6 )
              !ust_rm(i,j,k) = ust_ker * (1.0/Hz(i,j,k)) * (cff7-cff6)!(cff2-cff1)
               
              ust_rm(i,j,k) = act(i,j)*kw(i,j)
     &                       * 1.0/max(1.0-exp(-4.*kD(i,j)),eps)
     &                       * (1./Hz(i,j,k)) * (cff2 - cff1)
!    &                       * inv_dz(i,j,k) *( cff7-cff6 )




            !Define ust_r, vst_r for universal kvf calculation
            ! these are stokes at rho-points
              !print *, 'ust_r(i,j,k) =', ust_r(i,j,k) 
              ust_r(i,j,k) = ust_rm(i,j,k) * wdrx(i,j)
              vst_r(i,j,k) = ust_rm(i,j,k) * wdre(i,j)
              cff1=cff2
              !cff6=cff7

            !Heisenbug area
             !if (mynode.eq.0) then
             !   if (k.eq.N-1) then
             !      if (i.eq.5) then
             !         if (j.eq.5) then
                          ! printing ust_ker also removes heisenbug
             !             print *, 'ust_ker=', ust_ker
             !            print *, 'ust_rm BEFORE=', ust_rm(i,j,k)
             !            print *, 'ust_r BEFORE=', ust_r(i,j,k)
             !         endif
             !      endif
             !   endif
            ! endif

           enddo !<--- k
# if !defined MONOWAVE
           ! Broad spectra
           keff(i,j)=max(2.*pi/20./h(i,j),0.5*sqrt(ustr0(i,j)**2
     &               +vstr0(i,j)**2)/
     &          max(.000001,Dstp(i,j)*sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2)))
           usec=ustr0(i,j)*tanh(2*keff(i,j)*Dstp(i,j))
           vsec=vstr0(i,j)*tanh(2*keff(i,j)*Dstp(i,j))
! rescale kD/kw*keff
          inv_exef(i,j)=1./max(1.0-exp(-4.*kD(i,j)/kw(i,j)*keff(i,j)),eps)
! act--> .5.*usec/keff**2         
          ust_ker=0.5*usec/keff(i,j)*inv_exef(i,j)
          vst_ker=0.5*vsec/keff(i,j)*inv_exef(i,j)

          FACo=0.
          FAC=2.
          COUNTR=0.
          USTERR=.01
          ustm0=ust_rm(i,j,N)
          ! old code has commented out line, think its a mistake
          !ustm0w=sqrt(ustr0(i,j)**2+vstr0(i,j))
          ustm0w=sqrt(ustr0(i,j)**2+vstr0(i,j)**2)
          do while (abs(1-FAC).gt.USTERR.and.
     &         abs(FAC-FACo)/FAC.gt.1e-6.and.COUNTR.lt.15.)
          ! reset cff1
            cff1 = exp( 2.*keff(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*keff(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
            USEZ=0.
            VSEZ=0.
            COUNTR=COUNTR+1.
            if (FACo.ne.0.) then
               ust_rm(i,j,:)=ust_rm(i,j,:)/FAC
            endif
            FACo=FAC
            do k=1,N,+1              ! <-- irreversible
              if (kw(i,j)*Dstp(i,j).gt.2.*pi/20.
     &            .or.abs(ustm0-ustm0w)/ustm0w.gt..1) then ! should avoid if inside loop
!  abs(ustm0-ustm0w)/ustm0w> .1 helps in conditions of mixed wind sea and swell

                cff2 = exp( 2.*keff(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                -exp(-2.*keff(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
                vst_r(i,j,k)=ust_rm(i,j,k)*wdre(i,j)*
     &                       sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &                   (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &                       vst_ker*(1./Hz(i,j,k))*( cff2-cff1 )
                ust_r(i,j,k)=ust_rm(i,j,k)*wdrx(i,j)*
     &                       sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &                   (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &                       ust_ker*(1./Hz(i,j,k))*( cff2-cff1 )
                cff1=cff2
              else
                vst_r(i,j,k)= ust_rm(i,j,k)*wdre(i,j) !
                ust_r(i,j,k)= ust_rm(i,j,k)*wdrx(i,j)  ! monochromatic in shallow water
              endif
              USEZ=USEZ+ust_r(i,j,k)*Hz(i,j,k)
              VSEZ=VSEZ+vst_r(i,j,k)*Hz(i,j,k)
            enddo
            if (sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2).ne.0) then ! ensure consistency betwew ustk and the depth integraded
              FAC=sqrt(USEZ**2+VSEZ**2)/(
     &            sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2)*Dstp(i,j))
            else
              FAC=1.
            endif
          end do ! while
          fac2d(i,j)=FAC
          if (kw(i,j)*Dstp(i,j).lt.2.*pi/20.
     &        .or.abs(ustm0-ustm0w)/ustm0w.lt..1) then
            keff(i,j)=kw(i,j)  ! shallow water
          endif
          !!!! FILLER CODE !!!!
          ! do k=1,N,+1
          !    ust_rm(i,j,k)=0.0
          !    ust_r(i,j,k)=0.0
          !    vst_r(i,j,k)=0.0
          ! enddo
# endif /* !MONOWAVE */

         enddo ! i
       enddo   ! j

! ======================================================
! Horizontal interpolation of ust_r to u- & v- points.
! *** Notice that mask should NOT be applied here.
!
!       Finalize 3-D ust, vst arrays
! ==================================================
       !print *, 'Bounds of ust_r in first dimension: ', LBOUND(ust_r, 1), UBOUND(ust_r, 1) 
       !print *, 'PRE FINALIZE ust_r(4,4,N)', ust_r(4,4,N)
       !print *, 'PRE FINALIZE ust_rm(4,4,N)', ust_rm(4,4,N)
       !print *, 'PRE FINALIZE wdrx(4,4)', wdrx(4,4)
       do k=1,N
         do j=0,ny+1
            do i=0,nx+1
         !do j=jstrR,jendR
         ! do i=istr,iendR
             ust(i,j,k) = 0.5 * (ust_r(i-1,j,k) + ust_r(i,j,k))
# if !defined MONOWAVE 
     &                     /fac2d(i,j)
# endif            
#ifdef MASKING
     &              * umask(i,j)
#endif      
          enddo
        enddo

        do j=0,ny+1
           do i=0,nx+1
        !do j=jstr,jendR
        !  do i=istrR,iendR
            vst(i,j,k) = 0.5 * (vst_r(i,j-1,k) + vst_r(i,j,k))
# if !defined MONOWAVE 
     &                     /fac2d(i,j)
# endif 
#ifdef MASKING
     &              * vmask(i,j)
#endif      
         enddo
        enddo
      enddo
# endif /* WAVEOFFLINE */
!------------------------------------------------------------
      ! Stokes drift calculated now
!-----------------------------------------------------------

# if defined BREAKING && defined LMD_MIXING && !defined AKB0
!====================================================
!       Non-conservative effects
!====================================================
!
! Akw: primary wave-induced additional diffusivity for tracers at w-point
! =======================================================================
!
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstrR,jendR
      !  do i=istrR,iendR
#ifndef MONOWAVE
        keffkw=2*pi/max(1.,llmm(i,j))
#else
        keffkw=kw(i,j)
#endif
        cff2 = 0.25*wh(i,j)/max(1.-exp(-2.*keffkw*Dstp(i,j)),eps)
         do k=0,N
            cff1 = cff2*( exp( keffkw*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                   -exp(-keffkw*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
            cff3 = cff1**2
            ! Akw disabled - represents kinematic spreading, not physical mixing
            ! See Uchiyama 2010 Eq. 6 discussion
            Akw(i,j,k) = 0.0
            E_pre(i,j,k) = cff3    ! store this for the next time step.
          enddo
         !if (mynode==0 .and. j==1 .and. i==1) then
         !  print*,'Akw(i,j,N)=',Akw(i,j,N) ! Debug
         !endif
         enddo
      enddo
!----------------------------------------
!

! Akb: Eddy viscosity due to depth-induced wave breaking
! ======================================================
!
#  define kb0 wrk4
#  define vs  cff1
#  define kvb cff2
      !if (mynode==0) then
      ! print*,'Here'
      !endif
      do j=0,ny+1
         do i=0,nx+1
      !do j=jstrR,jendR
      !  do i=istrR,iendR
# ifdef Akb_U2010

          inv_k =min(max(a_kv*wh(i,j),eps),Dstp(i,j))
          kbrk  =1.0/inv_k         ! vertical scale for breaking

          fb1   =0.0               ! fb1 = \int gb dz at z = -h
          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))


          vs    =ebrk(i,j)**c1o3   ! representative velocity scale
          kvb   =bconst*vs*wh(i,j) ! depth-averaged Akb
          do k=1,N,+1
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))

            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            kb0(k) = kvb*fb*Dstp(i,j)
            fb1=fb2     ! recursive procedure
          enddo
          do k=1,N-1
            Akb(i,j,k)=0.5*(kb0(k)+kb0(k+1))
          enddo
          Akb(i,j,0)=0.D0
          Akb(i,j,N)=max(1.5*kb0(N)-0.5*kb0(N-1),0.0)
            fb2 =2.*inv_k*fn1/(1.+fn1)
# else   /* Akb_U2010 */

          ! energy dissipation profile based on scaling by Sutherland
          ! and Melville 2015 according to Romero et al. 2021
          ! Implemented by L. Romero
          hs=max(wh(i,j)*sqrt(2.),.1)
          if (wh(i,j)*sqrt(2.).gt..1) then
            ! bottom depth
           zzb=(z_w(i,j,0)-z_w(i,j,N)) !
           L=0.
           do k=1,N!,+1              ! <-- irreversible
            zz=(z_w(i,j,k)-z_w(i,j,N))/hs
            zzm=(z_w(i,j,k-1)-z_w(i,j,N))/hs
            hsfs=max(0.,(sign(1.,zz*hs+a_kv*hs)+1.)/2.)  ! surface
            if (zz.gt.-3.) then
                    if (mynode==0 .and. j==1 .and. i==1 .and. k==10) then
                            print*,'In the loop'
                    endif
            kb0(k)=1./4.*sqrt(3.)*(
     &             log(2.+sqrt(3.)-5.*min(0.,zz))
     &            -log(2.-sqrt(3.)-5.*min(0.,zz))
     &            -log(2.+sqrt(3.)-5.*min(0.,zzm))
     &            +log(2.-sqrt(3.)-5.*min(0.,zzm)))*hs/Hz(i,j,k)
            kb0(k)=kb0(k)*.9508  ! normalization to integrate to 1 for z'>-3
            else
            kb0(k)=0.0
            endif
            Akb(i,j,k)=(ebrk(i,j)/hs*kb0(k))  ! dissipation
            L= max(0.,vonKar*(a_kv*hs*hsfs+(1.-hsfs)*(-zz*hs)))  ! Terray 1999
            Akb(i,j,k)=(Akb(i,j,k)*L)**c1o3*L
            Akb(i,j,k)=Akb(i,j,k)*(1.-tanh(1.*zz)**4.) ! tapering function
           !if (mynode==0 .and. j==1 .and. i==1 .and. k==10) then
                   !print*,'in the if non zero' ! Debug
                   !print*,'Akb(i,j,10)=',Akb(i,j,k)
                   !print*,'Akb(i,j,10)=',ebrk(i,j)/hs*kb0(k)
                   !print*,'ebrk(i,j)=',ebrk(i,j)
                   !print*,'kb0(k)=',kb0(k)
                   !print*,'zz=',zz
                   !print*,'zz=',(z_w(i,j,k)-z_w(i,j,N))/hs
                   !print*,'L=',L
                   !print*,'c1o3=',c1o3
                   !print*,'hs=',hs
                   !print*,'zzb=',zzb
                   !print*,'zzm=',zzm
                   !print*,'a_kv=',a_kv
                   !print*,'vonKar=',vonKar
                   !print*,'hsfs=',hsfs
            !endif

           enddo

          else
            !if (mynode==0 .and. j==1 .and. i==1) then
            ! print*,'in the else zero' ! Debug
            !endif
            Akb(i,j,:)=0.
          endif

# endif /* Akb_U2010 */
#  ifdef MASKING
          do k=0,N
            Akb(i,j,k)=Akb(i,j,k)*rmask(i,j)
          enddo
#  endif
        !if (mynode==0 .and. j==1 .and. i==1) then
        !  print*,'Akb(i,j,N)=',Akb(i,j,N) ! Debug
        !endif
        enddo
      enddo
# endif    /* BREAKING && LMD_MIXING && !AKB0 */

# if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
!
! 3D breaking acceleration term defined as a body force.
! =====================================================
!

      do j=0,ny+1
         do i=0,nx+1
      !do j=jstr-1,jend+1
      !  do i=istr-1,iend+1
          inv_k =min(max(a_brk*wh(i,j),eps),Dstp(i,j))
          kbrk  =1.0/inv_k         ! vertical scale for breaking
          fb1   =0.0               ! fb1 = \int gb dz at z = -h
          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))

          do k=1,N,+1              ! <-- irreversible
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            brk_r(i,j,k) = fb*brk(i,j)
            fb1=fb2     ! recursive procedure
          enddo

          enddo
        enddo


      do k=1,N
        do j=0,ny+1
           do i=0,nx+1
        !do j=jstrR,jendR
        !  do i=istr,iendR
            brk3dx(i,j,k)=0.5*( brk_r(i-1,j,k)*wdrx(i-1,j)
     &                           +brk_r(i,j,k)*wdrx(i,j) )
#  ifdef MASKING
     &                                         *umask(i,j)
#  endif
          enddo
        enddo

        do j=0,ny+1
           do i=0,nx+1
        !do j=jstr,jendR
         ! do i=istrR,iendR
            brk3de(i,j,k)=0.5*( brk_r(i,j-1,k)*wdre(i,j-1)
     &                           +brk_r(i,j,k)*wdre(i,j) )
#  ifdef MASKING
     &                                         *vmask(i,j)
#  endif
          enddo
        enddo
      enddo  
# endif        /* BREAKING && SURFACE_BREAK && BODY_BREAK */

# if defined WAVE_FRICTION && defined BODY_FRICTION
!
! 3D streaming acceleration term defined as a body force.
! =======================================================
! Vertical profile using FF_FUNC3 shape function (src/wec_frc.F:2640-2739).
! abot = near-bed orbital displacement amplitude
! delta = WBBL thickness (U2010 Eq. 56)
! frc_r(k) = frc * fb(k), projected onto xi/eta → frc3dx, frc3de
!
      do j=0,ny+1
        do i=0,nx+1
         if (wh(i,j) > eps) then
#  ifndef UORBOFF
          abot  =wh(i,j)/max(2.D0*sinh(min(kD(i,j),khmax)),eps)
#  else
          abot  =Uorb_mag(i,j)*inv_f(i,j)
#  endif
          delta =0.160455*((abot/max(Zob,eps))**0.82)*max(Zob,eps)  ! exponent 0.82 (Grant-Madsen), was 0.83 typo
          inv_k =min(max(a_frc*delta,eps),Dstp(i,j))
          kfrc  =1.0/inv_k

          fb1   =1.0-exp(-2.*kfrc*Dstp(i,j))
          fb0   =1.D0/(exp(-2.*kfrc*Dstp(i,j))-1.0)
          do k=1,N,+1
            fb2 =exp( kfrc*(zeta(i,j,knew)-z_w(i,j,k)-Dstp(i,j)))
     &          -exp(-kfrc*(zeta(i,j,knew)-z_w(i,j,k)+Dstp(i,j)))
            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            frc_r(i,j,k) = fb*frc(i,j)
            fb1=fb2
          enddo
         else
          do k=1,N
            frc_r(i,j,k) = 0.0
          enddo
         endif
        enddo
      enddo

      do k=1,N
        do j=0,ny+1
          do i=0,nx+1
            frc3dx(i,j,k)=0.5*( frc_r(i-1,j,k)*wdrx(i-1,j)
     &                           +frc_r(i,j,k)*wdrx(i,j) )
#  ifdef MASKING
     &                                         *umask(i,j)
#  endif
          enddo
        enddo

        do j=0,ny+1
          do i=0,nx+1
            frc3de(i,j,k)=0.5*( frc_r(i,j-1,k)*wdre(i,j-1)
     &                           +frc_r(i,j,k)*wdre(i,j) )
#  ifdef MASKING
     &                                         *vmask(i,j)
#  endif
          enddo
        enddo
      enddo
# endif        /* WAVE_FRICTION && BODY_FRICTION */

     !COMPUTE wstokes
      call wstokes_wec_tile(istr,iend,jstr,jend)
# endif       /* SOLVE3D */
      !print *, 'SUCCESS SUCCESS SUCCESS wec_forces_tile'
      !print *, 'ust(4,4,N)=', ust(4,4,N)
      !print *, 'ust_r(4,4,n)=', ust_r(4,4,N)

#ifdef EXCHANGE
      ! MPI EXCHANGES
      call exchange_xxx(sup)
      call exchange_xxx(ust2d,brk2dx,frc2dx)
      call exchange_xxx(vst2d,brk2de,frc2de)
#ifdef BBL_S95
      call exchange_xxx(r_D) !WHY IS THIS EXCHANGED??
      call exchange_xxx(Uorb_mag)
      call exchange_xxx(tauw_orb)
#endif
#ifdef SOLVE_3D
      call exchange_xxx(ust_r)
      call exchange_xxx(Kapsrf,calP)
      call exchange_xxx(ust,vst,kvf)
      call exchange_xxx(wst)
      call exchange_xxx(Akw,Akb)
#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      call exchange_xxx(brk3dx,brk3de)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
      call exchange_xxx(frc3dx,frc3de)
#  endif
#endif
#endif



!# ifdef EXCHANGE
!      call exchange2d_tile(istr,iend,jstr,jend, sup)
!      call exchange2d_3_tile(istr,iend,jstr,jend, ust2d,brk2dx,frc2dx)
!      call exchange2d_3_tile(istr,iend,jstr,jend, vst2d,brk2de,frc2de)
!#  if defined BBL_S95
!      call exchange2d_tile (istr,iend,jstr,jend, r_D)
!      call exchange2d_tile (istr,iend,jstr,jend, Uorb_mag)
!      call exchange2d_tile (istr,iend,jstr,jend, tauw_orb)
!#  endif
!#  ifdef SOLVE3D
!      call exchange2d_2_tile(istr,iend,jstr,jend,Kapsrf,calP)
!      call exchange_3_tile(istr,iend,jstr,jend, ust,N, vst,N, kvf,N)
!      call exchange_2_tile(istr,iend,jstr,jend,Akw,N+1,Akb,N+1)
!#   if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
!      call exchange_2_tile (istr,iend,jstr,jend, brk3dx,N, brk3de,N)
!#   endif
!#  endif
!# endif
      end subroutine wec_forces_tile

!-----------------------------------------------------------


      ! ##################
      ! NETCDF SUBROUTINES
      ! ##################


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                 WRITING OUTPUT SUBROUTINES
      !                 these follow functionality in
      !                 basic_output.F and bgc.F of June 2024 stable repo

      ! subroutines here
      !        wrt_wec
      !        create_file_wec_vars
      !        def_vars_wec_vars
      !        calc_avg_wec
      !        wrt_wec_his
      !        wrt_wec_avg

      !         basic flow here:
      !         wrt_wec is master function
      !         this calls wrt_wec_his and/or wrt_wec_avg
      !         wrt_wec_his/avg --> create file, define
      !         variables,calc_avg, and then write
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine wrt_wec(init)  ![
      implicit none

      ! import/export
      logical,optional       :: init

      ! initialisation writting
      if (wrt_his.and.present(init)) then
         call wrt_wec_his! write his data
      else
      ! local
         if (wrt_his)     call wrt_wec_his ! write his data
         if (wrt_avg)     call wrt_wec_avg ! write avg data
      endif
      
      end subroutine wrt_wec  !]

      subroutine create_file_wec_vars(fname,avg) ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ncid,ierr,varid
      
      if (avg) then
         call create_file('_wec_avg',fname)
      else 
         call create_file('_wec_his',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      if (avg) then
        call def_vars_wec_vars(ncid, .true.  )
        ierr=nf90_put_att(ncid,nf90_global,'type','WEC averages file')
      else
        call def_vars_wec_vars(ncid, .false. )
        ierr=nf90_put_att(ncid,nf90_global,'type','WEC history file')
      endif

      ierr=nf90_close(ncid)

      end subroutine create_file_wec_vars !]
      !------------------------------------



      subroutine def_vars_wec_vars(ncid,avg)
      ! define output variable & attributes in netcdf results file
      ! for instantaneous or averaged variables
    
      implicit none
      
      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid
      character(len=64) :: text_lname
      character(len=7)  :: dn_aux='auxil'


! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
     ! varid = nccreate(ncid,'time_step',(/dn_aux,dn_tm/),(/iaux,0/),nf90_int)
     ! ierr=nf90_put_att (ncid, varid, 'long_name',
     !&       'time step and record numbers from initialization')




      ! ust2d
      if (wrt_ust2d) then
        varid = nccreate(ncid,ust2d_name,(/dn_xu,dn_yr,
     &                         dn_tm/),(/xi_u,eta_rho,0/))
        if (.not. avg) text_lname=ust2d_long_name
        if (      avg) text_lname='averaged '/ /ust2d_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',ust2d_units)
      endif

      ! vst2d
      if (wrt_vst2d) then
        varid = nccreate(ncid,vst2d_name,(/dn_xr,dn_yv,
     &                         dn_tm/),(/xi_rho,eta_v,0/))
        if (.not. avg) text_lname=vst2d_long_name
        if (      avg) text_lname='averaged '/ /vst2d_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vst2d_units)
      endif


      ! ust
      if (wrt_ust) then
        varid = nccreate(ncid,ust_name,(/dn_xu,dn_yr,
     &                       dn_zr,dn_tm/),(/xi_u,eta_rho,N,0/))
        if (.not. avg) text_lname=ust_long_name
        if (      avg) text_lname='averaged '/ /ust_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',ust_units)
      endif

      ! vst
      if (wrt_vst) then
        varid = nccreate(ncid,vst_name,(/dn_xr,dn_yv,
     &                        dn_zr, dn_tm/),(/xi_rho,eta_v,N,0/))
        if (.not. avg) text_lname=vst_long_name
        if (      avg) text_lname='averaged '/ /vst_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vst_units)
      endif

      ! sup
      if (wrt_sup) then
        varid = nccreate(ncid,sup_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=sup_long_name
        if (      avg) text_lname='averaged '/ /sup_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',sup_units)
      endif
      
      if (wrt_sup_zeta) then
        varid = nccreate(ncid,sup_zeta_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=sup_zeta_long_name
        if (      avg) text_lname='averaged '/ /sup_zeta_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',sup_units)
      endif

      !frc2dx
      if (wrt_frc2dx) then
        varid = nccreate(ncid,frc2dx_name,(/dn_xu,dn_yr,
     &                         dn_tm/),(/xi_u,eta_rho,0/))
        if (.not. avg) text_lname=frc2dx_long_name
        if (      avg) text_lname='averaged '/ /frc2dx_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',frc2dx_units)
      endif


      !frc2de
      if (wrt_frc2de) then
        varid = nccreate(ncid,frc2de_name,(/dn_xr,dn_yv,
     &                         dn_tm/),(/xi_rho,eta_v,0/))
        if (.not. avg) text_lname=frc2de_long_name
        if (      avg) text_lname='averaged '/ /frc2de_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',frc2de_units)
      endif

      !brk2dx
      if (wrt_brk2dx) then
        varid = nccreate(ncid,brk2dx_name,(/dn_xu,dn_yr,
     &                         dn_tm/),(/xi_u,eta_rho,0/))
        if (.not. avg) text_lname=brk2dx_long_name
        if (      avg) text_lname='averaged '/ /brk2dx_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',brk2dx_units)
      endif


      !brk2de
      if (wrt_brk2de) then
        varid = nccreate(ncid,brk2de_name,(/dn_xr,dn_yv,
     &                         dn_tm/),(/xi_rho,eta_v,0/))
        if (.not. avg) text_lname=brk2de_long_name
        if (      avg) text_lname='averaged '/ /brk2de_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',brk2de_units)
      endif

#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      !brk3dx
      if (wrt_brk3dx) then
        varid = nccreate(ncid,brk3dx_name,(/dn_xu,dn_yr,
     &                        dn_zr, dn_tm/),(/xi_u,eta_rho,N,0/))
        if (.not. avg) text_lname=brk3dx_long_name
        if (      avg) text_lname='averaged '/ /brk3dx_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',brk3dx_units)
      endif

      !brk3de
      if (wrt_brk3de) then
        varid = nccreate(ncid,brk3de_name,(/dn_xr,dn_yv,
     &                        dn_zr, dn_tm/),(/xi_rho,eta_v,N,0/))
        if (.not. avg) text_lname=brk3de_long_name
        if (      avg) text_lname='averaged '/ /brk3de_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',brk3de_units)
      endif
#  endif

#  if defined WAVE_FRICTION && defined BODY_FRICTION
      !frc3dx
      if (wrt_frc3dx) then
        varid = nccreate(ncid,frc3dx_name,(/dn_xu,dn_yr,
     &                        dn_zr, dn_tm/),(/xi_u,eta_rho,N,0/))
        if (.not. avg) text_lname=frc3dx_long_name
        if (      avg) text_lname='averaged '/ /frc3dx_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',frc3dx_units)
      endif

      !frc3de
      if (wrt_frc3de) then
        varid = nccreate(ncid,frc3de_name,(/dn_xr,dn_yv,
     &                        dn_zr, dn_tm/),(/xi_rho,eta_v,N,0/))
        if (.not. avg) text_lname=frc3de_long_name
        if (      avg) text_lname='averaged '/ /frc3de_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',frc3de_units)
      endif
#  endif


      ! wst
      if (wrt_wst) then
        varid = nccreate(ncid,wst_name,(/dn_xr,dn_yr,
     &                        dn_zr, dn_tm/),(/xi_rho,eta_rho,N,0/))
        if (.not. avg) text_lname=wst_long_name
        if (      avg) text_lname='averaged '/ /wst_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',wst_units)
      endif

      ! Akw 
      if (wrt_Akw) then
        varid = nccreate(ncid,Akw_name,(/dn_xr,dn_yr,
     &                        dn_zw, dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        if (.not. avg) text_lname=Akw_long_name
        if (      avg) text_lname='averaged '/ /Akw_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',Akw_units)
      endif


      ! Akb 
      if (wrt_Akb) then
        varid = nccreate(ncid,Akb_name,(/dn_xr,dn_yr,
     &                        dn_zw, dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        if (.not. avg) text_lname=Akb_long_name
        if (      avg) text_lname='averaged '/ /Akb_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',Akb_units)
      endif

      ! kvf 
      if (wrt_kvf) then
        varid = nccreate(ncid,kvf_name,(/dn_xr,dn_yr,
     &                        dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        if (.not. avg) text_lname=kvf_long_name
        if (      avg) text_lname='averaged '/ /kvf_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',kvf_units)
      endif

      ! calP 
      if (wrt_calP) then
        varid = nccreate(ncid,calP_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=calP_long_name
        if (      avg) text_lname='averaged '/ /calP_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',calP_units)
      endif


      ! Kapsrf 
      if (wrt_Kapsrf) then
        varid = nccreate(ncid,Kapsrf_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=Kapsrf_long_name
        if (      avg) text_lname='averaged '/ /Kapsrf_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',Kapsrf_units)
      endif

      ! wveb 
      if (wrt_wveb) then
        varid = nccreate(ncid,wveb_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=wveb_long_name
        if (      avg) text_lname='averaged '/ /wveb_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',wveb_units)
      endif

      ! wved 
      if (wrt_wved) then
        varid = nccreate(ncid,wved_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=wved_long_name
        if (      avg) text_lname='averaged '/ /wved_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',wved_units)
      endif

      ! hvf_int_xi (depth-integrated horizontal vortex force, xi)
      if (wrt_hvf_int) then
        varid = nccreate(ncid,hvf_int_xi_name,(/dn_xu,dn_yr,
     &                         dn_tm/),(/xi_u,eta_rho,0/))
        if (.not. avg) text_lname=hvf_int_xi_long_name
        if (      avg) text_lname='averaged '//hvf_int_xi_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',hvf_int_xi_units)
      endif

      ! hvf_int_eta (depth-integrated horizontal vortex force, eta)
      if (wrt_hvf_int) then
        varid = nccreate(ncid,hvf_int_eta_name,(/dn_xr,dn_yv,
     &                         dn_tm/),(/xi_rho,eta_v,0/))
        if (.not. avg) text_lname=hvf_int_eta_long_name
        if (      avg) text_lname='averaged '//hvf_int_eta_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',hvf_int_eta_units)
      endif

      ! r_D (bottom drag coefficient)
      if (wrt_rD) then
        varid = nccreate(ncid,rD_name,(/dn_xr,dn_yr,
     &                         dn_tm/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname=rD_long_name
        if (      avg) text_lname='averaged '//rD_long_name
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',rD_units)
      endif


      end subroutine def_vars_wec_vars
      !--------------------------------------------------


      subroutine calc_avg_wec
      !use dimensions !i0,i1,etc
      implicit none

      ! local
      real    :: coef
      integer :: itrc, itavg, k

      navg_wec = navg_wec +1
      coef = 1./navg_wec

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'wec :: started averaging. ',
     &   'output_period_avg (s) =', output_period_avg
      endif

      t_avg_wec = t_avg_wec*(1-coef) + time*coef
      ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin

      !-----------------------------
      !         2-D variables
      !----------------------------
      if (wrt_ust2d) then 
         ust2d_avg(1:i1,j0:j1) = ust2d_avg(1:i1,j0:j1) * (1-coef) 
     &                           + ust2d(1:i1,j0:j1) * coef
      endif

      if (wrt_vst2d) then 
         vst2d_avg(i0:i1,1:j1) = vst2d_avg(i0:i1,1:j1) * (1-coef) 
     &                           + vst2d(i0:i1,1:j1) * coef
      endif

      if (wrt_sup) then 
         sup_avg(i0:i1,j0:j1) = sup_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + sup(i0:i1,j0:j1) * coef
      endif

      if (wrt_sup_zeta) then 
         sup_zeta_avg(i0:i1,j0:j1) = sup_zeta_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + sup_zeta(i0:i1,j0:j1) * coef
      endif



      if (wrt_calP) then 
         calP_avg(i0:i1,j0:j1) = calP_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + calP(i0:i1,j0:j1) * coef
      endif

      if (wrt_Kapsrf) then 
         Kapsrf_avg(i0:i1,j0:j1) = Kapsrf_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + Kapsrf(i0:i1,j0:j1) * coef
      endif


      if (wrt_wveb) then 
         wveb_avg(i0:i1,j0:j1) = wveb_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + wveb(i0:i1,j0:j1) * coef
      endif

      if (wrt_wved) then 
         wved_avg(i0:i1,j0:j1) = wved_avg(i0:i1,j0:j1) * (1-coef) 
     &                           + wved(i0:i1,j0:j1) * coef
      endif


      if (wrt_frc2dx) then 
         frc2dx_avg(i0:i1,j0:j1) = frc2dx_avg(1:i1,j0:j1) * (1-coef) 
     &                           + frc2dx(1:i1,j0:j1) * coef
      endif
 
      if (wrt_frc2de) then 
         frc2de_avg(i0:i1,j0:j1) = frc2de_avg(i0:i1,1:j1) * (1-coef) 
     &                           + frc2de(i0:i1,1:j1) * coef
      endif


      if (wrt_brk2dx) then 
         brk2dx_avg(i0:i1,j0:j1) = brk2dx_avg(1:i1,j0:j1) * (1-coef) 
     &                           + brk2dx(1:i1,j0:j1) * coef
      endif
 
      if (wrt_brk2de) then 
         brk2de_avg(i0:i1,j0:j1) = brk2de_avg(i0:i1,1:j1) * (1-coef) 
     &                           + brk2de(i0:i1,1:j1) * coef
      endif

      if (wrt_hvf_int) then
         hvf_int_xi_avg(1:i1,j0:j1) =
     &        hvf_int_xi_avg(1:i1,j0:j1) * (1-coef)
     &      + hvf_int_xi(1:i1,j0:j1) * coef
      endif

      if (wrt_hvf_int) then
         hvf_int_eta_avg(i0:i1,1:j1) =
     &        hvf_int_eta_avg(i0:i1,1:j1) * (1-coef)
     &      + hvf_int_eta(i0:i1,1:j1) * coef
      endif


      !-----------------------------
      !         3-D variables
      !----------------------------
      if (wrt_ust) then
         ust_avg(1:i1,j0:j1,:) = ust_avg(1:i1,j0:j1,:) * (1-coef)
     &                           + ust(1:i1,j0:j1,:) * coef
      endif

      if (wrt_vst) then
         vst_avg(i0:i1,1:j1,:) = vst_avg(i0:i1,1:j1,:) * (1-coef)
     &                           + vst(i0:i1,1:j1,:) * coef
      endif

      if (wrt_wst) then
         wst_avg(i0:i1,j0:j1,:) = wst_avg(i0:i1,j0:j1,:) * (1-coef)
     &                           + wst(i0:i1,j0:j1,:) * coef
      endif

      if (wrt_Akw) then 
         Akw_avg(i0:i1,j0:j1,:) = Akw_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + Akw(i0:i1,j0:j1,:) * coef
      endif

      if (wrt_Akb) then 
         Akb_avg(i0:i1,j0:j1,:) = Akb_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + Akb(i0:i1,j0:j1,:) * coef
      endif

      if (wrt_kvf) then
         kvf_avg(i0:i1,j0:j1,:) = kvf_avg(i0:i1,j0:j1,:) * (1-coef)
     &                           + kvf(i0:i1,j0:j1,:) * coef
      endif

#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
      if (wrt_brk3dx) then
         brk3dx_avg(1:i1,j0:j1,:) = brk3dx_avg(1:i1,j0:j1,:)*(1-coef)
     &                             + brk3dx(1:i1,j0:j1,:) * coef
      endif

      if (wrt_brk3de) then
         brk3de_avg(i0:i1,1:j1,:) = brk3de_avg(i0:i1,1:j1,:)*(1-coef)
     &                             + brk3de(i0:i1,1:j1,:) * coef
      endif
#  endif

#  if defined WAVE_FRICTION && defined BODY_FRICTION
      if (wrt_frc3dx) then
         frc3dx_avg(1:i1,j0:j1,:) = frc3dx_avg(1:i1,j0:j1,:)*(1-coef)
     &                             + frc3dx(1:i1,j0:j1,:) * coef
      endif

      if (wrt_frc3de) then
         frc3de_avg(i0:i1,1:j1,:) = frc3de_avg(i0:i1,1:j1,:)*(1-coef)
     &                             + frc3de(i0:i1,1:j1,:) * coef
      endif
#  endif


      end subroutine calc_avg_wec


      subroutine wrt_wec_his

      ! Adpated from Devin and Jeroen's wrt_surface_flux_his in
      ! surf_flux.F

      !Made compataible with new I/O Oct 2022
      use coupling, only: r_D
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_his                              ! current file output record
      integer,save           :: total_rec_his=0                      ! total his output records so far
      real,save              :: output_time_his=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_his
      integer                :: tile, ierr

      if (wrt_his) then

      if (.not. first_step) output_time_his = output_time_his + dt   ! only start count after first timestep
      first_step=.false.                                             ! as first step the wec values are for t=0
      if (output_time_his>=output_period_his .or.
     &    output_time_his==0                      ) then             ! catch to get t=0 output

        if (mod(total_rec_his,nrpf_his)==0) then
          call create_file_wec_vars(fname_his,.false.)
          rec_his = 0
        endif
        total_rec_his = total_rec_his +1
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))
        ! start --> vector --> (/1,1,1,rec_his/) for 3-D variable 
        ! start --> vector --> (/1,1,rec_his/) for 2-D variable
        ! ---------------------
        ! 2-D Variables
        !--------------------
        start=1; start(3)=rec_his
        if (wrt_ust2d) call ncwrite(ncid,ust2d_name,
     &                           ust2d(1:i1,j0:j1),start)
        if (wrt_vst2d) call ncwrite(ncid,vst2d_name,
     &                           vst2d(i0:i1,1:j1),start)
        if (wrt_sup) call ncwrite(ncid,sup_name,
     &                           sup(i0:i1,j0:j1),start)
       if (wrt_sup_zeta) call ncwrite(ncid,sup_zeta_name,
     &                           sup_zeta(i0:i1,j0:j1),start)
        
        if (wrt_calP) call ncwrite(ncid,calP_name,
     &                           calP(i0:i1,j0:j1),start)
        if (wrt_Kapsrf) call ncwrite(ncid,Kapsrf_name,
     &                           Kapsrf(i0:i1,j0:j1),start)
        if (wrt_wveb) call ncwrite(ncid,wveb_name,
     &                           wveb(i0:i1,j0:j1),start)
        if (wrt_wved) call ncwrite(ncid,wved_name,
     &                           wved(i0:i1,j0:j1),start)
        if (wrt_frc2dx) call ncwrite(ncid,frc2dx_name,
     &                           frc2dx(1:i1,j0:j1),start)
        if (wrt_frc2de) call ncwrite(ncid,frc2de_name,
     &                           frc2de(i0:i1,1:j1),start)

        if (wrt_brk2dx) call ncwrite(ncid,brk2dx_name,
     &                           brk2dx(1:i1,j0:j1),start)
        if (wrt_brk2de) call ncwrite(ncid,brk2de_name,
     &                           brk2de(i0:i1,1:j1),start)

        if (wrt_hvf_int) call ncwrite(ncid,hvf_int_xi_name,
     &                           hvf_int_xi(1:i1,j0:j1),start)
        if (wrt_hvf_int) call ncwrite(ncid,hvf_int_eta_name,
     &                           hvf_int_eta(i0:i1,1:j1),start)
        if (wrt_rD) call ncwrite(ncid,rD_name,
     &                           r_D(i0:i1,j0:j1),start)


        ! ---------------------
        ! 3-D Variables
        !--------------------
        start=1; start(4)=rec_his
        if (wrt_ust) call ncwrite(ncid,ust_name,
     &                           ust(1:i1,j0:j1,:),start)
        if (wrt_vst) call ncwrite(ncid,vst_name,
     &                           vst(i0:i1,1:j1,:),start)
        if (wrt_wst) call ncwrite(ncid,wst_name,
     &                           wst(i0:i1,j0:j1,:),start)
        if (wrt_Akw) call ncwrite(ncid,Akw_name,
     &                           Akw(i0:i1,j0:j1,:),start)
        if (wrt_Akb) call ncwrite(ncid,Akb_name,
     &                           Akb(i0:i1,j0:j1,:),start)

        if (wrt_kvf) call ncwrite(ncid,kvf_name,
     &                           kvf(i0:i1,j0:j1,:),start)

#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
        if (wrt_brk3dx) call ncwrite(ncid,brk3dx_name,
     &                           brk3dx(1:i1,j0:j1,:),start)
        if (wrt_brk3de) call ncwrite(ncid,brk3de_name,
     &                           brk3de(i0:i1,1:j1,:),start)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
        if (wrt_frc3dx) call ncwrite(ncid,frc3dx_name,
     &                           frc3dx(1:i1,j0:j1,:),start)
        if (wrt_frc3de) call ncwrite(ncid,frc3de_name,
     &                           frc3de(i0:i1,1:j1,:),start)
#  endif

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wec :: wrote history, tdays =', tdays,
     &     'step =', iic-1, 'rec =', rec_his, '/', total_rec_his     ! -1 since iic=t(n+1)
        endif

        output_time_his=0
      endif
      endif  !wrt_file_his 

      end subroutine wrt_wec_his


      subroutine wrt_wec_avg
      ! write averaged variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      use coupling, only: r_D
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg                              ! current file output record
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      integer                :: itrc, ierr, itavg

      if (wrt_avg) then

      if (.not. first_step) then
        call calc_avg_wec                                         ! don't include t=0 in averaging
        output_time_avg = output_time_avg + dt                       ! only start count after first timestep
      endif
      first_step=.false.                                             ! as first step the bgc values are for t=0

      if (output_time_avg>=output_period_avg) then                   ! catch to get t=0 output

        if (mod(total_rec_avg,nrpf_avg)==0) then
          call create_file_wec_vars(fname_avg,.true.)
          rec_avg = 0
        endif
        total_rec_avg = total_rec_avg +1
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_wec/),(/rec_avg/))
      
        ! ---------------------
        ! 2-D Variables
        !--------------------
        start=1; start(3)=rec_avg
        if (wrt_ust2d) call ncwrite(ncid,ust2d_name,
     &                           ust2d_avg(1:i1,j0:j1),start)
        if (wrt_vst2d) call ncwrite(ncid,vst2d_name,
     &                           vst2d_avg(i0:i1,1:j1),start)
        if (wrt_sup) call ncwrite(ncid,sup_name,
     &                           sup_avg(i0:i1,j0:j1),start)
        if (wrt_sup_zeta) call ncwrite(ncid,sup_zeta_name,
     &                           sup_zeta_avg(i0:i1,j0:j1),start)
        if (wrt_calP) call ncwrite(ncid,calP_name,
     &                           calP_avg(i0:i1,j0:j1),start)
        if (wrt_Kapsrf) call ncwrite(ncid,Kapsrf_name,
     &                           Kapsrf_avg(i0:i1,j0:j1),start)
        if (wrt_wveb) call ncwrite(ncid,wveb_name,
     &                           wveb_avg(i0:i1,j0:j1),start)
        if (wrt_wved) call ncwrite(ncid,wved_name,
     &                           wved_avg(i0:i1,j0:j1),start)
        if (wrt_frc2dx) call ncwrite(ncid,frc2dx_name,
     &                           frc2dx_avg(1:i1,j0:j1),start)
        if (wrt_frc2de) call ncwrite(ncid,frc2de_name,
     &                           frc2de_avg(i0:i1,1:j1),start)
        if (wrt_brk2dx) call ncwrite(ncid,brk2dx_name,
     &                           brk2dx_avg(1:i1,j0:j1),start)
        if (wrt_brk2de) call ncwrite(ncid,brk2de_name,
     &                           brk2de_avg(i0:i1,1:j1),start)

        if (wrt_hvf_int) call ncwrite(ncid,hvf_int_xi_name,
     &                           hvf_int_xi_avg(1:i1,j0:j1),start)
        if (wrt_hvf_int) call ncwrite(ncid,hvf_int_eta_name,
     &                           hvf_int_eta_avg(i0:i1,1:j1),start)
        if (wrt_rD) call ncwrite(ncid,rD_name,
     &                           r_D(i0:i1,j0:j1),start)


        ! ---------------------
        ! 3-D Variables
        !--------------------
        start=1; start(4)=rec_avg
        if (wrt_ust) call ncwrite(ncid,ust_name,
     &                           ust_avg(1:i1,j0:j1,:),start)
        if (wrt_vst) call ncwrite(ncid,vst_name,
     &                           vst_avg(i0:i1,1:j1,:),start) 
        if (wrt_wst) call ncwrite(ncid,wst_name,
     &                           wst_avg(i0:i1,j0:j1,:),start) 
        if (wrt_Akw) call ncwrite(ncid,Akw_name,
     &                           Akw_avg(i0:i1,j0:j1,:),start) 

        if (wrt_Akb) call ncwrite(ncid,Akb_name,
     &                           Akb_avg(i0:i1,j0:j1,:),start) 
        if (wrt_kvf) call ncwrite(ncid,kvf_name,
     &                           kvf_avg(i0:i1,j0:j1,:),start)

#  if defined BREAKING && defined SURFACE_BREAK && defined BODY_BREAK
        if (wrt_brk3dx) call ncwrite(ncid,brk3dx_name,
     &                           brk3dx_avg(1:i1,j0:j1,:),start)
        if (wrt_brk3de) call ncwrite(ncid,brk3de_name,
     &                           brk3de_avg(i0:i1,1:j1,:),start)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
        if (wrt_frc3dx) call ncwrite(ncid,frc3dx_name,
     &                           frc3dx_avg(1:i1,j0:j1,:),start)
        if (wrt_frc3de) call ncwrite(ncid,frc3de_name,
     &                           frc3de_avg(i0:i1,1:j1,:),start)
#  endif

        navg_wec=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'wec_vars :: wrote averages, tdays =', tdays,
     &        'step =', iic-1, 'rec =', rec_avg, '/', total_rec_avg  ! -1 since iic=t(n+1)
     &         MYID
        endif
      endif

      endif  ! <-- wrt_file_avg


      end subroutine wrt_wec_avg
      !------------------------------------------------



#endif /* WEC for whole module */
        
      end module wec_frc
