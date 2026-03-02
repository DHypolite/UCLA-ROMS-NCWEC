      module mag
      ! -----------------------------------------------------------
      !                 ROMS module for:
      !                    MACROALGAE (=seaweed=kelp)
      !         
      !                 Originally coded following MAG model in Friedeer
      !                 et al. 2022 (FMARS), however, recently cleared
      !                 and made simpler to simulate seaweed CDR with
      !                 prescribed biomass and surface injections of
      !                 seaweed POC
      !                
      !                 Future revisions of this module will build back
      !                 in complexity, i.e., growth dynamics
      !            
      !             Coded by D. Dauhajre with help from P.
      !             Damien, J. Molemaker, and D. Dollery
      ! -----------------------------------------------------------

#include "cppdefs.opt"
! CPP flag MAG covering entire module so that the module is empty if
! MAG is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined MAG /* MAG for whole module! */


      use dimensions
      use nc_read_write
      use roms_read_write
      use netcdf
      use param
      use grid
      use ocean_vars
      use scalars
      use mixing
      use tracers
      use mpi_exchanges      

      implicit none


!  standard user-inputs (switches for forcing and writing)
#include "mag.opt"

!    constants that are not generally changed (e.g., drag coefficient)
#include "mag_params.opt"
      

      ! Make all variable private to this module unless public specified
      private
     
      ! save -> Make all public variables global
      ! (might not be needed depending on compiler)
      save

      ! -------------------------------------------------------------------
      !                     MAG VARIABLES
      !        ***********************************************
      ! sflx_pockelp_mass  | surface flux of kelp biomass (for sinking) | [g C / m^2 / day] 
      ! algae_b2d          | vert. integrated biomass | [g-dry / m^2]
      ! algae_b            | 3D biomass               | [g-dry / m^3]
      ! algae_lad          | leaf-area density        | [1/m]
      ! algae_shape        | vert. biomass shape func.| [1/m]
      ! algae_length       | length(height) of algae  | [m]
      ! algae_uptno3       | uptake rate by algae of no3 | [mmol N / m^3 / sec]
      ! algae_Nf           | fixed nitrogen of seaweed | [mg N / m^3]
      !       ************************************************
      ! ---------------------------------------------------------------------
      real,allocatable,dimension(:,:)          :: sflx_pockelp_mass
      real,allocatable,dimension(:,:)          :: algae_b2d
      real,allocatable,dimension(:,:,:)        :: algae_b
      real,allocatable,dimension(:,:,:),public :: algae_lad
      real,allocatable,dimension(:,:,:)        :: algae_shape
      real,allocatable,dimension(:,:)          :: algae_length
      real,allocatable,dimension(:,:,:),public :: algae_uptno3 !computed in bgc_ecosys_bec2.F
      real,allocatable,dimension(:,:,:),public :: algae_Nf !modifies PAR in bgc_ecosys_bec2.F (and eventually evolved with complex mag)

      !Averages for writing (same units as above)
      real,allocatable,dimension(:,:)          :: algae_b2d_avg
      real,allocatable,dimension(:,:,:)        :: algae_b_avg
      real,allocatable,dimension(:,:,:)        :: algae_lad_avg
      real,allocatable,dimension(:,:,:)        :: algae_shape_avg
      real,allocatable,dimension(:,:,:)        :: algae_uptno3_avg

      !surface flux of kelp biomass in units of [mmol C / m^2 / sec] 
      real,allocatable,dimension(:,:),public :: pockelp_prod_surf !used in bgc_ecosys_bec2.F

      ! sigma level defining bottom of algae (in (x,y))
      integer, allocatable,dimension(:,:) :: kbot 




      !Variables that are used locally in compute_horiz_rhs_uv_terms.h
      !for momentum sink
#ifdef ALGAE_DRAG
      ! Volume normalized drag term that is subtracted from ru (u-pts)
      real, allocatable, dimension (:,:,:),public :: vegdrag_u
      !Volume normalized drag term that is subtracted from rv (v-pts)
      real, allocatable, dimension (:,:,:),public :: vegdrag_v
      !Scratch variables to make code easier to read
      !         Eventually get rid of these for efficiency
      real, public :: vegcoeffu
      real, public  :: vegcoeffv
      real, public  :: vegvolu
      real, public  :: vegvolv
      real, public  :: v_ri1
      real, public  :: v_ri2
      real, public  :: u_rj1
      real, public :: u_rj2
      real, public :: velocmagu
      real, public :: velocmagv
#endif



      ! -------------------------------------------------------------------
      !
      !                      DECALRE THINGS FOR NETCDF OUTPUT
      !
      ! ---------------------------------------------------------------------
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Variable names and attributes for writing to netcdf
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !biomass
      character(len=1) :: b_his_name   = 'algae_b'
      character(len=7)  :: b_long_name  = 'biomass'
      character(len=10)  :: b_units = 'g-dry/m^3'

      !leaf-area density
      character(len=3) :: lad_his_name = 'algae_lad'
      character(len=17)   :: lad_long_name = 'leaf-area density'
      character(len=10) :: lad_units     = '1/m'

      !Vertical distribution function
      character(len=2) :: shape_his_name = 'algae_shape'
      character(len=20)   :: shape_long_name = 'biomass distribution'
      character(len=10) :: shape_units     = 'none'


      !Nitrate uptake
      character(len=10) :: uptno3_his_name = 'algae_uptno3'
      character(len=20)   :: uptno3_long_name = 'nitrate uptake'
      character(len=20) :: uptno3_units     = 'mmol N /( m^3 * sec)'


      ! Netcdf outputting for use in wrt_mag subroutine(s)
      ! -----------------
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_mag=0
      integer,save :: navg_mag = 0               ! number of samples in average
      integer :: record_avg = nrpf_avg           ! Triggers making of initial file
      integer :: record_his = nrpf_his
      real,save :: output_time_his = 0
      real,save :: output_time_avg = 0
     



      ! ---------------------
      ! Public functios
      ! ----------------------
      public set_mag_frc !called from set_forces.F
      public init_arrays_mag ! called from init_arrays.F
      public wrt_mag !called from main.F

      ! =======================
      ! END OF MODULE PRE-AMBLE
      ! =======================

      contains

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


      !Initialize MAG arrays (called from init_arrays.F)
      subroutine init_arrays_mag![
      implicit none

      !ncforce for prescribed forcing
      if (frc_biomass) allocate(nc_biomass%vdata(GLOBAL_2D_ARRAY,2) )
      if (frc_sflxpockelp) allocate(nc_sflxkelp%vdata(GLOBAL_2D_ARRAY,2) )
      
      !Read in forcing
      allocate( algae_b2d(GLOBAL_2D_ARRAY) )
      algae_b2d = 0.
      allocate( sflx_pockelp_mass(GLOBAL_2D_ARRAY) )
      sflx_pockelp_mass=0.


      !Standard MAG (algae) variables
      allocate( algae_b(GLOBAL_2D_ARRAY,N) )
      algae_b=0.
      allocate( algae_lad(GLOBAL_2D_ARRAY,N) )
      algae_lad=0.
      allocate( algae_length(GLOBAL_2D_ARRAY) )
      algae_length=0.
      allocate( algae_shape(GLOBAL_2D_ARRAY,N))
      algae_shape=0.
      allocate( algae_Nf(GLOBAL_2D_ARRAY,N))
      algae_Nf=0.

      !Averages for writing
      allocate( algae_b2d_avg(GLOBAL_2D_ARRAY) )
      algae_b2d_avg=0.
      allocate( algae_b_avg(GLOBAL_2D_ARRAY,N) )
      algae_b_avg=0.
      allocate( algae_lad_avg(GLOBAL_2D_ARRAY,N) )
      algae_lad_avg=0.
      allocate( algae_shape_avg(GLOBAL_2D_ARRAY,N)
      algae_shape_avg=0.

      !Nutrient uptake variables
      allocate( algae_uptno3(GLOBAL_2D_ARRAY,N) )
      algae_uptno3 = 0.
      allocate( algae_uptno3_avg(GLOBAL_2D_ARRAY,N) )
      algae_uptno3_avg = 0.

      !Converted POC kelp surface flux to BEC-compatible units
      allocate(pockelp_prod_surf(GLOBAL_2D_ARRAY))
      pockelp_prod_surf=0.


      !kbot
      allocate(kbot(GLOBAL_2D_ARRAY))
      kbot=0.

#ifdef ALGAE_DRAG
      allocate( vegdrag_u(GLOBAL_2D_ARRAY,N) )
      allocate( vegdrag_v(GLOBAL_2D_ARRAY,N) )
      vegdrag_u=0.
      vegdrag_v=0.
#endif


     


      end subroutine init_arrays_mag!]


      !----------------------------------------
      subroutine set_mag_frc ![
      ! Read in MAG forcing
      ! This is used to read in a biomass as a function of time 
      ! and/or a surface flux of POCkelp (to simulate biomass sinking)
      
      implicit none

      !Read in forcing from netcdf file
      if (frc_biomass) call set_frc_data(nc_biomass, algae_b2d, 'r')
      if (frc_sflxpockelp) call set_frc_data(nc_sflxkelp,sflx_pockelp_mass, 'r')

      !OR DEFINE FORCING ANALYTICALL IN SOME BLOCK OF CODE HERE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ----------------------------------------------------------------------
      !Convert MAG forcing to relevant quantities/units for ROMS-BEC coupling
      !-----------------------------------------------------------------------
      call conv_mag_frc

      end subroutine set_mag_frc!]

      !----------------------------------------

      
      subroutine conv_mag_frc(i,j)
      ! =========================================
      !  Convert mag forcing (biomass, sflx_pockelp)
      ! to what is needed for ROMS-BEC coupling 


      !                What this subroutine does
      !    1. converts 2D biomass to 3D biomass (computes algae_b)
      !          - computes algae_shape (shape function)
      !          - algae_b = algae_b2d * algae_shape
      !    2. computes leaf-area density
      !          - computes algae_lad, which needs algae_length
      !    3. converts surface flux of sinking biomass to BEC-compatible
      !       units
      ! =========================================
      implicit none
      ! inputs
      integer,intent(in) :: i,j
      integer :: k
      real :: conv


      !use i0,j0,i1,j1 as defined in init_dimensions in dimensions.F
      do j=j0,j1
         do i=i0,i1
            !Only compute kelp characteristics if there is biomass there
            if (algae_b2d(i,j).neq.0) then     
               
               !----------------------------------------------------
               !Define cultivation sigma level based on F_z_cult
               !(defined by user in mag.opt)
               !----------------------------------------------------
               call find_k(i,j,algae_zbot*-1,kbot(i,j))
               
               !----------------------
               !Calculate algae_length (based on z-integrated biomass)
               !-----------------------
               call calc_algaelength(i,j)

               ! -----------------------------------------------------
               ! Convert vert. integral of biomass to a vert. profile
               ! ----------------------------------------------------
               ! computes algae_shape(i,j)
               call calc_shape(i,j)
               !Now distribute biomass
               do k=1,N
                  algae_b(i,j,k) = algae_b2d(i,j) * algae_shape(i,j,k)
               enddo
               

               ! ---------------------------------------
               ! Convert biomass(z) to leaf-area density
               ! ---------------------------------------
               call calc_LAD(i,j)

               !---------------------------------
               ! Convert biomass to fixed nitrogen
               ! Nf used for PAR modulation
               !----------------------------------
               algae_Nf(i,j,:) = Qmin * algae_b(i,j,:)

            endif

            ! --------------------------------------
            !         Convert pockelp surface flux to 
            !         units of [mmol C / m^3 / sec]
            !         from     [g C / m^2 / day]
            ! --------------------------------------
            !Convert g C / m^2 / day to mmol C / m^2 / sec
                   ! g C --> mmol C   ! 1/day --> 1/sec
            conv = (1000. / 12.01) * (1./86400.)
            !divid by Hz at surface to make it a source term at surface and
            !give units of per m^3
            pockelp_prod_surf(i,j) = sflux_pockelp_mass(i,j)*conv / Hz(i,j,N)

        enddo !< -- i
      enddo !<-- j

      end subroutine conv_mag_frc 



      !======================================================
      ! 
      !                   Algae characteristic functions
      !
      !======================================================

      subroutine calc_shape(i,j)
      
      !------------------------------------
      !Calculate vert distribution function
      !------------------------------------
      
      implicit none
      integer :: k0,k_top
      real :: temp_shape(N)
      real :: temp_shape_int

      k0=kbot(i,j)
      !Get sigma level of top of plant
      call find_k(i,j,-algae_zbot + algae_length(i,j),k_top)
      
      !Make a dummy function that is constant in k
      temp_shape=0.
      do k=k0,k_top
         temp_shape(k)=1.
      enddo
      !Integrate temp_shape from kbot(kcult) up
      call ktrapz_kcult(i,j,temp_shape,temp_shape_int)
      !Normalize
      algae_shape(i,j) = temp_shape(:) / temp_shape_int
      end subroutine calc_shape

      subroutine calc_algaelength(i,j)
      !==============================================
      ! Calculate kelp length based on vertical integral of biomass
      !==============================================
      implicit none
      !inputs
      integer, intent(in):: i,j  
   
      algae_length(i,j) = (Hmax * algae_b2d(i,j) ) 
     &                   / (Kh + algae_b2d(i,j) )

      
      end subroutine calc_algaelength

      subroutine calc_LAD(i,j)
      !==============================================
      ! Calculate leaf-area density 
      !==============================================
      implicit none
      !inputs
      integer, intent(in):: i,j  

      !FILLER FUNCTION THAT NEEDS TO BE EXPANDED 

      !Convert  based on simple relationship that does not obey canopy
      ! LAD = Sa2b * B / 0.094, where dry_wet =0.094 and Sa2b is a
      ! surface-area to biomass relationship
      algae_lad(i,j,:) = (algae_b(i,j,:) / dry_wet) * Biomass_surfacearea_watercolumn

      end subroutine calc_LAD

      !======================================================
      ! 
      !                   UTILITY FUNCTIONS
      !
      !======================================================
      subroutine find_k(i,j,z_look,k_out)
      !================================
      !Find sigma level corresponding
      ! to a specified depth (z_look)
      !================================
      implicit none
      real, intent(in) :: z_look !should be negative
      integer, intent(in) :: i,j
      integer, intent(out) :: k_out
      !real :: diff(N)
      !diff = abs(z_r(i,j,:) - z_look)
      !k_out = minloc(diff,1)
      k_out=1
      do while (z_r(i,j,k_out)<z_look)
         k_out=k_out+1
      enddo

      end subroutine find_k


      subroutine ktrapz_kcult(i,j,var_in,var_out)
      ! Compute vertical integral from cultivation
      ! depth to surface
      implicit none
      integer, intent(in) :: i,j
      real, dimension(:),intent(in) :: var_in(N)
      real, intent(out) :: var_out
      integer :: k,k0
      k0=kbot(i,j) 
      var_out = 0
      call find_k(i,j,F_z_cult*-1,k0)
      do k=k0,N
         var_out=var_out
     &          + 0.5 * (var_in(k-1) + var_in(k))
     &          * (z_w(i,j,k) -z_w(i,j,k-1))
      enddo

      end subroutine ktrapz_kcult
      
      !======================================================





      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                 WRITING OUTPUT SUBROUTINES
      !                 these follow functionality in
      !                 surf_flux.F

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! NEW I/O Sept 2024
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !---------------------------------------------------------------------------
      subroutine wrt_mag  ![
      ! write surface flux variables to netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      character(len=99),save :: fname_his,fname_avg
      integer                :: ierr
      logical,save           :: first_step_avg=.true.
      logical,save           :: first_step_his=.true. 


      ! WRITE AVERAGE
      if (wrt_file_avg) then
           call calc_avg_mag                          ! don't include t=0 in averaging
           output_time_avg = output_time_avg + dt     ! only start count after first timestep
        if (output_time_avg>=output_period_avg) then  ! time for an output
          if (mod(record_avg,nrpf_avg)==0) then
             call create_mag_file(fname_avg,.true.)
             record_avg = 0
          endif
          record_avg = record_avg + 1
          output_time_avg = 0
          navg_mag = 0
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_avg/))
          if (mynode==0) then
             write(*,'(7x,A,1x,A)')
     &       'wrt_mag :: average'
          endif
          start=1; start(4)=record
          if (wrt_avg_b) call ncwrite(ncid,b_his_name,algae_b_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_shape) call ncwrite(ncid,shape_his_name,algae_shape_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_lad) call ncwrite(ncid,lad_his_name,algae_lad_avg(i0:i1,j0:j1,:),start)
          if (wrt_avg_uptno3) callncwrite(ncid,uptno3_his_name,algae_uptno3_avg(i0:i1,j0:j1,:),start)
          if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'mag :: wrote average, tdays =', tdays,
     &       'step =', iic, 'rec =', record_avg
          endif
      endif !-->output_time avg if
      endif !-->wrt_file_Avg
      
      ! WRITE HISTORY
      if (wrt_file_his) then
        if (.not. first_step_his) then
           output_time_his = output_time_his + dt
        endif
        first_step_his=.false.
        if (output_time_his>=output_period_his .or.
     &      output_time_his==0                 ) then  ! time for an output
          if (mod(record_his,nrpf_his)==0) then
             call create_mag_file(fname_his,.false.)
             record_his = 0
          endif
          record_his = record_his + 1
          output_time_his=0
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_his/))
          start=1; start(4)=record
          if (wrt_b) call ncwrite(ncid,b_his_name,algae_b(i0:i1,j0:j1,:),start)
          if (wrt_shape) call ncwrite(ncid,shape_his_name,algae_shape(i0:i1,j0:j1,:),start)
          if (wrt_lad) call ncwrite(ncid,lad_his_name,algae_lad(i0:i1,j0:j1,:),start)
          if (wrt_uptno3) callncwrite(ncid,uptno3_his_name,algae_uptno3(i0:i1,j0:j1,:),start)
 
          if (mynode == 0) then
             write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'mag :: wrote history, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record_his
          endif
      endif !-->wrt_file_his        
      end subroutine wrt_mag !]
      !---------------------------------------------------------------------------



      ! ----------------------------------------------------------------------
      subroutine create_mag_file(fname,avg)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ncid,ierr,varid

      if (avg) then
        call create_file('_mag_avg',fname)
      else
        call create_file('_mag',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call def_vars_mag( ncid,avg )

      ierr = nf90_enddef(ncid)

      end subroutine create_mag_file !]
      ! ----------------------------------------------------------------------


      ! ----------------------------------------------------------------------
      subroutine def_vars_mag( ncid, avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for actual or averaged variables
      implicit none

      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid, itrc
      character(len=64) :: text_lname



      !Biomass
      if ( ( wrt_b .and. .not. avg) .or. (wrt_avg_b.and. avg ) then
        if (.not. avg) text_lname=b_long_name
        if (      avg) text_lname='averaged '/ b_long_name
        varid =
     &  nccreate(ncid,b_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',b_units)
      endif

      !Leaf-area density
      if ( ( wrt_lad .and. .not. avg) .or. (wrt_avg_lad .and. avg ) then
        if (.not. avg) text_lname=lad_long_name
        if (      avg) text_lname='averaged '/ lad_long_name
        varid =
     &  nccreate(ncid,lad_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',lad_units)
      endif

      !Biomass shape function
      if ( ( wrt_shape .and. .not. avg) .or. (wrt_avg_shape .and. avg ) then
        if (.not. avg) text_lname=shape_long_name
        if (      avg) text_lname='averaged '/ shape_long_name
        varid =
     &  nccreate(ncid,shape_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',shape_units)
      endif


      !Nitrate Uptake
      if ( ( wrt_uptno3 .and. .not. avg) .or. (wrt_avg_uptno3 .and. avg ) then
        if (.not. avg) text_lname=uptno3_long_name
        if (      avg) text_lname='averaged '/ uptno3_long_name
        varid =
     &  nccreate(ncid,uptno3_his_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',uptno3_units)
      endif

      end subroutine def_vars_bgc  !]
      ! ----------------------------------------------------------------------

     ! ----------------------------------------------------------------------
      subroutine calc_avg_mag ![
      implicit none

      ! local
      real    :: coef

      if (wrt_file_avg) then 

      navg_mag = navg_mag +1
      coef = 1./navg_mag

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'mag :: started averaging. ',
     &   'output_period_avg (s) =', output_period_avg
      endif

      t_avg_mag = t_avg_mag*(1-coef) + time*coef
      if (wrt_avg_b) then 
         algae_b_avg(i0:i1,j0:j1,:) = algae_b_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_b(i0:i1,j0:j1,:) * coef
      endif

      if (wrt_avg_lad) then 
         algae_lad_avg(i0:i1,j0:j1,:) = algae_lad_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_lad(i0:i1,j0:j1,:) * coef
      endif
      if (wrt_avg_shape) then 
         algae_shape_avg(i0:i1,j0:j1,:) = algae_shape_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_shape(i0:i1,j0:j1,:) * coef
      endif
      if (wrt_avg_uptno3) then 
         algae_uptno3_avg(i0:i1,j0:j1,:) = algae_uptno3_avg(i0:i1,j0:j1,:) * (1-coef) 
     &                           + algae_uptno3(i0:i1,j0:j1,:) * coef
      endif
      endif !wrt_file_avg

      end subroutine calc_avg_mag !]
      !-----------------------------------------------------------------


#endif /* MAG for whole module */
      end module mag  
