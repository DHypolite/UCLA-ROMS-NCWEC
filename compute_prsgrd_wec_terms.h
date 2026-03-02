! Compute kvf, calP and Kapsrf for WEC with u(nrhs)
! nrhs = nstp at predictor, nrhs=3 at corrector

! calP & Kapsrf: surface pressure & Bernoulli head at rho-point
! =============================================================
       do j=0,ny
          do i=0,nx
       !do j=j0,j1
       !  do i=i0,i1
       !do j=jstr,jend        ! kv (k dot v) at rho-point
       !   do i=istr,iend      ! and its 1st & 2nd derivertives
            do k=1,N          ! at rho-point
               kv(k) =0.5*kw(i,j)*(
     &           wdrx(i,j)* ( u(i,j,k,nrhs)+u(i+1,j,k,nrhs) )
     &           + wdre(i,j)* ( v(i,j,k,nrhs)+v(i,j+1,k,nrhs) ) )
            enddo

# if defined MONOWAVE
	    kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
# else
	    ! should this be done in MONOWAVE too????
	    kvsurf = 0.5*(ust_r(i,j,N)*( u(i,j,N,nrhs)+u(i+1,j,N,nrhs))
     &              + vst_r(i,j,N)*(v(i,j,N,nrhs)+v(i,j+1,N,nrhs) ))
# endif
            do k=1,N-1
               dkvdz(k) =2.0*(kv(k+1)-kv(k))/(Hz(i,j,k+1)+Hz(i,j,k))
            enddo
            dkvdz(0) = dkvdz(1)    !2.*dkvdz(1)-dkvdz(2) ! severe!
            dkvdz(N) = dkvdz(N-1)  !2.*dkvdz(N-1)-dkvdz(N-2)
            do k=1,N
              d2kv(k) =dkvdz(k)-dkvdz(k-1) ! d^2kv/dz^2 x Hz
            enddo
            cff3VF = 0.0  
            do k=1,N
               dd   = z_r(i,j,k)-z_w(i,j,N)
               cff3VF = cff3VF + d2kv(k)*(
     &             exp( 2.*kw(i,j)*(dd-Dstp(i,j)))
     &           + exp(-2.*kw(i,j)*(dd+Dstp(i,j))) )
            enddo
            !cff1VF =-2.0*exp(-2.*kD(i,j))*inv_ex(i,j)*dkvdz(N)
	    !Removed inv_ex array, now code directly into here
	    cff1VF = -2.0 * exp(-2.*kD(i,j)) * 1.0/max(1.0-exp(-4.*kD(i,j)),1e-10) * dkvdz(N) 
            cff2VF = dkvdz(0)/max(tanh(2.*kD(i,j)),1e-10)
            !cff3VF = cff3VF*inv_ex(i,j)
	    cff3VF = cff3VF * 1.0/max(1.0-exp(-4.*kD(i,j)),1e-10)
# if defined MONOWAVE
            cff4VF =-2.0*kw(i,j)*kvsurf
# else 
!	    ! will need to make these public in WEC!!!
!	    cff4VF = -2.0*(sinh(min(khmax,keff(i,j)*Dstp(i,j))))**2/
!     &                 cosh(min(2.*khmax,2.*keff(i,j)*h(i,j)))*kvsurf
! Two scale approximation Romero et al.
          cff4VF = (sinh(min(khmax,keff(i,j)*Dstp(i,j))))**2/
     &                 cosh(min(2.*khmax,2.*keff(i,j)*h(i,j)))*tanh(kD(i,j))+
     &            (sinh(min(khmax,kw(i,j)*Dstp(i,j))))**2/
     &                 cosh(min(2.*khmax,2.*kw(i,j)*h(i,j)))*(1-tanh(kD(i,j)))
          cff4VF = -2.0*kvsurf*cff4VF

# endif
# if defined SUP_OFF
            actp=actf(i,j)
# else
            actp=act(i,j)
# endif
# if defined MONOWAVE
            calP(i,j) = actp*tanh(kD(i,j))
     &                    *( cff1VF+cff2VF+cff3VF+cff4VF )
# else
            calP(i,j) = (actp*tanh(kD(i,j))
     &                    *( cff1VF+cff3VF) +act(i,j)*tanh(kD(i,j))*cff2VF + cff4VF) ! cff2 from total act, and cff4 computed from surface Ustk with u (wrongly named kvsurf)
# endif
#  ifdef MASKING
     &                                 *rmask(i,j)
#  endif
             cff5VF=0.0
             do k=1,N
                cff5VF = cff5VF +  d2kv(k)*
     &           ( exp( 2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)-Dstp(i,j)))
     &            -exp(-2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)+Dstp(i,j))) )
             enddo
             !Kapsrf(i,j) = cff5VF*actp*inv_ex(i,j)
             !Kapsrf(i,j) = cff5VF*actp* 1.0/max(1.0-exp(-4.*kD(i,j)),1e-10)
             Kapsrf(i,j) = cff5VF*act(i,j)* 1.0/max(1.0-exp(-4.*kD(i,j)),1e-10)


#  ifdef MASKING
     &                           *rmask(i,j)
#  endif



!
! kvf : vertical vortex force term (K term) at rho-point
! ======================================================
!
             do k=1,N-1       
	     !This is now done in the same manner for MONOWAVE or broad spectra
             !since we now define a ust_r, vst_r for MONOWAVE
             !kvr --> K term at horiz rho- and vert-w point
                !print *, 'ust_r = ', ust_r(i,j,k)
                kvr(k) = 0.25*(
     &                (ust_r(i,j,k)+ust_r(i,j,k+1))*
     &                (u(i,j,k+1,nrhs)-u(i,j,k,nrhs)
     &                        +u(i+1,j,k+1,nrhs)-u(i+1,j,k,nrhs))
     &              + (vst_r(i,j,k)+vst_r(i,j,k+1))*
     &                (v(i,j,k+1,nrhs)-v(i,j,k,nrhs)
     &                        +v(i,j+1,k+1,nrhs)-v(i,j+1,k,nrhs)))
     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
 

	     enddo

             ! Apply top and bottom B.C.s
	     ! for bottom --> kvr(k=1/2) = 0 --> ust=0 at sea-floor
             kvr(0)= 0.D0   
	     !Extrapolation from interior points for surface (k=N+1/2) 
             kvr(N)=2.*kvr(N-1)-kvr(N-2) 

	     !Move kvr to vert. rho-pt and place in public kvf shared array
             do k=1,N                            
	      
                kvf(i,j,k)=0.5*( kvr(k)+kvr(k-1) ) 
#  ifdef MASKING
     &                             *rmask(i,j)
#  endif
             enddo
          enddo !<---i
       enddo           ! <---- j



