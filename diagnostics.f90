MODULE diagnostics

USE parameters
USE mpi
USE fft
USE files
USE derivatives
USE special

IMPLICIT NONE

CONTAINS



     subroutine diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)   !Would be optimal to import 3 scratch for omega (zxk,zyk,zzk) and 3 others for grad b (txk,tyk,tzk)

       !In PV subroutine it would be great to send 3 more scratch arrays for grad(b)

       double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk
       double complex, dimension(iktx,ikty,n3h1) :: psik
       double complex, dimension(iktx,ikty,n3h1) :: u_rot,v_rot,b_rot,wak
       
       double complex, dimension(iktx,ikty,n3h1) :: zxk,zyk,zzk     !From scratch arrays?
       double precision,    dimension(n1d,n2d,n3h1) :: zxr,zyr,zzr

       real, dimension(n3h0) :: ks,ku,ps,ps_quad   !Staggered (s) and unstaggered (u) energy (k: kinetic, p: potential)                                                    
       real :: ktot_p,ktot,ptot_p,ptot,ptot_quad,ptot_quad_p

       real, dimension(n3h0) :: ks_rot,pu_rot   !Staggered (s) and unstaggered (u) energy (k: kinetic, p: potential)                                                    
       real :: ktot_rot_p,ktot_rot,ptot_rot_p,ptot_rot

       equivalence(zxk,zxr)
       equivalence(zyk,zyr)
       equivalence(zzk,zzr)

      !Compute velocity/buoyancy fields                                                                                                                                                                                           
      do izh1=1,n3h1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  u_rot(ikx,iky,izh1) =  - i*ky*psik(ikx,iky,izh1)
                  v_rot(ikx,iky,izh1) =    i*kx*psik(ikx,iky,izh1)
                  b_rot(ikx,iky,izh1) =    ( psik(ikx,iky,izh1+1) - psik(ikx,iky,izh1) )/(r_1(izh2)*dz)    !1/r_1 d psi/dz                                                                                                                
               else
                  u_rot(ikx,iky,izh1) =  (0.D0,0.D0)
                  v_rot(ikx,iky,izh1) =  (0.D0,0.D0)
                  b_rot(ikx,iky,izh1) =  (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

       
       !Explicit treatment of boundaries:  in QG, uz=vz=0 => b=0 at the top
       if(mype==(npe-1)) then
          do iky=1,ikty
             do ikx=1,iktx

                b_rot(ikx,iky,iztop1) = (0.D0,0.D0)

             end do
          end do
       end if



       !Compute k_H spectra at various heights!
       !**************************************!

       do iz=1,num_spec
          if(out_hspec ==1 .and. mod(iter,freq_hspec) ==0) call hspec(uk,vk,wk,bk,u_rot,v_rot,b_rot,iz)
       end do

       !Compute vorticity for other routines!
       !************************************!
      
!       if((out_ens ==1 .and. mod(iter,freq_ens) ==0) .or. (out_pv ==1 .and. mod(iter,freq_pv) ==0) .or. (iter==0 .and. ( out_ens==1 .or. out_pv==1)  ) ) then 
!          call vort(uk,vk,wk,zxk,zyk,zzk)
!          if(out_ens ==1 .and. ( mod(iter,freq_ens)==0 .or. iter==0 ))  call enstrophy(zxk,zyk,zzk)
!       end if

       !Conditions of integrability!
       !***************************!

       if(out_cond ==1 .and. (mod(iter,freq_cond ) ==0 .or. iter==0 ))  call cond_integrability(uk,vk,wk,bk)
       if(out_grow ==1 .and. (mod(iter,freq_grow ) ==0 .or. iter==0 ))  call compute_growth(uk,vk,bk)
       if(out_condwz==1 .and. (mod(iter,freq_condwz) ==0))                call cond_wz(wak)

       !Vertical scale of buoyancy and vertical velocity!
       !************************************************!

       if(out_vbuoy ==1 .and. (mod(iter,freq_vbuoy)  ==0))         call buoy_vert_scale(wk,bk)          !For full w and staggered buoyancy 
       if(out_vbuoyr==1 .and. (mod(iter,freq_vbuoyr) ==0))  then
          call generate_halo_q(b_rot)
          call buoy_vert_scale_rot(wak,b_rot)  !For ageo w and unstaggered QG buoyancy
       end if

       !Integrate over k_H to plot fields as a function of z!
       !****************************************************!

       ! 1/(2pi^)^2 int(int( e^2(x,y,z)/2 dx)dy = 1/2 sum over kx ky |eh|^2(z)                                                                                              
       
       ks=0.
       ku=0.
       ps=0.
       ps_quad=0.

       ktot_p=0.
       ktot=0.
       ptot_p=0.
       ptot=0.
       ptot_quad_p=0.
       ptot_quad=0.


       ks_rot=0.
       pu_rot=0.       

       ktot_rot_p=0.
       ktot_rot=0.
       ptot_rot_p=0.
       ptot_rot=0.

       !With dealiasing, sum_k 1/2 |u(kx,ky,z)|^2 = sum_k L |u|^2 - 0.5 |u(0,0,z)|^2                                                                             
       
       do iz=1,n3h0
          izh1=iz+1
          izh2=iz+2
          
          do iky=1,ikty
             do ikx=1,iktx
                
                if(L(ikx,iky)==1) then
                   ks(iz) = ks(iz) + real( uk(ikx,iky,izh2)*CONJG( uk(ikx,iky,izh2) ) + vk(ikx,iky,izh2)*CONJG( vk(ikx,iky,izh2) ) )
                   ku(iz) = ku(iz) + real( wk(ikx,iky,izh2)*CONJG( wk(ikx,iky,izh2) ) )*Ar2
                   ps_quad(iz) = ps_quad(iz) + real( bk(ikx,iky,izh2)*CONJG( bk(ikx,iky,izh2) ) )*(Bu*r_1s(izh2)/r_2s(izh2))         

                   ks_rot(iz) = ks_rot(iz) + real( u_rot(ikx,iky,izh1)*CONJG( u_rot(ikx,iky,izh1) ) + v_rot(ikx,iky,izh1)*CONJG( v_rot(ikx,iky,izh1) ) )
                   pu_rot(iz) = pu_rot(iz) + real( b_rot(ikx,iky,izh1)*CONJG( b_rot(ikx,iky,izh1) ) )*(Bu*r_1(izh2)/r_2(izh2))

                end if
                
             enddo
          enddo

          ps(iz) = -zash0(iz)*real(bk(1,1,izh2))/Ro          

          !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                                    
          
          ks(iz) = ks(iz) - 0.5*real( uk(1,1,izh2)*CONJG( uk(1,1,izh2) ) + vk(1,1,izh2)*CONJG( vk(1,1,izh2) ) )
          ku(iz) = ku(iz) - 0.5*real( wk(1,1,izh2)*CONJG( wk(1,1,izh2) ) )*Ar2
          ps_quad(iz) = ps_quad(iz) - 0.5*real( bk(1,1,izh2)*CONJG( bk(1,1,izh2) ) )*(Bu*r_1s(izh2)/r_2s(izh2))    

          ks_rot(iz) = ks_rot(iz) - 0.5*real( u_rot(1,1,izh1)*CONJG( u_rot(1,1,izh1) ) + v_rot(1,1,izh1)*CONJG( v_rot(1,1,izh1) ) )
          pu_rot(iz) = pu_rot(iz) - 0.5*real( b_rot(1,1,izh1)*CONJG( b_rot(1,1,izh1) ) )*(Bu*r_1(izh2)/r_2(izh2))

          
       end do

       

       !If desired, we can plot energy as a function of z and also the vertical energy spectrum                                                                                                           

       if(out_ez ==1   .and. ( mod(iter,freq_ez)==0 .or. iter == 0)) call plot_ez(ks,ku,ps)
       if(out_rotz ==1 .and. ( mod(iter,freq_rotz)==0 .or. iter==0)) call plot_rotz(ks_rot,pu_rot)              
!       if(out_specz ==1 .and. ( mod(iter,freq_specz)==0 .or. iter == 0)) call specz(ks,ku,pu)

       !Compute the total energy by integrating over z!
       !**********************************************!

       !1/L int(e(z) dz) becomes just a dummy sum over all grid points without interpolation needed...                                                           
       !First some locally to each processor                                                                                                                            
       
       do izh0=1,n3h0
          izh2=izh0+2
          
          ktot_p     =ktot_p      + rho_s(izh2)*ks(izh0) + rho_u(izh2)*ku(izh0)
          ptot_p     =ptot_p      + rho_s(izh2)*ps(izh0)
          ptot_quad_p=ptot_quad_p + rho_s(izh2)*ps_quad(izh0) 


          ktot_rot_p=ktot_rot_p + rho_s(izh2)*ks_rot(izh0) 
          ptot_rot_p=ptot_rot_p + rho_u(izh2)*pu_rot(izh0)      

          
       end do

       ktot_p     =ktot_p     /n3
       ptot_p     =ptot_p     /n3
       ptot_quad_p=ptot_quad_p/n3

       ktot_rot_p=ktot_rot_p/n3
       ptot_rot_p=ptot_rot_p/n3

       !Sum results from each processor                                                                                                            
       
       call mpi_reduce(ktot_p     ,ktot     , 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(ptot_p     ,ptot     , 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(ptot_quad_p,ptot_quad, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

       call mpi_reduce(ktot_rot_p,ktot_rot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(ptot_rot_p,ptot_rot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

       if(mype==0) write(unit=unit_energy ,fmt=*) time,ktot,ptot,ptot_quad
       if(mype==0) write(unit=unit_energyr,fmt=*) time,ktot_rot,ptot_rot
       
       
     end subroutine diag_zentrum



     subroutine normalize_trop(uk,vk,wk,bk,psik,qk,wak)   !This subroutine normalizes fields with the RMS horizontal velocity at the tropopause

       !In PV subroutine it would be great to send 3 more scratch arrays for grad(b)

       double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk
       double complex, dimension(iktx,ikty,n3h1) :: wak,psik,qk

       real :: horu    !Initial U_rms
       real :: norm    !Norm with which we normalize to get the new rms velo to be URMS.
       
       integer :: which_mype  !Which mype contains the level we normalize 

       integer :: iz0,proc

       !Find out which mype contains the grid point on which to compute U_RMS: all processors must know which it is for future broadcast
       do proc=1,npe
          if( trop_height > (proc-1)*n3h0 .AND. trop_height <= proc*n3h0 )  which_mype=proc-1       
       end do


       if(mype==which_mype) then

          iz0=trop_height - mype*n3h0 + 2  !2-halo velocity fields 

          !Compute U_RMS

          
          horu=0.
          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   horu = horu + real( uk(ikx,iky,iz0)*CONJG( uk(ikx,iky,iz0) ) + vk(ikx,iky,iz0)*CONJG( vk(ikx,iky,iz0) ) )
                end if

             enddo
          enddo

          !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                                                                                                     
          horu = horu - 0.5*real( uk(1,1,iz0)*CONJG( uk(1,1,iz0) ) + vk(1,1,iz0)*CONJG( vk(1,1,iz0) ) )

          horu = sqrt(horu)

          norm = horu/URMS

       end if

       !Broadcast the norm to all mypes
       call mpi_bcast(norm,1,MPI_REAL,which_mype,MPI_COMM_WORLD,ierror)

       !Not just normalize everybody
          uk=uk/norm
          vk=vk/norm
          wk=wk/norm
          bk=bk/norm
      psik=psik/norm
          qk=qk/norm
        wak=wak/norm

      end subroutine normalize_trop

     subroutine hspec(uk,vk,wk,bk,u_rot,v_rot,b_rot,level)
       
       double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk
       double complex, dimension(iktx,ikty,n3h1) :: u_rot,v_rot,b_rot
       
       integer :: level,proc,izh1s,izh2s  
       
       real, dimension(0:ktx,2) :: spz        !spectrum for kinetic (1) and potential (2)                                                                            
       real, dimension(0:ktx)   :: num_modes  !number of modes for a given kh                                                                                            
       integer :: mode                        !Mode no                                                                                                                      
       integer :: unit,unit_ro
       
       double precision :: kh
       
       !To compute Ro_macro: compute the average U_peak=sqrt(u^2+v^2) and L_peak, the scale at which energy peaks for each height (bot mid top) 

       real, dimension(0:ktx) :: U_peak_mode        !u^2 + v^2 at each kh mode.
       real                   :: U_peak             !U_peak_mode summed over all modes

       real :: L_peak                               !Actual energy bearing scale
       real :: L_peak_mem                           !Actual energy bearing scale (recalls the max result before regression, incase regression failed)
       real :: max_L_peak                           !Value of the spectrum at a given height (temp variable)

       double precision :: vx(npt),vy(npt),poly(3)

       !Use the concerned mype                                                                                                                                   
 
       if(height(level)<1 .or. height(level)>n3-1) write(*,*) "Problem with level in horizontal_spectrum"  !don't include n3, otherwise need a special treatment of boundary conditions...
   
       proc = (height(level)-1)/n3h0             !which processor                                                                                                                  
       izh1s = height(level) - proc*n3h0 + 1      !position in the processor (valid for n3h2 fields only)  
       izh2s = height(level) - proc*n3h0 + 2      !position in the processor (valid for n3h1 fields only)  
       
       if(mype==proc) then   
     
          !Which file to write on 
          if(level==1) then
             unit   =unit_h0
             unit_ro=unit_ro0
          elseif(level==2) then
             unit   =unit_h1
             unit_ro=unit_ro1
          elseif(level==3) then
             unit   =unit_h2
             unit_ro=unit_ro2
          elseif(level==4) then
             unit   =unit_h3
             unit_ro=unit_ro3
          elseif(level==5) then
             unit   =unit_h4
             unit_ro=unit_ro4
          elseif(level==6) then
             unit   =unit_h5
             unit_ro=unit_ro5
          elseif(level==7) then
             unit   =unit_h6
             unit_ro=unit_ro6
          elseif(level==8) then
             unit   =unit_h7
             unit_ro=unit_ro7
          elseif(level==9) then
             unit   =unit_h8
             unit_ro=unit_ro8
          elseif(level==10) then
             unit   =unit_h9
             unit_ro=unit_ro9
          else
             write(*,*) "Problem with level... can't find the file..."
          end if

!          if(level==top_height) then
!             unit=unit_htop
!             unit_ro=unit_ro_top
!          else if(level==mid_height) then
!             unit=unit_hmid
!             unit_ro=unit_ro_mid
!          else if(level==bot_height) then
!             unit=unit_hbot
!             unit_ro=unit_ro_bot
!          else
!             write(*,*) "Problem with level... can't find the file..."
!          end if
          
          spz=0.
          num_modes=0.
          U_peak=0.

          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2  = kx*kx+ky*ky
                kh   = sqrt(1.D0*kh2)
                
                mode = ifix(real(kh*L1/twopi+0.5))
                
                if (L(ikx,iky).eq.1) then
                   
                   spz(mode,1)   = spz(mode,1) + real( uk(ikx,iky,izh2s)*CONJG(uk(ikx,iky,izh2s)) )   !u on staggered grid                                        
                   spz(mode,1)   = spz(mode,1) + real( vk(ikx,iky,izh2s)*CONJG(vk(ikx,iky,izh2s)) )   !v on staggered grid 

                   U_peak = U_peak + spz(mode,1)
                  
                   spz(mode,1)   = spz(mode,1) + 0.5*Ar2*( real( wk(ikx,iky,izh2s)*CONJG(wk(ikx,iky,izh2s)) )  + real( wk(ikx,iky,izh2s-1)*CONJG(wk(ikx,iky,izh2s-1)) )  )  !w interpol
!                   spz(mode,2)   = spz(mode,2) + 0.5*( real( tk(ikx,iky,izh2s)*CONJG(tk(ikx,iky,izh2s)) )  + real( tk(ikx,iky,izh2s-1)*CONJG(tk(ikx,iky,izh2s-1)) )  )*(Bu*r_1s(izh2s)/r_2s(izh2s))  !t interpol
                    spz(mode,2)   = spz(mode,2) + real( bk(ikx,iky,izh2s)*CONJG(bk(ikx,iky,izh2s)) )*(Bu*r_1s(izh2s)/r_2s(izh2s))  !b is not interpolated (avail at stag points
                   num_modes(mode) = num_modes(mode) + 2
                   
                endif
             enddo
          enddo



          !Finish computing U_peak and L_peak, then print Ro_macro
          
          U_peak = U_peak - 0.5*( real( uk(1,1,izh2s)*CONJG(uk(1,1,izh2s))) + real( vk(1,1,izh2s)*CONJG(vk(1,1,izh2s)) ))   !To be very precise, substract half of kh=0 mode...
          U_peak = sqrt(U_peak)  !Since sum_kh L |u|^2 - 0.5* |u(kh=0)|^2 = 0.5*sum |u|^2
          
          !Method 1: easy but leading to discontinuity in L_peak
          
          L_peak=0.
          max_L_peak=0.
          
          do mode=0,ktx-1
             if(spz(mode,1)>max_L_peak) then
                max_L_peak=spz(mode,1)
                L_peak=mode
             end if
          end do
          
          
          if(parabolic_L_peak==1) then
             !Fits parabola with npt points.
             L_peak_mem=L_peak
             if(L_peak<((npt-1)/2+1)) then
                do mode=1,npt
                   vx(mode)=mode
                end do
             elseif(L_peak>=((npt-1)/2+1) ) then
                do mode=1,npt
                   vx(mode)=L_peak-((npt-1)/2+1)+mode
                end do
             end if
             do mode=1,npt
                vy(mode)=spz(vx(mode),1)
             end do
             poly =  polyfit(vx,vy,2)
             L_peak = -poly(2)/(2.*poly(3))  !Maximum of the parabola a + bx + cx^2 is x_max = -b/2c                                                                                                                           
             
             !Make sure the value is sensible, that is, it is contained within L_peak_mem +/- 1
             if(L_peak>L_peak_mem+((npt-1)/2) .or. L_peak<L_peak_mem-((npt-1)/2)) L_peak=L_peak_mem 
             
             
          end if
          
          
          L_peak=L1/L_peak   !L_peak = 2pi/mode for which spectrum(mode) is max)
          !Print Ro_macro = Ro*U_peak/L_peak
          write(unit_ro,fmt=*) time,Ro*U_peak/L_peak,U_peak,L_peak

          
          spz=0.5*rho_s(izh2s)*spz
          
          do mode=0,ktx-1
             if (num_modes(mode).ne.0) then         !mode, kin energy, pot energy, num of modes                                                               
                
                write(unit,fmt=*) float(mode),spz(mode,1),spz(mode,2),num_modes(mode)    !There was a problem here.[r_1/r_2 multiplied a second time to spz...]
             endif

          enddo
          write(unit,*) '           '
          call flush(unit)

   if(out_hg == 1) then
       !Now, plot the rot/div spectrum
       !******************************

       !Which file to write on                                                                                                                                             
       
          if(level==1) then
             unit   =unit_hg0
          elseif(level==2) then
             unit   =unit_hg1
          elseif(level==3) then
             unit   =unit_hg2
          elseif(level==4) then
             unit   =unit_hg3
          elseif(level==5) then
             unit   =unit_hg4
          elseif(level==6) then
             unit   =unit_hg5
          elseif(level==7) then
             unit   =unit_hg6
          elseif(level==8) then
             unit   =unit_hg7
          elseif(level==9) then
             unit   =unit_hg8
          elseif(level==10) then
             unit   =unit_hg9
          else
             write(*,*) "Problem with level... can't find the file..."
          end if

!       if(level==top_height) then
!          unit=unit_htopg
!       else if(level==mid_height) then
!          unit=unit_hmidg
!       else if(level==bot_height) then
!          unit=unit_hbotg
!       else
!          write(*,*) "Problem with level... can't find the file..."
!       end if
       

     spz=0.
     num_modes=0.

     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           kh2  = kx*kx+ky*ky
           kh   = sqrt(1.D0*kh2)

           mode = ifix(real(kh*L1/twopi+0.5))

           if (L(ikx,iky).eq.1) then

              !Assuming energy_rot goes like u_r^2 + v_r^2 + b_r^2
              spz(mode,1)   = spz(mode,1) + real( u_rot(ikx,iky,izh1s)*CONJG(u_rot(ikx,iky,izh1s)) )   !u on staggered grid                                        
              spz(mode,1)   = spz(mode,1) + real( v_rot(ikx,iky,izh1s)*CONJG(v_rot(ikx,iky,izh1s)) )   !v on staggered grid    
!              spz(mode,1)   = spz(mode,1) + 0.5*( real( b_rot(ikx,iky,izh1s)*CONJG(b_rot(ikx,iky,izh1s)) )  + real( b_rot(ikx,iky,izh1s-1)*CONJG(b_rot(ikx,iky,izh1s-1)) )  )*(Bu*r_1s(izh2s)/r_2s(izh2s)) 
              !Bu was forgotten in the line above

              !Assuming energy_div goes like u_d^2 + v_d^2 + w_d^2 + b_d^2
              spz(mode,2)   = spz(mode,2) + real( (uk(ikx,iky,izh2s) - u_rot(ikx,iky,izh1s))*CONJG( uk(ikx,iky,izh2s) - u_rot(ikx,iky,izh1s)   ) )   !u on staggered grid 
              spz(mode,2)   = spz(mode,2) + real( (vk(ikx,iky,izh2s) - v_rot(ikx,iky,izh1s))*CONJG( vk(ikx,iky,izh2s) - v_rot(ikx,iky,izh1s)   ) )   !v on staggered grid 
              spz(mode,2)   = spz(mode,2) + 0.5*( real( wk(ikx,iky,izh2s)*CONJG(wk(ikx,iky,izh2s)) )  + real( wk(ikx,iky,izh2s-1)*CONJG(wk(ikx,iky,izh2s-1)) )  )*Ar2  !w_div=wk  
!              spz(mode,2)   = spz(mode,2) + 0.5*( real(( tk(ikx,iky,izh2s) - b_rot(ikx,iky,izh1s)  )*CONJG( tk(ikx,iky,izh2s) - b_rot(ikx,iky,izh1s) ) )  + real((tk(ikx,iky,izh2s-1) - b_rot(ikx,iky,izh1s-1) )*CONJG( tk(ikx,iky,izh2s-1) - b_rot(ikx,iky,izh1s-1) ) )  )*(Bu*r_1s(izh2s)/r_2s(izh2s))  !t interpolated  


              num_modes(mode) = num_modes(mode) + 2

           endif
        enddo
     enddo

     spz=0.5*rho_s(izh2s)*spz

     do mode=0,ktx-1
        if (num_modes(mode).ne.0) then         !mode, rot kinetic energy, div kinetic energy
           write(unit,fmt=*) float(mode),spz(mode,1),spz(mode,2),num_modes(mode)
        endif

     enddo
     write(unit,*) '           '
     call flush(unit)
  end if! if(out_hg ==1)
  end if

end subroutine hspec









     subroutine hspec_waves(BRk,BIk,CRk,CIk,level)
       
       double complex, dimension(iktx,ikty,n3h0) :: BRk,BIk
       double complex, dimension(iktx,ikty,n3h0) :: CRk,CIk
       
       integer :: level,proc,izh0s,izh1s,izh2s  
       
       real, dimension(0:ktx,2) :: spz        !spectrum for kinetic (1) and potential (2)                                                                            
       real, dimension(0:ktx)   :: num_modes  !number of modes for a given kh                                                                                            
       integer :: mode                        !Mode no                                                                                                                      
       integer :: unit
       
       double precision :: kh
       
       !Use the concerned mype                                                                                                                                   
 
       if(height(level)<1 .or. height(level)>n3) write(*,*) "Problem with level in horizontal_spectrum"  !don't include n3, otherwise need a special treatment of boundary conditions...
   
       proc  = (height(level)-1)/n3h0             !which processor     
       izh0s = height(level) - proc*n3h0          !position in the processor (valid for n3h0 fields only)  
       izh1s = height(level) - proc*n3h0 + 1      !position in the processor (valid for n3h1 fields only)  
       izh2s = height(level) - proc*n3h0 + 2      !position in the processor (valid for n3h2 fields only)  

       
       
       if(mype==proc) then   
     
          !Which file to write on 
          if(level==1) then
             unit   =unit_h0w
          elseif(level==2) then
             unit   =unit_h1w
          elseif(level==3) then
             unit   =unit_h2w
          elseif(level==4) then
             unit   =unit_h3w
          elseif(level==5) then
             unit   =unit_h4w
          elseif(level==6) then
             unit   =unit_h5w
          elseif(level==7) then
             unit   =unit_h6w
          elseif(level==8) then
             unit   =unit_h7w
          elseif(level==9) then
             unit   =unit_h8w
          elseif(level==10) then
             unit   =unit_h9w
          else
             write(*,*) "Problem with level... can't find the file..."
          end if

          spz=0.
          num_modes=0.

          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2  = kx*kx+ky*ky
                kh   = sqrt(1.D0*kh2)
                
                mode = ifix(real(kh*L1/twopi+0.5))
                
                if (L(ikx,iky).eq.1) then
                   
                   spz(mode,1)   = spz(mode,1) + real( (BRk(ikx,iky,izh0s)+i*BIk(ikx,iky,izh0s))  *CONJG(BRk(ikx,iky,izh0s)+i*BIk(ikx,iky,izh0s)) )   !KE = 0.5 |LA|^2     

                   if(height(level) .eq. 1) then   !Az = 0 at the bot, so we take 1/2 the value of the grid point below.
                      spz(mode,2)   = spz(mode,2) + 0.25*(  1/(Bu*r_2(izh2s)))*kh2*real( CRk(ikx,iky,izh0s)*CONJG( CRk(ikx,iky,izh0s) ) + CIk(ikx,iky,izh0s)*CONJG( CIk(ikx,iky,izh0s)) )
                   else   !We interpolate to get spectra at staggered grid points.
                      spz(mode,2)   = spz(mode,2) + 0.25*( (1/(Bu*r_2(izh2s-1)))*kh2*real( CRk(ikx,iky,izh0s-1)*CONJG( CRk(ikx,iky,izh0s-1) ) + CIk(ikx,iky,izh0s-1)*CONJG( CIk(ikx,iky,izh0s-1) )) + (1/(Bu*r_2(izh2s)))*kh2*real( CRk(ikx,iky,izh0s)*CONJG( CRk(ikx,iky,izh0s) ) + CIk(ikx,iky,izh0s)*CONJG( CIk(ikx,iky,izh0s) ))   )
                   end if

                   num_modes(mode) = num_modes(mode) + 2
                   
                endif
             enddo
          enddo

          spz=0.5*spz

          do mode=0,ktx-1
             if (num_modes(mode).ne.0) then         !mode, kin energy, pot energy, num of modes                             
                write(unit,fmt=*) float(mode),spz(mode,1),spz(mode,2),num_modes(mode)    
             endif
          enddo
          write(unit,*) '           '
          call flush(unit)

       end if
     end subroutine hspec_waves








     subroutine wave_energy(ARk,AIk,BRk,BIk,CRk,CIk)
       
       !Computes the volume integrated waves energy: WKE = int( 0.5* |LA|^2 dV ) = Uw_scale^2 int( 0.5* |LA|^2 dV ) in nondim form.
       !Computes the volume integrated waves potential energy: WPE = int( 0.25* (f/N)^2 |grad(Az)|^2 dV ). I defined a field C=Az to simplifiy.
       !Nondim: WPE = int( 0.25* Uw_scalle^2/Bu 1/N'^2(z) |grad(Az)|^2 dV )   
       !Then, WPE = int( 0.25 (f/N)^2 sum_k kh^2 [ |CRk|^2 + |CIk|^2 ] dz

       !In the real space, WKE = int( 0.5* |B|^2 dV ) = int( 0.5*( Br^2 + Bi^2 ) dV ) which is analoguous to KE ~ int ( 0.5* (u^2 + v^2) dV)
       !The integral in real space can be converted to a sum on horizontal wavenumbers at each vertical level:
       ! 
       ! KE(z) = sum over all k's such that L=1 of |uk|^2 minus half the kh=0 mode.

       double complex, dimension(iktx,ikty,n3h0) :: ARk,AIk
       double complex, dimension(iktx,ikty,n3h0) :: BRk,BIk
       double complex, dimension(iktx,ikty,n3h0) :: CRk,CIk

       real, dimension(n3h0) :: k_p
       real, dimension(n3h0) :: p_p
       real, dimension(n3h0) :: c_p      !Correction to K: 1/16 nabla A
       
       real :: ktot_p,ktot
       real :: ptot_p,ptot
       real :: ctot_p,ctot

       real :: k0tot_p,k0tot  !Kinetic energy in the kh=0 mode

       ktot_p = 0.
       ptot_p = 0.
       ctot_p = 0.

       k0tot_p = 0.

       k_p = 0.
       p_p = 0.
       c_p = 0.

       do izh0=1,n3h0
          izh2=izh0+2

          do ikx=1,iktx
             kx=kxa(ikx)
             do iky=1,ikty
                ky=kya(iky)
                kh2=kx*kx+ky*ky

                if(L(ikx,iky)==1) then
                   k_p(izh0) = k_p(izh0) + real( BRk(ikx,iky,izh0)*CONJG( BRk(ikx,iky,izh0) ) + BIk(ikx,iky,izh0)*CONJG( BIk(ikx,iky,izh0) ) )
                   p_p(izh0) = p_p(izh0) + (0.5/(r_2(izh2)*Bu))*kh2*real( CRk(ikx,iky,izh0)*CONJG( CRk(ikx,iky,izh0) ) + CIk(ikx,iky,izh0)*CONJG( CIk(ikx,iky,izh0) ) )
                   c_p(izh0) = c_p(izh0) + (1./8.)*(1./(Bu*Bu))*kh2*kh2*real( ARk(ikx,iky,izh0)*CONJG( ARk(ikx,iky,izh0) ) + AIk(ikx,iky,izh0)*CONJG( AIk(ikx,iky,izh0) ) )
                end if


             enddo
          enddo

          !With dealiasing, sum_k 1/2 |u(kx,ky,z)|^2 = sum_k L |u|^2 - 0.5 |u(0,0,z)|^2     
          k_p(izh0) = k_p(izh0) - 0.5*real( BRk(1,1,izh0)*CONJG( BRk(1,1,izh0) ) + BIk(1,1,izh0)*CONJG( BIk(1,1,izh0) ) )

          !Sum local to the processor
          ktot_p = ktot_p + k_p(izh0)
          ptot_p = ptot_p + p_p(izh0)
          ctot_p = ctot_p + c_p(izh0)

          k0tot_p = k0tot_p + 0.5*real( BRk(1,1,izh0)*CONJG( BRk(1,1,izh0) ) + BIk(1,1,izh0)*CONJG( BIk(1,1,izh0) ) )

       end do

       !--- Vertical slices ---!

       !Get the dimensional values of WKE, WPE, and inverse wave Richardson number (this last is actually nondimensional...)
       k_p = Uw_scale*Uw_scale*k_p
       p_p = Uw_scale*Uw_scale*p_p
       c_p = Uw_scale*Uw_scale*c_p

       !If desired, we can plot these quantities as a function of z
       if(out_wz ==1   .and. ( mod(iter,freq_wz)==0 .or. iter == 0)) call plot_wz(k_p,p_p,c_p)


       !--- Volume integrals ---!

       !Sum results from each processor                                                                                                                                                 
       call mpi_reduce(ktot_p,ktot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(ptot_p,ptot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(ctot_p,ctot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)

       call mpi_reduce(k0tot_p,k0tot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)

       ktot = Uw_scale*Uw_scale*ktot/n3
       ptot = Uw_scale*Uw_scale*ptot/n3
       ctot = Uw_scale*Uw_scale*ctot/n3

       k0tot = Uw_scale*Uw_scale*k0tot/n3

       if(mype==0) write(unit_we,fmt=*) time*(L_scale/U_scale)/(3600*24),ktot,k0tot
       if(mype==0) write(unit_ce,fmt=*) time*(L_scale/U_scale)/(3600*24),ptot,ctot

     end subroutine wave_energy





     !Calculate the vertically-averaged wave energy (only total wave energy for now)
     subroutine we_vave(BRk,BIk,BRr,BIr)

       double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk 
       double precision, dimension(n1d,n2d,n3h0)   :: BRr, BIr

       double complex,   dimension(iktx,ikty,n3h0) :: BRmem, BImem     !Save temporarily in Xmem to avoid fft back

       real, dimension(n1,n2) :: WE       !Total wave energy density 0.5 |LA+|^2
       real, dimension(n1,n2) :: WE_ave   !Total wave energy density 0.5 |LA+|^2, vertically averaged

       character(len = 32) :: fname                !future file name                                                                                                           


       !Save onto Bmem to avoid fft back
       BRmem=BRk
       BImem=BIk

       call fft_c2r(BRk,BRr,n3h0)
       call fft_c2r(BIk,BIr,n3h0)

       !Initialize WE to zero
       WE = 0.
       WE_ave=0.

       !Compute the total wave energy:
       do izh0=1,n3h0
          do ix=1,n1
             do iy=1,n2
                WE(ix,iy) = WE(ix,iy) + 0.5*(BRr(ix,iy,izh0)*BRr(ix,iy,izh0) + BIr(ix,iy,izh0)*BIr(ix,iy,izh0))
             end do
          end do
       end do

       !Recover the k-space versions of B
       BRk=BRmem
       BIk=BImem


       !Perform the vertical integral:
       call mpi_reduce(WE,WE_ave,n1*n2,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)


       !Now in mype==0, print slice
       if(mype==0) then

          !Normalize and make dimensional
          WE_ave = WE_ave*Uw_scale*Uw_scale/n3

          write (fname, "(A7,I3,A4)") "WE_vave",count_vave,".dat"
          open (unit=unit_wvave,file=fname,action="write",status="replace")
          
          do iy=1,n2
             write(unit=unit_wvave,fmt=333) (WE_ave(ix,iy),ix=1,n1)
             write(unit=unit_wvave,fmt=*) '           '
          enddo
333       format(1x,E12.5,1x)
          
          close (unit=unit_wvave)

          count_vave=count_vave+1
       end if


     end subroutine we_vave


SUBROUTINE plot_wz(ks,ku,ps)    !Exact copy of plot_ez (I just changed the name of the file and the mpi tags)

   real, dimension(n3h0) :: ks,ku,ps         !Staggered (s) and unstaggered (u) energy         
   real, dimension(n3h0) :: ks_r,ku_r,ps_r   !Copies of eu and es to keep es and eu valid for energy subroutine.

   real,dimension(n3)   :: kuz,ksz,psz             !energy(z) for entire domain
   
   integer :: processor,izp,nrec
   
   ksz=0.
   kuz=0.
   psz=0.

    !Send err_p from other processors to mype = 0                                                                                                                    
    if(mype>0) call mpi_send(ks,n3h0,MPI_REAL,0,tag_kzs2,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(ku,n3h0,MPI_REAL,0,tag_kzu2,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(ps,n3h0,MPI_REAL,0,tag_pzs2,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part onto err                  
       do iz=1,n3h0
             ksz(iz) = ks(iz)
             kuz(iz) = ku(iz)
             psz(iz) = ps(iz)
       end do


       !Receive other parts from other processors    
       do nrec=1,npe-1

          call mpi_recv(ks_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_kzs2,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                       
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             ksz(izp) = ksz(izp) + ks_r(iz)

          end do

          !Kinetic unstag energy part
          call mpi_recv(ku_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_kzu2,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             kuz(izp) = kuz(izp) + ku_r(iz)

          end do

          
          !Potential unstag energy part
          call mpi_recv(ps_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_pzs2,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             psz(izp) = psz(izp) + ps_r(iz)

          end do

       end do

       !PLOT!
     do iz=1,n3
        write(unit_wz,fmt=*) zas(iz)*H_scale,ksz(iz),kuz(iz),psz(iz)      !iz,URMS, WRMS, BRMS
     enddo
     write(unit_wz,*) '           '
     call flush(unit_wz)

   end if


 END SUBROUTINE plot_wz







     subroutine wave_energy_original(BRk,BIk,CRk,CIk)
       
       !Computes the volume integrated waves kinetic energy: WKE = int( 0.5* |LA|^2 dV ) = Uw_scale^2 int( 0.5* |LA|^2 dV ) in nondim form.
       !Computes the volume integrated waves potential energy: WPE = int( 0.25* (f/N)^2 |grad(Az)|^2 dV ). I defined a field C=Az to simplifiy.
       !Nondim: WPE = int( 0.25* Uw_scalle^2/Bu 1/N'^2(z) |grad(Az)|^2 dV )
       
       !Then, WPE = int( 0.25 (f/N)^2 sum_k kh^2 [ |CRk|^2 + |CIk|^2 ] dz


       !In the real space, WKE = int( 0.5* |B|^2 dV ) = int( 0.5*( Br^2 + Bi^2 ) dV ) which is analoguous to KE ~ int ( 0.5* (u^2 + v^2) dV)
       !The integral in real space can be converted to a sum on horizontal wavenumbers at each vertical level:
       ! 
       ! KE(z) = sum over all k's such that L=1 of |uk|^2 minus half the kh=0 mode.

       
       double complex, dimension(iktx,ikty,n3h0) :: BRk,BIk
       double complex, dimension(iktx,ikty,n3h0) :: CRk,CIk

       real :: k_p, ktot
       real :: p_p, ptot

       k_p = 0.
       p_p = 0.

       do izh0=1,n3h0
          izh2=izh0+2
          do ikx=1,iktx
             kx=kxa(ikx)
             do iky=1,ikty
                ky=kya(iky)
                kh2=kx*kx+ky*ky

                if(L(ikx,iky)==1) then
                   k_p = k_p + real( BRk(ikx,iky,izh0)*CONJG( BRk(ikx,iky,izh0) ) + BIk(ikx,iky,izh0)*CONJG( BIk(ikx,iky,izh0) ) )
                   p_p = p_p + 0.5*(1/(Bu*r_2(izh2)))*kh2*real( CRk(ikx,iky,izh0)*CONJG( CRk(ikx,iky,izh0) ) + CIk(ikx,iky,izh0)*CONJG( CIk(ikx,iky,izh0) ) )
                end if


             enddo
          enddo

          !With dealiasing, sum_k 1/2 |u(kx,ky,z)|^2 = sum_k L |u|^2 - 0.5 |u(0,0,z)|^2     
          k_p = k_p - 0.5*real( BRk(1,1,izh0)*CONJG( BRk(1,1,izh0) ) + BIk(1,1,izh0)*CONJG( BIk(1,1,izh0) ) )

       end do

       k_p = Uw_scale*Uw_scale*k_p/n3
       p_p = Uw_scale*Uw_scale*p_p/n3

       !Sum results from each processor                                                                                                                                                 
       call mpi_reduce(k_p,ktot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(p_p,ptot,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)

       if(mype==0) write(unit_we,fmt=*) time,ktot,ptot

     end subroutine wave_energy_original


     subroutine we_conversion(ARk, AIk, nBRk, nBIk, rBRk, rBIk, nBRr, nBIr, rBRr, rBIr)

       !Computes the wave potential energy conversion terms:
       !\Gamma_a = 0.25 int( nabla^2 A* J(psi,LA) + nabla^2 A J(psi,LA)* )dV = 0.5 int ( DR JR + DI JI ) dV where DR = Re(nabla^2 A) and JI = Im( J(psi,LA)   ) etc.
       !\Gamma_r = 0.25 int( nabla^2 A* R         + nabla^2 A R*         )dV = 0.5 int ( DR RR + DI RI ) dV where DR = Re(nabla^2 A) and RI = Im( i zeta LA/2 ) etc.

       double complex,   dimension(iktx,ikty,n3h0) :: ARk, AIk

       double complex,   dimension(iktx,ikty,n3h0) :: nBRk, nBIk, rBRk, rBIk
       double precision, dimension(n1d,n2d,n3h0)   :: nBRr, nBIr, rBRr, rBIr

       !Nabla^2 A --> -kh2
       double complex,   dimension(iktx,ikty,n3h0) :: nARk, nAIk
       double precision, dimension(n1d,n2d,n3h0)   :: nARr, nAIr

       !Store REF and ADV terms temporarily (to avoid fft back)
       double complex,   dimension(iktx,ikty,n3h0) :: tRk, tIk

       !Conversion terms for advection and refraction.
       real ::  ca_p,cr_p,ca_tot,cr_tot 

       equivalence(nARk,nARr)
       equivalence(nAIk,nAIr)

       ca_p = 0.
       cr_p = 0.
       ca_tot = 0.
       cr_tot = 0.

       !Start with the advective conversion term
       tRk = nBRk
       tIk = nBIk

       !Compute nabla^2 A == nA
       do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               
               nARk(ikx,iky,izh0) =  - kh2*ARk(ikx,iky,izh0)
               nAIk(ikx,iky,izh0) =  - kh2*AIk(ikx,iky,izh0)
               
            enddo
         enddo
      enddo
      
      !FFT advection and nA to real-space to compute the conversion term
      call fft_c2r(nARk,nARr,n3h0)
      call fft_c2r(nAIk,nAIr,n3h0)

      call fft_c2r(nBRk,nBRr,n3h0)
      call fft_c2r(nBIk,nBIr,n3h0)

      !Compute the local integral for advection conversion
      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   ca_p = ca_p + nARr(ix,iy,izh0)*nBRr(ix,iy,izh0) + nAIr(ix,iy,izh0)*nBIr(ix,iy,izh0)

                end if
             end do
          end do
       end do

       !Recover k-space advective terms
       nBRk = tRk
       nBIk = tIk
   

       !Next: refractive terms
       tRk = rBRk
       tIk = rBIk


       !FFT refraction to real-space to compute the conversion term
       call fft_c2r(rBRk,rBRr,n3h0)
       call fft_c2r(rBIk,rBIr,n3h0)

       !Compute the local integral for refraction conversion
       do izh0=1,n3h0
          do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then
                   
                   cr_p = cr_p + 0.5*( -nARr(ix,iy,izh0)*rBIr(ix,iy,izh0) + nAIr(ix,iy,izh0)*rBRr(ix,iy,izh0) )    !Because I define rB = zeta LA, and refraction = (i/2) zeta LA 
                   
                end if
             end do
          end do
       end do

       !Recover k-space advective terms
       rBRk = tRk
       rBIk = tIk

       !Sum over all processors
       call mpi_reduce(ca_p,ca_tot, 1,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(cr_p,cr_tot, 1,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierror)

       !Normalize
       ca_tot = ca_tot*Uw_scale*Uw_scale/(2.*n1*n2*n3*Bu)
       cr_tot = cr_tot*Uw_scale*Uw_scale/(2.*n1*n2*n3*Bu)

       if(mype==0) write(unit=unit_conv ,fmt=*) time,ca_tot,cr_tot

     end subroutine we_conversion



    subroutine gamma_conversion(ARk, AIk, BRk, BIk, nBRk, nBIk, rBRk, rBIk, nBRr, nBIr, rBRr, rBIr)

       !Computes the wave energy conversion terms: SIGNS ARE INVERTED FOR GAMMA_R AND GAMMA_A COMPARED WITH WE_ AND WPE_CONVERSION! ALSO EVERYBODY IS DIMENSIONAL.                   
       !\Gamma_a = - 0.25 int( nabla^2 A* J(psi,LA) + nabla^2 A J(psi,LA)* )dV = - 0.5 int ( DR JR + DI JI ) dV where DR = Re(nabla^2 A) and JI = Im( J(psi,LA)   ) etc.             
       !\Gamma_r = - 0.25 int( nabla^2 A* R         + nabla^2 A R*         )dV = - 0.5 int ( DR RR + DI RI ) dV where DR = Re(nabla^2 A) and RI = Im( i zeta LA/2 ) etc.            

       double complex,   dimension(iktx,ikty,n3h0) :: ARk, AIk
       double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk

       double complex,   dimension(iktx,ikty,n3h0) :: nBRk, nBIk, rBRk, rBIk
       double precision, dimension(n1d,n2d,n3h0)   :: nBRr, nBIr, rBRr, rBIr

       !Nabla^2 A --> -kh2                                                                                                                                               
       double complex,   dimension(iktx,ikty,n3h0) :: nARk, nAIk
       double precision, dimension(n1d,n2d,n3h0)   :: nARr, nAIr

       !Store REF and ADV terms temporarily (to avoid fft back)                                                                                                          
       double complex,   dimension(iktx,ikty,n3h0) :: tRk, tIk
       double precision, dimension(n1d,n2d,n3h0)   :: tRr, tIr

       !Conversion terms for advection, refraction, forcing and dissipation                                                                                               
       real ::  ca_p, ca_tot
       real ::  cr_p, cr_tot
       real ::  cd_p, cd_tot

       equivalence(nARk,nARr)
       equivalence(nAIk,nAIr)

       equivalence(tRk,tRr)
       equivalence(tIk,tIr)

       ca_p = 0.
       cr_p = 0.
       cd_p = 0.

       ca_tot = 0.
       cr_tot = 0.
       cd_tot = 0.


       !Compute nabla^2 A == nA                                                                                                                                        
       do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky

               nARk(ikx,iky,izh0) =  - kh2*ARk(ikx,iky,izh0)
               nAIk(ikx,iky,izh0) =  - kh2*AIk(ikx,iky,izh0)

            enddo
         enddo
      enddo

      !FFT nA to real-space to compute the conversion terms                                                                                                            
      call fft_c2r(nARk,nARr,n3h0)
      call fft_c2r(nAIk,nAIr,n3h0)


       !-----------------!                                                                                                                                            
       !--- Advection ---!                                                                                                                                          
       !-----------------!                                                                                                                                               

       !Start with the advective conversion term                                                                                                                          
       tRk = nBRk
       tIk = nBIk

      !FFT advection to real-space to compute the conversion term                                                                                                         
      call fft_c2r(nBRk,nBRr,n3h0)
      call fft_c2r(nBIk,nBIr,n3h0)

      !Compute the local integral for advection conversion                                                                                                                    
      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   ca_p = ca_p - 0.5*(nARr(ix,iy,izh0)*nBRr(ix,iy,izh0) + nAIr(ix,iy,izh0)*nBIr(ix,iy,izh0))

                end if
             end do
          end do
       end do

       !Recover k-space advective terms                                                                                                                                     
       nBRk = tRk
       nBIk = tIk


       !------------------!                                                                                                                                                  
       !--- Refraction ---!                                                                                                                                                  
       !------------------!                                                                                                                                                  

       !Next: refractive terms                                                                                                                                                        
       tRk = rBRk
       tIk = rBIk

       !FFT refraction to real-space to compute the conversion term                                                                                                             
       call fft_c2r(rBRk,rBRr,n3h0)
       call fft_c2r(rBIk,rBIr,n3h0)

       !Compute the local integral for refraction conversion                                                                                                                          
       do izh0=1,n3h0
          do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   cr_p = cr_p + 0.25*(nARr(ix,iy,izh0)*rBIr(ix,iy,izh0) - nAIr(ix,iy,izh0)*rBRr(ix,iy,izh0) )    !Because I define rB = zeta LA, and refraction = (i/2) zeta LA      

                end if
             end do
          end do
       end do

       !Recover k-space advective terms                                                                                                                                               
       rBRk = tRk
       rBIk = tIk


       !-------------------!                                                                                                                                            
       !--- Dissipation ---!                                                                                                                
       !-------------------!                                                                                                                                                  

       !Compute dissipation - nuhX*nabla^{2*ilapX} LA and store in the temporary array                                                                                             
       do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky

               tRk(ikx,iky,izh0) =  - 0.5*( nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w)) )*BRk(ikx,iky,izh0)
               tIk(ikx,iky,izh0) =  - 0.5*( nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w)) )*BIk(ikx,iky,izh0)

            enddo
         enddo
      enddo

      !FFT dissipation to real-space to compute the conversion term                                                                                                           
      call fft_c2r(tRk,tRr,n3h0)
      call fft_c2r(tIk,tIr,n3h0)

      !Compute the local integral for advection conversion                                                                                                                 
                                                                                                                                                                      
      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   cd_p = cd_p + nARr(ix,iy,izh0)*tRr(ix,iy,izh0) + nAIr(ix,iy,izh0)*tIr(ix,iy,izh0)

                end if
             end do
          end do
       end do


       !-------------------!                                                                                                                                             
       !--- Write files ---!                                                                                                                                                    
       !-------------------!  


       !Sum over all processors                                                                                                                                              
       call mpi_reduce(ca_p,ca_tot, 1,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(cr_p,cr_tot, 1,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierror)
       call mpi_reduce(cd_p,cd_tot, 1,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierror)

       !Normalize                                                                                                                                                                     
       ca_tot = ca_tot*Uw_scale*Uw_scale*U_scale/(n1*n2*n3*Bu*L_scale)
       cr_tot = cr_tot*Uw_scale*Uw_scale*U_scale/(n1*n2*n3*Bu*L_scale)
       cd_tot = cd_tot*Uw_scale*Uw_scale*U_scale/(n1*n2*n3*Bu*L_scale)

       if(mype==0) write(unit=unit_gamma ,fmt=*) time*(L_scale/U_scale)/(3600*24),ca_tot,cr_tot,cd_tot


     end subroutine gamma_conversion





 SUBROUTINE enstrophy(zxk,zyk,zzk) 

   !This subroutine computes enstrophy cheaply without any form of interpolation.                                                                                    
   !It uses the mid-point (or rectangle rule) for the staggered fields (omega_3)                                                                                        
   !And the trapezoidal rule for omega_1,2. There is therefore no interpolation needed.                                                                                          
   !Ens = 0.5 * int( omega^2 ) dV

   !DIMENSIONAL FORM
   !We also use the occasion to compute "small-scale" Ro and Fr numbers, defined as
   !Ro_micro = sqrt{<omega_z^2>}/f                                                                                                                                                                                         
   !Fr_micro = sqrt{<omega_H^2>}/N 
   !so that
   !Ens = 0.5*[ (f*Ro)^2 + (N*Fr)^2 ]

   !NONDIMENSIONAL FORM:
   !We also use the occasion to compute "small-scale" Ro and Fr numbers, defined as
   !Ro_micro = Ro_macro*sqrt{<omega_z^2>}                                                                                                                                                                                         
   !Fr_micro = Fr_macro*sqrt{<omega_H^2>} 
   !so that
   !Ens = 0.5*(U/H)^2 *[<omega_H^2> + Ar2*<omega_z^2> ] 
             !\___/!
   !      we omit this prop constant for now !


   double complex, dimension(iktx,ikty,n3h1) :: zxk,zyk,zzk

   real, dimension(n3h0) :: ws,wu   !Staggered (s) and unstaggered (u) omega^2
   real :: Ro_micro,Fr_micro,Ens,ws_p,wu_p,ws_tot,wu_tot

   !Let's first sum over k_h to get a function of z alone.                                                                                                                  
   ! 1/(2pi^)^2 int(int( w^2(x,y,z)/2 dx)dy = 1/2 sum over kx ky |wh|^2(z)                                                                                                 

   ws=0.
   wu=0.
   ws_p=0.
   wu_p=0.

    do iz=1,n3h0
       izh1=iz+1

       do iky=1,ikty
          ky = kya(iky)
          do ikx=1,iktx
             kx = kxa(ikx)
             
             if(L(ikx,iky)==1) then

                wu(iz) = wu(iz) + 2.*real(   zxk(ikx,iky,izh1)*CONJG( zxk(ikx,iky,izh1) )   +    zyk(ikx,iky,izh1)*CONJG( zyk(ikx,iky,izh1) )       )
                ws(iz) = ws(iz) + 2.*real(   zzk(ikx,iky,izh1)*CONJG( zzk(ikx,iky,izh1) )   )

!     wu(iz) = wu(iz) + real( ( i*ky*wk(ikx,iky,izh2) - (vk(ikx,iky,izh2+1)-vk(ikx,iky,izh2))/dz  )*CONJG( i*ky*wk(ikx,iky,izh2) - (vk(ikx,iky,izh2+1)-vk(ikx,iky,izh2))/dz  )  )
!     wu(iz) = wu(iz) + real( ( (uk(ikx,iky,izh2+1)-uk(ikx,iky,izh2))/dz - i*kx*wk(ikx,iky,izh2)  )*CONJG( (uk(ikx,iky,izh2+1)-uk(ikx,iky,izh2))/dz - i*kx*wk(ikx,iky,izh2)  )  )
!     ws(iz) = ws(iz) + real( ( i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)  )*CONJG( i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)  )  )

             end if

          enddo
       enddo

       !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.  

       wu(iz) = wu(iz) - real(   zxk(1,1,izh1)*CONJG( zxk(1,1,izh1) )   +    zyk(1,1,izh1)*CONJG( zyk(1,1,izh1) )       )
!       wu(iz) = wu(iz) - 0.5*real( ( - (vk(1,1,izh2)-vk(1,1,izh2-1))/dz  )*CONJG(- (vk(1,1,izh2)-vk(1,1,izh2-1))/dz  )  )
!       wu(iz) = wu(iz) - 0.5*real( (   (uk(1,1,izh2)-uk(1,1,izh2-1))/dz  )*CONJG(  (uk(1,1,izh2)-uk(1,1,izh2-1))/dz  )  )

    end do

    !Explicit treatment of boundary conditions
    if(mype==(npe-1))  wu(n3h0) = (0.D0,0.D0)

    !Plot as a function of z if desired!
    !**********************************!

    call plot_ensz(ws,wu)

    !1/L int(e(z) dz) becomes just a dummy sum over all grid points without interpolation needed...                                                                        
    !First some locally to each processor                                                                                                                                 

    do izh0=1,n3h0

       ws_p=ws_p + ws(izh0)
       wu_p=wu_p + wu(izh0)

    end do

    ws_p=ws_p/N3 
    wu_p=wu_p/N3

    !Sum results from each processor                                                                                                                                                                                      

    call mpi_reduce(ws_p,ws_tot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(wu_p,wu_tot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

    !This is now just <omega_X^2>. Now get the desired dimensional values.
                                                                                                          
    Ro_micro = Ro*sqrt(ws_tot)
    Fr_micro = Fr*sqrt(wu_tot)
    Ens      = 0.5*(ws_tot + Ar2*wu_tot)  !There is a (U/H)^2 constant missing for right scaling.

    if(mype==0) write(unit=unit_ens,fmt=*) time,Ens,Ro_micro,Fr_micro

  END SUBROUTINE enstrophy







!************************************************************************!
!!!! Subroutines to plot z-dependent horizontally averaged quantities !!!!
!************************************************************************!







 SUBROUTINE plot_ez(ks,ku,ps)    !Now plots z, HOR KIN, WRMS and BRMS 

   real, dimension(n3h0) :: ks,ku,ps         !Staggered (s) and unstaggered (u) energy         
   real, dimension(n3h0) :: ks_r,ku_r,ps_r   !Copies of eu and es to keep es and eu valid for energy subroutine.

   real,dimension(n3)   :: kuz,ksz,psz             !energy(z) for entire domain
   
   integer :: processor,izp,nrec
   
   ksz=0.
   kuz=0.
   psz=0.

    !Send err_p from other processors to mype = 0                                                                                                                    
    if(mype>0) call mpi_send(ks,n3h0,MPI_REAL,0,tag_kzs,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(ku,n3h0,MPI_REAL,0,tag_kzu,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(ps,n3h0,MPI_REAL,0,tag_pzs,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part onto err                  
       do iz=1,n3h0
             ksz(iz) = ks(iz)
             kuz(iz) = ku(iz)
             psz(iz) = ps(iz)
       end do


       !Receive other parts from other processors    
       do nrec=1,npe-1

          call mpi_recv(ks_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_kzs,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                       
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             ksz(izp) = ksz(izp) + ks_r(iz)

          end do

          !Kinetic unstag energy part
          call mpi_recv(ku_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_kzu,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             kuz(izp) = kuz(izp) + ku_r(iz)

          end do

          
          !Potential unstag energy part
          call mpi_recv(ps_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_pzs,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             psz(izp) = psz(izp) + ps_r(iz)

          end do

       end do

       !PLOT!
     do iz=1,n3
        write(unit_ez,fmt=*) iz,ksz(iz),kuz(iz),psz(iz)      !iz,URMS, WRMS, BRMS
     enddo
     write(unit_ez,*) '           '
     call flush(unit_ez)

   end if


 END SUBROUTINE plot_ez





 SUBROUTINE plot_rotz(ks,pu) 

   real, dimension(n3h0) :: ks,pu         !Staggered (s) and unstaggered (u) energy         
   real, dimension(n3h0) :: ks_r,pu_r   !Copies of eu and es to keep es and eu valid for energy subroutine.

   real,dimension(n3)   :: kz,pz             !energy(z) for entire domain
   
   integer :: processor,izp,nrec
   
   kz=0.
   pz=0.

    !Send err_p from other processors to mype = 0                                                                                                                    
    if(mype>0) call mpi_send(ks,n3h0,MPI_REAL,0,tag_rzs,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(pu,n3h0,MPI_REAL,0,tag_rzu,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part onto err                  
       kz(1) = ks(1) 
       pz(1) =       + 0.5*pu(1)
       do iz=2,n3h0
             kz(iz) = ks(iz) 
             pz(iz) =        + 0.5*(pu(iz) + pu(iz-1))
       end do
       pz(n3h0+1) = 0.5*pu(n3h0)   

       !Receive other parts from other processors    
       do nrec=1,npe-1

          call mpi_recv(ks_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_rzs,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                       
          processor=status(MPI_SOURCE)

          !Add the awaiting half level
          do iz=1,n3h0
             izp=processor*n3h0+iz

             kz(izp) = kz(izp) + ks_r(iz)

          end do

          !Potential unstag energy part
          call mpi_recv(pu_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_rzu,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          !Add the awaiting half level  !except for last mype
          if(processor<(npe-1))   pz((processor+1)*n3h0+1)=pz((processor+1)*n3h0+1) + 0.5*pu_r(n3h0)
      
          !iz=1 case
          izp=processor*n3h0+1
          pz(izp) = pz(izp) + 0.5* pu_r(1) 

          do iz=2,n3h0
             izp=processor*n3h0+iz

             pz(izp) = pz(izp) + 0.5*( pu_r(iz) + pu_r(iz-1) )

          end do

       end do

       !PLOT!

     do iz=1,n3
        write(unit_rotz,fmt=*) zas(iz),rho_st(iz)*kz(iz),rho_st(iz)*pz(iz)!,rho_st(iz)*(kz(iz)+pz(iz))
     enddo
     write(unit_rotz,*) '           '
     call flush(unit_rotz)

   end if


 END SUBROUTINE plot_rotz








 SUBROUTINE plot_ensz(ks,ku)

   real, dimension(n3h0) :: ks,ku            !Staggered (s) and unstaggered (u) parts of enstrophy                                                                                   
   real, dimension(n3h0) :: ks_r,ku_r        !Copies of wu and ws to keep ws and wu valid for enstrophy subroutine.                                                                  

   real,dimension(n3)   :: kz           !enstrophy(z) for entire domain                                                                                                              

   integer :: processor,izp,nrec

   kz=0.

    !Send err_p from other processors to mype = 0                                                                                                                                   \
                                                                                                                                                                                     
    if(mype>0) call mpi_send(ks,n3h0,MPI_REAL,0,tag_ezs,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(ku,n3h0,MPI_REAL,0,tag_ezu,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part onto err                                                                                                                                                  \
                                                                                                                                                                                     
       kz(1) = ks(1) + 0.5*ku(1)
       do iz=2,n3h0
             kz(iz) = ks(iz) + 0.5*(ku(iz) + ku(iz-1))
       end do
       kz(n3h0+1) = 0.5*ku(n3h0)


       !Receive other parts from other processors                                                                                                                                    
       do nrec=1,npe-1

          call mpi_recv(ks_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_ezs,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                              
          processor=status(MPI_SOURCE)

          !Add the awaiting half level                                                                                                                                               
          do iz=1,n3h0
             izp=processor*n3h0+iz

             kz(izp) = kz(izp) + ks_r(iz)

          end do

          !Kinetic unstag energy part                                                                                                                                                
          call mpi_recv(ku_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_ezu,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                              
          processor=status(MPI_SOURCE)

          !Add the awaiting half level  !except for last mype                                                                                                                        
          if(processor<(npe-1))   kz((processor+1)*n3h0+1)=kz((processor+1)*n3h0+1) + 0.5*ku_r(n3h0)


          !iz=1 case                                                                                                                                                                 
          izp=processor*n3h0+1
          kz(izp) = kz(izp) + 0.5* ku_r(1)

          do iz=2,n3h0
             izp=processor*n3h0+iz

             kz(izp) = kz(izp) + 0.5*( ku_r(iz) + ku_r(iz-1) )

          end do

       end do

       !PLOT!                                                                                                                                                                        
     do iz=1,n3
        write(unit_ensz,fmt=*) zas(iz),rho_st(iz)*kz(iz)
     enddo
     write(unit_ensz,*) '           '
     call flush(unit_ensz)

   end if


 END SUBROUTINE plot_ensz





  !****************!
  !!!!! SLICES !!!!!
  !****************!

  subroutine slices(ARk,AIK,ARr,AIr,BRk,BIk,BRr,BIr,CRk,CIk,CRr,CIr,dBRk,dBIk,dBRr,dBIr,id_field)

    double complex, dimension(iktx,ikty,n3h0) :: ARk,AIk
    double precision,    dimension(n1d,n2d,n3h0) :: ARr, AIr

    double complex, dimension(iktx,ikty,n3h0) :: BRk,BIk
    double precision,    dimension(n1d,n2d,n3h0) :: BRr, BIr

    double complex, dimension(iktx,ikty,n3h0) :: dBRk,dBIk
    double precision,    dimension(n1d,n2d,n3h0) :: dBRr, dBIr


    double complex, dimension(iktx,ikty,n3h0) :: CRk,CIk
    double precision,    dimension(n1d,n2d,n3h0) :: CRr, CIr

    !Temp arrays for convenience
    double complex, dimension(iktx,ikty,n3h0) :: Rmemk, Imemk
    double complex, dimension(iktx,ikty,n3h0) :: Rmemk2, Imemk2
    double precision,    dimension(n1d,n2d,n3h0) :: Rmem,Imem,Rmem2,Imem2

    double precision,    dimension(n1d,n2d,n3h0+2*hlvl(id_field)) :: field

    real, dimension(n1,n3h0) :: XZ_slice_p        !Scratch array for xz slices (divided amongst processors)                                                                                                                               
    real, dimension(n1,n3)   :: XZ_slice          !Scratch array for xz slices                                                                                                                                                            

    integer :: unit
    integer :: id_field
    character(len = 32) :: fname                !future file name                                                                                                                                                                         

    integer :: nrec
    integer :: processor

    equivalence(Rmem,Rmemk)
    equivalence(Imem,Imemk)
    equivalence(Rmem2,Rmemk2)
    equivalence(Imem2,Imemk2)

    !For z-slices, pick xy section!
    do ix=1,n1
!       yval(ix)=n2/4
       if(x_equal_minus_y_transect==1) then
          yval(ix)=n2-(ix-1)+y_trans   !transect = y - x + y_trans 
          if( yval(ix) > n2 ) then
             yval(ix) = yval(ix) - n2
          else if( yval(ix) < 1 ) then
             yval(ix) = yval(ix) + n2
          end if
       else
          yval(ix)=ix             !transect x = +y
       end if
    end do

    !LAR with nabla A removed
    if(id_field==1) then
       Rmemk = BRk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Akv/ Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BRk(ikx,iky,izh0) = BRk(ikx,iky,izh0) + 0.25*kh2*ARk(ikx,iky,izh0)/Bu
                   else
                      BRk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BRk,BRr,n3h0)
       field = Uw_scale*BRr
    !LAI with nabla A removed 
    elseif(id_field==2) then
       Imemk = BIk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Ak / Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BIk(ikx,iky,izh0) = BIk(ikx,iky,izh0) + 0.25*kh2*AIk(ikx,iky,izh0)/Bu
                   else
                      BIk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BIk,BIr,n3h0)
       field = Uw_scale*BIr

    !Re(L+A) for comparing to netCDF
    elseif(id_field==3) then
       Rmemk = BRk
       call fft_c2r(BRk,BRr,n3h0)
       field = Uw_scale*BRr       
    !Im(L+A) for comparing to netCDF
    elseif(id_field==4) then
       Imemk = BIk
       call fft_c2r(BIk,BIr,n3h0)
       field = Uw_scale*BIr

       !LAR_x
    elseif(id_field==5) then
       Rmemk = BRk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Akv/ Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BRk(ikx,iky,izh0) = i*kx*(BRk(ikx,iky,izh0) + 0.25*kh2*ARk(ikx,iky,izh0)/Bu)
                   else
                      BRk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BRk,BRr,n3h0)
       field = Uw_scale*BRr/L_scale

       !LAI_x
    elseif(id_field==6) then
       Imemk = BIk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Ak / Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BIk(ikx,iky,izh0) = i*kx*(BIk(ikx,iky,izh0) + 0.25*kh2*AIk(ikx,iky,izh0)/Bu)
                   else
                      BIk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BIk,BIr,n3h0)
       field = Uw_scale*BIr/L_scale



       !LAR_y
    elseif(id_field==7) then
       Rmemk = BRk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Akv/ Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BRk(ikx,iky,izh0) = i*ky*(BRk(ikx,iky,izh0) + 0.25*kh2*ARk(ikx,iky,izh0)/Bu)
                   else
                      BRk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BRk,BRr,n3h0)
       field = Uw_scale*BRr/L_scale

       !LAI_y
    elseif(id_field==8) then
       Imemk = BIk
       if(ybj_plus==1) then  !In YBJ+, B = L^+ A = LA + 0.25 nabla(A), or in nondimensional k-space, Bk = LAk - 0.25 k_h^2 Ak /Bu. Here we transform Bk into LAk by adding 0.25 kh2 Ak / Bu
          do izh0=1,n3h0
             do iky=1,ikty
                ky = kya(iky)
                do ikx=1,iktx
                   kx = kxa(ikx)
                   kh2=kx*kx+ky*ky
                   if(L(ikx,iky)==1) then
                      BIk(ikx,iky,izh0) = i*ky*(BIk(ikx,iky,izh0) + 0.25*kh2*AIk(ikx,iky,izh0)/Bu)
                   else
                      BIk(ikx,iky,izh0) = (0.D0,0.D0)
                   end if
                enddo
             enddo
          end do
       end if
       call fft_c2r(BIk,BIr,n3h0)
       field = Uw_scale*BIr/L_scale

    end if


    !Print bottom slice
    if( bot_height > mype*n3h0 .AND. bot_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehbot",id_field,count_slice(id_field),".dat"
       write (fname, "(A9,I1,I3,A4)") "slicehbot",id_field,count_slice(id_field),".dat"
       open (unit=unit_slices,file=fname,action="write",status="replace")
       
       iz=bot_height - mype*n3h0 + hlvl(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices,fmt=*) '           '
       enddo
333    format(1x,E12.5,1x)
       
       close (unit=unit_slices)
       
    end if
    
    
    !Print mid-height slice
    if( mid_height > mype*n3h0 .AND. mid_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehmid",id_field,count_slice(id_field),".dat"
       write (fname, "(A9,I1,I3,A4)") "slicehmid",id_field,count_slice(id_field),".dat"
       open (unit=unit_slices,file=fname,action="write",status="replace")
       
       iz=mid_height - mype*n3h0 + hlvl(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices)
       
    end if
    
    !Print top slice
    if( top_height > mype*n3h0 .AND. top_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehtop",id_field,count_slice(id_field),".dat"
       write (fname, "(A9,I1,I3,A4)") "slicehtop",id_field,count_slice(id_field),".dat"
       open (unit=unit_slices,file=fname,action="write",status="replace")
       
       iz=top_height - mype*n3h0 + hlvl(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices)
       
    end if




    !Print vertical slice

       if(mype==0) then
!          write (fname, "(A6,I1,I1,A4)") "slicev",id_field,count_slice(id_field),".dat"
          write (fname, "(A6,I1,I3,A4)") "slicev",id_field,count_slice(id_field),".dat"
          open (unit=unit_slicesv,file=fname,action="write",status="replace")

          !Copy ur slice on XY_slice (NOTICE IT'S NOT ON XY_slice_p)                                                                                                                                                                      
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0+hlvl(id_field)
                XZ_slice(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do

          !Receive from other processors                                                                                                                                                                                                  
          do nrec=1,npe-1
             call mpi_recv(XZ_slice_p,n1*n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_slice_xz(id_field),MPI_COMM_WORLD,status,ierror)
             processor=status(MPI_SOURCE)
             !Copy onto scratch array                                                                                                                                                                                                     
             do ix=1,n1
                do iz=1,n3h0
                   XZ_slice(ix,iz+n3h0*processor) = XZ_slice_p(ix,iz)
                end do
             end do
          end do

          !Now print the complete slice onto file                                                                                                                                                                                         
!          do ix=1,n1
!            do iz=1,n3
!               write(unit=count_slice(id_field),fmt=*) real(xa(ix)),real(za(iz)),XZ_slice(ix,iz)
!             end do
!          end do
          do iz=1,n3
             write(unit=unit_slicesv,fmt=333) (XZ_slice(ix,iz),ix=1,n1)
             write(unit=unit_slicesv,fmt=*) '           '
          enddo
          close (unit=unit_slicesv)

          
       end if

       !All other processors (mype>0) send their XZ_field_p to mype 0                                                                                                                                                                     
       if(mype/=0) then
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0 + hlvl(id_field)
                XZ_slice_p(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do
          
          !Now send these chunks to mype 0                                                                                                                                                                                                
          call mpi_send(XZ_slice_p,n1*n3h0,MPI_REAL,0,tag_slice_xz(id_field),MPI_COMM_WORLD,ierror)
          
       end if





          if(id_field==1)    then 
             BRk=Rmemk
          elseif(id_field==2)    then 
             BIk=Imemk
          elseif(id_field==3)    then
             BRk=Rmemk
          elseif(id_field==4)    then
             BIk=Imemk
          elseif(id_field==5)    then
             BRk=Rmemk
          elseif(id_field==6)    then
             BIk=Imemk
          elseif(id_field==7)    then
             BRk=Rmemk
          elseif(id_field==8)    then
             BIk=Imemk
          end if

          count_slice(id_field)=count_slice(id_field)+1


     end subroutine slices





  subroutine slices2(uk,vk,wak,bk,psik,ur,vr,war,br,psir,id_field)

    double complex, dimension(iktx,ikty,n3h2) :: uk,vk,bk
    double complex, dimension(iktx,ikty,n3h1) :: zzk,psik,wak
    
    double precision,    dimension(n1d,n2d,n3h2) :: ur,vr,br
    double precision,    dimension(n1d,n2d,n3h1) :: zzr,psir,war

    double complex, dimension(iktx,ikty,n3h2) :: bmem
    double complex, dimension(iktx,ikty,n3h1) :: qmem
    
    !Suboptimal storage use: only to avoid fft back in field=9: the Jacobian of psi with zeta
    double complex, dimension(iktx,ikty,n3h2) :: umem
    double complex, dimension(iktx,ikty,n3h2) :: vmem

    double complex, dimension(iktx,ikty,n3h1)    :: zxk,zyk
    double precision,    dimension(n1d,n2d,n3h1) :: zxr,zyr



    double precision,    dimension(n1d,n2d,n3h0+2*hlvl2(id_field)) :: field

    real, dimension(n1,n3h0) :: XZ_slice_p        !Scratch array for xz slices (divided amongst processors)                                                                                                                               
    real, dimension(n1,n3)   :: XZ_slice          !Scratch array for xz slices                                                                                                                                                            

    integer :: unit
    integer :: id_field
    character(len = 32) :: fname                !future file name                                                                                                                                                                         

    integer :: nrec
    integer :: processor

    equivalence(zxk,zxr)
    equivalence(zyk,zyr)
    equivalence(zzk,zzr)

    !For z-slices, pick xy section!                                                                                                                                                
    do ix=1,n1
!       yval(ix)=n2/4                                                                                                                                                              
       if(x_equal_minus_y_transect==1) then
          yval(ix)=n2-(ix-1)+y_trans   !transect = y - x + y_trans                                                                                                                 
          if( yval(ix) > n2 ) then
             yval(ix) = yval(ix) - n2
          else if( yval(ix) < 1 ) then
             yval(ix) = yval(ix) + n2
          end if
       else
          yval(ix)=ix             !transect x = +y                                                                                                                                 
       end if
    end do



    if(id_field==1)   then
       bmem=uk
       call fft_c2r(uk,ur,n3h2)
       field = U_scale*ur  
    else if(id_field==2) then
       bmem=vk
       call fft_c2r(vk,vr,n3h2)
       field = U_scale*vr  
    else if(id_field==3) then
       qmem=psik
       call fft_c2r(psik,psir,n3h1)
       field = U_scale*L_scale*psir    !In fact, w_DIM = UH/L w_NDIM (in the code) = UH/L Ro w_NDIM_1 (since w_NDIM_0 =0 in QG) => w_real_life = U Ar Ro w_1_computed
    else if(id_field==4) then             !zeta_x for comparison with J(psi,zeta)
       !Calculate the vorticity gradient                                                                                                                                                
       do izh1=1,n3h1
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2=kx*kx+ky*ky

                if(L(ikx,iky)==1) then
                   zxk(ikx,iky,izh1) = -i*kx*kh2*psik(ikx,iky,izh1)
                else
                   zxk(ikx,iky,izh1) = (0.D0,0.D0)
                end if

             enddo
          enddo
       end do
       call fft_c2r(zxk,zxr,n3h1)
       field = zxr*U_scale/(L_scale*L_scale)
    else if(id_field==5) then             !zeta_y for comparison with J(psi,zeta)
       !Calculate the vorticity gradient                                                                                                                                                
       do izh1=1,n3h1
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2=kx*kx+ky*ky

                if(L(ikx,iky)==1) then
                   zyk(ikx,iky,izh1) = -i*ky*kh2*psik(ikx,iky,izh1)
                else
                   zyk(ikx,iky,izh1) = (0.D0,0.D0)
                end if

             enddo
          enddo
       end do
       call fft_c2r(zyk,zyr,n3h1)
       field = zyr*U_scale/(L_scale*L_scale)
    elseif(id_field==6) then   !Fr^2/Ro * b_z/r_2 << 1 ? (Limiting QG assumption)                                                                                                                                                           
       bmem=bk
       call fft_c2r(bk,br,n3h2)

       if(where_bz==stag) then
          do izh1=1,n3h1
             izh2=izh1+1
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      field(ix,iy,izh1) = (br(ix,iy,izh2+1)-br(ix,iy,izh2-1))/(2.*dz)  
                   else
                      field(ix,iy,izh1) = 0.
                   end if
                end do
             end do
          end do
          
          if(mype==0) then
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      field(ix,iy,izbot1) =  (br(ix,iy,izbot2+1) - br(ix,iy,izbot2)   )/(2.D0*dz)
                   else
                      field(ix,iy,izbot1) = 0.
                   end if
                end do
             end do
          elseif(mype==(npe-1)) then
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      field(ix,iy,iztop1) =  (br(ix,iy,iztop2) - br(ix,iy,iztop2-1)   )/(2.D0*dz)
                   else
                      field(ix,iy,iztop1) = 0.
                   end if
                end do
             end do
          end if
       elseif(where_bz==unstag) then     !We then want bz/N^2 on UNSTAG points.     
          do izh1=1,n3h1
             izh2=izh1+1
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      field(ix,iy,izh1) =  (br(ix,iy,izh2+1)-br(ix,iy,izh2))/dz 
                   else
                      field(ix,iy,izh1) = 0.
                   end if
                end do
             end do
          end do
          
          if(mype==(npe-1)) then
             do ix=1,n1d
                do iy=1,n2d
                      field(ix,iy,iztop1) = 0.
                end do
             end do
          end if


       end if

       field = (U_scale*U_scale/(H_scale*H_scale*Ro))*field     !In fact, b_real_life = fUL/H b_computed = U^2/(Ro*H) b_computed                                                                                                            
    else if(id_field==7) then             !micro Rossby number (or in QG, relative vorticity)                                                                                                                                 
       !Compute the z-component of vorticity

       zzk = (0.D0,0.D0)
       do izh1=1,n3h1
          izh2=izh1+1
          
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                
                if(L(ikx,iky)==1) then
                   
                   zzk(ikx,iky,izh1) = i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)
                   
                else
                   
                   zzk(ikx,iky,izh1) = (0.D0,0.D0)
                   
                end if
                
             enddo
          enddo
       end do
       
       call fft_c2r(zzk,zzr,n3h1)
       field = Ro*zzr
    else if(id_field==8) then             !Calculate the normalized (shear componet of) strain: 1/2 (psi_xx - psi_yy) / f
       qmem=psik
       do izh1=1,n3h1
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)

                if(L(ikx,iky)==1) then

                   psik(ikx,iky,izh1) = -(kx*kx-ky*ky)*psik(ikx,iky,izh1)

                else

                   psik(ikx,iky,izh1) = (0.D0,0.D0)

                end if

             enddo
          enddo
       end do
       call fft_c2r(psik,psir,n3h1)
       field = Ro*psir

    else if(id_field==9) then    !Calculate J(psi,zeta) = U zeta_x + V zeta_y
       umem=uk
       vmem=vk

       call fft_c2r(uk,ur,n3h2)      
       call fft_c2r(vk,vr,n3h2)      


       !Calculate the vorticity gradient
       do izh1=1,n3h1
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2=kx*kx+ky*ky

                if(L(ikx,iky)==1) then

                   zxk(ikx,iky,izh1) = -i*kx*kh2*psik(ikx,iky,izh1)
                   zyk(ikx,iky,izh1) = -i*ky*kh2*psik(ikx,iky,izh1)

                else

                   zxk(ikx,iky,izh1) = (0.D0,0.D0)
                   zyk(ikx,iky,izh1) = (0.D0,0.D0)

                end if

             enddo
          enddo
       end do

       call fft_c2r(zxk,zxr,n3h1)
       call fft_c2r(zyk,zyr,n3h1)

       
       do izh1=1,n3h1
          izh2=izh1+1
          do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then
                   
                   field(ix,iy,izh1) = ur(ix,iy,izh2)*zxr(ix,iy,izh1) + vr(ix,iy,izh2)*zyr(ix,iy,izh1)  

                else

                   field(ix,iy,izh1) = 0.

                end if
             end do
          end do
       end do
       field = field*(U_scale*U_scale)/(L_scale*L_scale)
    end if


    !Print bottom slice
    if( bot_height > mype*n3h0 .AND. bot_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehbot",id_field,count_slice(id_field),".dat"
       write (fname, "(A10,I1,I3,A4)") "slice2hbot",id_field,count_slice2(id_field),".dat"
       open (unit=unit_slices2,file=fname,action="write",status="replace")
       
       iz=bot_height - mype*n3h0 + hlvl2(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices2,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices2,fmt=*) '           '
       enddo
333    format(1x,E12.5,1x)
       
       close (unit=unit_slices2)
       
    end if
    
    
    !Print mid-height slice
    if( mid_height > mype*n3h0 .AND. mid_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehmid",id_field,count_slice(id_field),".dat"
       write (fname, "(A10,I1,I3,A4)") "slice2hmid",id_field,count_slice2(id_field),".dat"
       open (unit=unit_slices2,file=fname,action="write",status="replace")
       
       iz=mid_height - mype*n3h0 + hlvl2(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices2,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices2,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices2)
       
    end if
    
    !Print top slice
    if( top_height > mype*n3h0 .AND. top_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehtop",id_field,count_slice(id_field),".dat"
       write (fname, "(A10,I1,I3,A4)") "slice2htop",id_field,count_slice2(id_field),".dat"
       open (unit=unit_slices2,file=fname,action="write",status="replace")
       
       iz=top_height - mype*n3h0 + hlvl2(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices2,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices2,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices2)
       
    end if




    !Print vertical slice

       if(mype==0) then
!          write (fname, "(A6,I1,I1,A4)") "slicev",id_field,count_slice(id_field),".dat"
          write (fname, "(A7,I1,I3,A4)") "slice2v",id_field,count_slice2(id_field),".dat"
          open (unit=unit_slices2v,file=fname,action="write",status="replace")

          !Copy ur slice on XY_slice (NOTICE IT'S NOT ON XY_slice_p)                                                                                                                                                                      
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0+hlvl2(id_field)
                XZ_slice(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do

          !Receive from other processors                                                                                                                                                                                                  
          do nrec=1,npe-1
             call mpi_recv(XZ_slice_p,n1*n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_slice_xz2(id_field),MPI_COMM_WORLD,status,ierror)
             processor=status(MPI_SOURCE)
             !Copy onto scratch array                                                                                                                                                                                                     
             do ix=1,n1
                do iz=1,n3h0
                   XZ_slice(ix,iz+n3h0*processor) = XZ_slice_p(ix,iz)
                end do
             end do
          end do

          !Now print the complete slice onto file                                                                                                                                                                                         
!          do ix=1,n1
!            do iz=1,n3
!               write(unit=count_slice(id_field),fmt=*) real(xa(ix)),real(za(iz)),XZ_slice(ix,iz)
!             end do
!          end do
          do iz=1,n3
             write(unit=unit_slices2v,fmt=333) (XZ_slice(ix,iz),ix=1,n1)
             write(unit=unit_slices2v,fmt=*) '           '
          enddo
          close (unit=unit_slices2v)

          
       end if

       !All other processors (mype>0) send their XZ_field_p to mype 0                                                                                                                                                                     
       if(mype/=0) then
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0 + hlvl2(id_field)
                XZ_slice_p(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do
          
          !Now send these chunks to mype 0                                                                                                                                                                                                
          call mpi_send(XZ_slice_p,n1*n3h0,MPI_REAL,0,tag_slice_xz2(id_field),MPI_COMM_WORLD,ierror)
          
       end if





          if(id_field==1)    uk=bmem
          if(id_field==2)    vk=bmem
          if(id_field==3)    psik=qmem
          if(id_field==6)    bk=bmem
          if(id_field==8)    psik=qmem
          if(id_field==9) then
             uk=umem
             vk=vmem
          end if

          count_slice2(id_field)=count_slice2(id_field)+1


        end subroutine slices2










  subroutine slices3(ARk,AIK,ARr,AIr,dBRk,dBIk,dBRr,dBIr,nBRk,nBIk,nBRr,nBIr,rBRk,rBIk,rBRr,rBIr,id_field)

    double complex, dimension(iktx,ikty,n3h0) :: ARk,AIk
    double precision,    dimension(n1d,n2d,n3h0) :: ARr, AIr

    double complex, dimension(iktx,ikty,n3h0) :: dBRk,dBIk
    double precision,    dimension(n1d,n2d,n3h0) :: dBRr, dBIr

    double complex, dimension(iktx,ikty,n3h0) :: nBRk,nBIk
    double precision,    dimension(n1d,n2d,n3h0) :: nBRr, nBIr

    double complex, dimension(iktx,ikty,n3h0) :: rBRk,rBIk
    double precision,    dimension(n1d,n2d,n3h0) :: rBRr, rBIr


    !Temp arrays for convenience
    double complex, dimension(iktx,ikty,n3h0) :: Rmemk, Imemk
    double complex, dimension(iktx,ikty,n3h0) :: Rmemk2, Imemk2
    double precision,    dimension(n1d,n2d,n3h0) :: Rmem,Imem,Rmem2,Imem2

    double precision,    dimension(n1d,n2d,n3h0+2*hlvl3(id_field)) :: field

    real, dimension(n1,n3h0) :: XZ_slice_p        !Scratch array for xz slices (divided amongst processors)                                                             
    real, dimension(n1,n3)   :: XZ_slice          !Scratch array for xz slices                                                                                             

    integer :: unit
    integer :: id_field
    character(len = 32) :: fname                !future file name                                                                         

    integer :: nrec
    integer :: processor

    equivalence(Rmem,Rmemk)
    equivalence(Imem,Imemk)
    equivalence(Rmem2,Rmemk2)
    equivalence(Imem2,Imemk2)

    !For z-slices, pick xy section!
    do ix=1,n1
!       yval(ix)=n2/4
       if(x_equal_minus_y_transect==1) then
          yval(ix)=n2-(ix-1)+y_trans   !transect = y - x + y_trans 
          if( yval(ix) > n2 ) then
             yval(ix) = yval(ix) - n2
          else if( yval(ix) < 1 ) then
             yval(ix) = yval(ix) + n2
          end if
       else
          yval(ix)=ix             !transect x = +y
       end if
    end do

    !d/dt LAR
    if(id_field==1) then
       Rmemk = dBRk
       call fft_c2r(dBRk,dBRr,n3h0)
       field = Uw_scale*dBRr*U_scale/L_scale
    !d/dt LAI
    elseif(id_field==2) then
       Rmemk = dBIk
       call fft_c2r(dBIk,dBIr,n3h0)
       field = Uw_scale*dBIr*U_scale/L_scale
    !J(psi,LAR)
    elseif(id_field==3) then
       Rmemk = nBRk
       call fft_c2r(nBRk,nBRr,n3h0)
       field = Uw_scale*nBRr*U_scale/L_scale
    !J(psi,LAI)
    elseif(id_field==4) then
       Rmemk = nBIk
       call fft_c2r(nBIk,nBIr,n3h0)
       field = Uw_scale*nBIr*U_scale/L_scale
    !Re(i/2 zeta*LA) = - 0.5 rBIk
    elseif(id_field==5) then
       Rmemk = rBIk
       call fft_c2r(rBIk,rBIr,n3h0)
       field = -0.5*Uw_scale*rBIr*U_scale/L_scale
    !Im(i/2 zeta*LA) = + 0.5 rBRk
    elseif(id_field==6) then
       Rmemk = rBRk
       call fft_c2r(rBRk,rBRr,n3h0)
       field =  0.5*Uw_scale*rBRr*U_scale/L_scale
    !Re(i/2 cor*nabla(A)) = -0.5 cor nAI
    elseif(id_field==7) then
       Rmemk = AIk
       do izh0=1,n3h0
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2=kx*kx+ky*ky
                
                AIk(ikx,iky,izh0) = -kh2*AIk(ikx,iky,izh0)

             end do
          end do
       end do

       call fft_c2r(AIk,AIr,n3h0)
       field = -0.5*cor*Uw_scale*AIr/Bu
    !Im(i/2 cor*nabla(A)) = +0.5 cor nAR
    elseif(id_field==8) then
       Rmemk = ARk
       do izh0=1,n3h0
          do iky=1,ikty
             ky = kya(iky)
             do ikx=1,iktx
                kx = kxa(ikx)
                kh2=kx*kx+ky*ky
                
                ARk(ikx,iky,izh0) = -kh2*ARk(ikx,iky,izh0)

             end do
          end do
       end do
       call fft_c2r(ARk,ARr,n3h0)
       field =  0.5*cor*Uw_scale*ARr/Bu
    end if


    !Print bottom slice
    if( bot_height > mype*n3h0 .AND. bot_height <= (mype+1)*n3h0 ) then
       
       write (fname, "(A10,I1,I3,A4)") "slicehbot3",id_field,count_slice3(id_field),".dat"
       open (unit=unit_slices3,file=fname,action="write",status="replace")
       
       iz=bot_height - mype*n3h0 + hlvl3(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices3,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices3,fmt=*) '           '
       enddo
333    format(1x,E12.5,1x)
       
       close (unit=unit_slices3)
       
    end if
    
    
    !Print mid-height slice
    if( mid_height > mype*n3h0 .AND. mid_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehmid",id_field,count_slice(id_field),".dat"
       write (fname, "(A10,I1,I3,A4)") "slicehmid3",id_field,count_slice3(id_field),".dat"
       open (unit=unit_slices3,file=fname,action="write",status="replace")
       
       iz=mid_height - mype*n3h0 + hlvl3(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices3,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices3,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices3)
       
    end if
    
    !Print top slice
    if( top_height > mype*n3h0 .AND. top_height <= (mype+1)*n3h0 ) then
       
       !          write (fname, "(A9,I1,I1,A4)") "slicehtop",id_field,count_slice(id_field),".dat"
       write (fname, "(A10,I1,I3,A4)") "slicehtop3",id_field,count_slice3(id_field),".dat"
       open (unit=unit_slices3,file=fname,action="write",status="replace")
       
       iz=top_height - mype*n3h0 + hlvl3(id_field)
       
       
       do iy=1,n2
          write(unit=unit_slices3,fmt=333) (real(field(ix,iy,iz)),ix=1,n1)
          write(unit=unit_slices3,fmt=*) '           '
       enddo
       
       
       close (unit=unit_slices3)
       
    end if




    !Print vertical slice

       if(mype==0) then
!          write (fname, "(A6,I1,I1,A4)") "slicev",id_field,count_slice(id_field),".dat"
          write (fname, "(A7,I1,I3,A4)") "slicev3",id_field,count_slice3(id_field),".dat"
          open (unit=unit_slices3v,file=fname,action="write",status="replace")

          !Copy ur slice on XY_slice (NOTICE IT'S NOT ON XY_slice_p)                                                                                                                                                                      
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0+hlvl3(id_field)
                XZ_slice(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do

          !Receive from other processors                                                                                                                                                                                                  
          do nrec=1,npe-1
             call mpi_recv(XZ_slice_p,n1*n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_slice_xz3(id_field),MPI_COMM_WORLD,status,ierror)
             processor=status(MPI_SOURCE)
             !Copy onto scratch array                                                                                                                                                                                                     
             do ix=1,n1
                do iz=1,n3h0
                   XZ_slice(ix,iz+n3h0*processor) = XZ_slice_p(ix,iz)
                end do
             end do
          end do

          !Now print the complete slice onto file                                                                                                                                                                                         
!          do ix=1,n1
!            do iz=1,n3
!               write(unit=count_slice(id_field),fmt=*) real(xa(ix)),real(za(iz)),XZ_slice(ix,iz)
!             end do
!          end do
          do iz=1,n3
             write(unit=unit_slices3v,fmt=333) (XZ_slice(ix,iz),ix=1,n1)
             write(unit=unit_slices3v,fmt=*) '           '
          enddo
          close (unit=unit_slices3v)

          
       end if

       !All other processors (mype>0) send their XZ_field_p to mype 0                                                                                                                                                                     
       if(mype/=0) then
          do ix=1,n1
             do izh0=1,n3h0
                iz  =izh0 + hlvl3(id_field)
                XZ_slice_p(ix,izh0) = field(ix,yval(ix),iz)
             end do
          end do
          
          !Now send these chunks to mype 0                                                                                                                                                                                                
          call mpi_send(XZ_slice_p,n1*n3h0,MPI_REAL,0,tag_slice_xz3(id_field),MPI_COMM_WORLD,ierror)
          
       end if





          if(id_field==1)    then 
             dBRk=Rmemk
          elseif(id_field==2)    then 
             dBIk=Rmemk
          elseif(id_field==3)    then 
             nBRk=Rmemk
          elseif(id_field==4)    then 
             nBIk=Rmemk
          elseif(id_field==5)    then 
             rBIk=Rmemk
          elseif(id_field==6)    then 
             rBRk=Rmemk
          elseif(id_field==7)    then 
             AIk=Rmemk
          elseif(id_field==8)    then 
             ARk=Rmemk
          end if

          count_slice3(id_field)=count_slice3(id_field)+1


        end subroutine slices3







  !**************!
  !!!! OTHERS !!!!
  !**************!

  SUBROUTINE validate_run

    integer :: dummy_count

    !Makes sure that the run is sensible

    !1. Check time step
       
!    write(*,*) "dt=",delt
!    if(base_state >= linear_prof) write(*,*) "GW=",sqrt(Ar2)*Fr/( sqrt(r_1st(n3/2+1)*r_2st(n3/2+1) ) )
!    if(base_state <3) write(*,*) "GW=",sqrt(Ar2)*Fr 
    

       !In nondim, stability for inertial waves becomes 1/Ro instead of f...
    if(delt >= Ro) then
       write(*,*) "Can't resolve inertial waves, Ro=",Ro
!       stop
!    elseif( (normalize==1 .and.   delt >= dz/sqrt(k_init+p_init)) .or. (norm_trop==1) .and. delt >= dz/URMS  ) then   !Approximate CFL depending on the normalization process. Watch out if no normalization! 
!       write(*,*) "CFL fails"
!       stop
    end if
       
 
    !2. Check grid and resolution

    if(mod(n3,npe) /=0) then
       write(*,*) "n3/npe not an integer"  
       stop
    elseif(n3/npe <2)then
       write(*,*) "n3/npe must be >= 2 for differentiation"
       stop
    end if




    !4 Accurate resolving of all the horizontal wavenumbers (see A Note on the Numerical Representation ... Tulloch & Smith 2009)

 !   if( (4./3.)*n3 < sqrt(N_2_stra)*sqrt(Ar2)*(1.*n1)/cor) then  !N3 is NOT > (N0/f) Ar N1.  (I choose N0 is N_s, the max N)  !Relax this to 1.5*N3 is NOT > ...
 !      write(*,*) "Insufficient vertical resolution: sig*Ar=",sqrt(N_2_stra*Ar2)/cor
       !stop
 !   end if


    !5. Monitoring parameters...

    write(unit_run,*) "Resolution = ",n1," x ",n2," x ",n3
    write(unit_run,*) "Processors = ",npe
    write(unit_run,*)
    write(unit_run,*) "Viscosity"
    write(unit_run,*)
    write(unit_run,*) "Laplacian order =",ilap1,ilap2,ilap1w,ilap2w
    write(unit_run,*) "Hor. coeffs=",coeff1,coeff2,coeff1w,coeff2w
    write(unit_run,*) "Ver. coeff=",coeffz
   


    write(unit_run,*) "Stratification=",stratification
    write(unit_run,*)
!    write(unit_run,*) "N_troposphere  = ",real(sqrt(N_2_trop))
!    write(unit_run,*) "N_stratosphere = ",real(sqrt(N_2_stra))
    write(unit_run,*)
 !   write(unit_run,*) "h = ",real(H_scale*H_N),"meters"
 !   write(unit_run,*) "h_NONDIM = 2pi/",fraction," or ",real(1.*n3/(1.*fraction)),"dz"
 !   write(unit_run,*) "z0=",z0
 !   write(unit_run,*) "where_bz=",where_bz
 !   write(unit_run,*) "norm_trop=",norm_trop
    

 !   write(unit_run,*)  "Vertial structure = ", init_vertical_structure

    write(unit_run,*)  
    write(unit_run,*)

    write(unit_run,*) "Coriolis parameter (f)= ",cor,"s^-1   or ",cor/0.0001,"f_earth"  
    write(unit_run,*) "U = ",U_scale,"ms^-1"
    write(unit_run,*) "L = ",L_scale,"m"
    write(unit_run,*) "H = ",H_scale,"m"
   
    write(unit_run,*)  
    write(unit_run,*)



 !   write(unit_run,*) "N/f = ",real(sqrt(N_2_trop)/cor)," in the troposphere"
 !   write(unit_run,*) "N/f = ",real(sqrt(N_2_stra)/cor)," in the stratosphere"
    write(unit_run,*) "Ar = H/L = 1/",L_scale/H_scale
    write(unit_run,*) "Ro = ",Ro,"   Fr = ",Fr,"   Bu = ",Bu    !Average Fr
    
 !   write(unit_run,*) "Transition wavelengths (DIM)"
 !   write(unit_run,*)  

 !   write(unit_run,*) "lambda_D = 2pi/kD = ",real(twopi*0.5*(sqrt(N_2_trop)+sqrt(N_2_stra))*10000/cor)            ,"meters"          !kD = f/ND with D is the depth of the domain = 10000m 
 !   write(unit_run,*) "lambda_J = 2pi/kJ = ",real(twopi*0.5*(sqrt(N_2_trop)+sqrt(N_2_stra))*real(H_scale*H_N)/cor),"meters"          !kD = f/ND with D is the depth of the domain = 10000m 

    write(unit_run,*)  
    write(unit_run,*)  

 !   write(unit_run,*) "Transition wavenumbers (NONDIM) - for SQG behavior, kD << kh << kJ" 
 !   write(unit_run,*) 

 !   write(unit_run,*) "kD = ",sqrt(Bu)/(0.5*twopi)   !fL/NH *1/D' where D' is the nondim dist between tropopause and bot or top.
 !   write(unit_run,*) "kJ = ",sqrt(Bu)/H_N
 !   write(unit_run,*) "kT = N/3 = ",n1/3          !This is the max well-resolved wavenumber

 !   write(unit_run,*) "Numerical Stability"
 !   write(unit_run,*)      
    write(unit_run,*) "dt  =",delt 
 !   write(unit_run,*) "CFL =",dz/sqrt(k_init+p_init)
  !  write(unit_run,*) "GW  =",U_scale/(L_scale*sqrt(N_2_stra))
  !  write(unit_run,*) "IW  =",Ro

    write(unit_run,*)      
    write(unit_run,*)      

    write(unit_run,*) "Outputs"
    write(unit_run,*)      

    
!    write(unit_run,*) "Location of horizontal spectra"
    
!    write(unit_run,*) 

 
 !   write(unit_run,*) "Distance in h units from z0 [ Dist = (z-z0)/h ]"

 !   do dummy_count=1,num_spec
 !          write(unit_run,*) "spec,dist=",dummy_count-1, (zas(height(dummy_count))-z0)/H_N 
 !   end do

 !   write(unit_run,*) 
!    write(unit_run,*) 

!    write(unit_run,*) "Location of horizontal slices"
    
!    write(unit_run,*) 
 
!    write(unit_run,*) "Distance in h units from z0 [ Dist = (z-z0)/h ]"

!    write(unit_run,*) "bot=", (zas(bot_height)-z0)/H_N 
!    write(unit_run,*) "mid=", (zas(mid_height)-z0)/H_N 
!    write(unit_run,*) "top=", (zas(top_height)-z0)/H_N 
    

    write(unit_run,*) 
    write(unit_run,*) 

    write(unit_run,*) "Period of total energy        output:",freq_etot*delt,  freq_etot
    write(unit_run,*) "Period of total wave energy   output:",freq_we*delt,  freq_we
    write(unit_run,*) "Period of flow z-profile      output:",freq_ez*delt,  freq_ez
    write(unit_run,*) "Period of wave z-profile      output:",freq_wz*delt,  freq_wz
    write(unit_run,*) "Period of horziontal spectrum output:",freq_hspec*delt, freq_hspec
    write(unit_run,*) "Period of            slab     output:",freq_slab*delt,  freq_slab
    write(unit_run,*) "Period of            slices   output:",freq_slice*delt, freq_slice

    
    

  END SUBROUTINE VALIDATE_RUN




!**************************************!
!!!!! ARCHIVES - VESTIGES EVOLUTIF !!!!!
!**************************************!






 SUBROUTINE continuity_anelastic(uk,vk,wk)

   !Computes the extent to which the continuity equation fails (l2norm of div(u) in h-space)                                                                                                              
   double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk
   double complex, dimension(iktx,ikty,n3h0) :: div    !Div not squared                                                                                                            
   double precision, dimension(n1d,n2d,n3h0) :: divr   !Div in real-space
                                                      
   double precision :: dtot_p,dtot                 !Average value of div
   real, dimension(4) :: dmax_p,dmax,dmax_recv       !Maximum absolute value of div and its 3-coordinate location in the domain. 



   integer :: nrec

   equivalence(divr,div)

   div=(0.D0,0.D0)

   dtot_p=0.D0
   dtot=0.D0


 do izh0=1,n3h0
    izh2=izh0+2
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          if (L(ikx,iky).eq.1) then
             div(ikx,iky,izh0) = i*kx*uk(ikx,iky,izh2) + i*ky*vk(ikx,iky,izh2) + 0.5*r_3(izh2)*wk(ikx,iky,izh2) + wk(ikx,iky,izh2)/dz + 0.5*r_3(izh2)*wk(ikx,iky,izh2-1) - wk(ikx,iky,izh2-1)/dz
          endif
       enddo
    enddo
 enddo

 !Explicit Boundary /treatment (required for top only)                                                                                                                                                                             

 if(mype==0) then

    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          if (L(ikx,iky).eq.1) then
             div(ikx,iky,1) = i*kx*uk(ikx,iky,izbot2) + i*ky*vk(ikx,iky,izbot2) + 0.5*r_3(izbot2)*wk(ikx,iky,izbot2) + wk(ikx,iky,izbot2)/dz
          endif
       enddo
    enddo



 elseif(mype==(npe-1)) then

    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          if (L(ikx,iky).eq.1) then
             div(ikx,iky,n3h0) = i*kx*uk(ikx,iky,iztop2) + i*ky*vk(ikx,iky,iztop2) + 0.5*r_3(iztop2)*wk(ikx,iky,iztop2-1) - wk(ikx,iky,iztop2-1)/dz
          endif
       enddo
    enddo

 end if

!  1/V triple integral ( DIV^2 ) dV = 1/N3 [ sum_i=1,N3 sum over all k's 2*L(kx,ky) |DIVK|^2 - |DIVK(kh=0)|^2 ]

!z-integral using the mid-point rule (just a dummy sum for staggered field)                                                                                                                                                             
 do izh0=1,n3h0
    do iky=1,ikty
       do ikx=1,iktx
          if (L(ikx,iky).eq.1) then
             dtot_p=dtot_p+2.D0*DBLE(div(ikx,iky,izh0)*CONJG(div(ikx,iky,izh0)))
          endif
       enddo
    enddo
 enddo

 do izh0=1,n3h0
    dtot_p=dtot_p - DBLE(div(1,1,izh0)*CONJG(div(1,1,izh0)))
 enddo

 dtot_p=dtot_p/(1.D0*n3)



 !Sum results from each processor           


 call mpi_reduce(dtot_p,dtot, 1,MPI_DOUBLE,   MPI_SUM,0,MPI_COMM_WORLD,ierror)



 !Normalize with enstrophy... (store enstrophy in dtot_p)                                                                                                                                                                             
! call enstrophy(uk,vk,wk,dtot_p)




 !Compute the domain maximum of divergence in real-space!
 !------------------------------------------------------!

 call fft_c2r(div,divr,n3h0)

 dmax=0.
 dmax_p=0.

 !For each processor, loop over all the values and find the max absolute value
 do izh0=1,n3h0
    do ix=1,n1d
       do iy=1,n2d
          if(ix<=n1) then

             if( dmax_p(1) < real(abs(divr(ix,iy,izh0))) )   then

                dmax_p(1) = real(abs(divr(ix,iy,izh0)))
                dmax_p(2) = ix
                dmax_p(3) = iy
                dmax_p(4) = izh0*(mype+1)
                
             end if


          end if
       end do
    end do
 end do

 !Now find the max amongst all processors
 if(mype>0) call mpi_send(dmax_p,4,MPI_REAL,0,tag_cont,MPI_COMM_WORLD,ierror)

 if(mype==0) then
    
    dmax = dmax_p   !initialize to min/max of mype 0                                                                                                                         
    
    do nrec=1,npe-1
       
       call mpi_recv(dmax_recv,4,MPI_REAL,MPI_ANY_SOURCE,tag_cont,MPI_COMM_WORLD,status,ierror)
       
       if(dmax(1) < dmax_recv(1)) dmax = dmax_recv
       
    end do
        
 end if

 !Finally, print the domain maximum and the average of divergence.
 if(mype==0)   write(unit=unit_divtot,fmt=*) time,sqrt(dtot),sqrt(dmax(1))
 if(mype==0)   write(unit=unit_dwhere,fmt=*) time,dmax(2),dmax(3),dmax(4)


END SUBROUTINE continuity_anelastic






 SUBROUTINE energy_linear(uk,vk,wk,bk,ktot,ptot)  

   !*** INCOMPLETE:
   !    still has r_1 and r_2...


   !This subroutine computes the quadratic form of energy (conserved quantity for constant-N Boussinesq and all QG versions) 
   !It uses the mid-point (or rectangle rule) for the staggered fields (u,v,b)                                                                   
   !And the trapezoidal rule for w. There is therefore no interpolation needed.                                                                             

   !It is now used only to normalize energy when initializing. C'est un vestige evolutif.
   !It has inspired the main subroutine diag_zentrum, and is now pretty much contained into it.

   double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk

   real, dimension(n3h0) :: ks,ku,ps   !Staggered (s) and unstaggered (u) energy (k: kinetic, p: potential)                                                      

   real :: ktot_p,ktot,ptot_p,ptot

   !Let's first sum over k_h to get a function of z alone.                                                                                                         
   ! 1/(2pi^)^2 int(int( e^2(x,y,z)/2 dx)dy = 1/2 sum over kx ky |eh|^2(z)                                                                                        

   ks=0.
   ku=0.
   ps=0.
   ktot_p=0.
   ktot=0.
   ptot_p=0.
   ptot=0.

   !With dealiasing, sum_k 1/2 |u(kx,ky,z)|^2 = sum_k L |u|^2 - 0.5 |u(0,0,z)|^2                                                                          

    do iz=1,n3h0
       izh2=iz+2

       do iky=1,ikty
          do ikx=1,iktx

             if(L(ikx,iky)==1) then
                ks(iz) = ks(iz) + real( uk(ikx,iky,izh2)*CONJG( uk(ikx,iky,izh2) ) + vk(ikx,iky,izh2)*CONJG( vk(ikx,iky,izh2) ) )
                ku(iz) = ku(iz) + real( wk(ikx,iky,izh2)*CONJG( wk(ikx,iky,izh2) ) )*Ar2
                ps(iz) = ps(iz) + real( bk(ikx,iky,izh2)*CONJG( bk(ikx,iky,izh2) ) )*(Bu*r_1(izh2)/r_2(izh2))
             end if

          enddo
       enddo

       !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                             

       ks(iz) = ks(iz) - 0.5*real( uk(1,1,izh2)*CONJG( uk(1,1,izh2) ) + vk(1,1,izh2)*CONJG( vk(1,1,izh2) ) )
       ku(iz) = ku(iz) - 0.5*real( wk(1,1,izh2)*CONJG( wk(1,1,izh2) ) )*Ar2
       ps(iz) = ps(iz) - 0.5*real( bk(1,1,izh2)*CONJG( bk(1,1,izh2) ) )*(Bu*r_1(izh2)/r_2(izh2))

    end do


    !If desired, we can plot energy as a function of z
!    if(out_ez ==1 .and. mod(iter,freq_ez)==0 .and. iter>0) call plot_ez(ks,ku,pu)   NOT WORKING WITH STAG B I GUESS...

    !1/L int(e(z) dz) becomes just a dummy sum over all grid points without interpolation needed...                                                            
    !First some locally to each processor                                                                                                                         

    do izh0=1,n3h0
       izh2=izh0+2

       ktot_p=ktot_p + rho_s(izh2)*ks(izh0) + rho_u(izh2)*ku(izh0)
       ptot_p=ptot_p                        + rho_s(izh2)*ps(izh0)

    end do

    ktot_p=ktot_p/n3
    ptot_p=ptot_p/n3


!    if(N_2<1e-15) ptot_p=0.

    !Sum results from each processor                                                                                                                         

    call mpi_reduce(ktot_p,ktot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call mpi_reduce(ptot_p,ptot, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(mype==0 .and. iter>0) write(unit=unit_energy,fmt=*) time,ktot,ptot,ktot+ptot

  END SUBROUTINE energy_linear







SUBROUTINE cond_integrability(uk,vk,wk,bk)

      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk

      double complex,   dimension(iktx,ikty,n3h0) :: c1k,c2k,c3k,c4k
      double precision, dimension(n1d,n2d,n3h0)   :: c1r,c2r,c3r,c4r

      double precision :: myz !Position of stag or unstag vertical level

      real :: m_cond(3,2,2),m_cond_recv(3,2,2),minmax(3,2,2)  !min(X,1,Y) and max(x,2,Y) for condition X. Y=1 means whole domain, Y=2 is for the jump region only.
      real :: fraction_domain_p(3,2),fraction_domain(3,2)  !Fraction(X,Y) of the domain where condition X is met. (divided _p or not amongst processors)

      integer :: nrec,dummy_count_1, dummy_count_2

      equivalence(c1r,c1k)
      equivalence(c2r,c2k)
      equivalence(c3r,c3k)
      equivalence(c4r,c4k)


      do izh0=1,n3h0
         izh2=izh0+2

         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)

               if(L(ikx,iky)==1) then

                  if(where_bz==unstag) then
                     c1k(ikx,iky,izh0) = (1/r_2(izh2))*(bk(ikx,iky,izh2+1)-bk(ikx,iky,izh2))/dz      
                  elseif(where_bz==stag) then
                     c1k(ikx,iky,izh0) = (1/r_2s(izh2) )*(bk(ikx,iky,izh2+1)-bk(ikx,iky,izh2-1))/(2.*dz) 
                  end if

                  c2k(ikx,iky,izh0) = i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)                !vertical vorticity
                  c3k(ikx,iky,izh0) = i*kx*uk(ikx,iky,izh2) - i*ky*vk(ikx,iky,izh2)                !u_x - v_y
                  c4k(ikx,iky,izh0) = i*kx*vk(ikx,iky,izh2) + i*ky*uk(ikx,iky,izh2)                !v_x + u_y

               else

                  c1k(ikx,iky,izh0) = (0.D0,0.D0)
                  c2k(ikx,iky,izh0) = (0.D0,0.D0)
                  c3k(ikx,iky,izh0) = (0.D0,0.D0)
                  c4k(ikx,iky,izh0) = (0.D0,0.D0)

               end if

            enddo
         enddo
      end do

      !Boundary conditions.
      if(where_bz==stag) then
         if(mype==0) then
            
            do iky=1,ikty
               do ikx=1,iktx
                  
                  if(L(ikx,iky)==1) then
                     
                     c1k(ikx,iky,1) = (1/r_2s(izbot2))*(bk(ikx,iky,izbot2+1) - bk(ikx,iky,izbot2)   )/(2.D0*dz)
                     
                  else
                     
                     c1k(ikx,iky,1) = (0.D0,0.D0)
                  endif
               enddo
            enddo

         elseif(mype==(npe-1)) then
            
            do iky=1,ikty
               do ikx=1,iktx
                  
                  if(L(ikx,iky)==1) then
                     
                     c1k(ikx,iky,n3h0) = (1/r_2s(iztop2))*(bk(ikx,iky,iztop2) - bk(ikx,iky,iztop2-1)   )/(2.D0*dz)
                     
                  else
                     
                     c1k(ikx,iky,n3h0) = (0.D0,0.D0)
                  endif
               enddo
            enddo

         end if
      elseif(where_bz==unstag) then
         if(mype==(npe-1)) then
               
            do iky=1,ikty
               do ikx=1,iktx
                     c1k(ikx,iky,n3h0) = (0.D0,0.D0)
               enddo
            enddo

         end if
      end if

      !Bring to r-space.

      call fft_c2r(c1k,c1r,n3h0)
      call fft_c2r(c2k,c2r,n3h0)
      call fft_c2r(c3k,c3r,n3h0)
      call fft_c2r(c4k,c4r,n3h0)


      !Compute the conditions:

      !1. Gravitational (N^2<0): 1 + Fr^2/Ro 1/r2 db/dz < 0
      !2. Centrifugal/inertial (A<0): 1 + Ro*Xi < 0  
      !3. Anticyclonic ageostrophic baroclinic (A-S<0): 1 + Ro*(Xi+S) < 0 where S = sqrt( (ux-vy)^2 + (vx+uy)^2)    (in principle this should be along isopycnals...)

      dummy_count_1 = 0
      dummy_count_2 = 0
      fraction_domain_p  = 0.

      m_cond(1,1,1)=1.
      m_cond(1,2,1)=1.
      m_cond(2,1,1)=1.
      m_cond(2,2,1)=1.
      m_cond(3,1,1)=1. 
      m_cond(3,2,1)=1. 

      m_cond(1,1,2)=1.
      m_cond(1,2,2)=1.
      m_cond(2,1,2)=1.
      m_cond(2,2,2)=1.
      m_cond(3,1,2)=1. 
      m_cond(3,2,2)=1. 


   do izh0=1,n3h0
         izh2=izh0+2
         do ix=1,n1d
            do iy=1,n2d
               if(ix<=n1) then
                  
                   c1r(ix,iy,izh0) =       1. + (Fr*Fr/Ro)*c1r(ix,iy,izh0)
                   c3r(ix,iy,izh0) =       1. + Ro*(c2r(ix,iy,izh0) - sqrt( c3r(ix,iy,izh0)*c3r(ix,iy,izh0) + c4r(ix,iy,izh0)*c4r(ix,iy,izh0)))
                   c2r(ix,iy,izh0) =       1. + Ro*c2r(ix,iy,izh0)


                   !Update min/max for the whole domain
                   if( c1r(ix,iy,izh0) < m_cond(1,1,1)) m_cond(1,1,1) = c1r(ix,iy,izh0)
                   if( c2r(ix,iy,izh0) < m_cond(2,1,1)) m_cond(2,1,1) = c2r(ix,iy,izh0)
                   if( c3r(ix,iy,izh0) < m_cond(3,1,1)) m_cond(3,1,1) = c3r(ix,iy,izh0)

                   if( c1r(ix,iy,izh0) > m_cond(1,2,1)) m_cond(1,2,1) = c1r(ix,iy,izh0)
                   if( c2r(ix,iy,izh0) > m_cond(2,2,1)) m_cond(2,2,1) = c2r(ix,iy,izh0)
                   if( c3r(ix,iy,izh0) > m_cond(3,2,1)) m_cond(3,2,1) = c3r(ix,iy,izh0)

                   !Update fraction of the whole domain where conditions are met.
                   if( c1r(ix,iy,izh0) < 0) fraction_domain_p(1,1)=fraction_domain_p(1,1)+1.
                   if( c2r(ix,iy,izh0) < 0) fraction_domain_p(2,1)=fraction_domain_p(2,1)+1.
                   if( c3r(ix,iy,izh0) < 0) fraction_domain_p(3,1)=fraction_domain_p(3,1)+1.                  

                   !Jump-region quantities!
                   !----------------------!

                   !Condition 1 may be computed on staggered or unstaggered grid points
                   if(where_bz==stag) then
                      myz=zash0(izh0)
                   elseif(where_bz==unstag) then
                      myz=zah0(izh0)
                   end if
                      
                   if(abs(z0-myz) <= jump_region_width*H_N) then
                   
                      !Update min/max for the jump region
                      if( c1r(ix,iy,izh0) < m_cond(1,1,2)) m_cond(1,1,2) = c1r(ix,iy,izh0)                      
                      if( c1r(ix,iy,izh0) > m_cond(1,2,2)) m_cond(1,2,2) = c1r(ix,iy,izh0)

                      !Update fraction of the jump region where conditions are met.
                      if( c1r(ix,iy,izh0) < 0) fraction_domain_p(1,2)=fraction_domain_p(1,2)+1.

                      dummy_count_1 = dummy_count_1 + 1 !How many points were in the jump region

                   end if

                   !Conditions 2-3 are computed on staggered grid points
                   myz = zash0(izh0)

                   if(abs(z0-myz) <= jump_region_width*H_N) then
                   
                      !Update min/max for the jump region
                      if( c2r(ix,iy,izh0) < m_cond(2,1,2)) m_cond(2,1,2) = c2r(ix,iy,izh0)
                      if( c3r(ix,iy,izh0) < m_cond(3,1,2)) m_cond(3,1,2) = c3r(ix,iy,izh0)
                      
                      if( c2r(ix,iy,izh0) > m_cond(2,2,2)) m_cond(2,2,2) = c2r(ix,iy,izh0)
                      if( c3r(ix,iy,izh0) > m_cond(3,2,2)) m_cond(3,2,2) = c3r(ix,iy,izh0)
                      
                      !Update fraction of the jump region where conditions are met.
                      if( c2r(ix,iy,izh0) < 0) fraction_domain_p(2,2)=fraction_domain_p(2,2)+1.
                      if( c3r(ix,iy,izh0) < 0) fraction_domain_p(3,2)=fraction_domain_p(3,2)+1.                                     

                      dummy_count_2 = dummy_count_2 + 1  !How many points were in the jump region

                   end if


                else

                   c1r(ix,iy,izh0) =       0.
                   c3r(ix,iy,izh0) =       0.
                   c2r(ix,iy,izh0) =       0.

                endif
             end do
          end do
       end do

       !Normalize fractions
       fraction_domain_p(1,1) = fraction_domain_p(1,1)/(n1*n2*n3) !divide by the number of grid cells
       fraction_domain_p(2,1) = fraction_domain_p(2,1)/(n1*n2*n3) !divide by the number of grid cells
       fraction_domain_p(3,1) = fraction_domain_p(3,1)/(n1*n2*n3) !divide by the number of grid cells

        if(dummy_count_1 > 0) then
          fraction_domain_p(1,2) = fraction_domain_p(1,2)/dummy_count_1
       end if
       if(dummy_count_2 > 0) then
          fraction_domain_p(2,2) = fraction_domain_p(2,2)/dummy_count_2
          fraction_domain_p(3,2) = fraction_domain_p(3,2)/dummy_count_2
       end if

       !Add fraction of domain and growth for each processor.
       call mpi_reduce(fraction_domain_p,fraction_domain, 6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

       if(mype>0) call mpi_send(m_cond,12,MPI_REAL,0,tag_cond,MPI_COMM_WORLD,ierror)

       if(mype==0) then

          minmax = m_cond   !initialize to min/max of mype 0

          do nrec=1,npe-1

             call mpi_recv(m_cond_recv,12,MPI_REAL,MPI_ANY_SOURCE,tag_cond,MPI_COMM_WORLD,status,ierror)

             if(m_cond_recv(1,1,1) < minmax(1,1,1))  minmax(1,1,1) = m_cond_recv(1,1,1)
             if(m_cond_recv(2,1,1) < minmax(2,1,1))  minmax(2,1,1) = m_cond_recv(2,1,1)
             if(m_cond_recv(3,1,1) < minmax(3,1,1))  minmax(3,1,1) = m_cond_recv(3,1,1)

             if(m_cond_recv(1,2,1) > minmax(1,2,1))  minmax(1,2,1) = m_cond_recv(1,2,1)
             if(m_cond_recv(2,2,1) > minmax(2,2,1))  minmax(2,2,1) = m_cond_recv(2,2,1)
             if(m_cond_recv(3,2,1) > minmax(3,2,1))  minmax(3,2,1) = m_cond_recv(3,2,1)

             if(m_cond_recv(1,1,2) < minmax(1,1,2))  minmax(1,1,2) = m_cond_recv(1,1,2)
             if(m_cond_recv(2,1,2) < minmax(2,1,2))  minmax(2,1,2) = m_cond_recv(2,1,2)
             if(m_cond_recv(3,1,2) < minmax(3,1,2))  minmax(3,1,2) = m_cond_recv(3,1,2)

             if(m_cond_recv(1,2,2) > minmax(1,2,2))  minmax(1,2,2) = m_cond_recv(1,2,2)
             if(m_cond_recv(2,2,2) > minmax(2,2,2))  minmax(2,2,2) = m_cond_recv(2,2,2)
             if(m_cond_recv(3,2,2) > minmax(3,2,2))  minmax(3,2,2) = m_cond_recv(3,2,2)



          end do

      !Plot results.  

          !minmax(condition,min or max,whole or jump domain) 
          !fraction_domain(condition,  whole or jump domain) 

          write(unit=unit_cond1,fmt=*) time,minmax(1,1,1),minmax(1,2,1),fraction_domain(1,1)
          write(unit=unit_cond2,fmt=*) time,minmax(2,1,1),minmax(2,2,1),fraction_domain(2,1)
          write(unit=unit_cond3,fmt=*) time,minmax(3,1,1),minmax(3,2,1),fraction_domain(3,1)

          write(unit=unit_cond1j,fmt=*) time,minmax(1,1,2),minmax(1,2,2),fraction_domain(1,2)
          write(unit=unit_cond2j,fmt=*) time,minmax(2,1,2),minmax(2,2,2),fraction_domain(2,2)
          write(unit=unit_cond3j,fmt=*) time,minmax(3,1,2),minmax(3,2,2),fraction_domain(3,2)

       end if

END SUBROUTINE cond_integrability

       




SUBROUTINE compute_growth(uk,vk,bk) !Compute the horizontally-averaged growth rates of static and inertial instability at every vertical level

      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,bk

      double complex,   dimension(iktx,ikty,n3h0) :: c1k,c2k
      double precision, dimension(n1d,n2d,n3h0)   :: c1r,c2r

      real, dimension (n3h0,2) :: growth_p
      real, dimension (n3,2) :: growth 

      integer :: nrec, izp
      integer :: processor
   
      equivalence(c1r,c1k)
      equivalence(c2r,c2k)

      do izh0=1,n3h0
         izh2=izh0+2

         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)

               if(L(ikx,iky)==1) then

                  if(where_bz==unstag) then
                     c1k(ikx,iky,izh0) = (1/r_2(izh2))*(bk(ikx,iky,izh2+1)-bk(ikx,iky,izh2))/dz      
                  elseif(where_bz==stag) then
                     c1k(ikx,iky,izh0) = (1/r_2s(izh2) )*(bk(ikx,iky,izh2+1)-bk(ikx,iky,izh2-1))/(2.*dz) 
                  end if

                  c2k(ikx,iky,izh0) = i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)                !vertical vorticity

               else

                  c1k(ikx,iky,izh0) = (0.D0,0.D0)
                  c2k(ikx,iky,izh0) = (0.D0,0.D0)

               end if

            enddo
         enddo
      end do

      !Boundary conditions.
      if(where_bz==stag) then
         if(mype==0) then
            
            do iky=1,ikty
               do ikx=1,iktx
                  
                  if(L(ikx,iky)==1) then
                     
                     c1k(ikx,iky,1) = (1/r_2s(izbot2))*(bk(ikx,iky,izbot2+1) - bk(ikx,iky,izbot2)   )/(2.D0*dz)
                     
                  else
                     
                     c1k(ikx,iky,1) = (0.D0,0.D0)
                  endif
               enddo
            enddo

         elseif(mype==(npe-1)) then
            
            do iky=1,ikty
               do ikx=1,iktx
                  
                  if(L(ikx,iky)==1) then
                     
                     c1k(ikx,iky,n3h0) = (1/r_2s(iztop2))*(bk(ikx,iky,iztop2) - bk(ikx,iky,iztop2-1)   )/(2.D0*dz)
                     
                  else
                     
                     c1k(ikx,iky,n3h0) = (0.D0,0.D0)
                  endif
               enddo
            enddo

         end if
      elseif(where_bz==unstag) then
         if(mype==(npe-1)) then
               
            do iky=1,ikty
               do ikx=1,iktx
                     c1k(ikx,iky,n3h0) = (0.D0,0.D0)
               enddo
            enddo

         end if
      end if

      !Bring to r-space.

      call fft_c2r(c1k,c1r,n3h0)
      call fft_c2r(c2k,c2r,n3h0)


      !Compute the conditions:

      !1. Gravitational (N^2<0): 1 + Fr^2/Ro 1/r2 db/dz < 0
      !2. Centrifugal/inertial (A<0): 1 + Ro*Xi < 0  
      
      growth   = 0.
      growth_p = 0.

      do izh0=1,n3h0
         izh2=izh0+2
         do ix=1,n1d
            do iy=1,n2d
               if(ix<=n1) then
                  
                   c1r(ix,iy,izh0) =       1. + (Fr*Fr/Ro)*c1r(ix,iy,izh0)
                   c2r(ix,iy,izh0) =       1. + Ro*c2r(ix,iy,izh0)

                   !Compute growth rate
                   if( c1r(ix,iy,izh0) < 0) then
                      if(where_bz==stag)   growth_p(izh0,1)=growth_p(izh0,1) + sqrt(-r_1s(izh2)*r_2s(izh2)*c1r(ix,iy,izh0)/(Ar2*Fr*Fr)) 
                      if(where_bz==unstag) growth_p(izh0,1)=growth_p(izh0,1) + sqrt(-r_1(izh2) *r_2(izh2) *c1r(ix,iy,izh0)/(Ar2*Fr*Fr)) 
                   end if

                   if( c2r(ix,iy,izh0) < 0) then
                      growth_p(izh0,2)=growth_p(izh0,2) + sqrt(-c2r(ix,iy,izh0)) 
                   end if

                   
                endif
             end do
          end do
       end do

       growth_p = growth_p/(n1*n2) !divide by the number of grid cells 


       if(mype>0) call mpi_send(growth_p,n3h0*2,MPI_REAL,0,tag_grow,MPI_COMM_WORLD,ierror)

       if(mype==0) then

          !Copy mype == 0's growth 
          do iz=1,n3h0
             izp=iz
             
             growth(izp,1) = growth(izp,1) + growth_p(iz,1)
             growth(izp,2) = growth(izp,2) + growth_p(iz,2)
             
             
          end do


          do nrec=1,npe-1

             call mpi_recv(growth_p,n3h0*2,MPI_REAL,MPI_ANY_SOURCE,tag_grow,MPI_COMM_WORLD,status,ierror)

             !Copy received share onto err                                                                                                                                                                                                      
             processor=status(MPI_SOURCE)

             !Add the awaiting half level                                                                                                                                                                                                       
             do iz=1,n3h0
                izp=processor*n3h0+iz
                   
                growth(izp,1) = growth(izp,1) + growth_p(iz,1)
                growth(izp,2) = growth(izp,2) + growth_p(iz,2)
                   

             end do

          end do

          !PLOT!                                                                                                                                                                                                                        
          do iz=1,n3
             write(unit_grow,fmt=*) iz,growth(iz,1),growth(iz,2)      !iz,sigma_static, sigma_inertial                                                                                                                  
          enddo
          write(unit_grow,*) '           '
          call flush(unit_grow)


       end if

     END SUBROUTINE compute_growth






SUBROUTINE cond_wz(wak)

      !This subroutine computes the condition on w_z, namely  Ro w_1_z / sqrt( |u_x|^2 + |v_x|^2 ) << 1 to check if QG holds.
      !In anelastic, we have w_1_z --> 1/rho_0 d\dz (rho_0 w)
      !To simplify, we assume u'~l'~1 from the initial condition so that Ro is representative of the "true" flow-imposed velocity-based (macro) Ro.


      double complex, dimension(iktx,ikty,n3h1) :: wak

      double complex, dimension(iktx,ikty,n3h0) :: c1k
      double precision,    dimension(n1d,n2d,n3h0) :: c1r

      real :: m_cond(2),m_cond_recv(2),minmax(2)       !min(x,1) and max(x,2) for condition x.  4th is QG limiting assumption |Fr^2b_z/Ror_2|<<1 or in dim terms, r1bz<<N^2                                                                 
      real :: fraction_domain_p, fraction_domain
      integer :: nrec

      equivalence(c1r,c1k)

      fraction_domain_p=0.
      fraction_domain  =0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2

         do iky=1,ikty
            do ikx=1,iktx

               if(L(ikx,iky)==1) then


                  c1k(ikx,iky,izh0) = (rho_u(izh2)*wak(ikx,iky,izh1)-rho_u(izh2-1)*wak(ikx,iky,izh1-1))/(rho_s(izh2)*dz)       !1/rho_s (rho_u w)_z

               else

                  c1k(ikx,iky,izh0) = (0.D0,0.D0)
   
               end if

            enddo
         enddo
      end do

      !Boundary conditions.                                                                                                                                                                                                       

      if(mype==0) then

         do iky=1,ikty
            do ikx=1,iktx

               if(L(ikx,iky)==1) then

                  c1k(ikx,iky,1) = wak(ikx,iky,izbot1)/dz    ! r_3 w + dw/dz = dw/dz since w=0 at z=0

                         else

                  c1k(ikx,iky,1) = (0.D0,0.D0)
               endif
            enddo
         enddo

      end if


      !Bring to r-space.                                                                                                                                                                                                             

      call fft_c2r(c1k,c1r,n3h0)

      !Compute the conditions:                                                                                                                                                                                                       
      m_cond(1)= Ro*abs(c1r(1,1,1))
      m_cond(2)= Ro*abs(c1r(1,1,1))


      do izh0=1,n3h0
         do ix=1,n1d
            do iy=1,n2d
               if(ix<=n1) then

                   c1r(ix,iy,izh0) =      Ro*abs(c1r(ix,iy,izh0))


                   !Update min/max                                                                                                                                                                                                    
                   if( c1r(ix,iy,izh0) < m_cond(1)) m_cond(1) = c1r(ix,iy,izh0)
                   if( c1r(ix,iy,izh0) > m_cond(2)) m_cond(2) = c1r(ix,iy,izh0)

                   !Update fraction of the domain where conditions are met.                                                                                                                                                              
                   if( c1r(ix,iy,izh0) > 1) fraction_domain_p=fraction_domain_p+1.

                else

                   c1r(ix,iy,izh0) =       0.

                endif
             end do
          end do
       end do

       fraction_domain_p = fraction_domain_p/(n1*n2*n3) !divide by the number of grid cells                                                                                                                                                 

       !Add fraction of domain for each processor.                                                                                                                                                                                          
       call mpi_reduce(fraction_domain_p,fraction_domain, 1,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierror)

       if(mype>0) call mpi_send(m_cond,2,MPI_REAL,0,tag_condwz,MPI_COMM_WORLD,ierror)

       if(mype==0) then

          minmax = m_cond   !initialize to min/max of mype 0                                                                                                                                                                              

          do nrec=1,npe-1

             call mpi_recv(m_cond_recv,2,MPI_REAL,MPI_ANY_SOURCE,tag_condwz,MPI_COMM_WORLD,status,ierror)

             if(m_cond_recv(1) < minmax(1))  minmax(1) = m_cond_recv(1)
             if(m_cond_recv(2) > minmax(2))  minmax(2) = m_cond_recv(2)

          end do

      !Plot results.                                                                                                                                                                                                                         
                                      
          write(unit=unit_condwz,fmt=*) time,minmax(1),minmax(2),fraction_domain
       end if

     END SUBROUTINE cond_wz






  subroutine buoy_vert_scale(wk,bk) !Vertical scale of w (unstaggered) and b (staggered) at UNstaggered grid points. 

       double complex, dimension(iktx,ikty,n3h2) :: wk,bk

       real, dimension(n3h0) :: bsum_p,bzsum_p,wsum_p,wzsum_p 
       real, dimension(n3h0) :: hb_p,hb_r,hw_p,hw_r
       real, dimension(n3)   :: hb,hw

       integer :: processor,izp,nrec


        bsum_p=0.
       bzsum_p=0.
        wsum_p=0.
       wzsum_p=0.

       do iz=1,n3h0
          izh2=iz+2

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                    bsum_p(iz) =  bsum_p(iz) + real( ( bk(ikx,iky,izh2) + bk(ikx,iky,izh2+1) )   *CONJG(  bk(ikx,iky,izh2) + bk(ikx,iky,izh2+1)  ) )/4.
                   bzsum_p(iz) = bzsum_p(iz) + real( ( bk(ikx,iky,izh2+1) - bk(ikx,iky,izh2) )   *CONJG(  bk(ikx,iky,izh2+1) - bk(ikx,iky,izh2)  ) )/(dz*dz)
                    wsum_p(iz) =  wsum_p(iz) + real(   wk(ikx,iky,izh2)*CONJG( wk(ikx,iky,izh2) ) )
                   wzsum_p(iz) = wzsum_p(iz) + real( ( wk(ikx,iky,izh2+1)-wk(ikx,iky,izh2-1) )   *CONJG(  wk(ikx,iky,izh2+1)-wk(ikx,iky,izh2-1)  ) )/(4.*dz*dz)
                end if

             enddo
          enddo

          !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                                              
           bsum_p(iz) =  bsum_p(iz) - 0.5*real( ( bk(1,1,izh2) + bk(1,1,izh2+1) )   *CONJG(  bk(1,1,izh2) + bk(1,1,izh2+1)  ) )/4.
          bzsum_p(iz) = bzsum_p(iz) - 0.5*real( ( bk(1,1,izh2+1) - bk(1,1,izh2) )   *CONJG(  bk(1,1,izh2+1) - bk(1,1,izh2)  ) )/(dz*dz)
           wsum_p(iz) =  wsum_p(iz) - 0.5*real( wk(1,1,izh2)*CONJG( wk(1,1,izh2) ) )
          wzsum_p(iz) = wzsum_p(iz) - 0.5*real( (wk(1,1,izh2+1)-wk(1,1,izh2-1))   *CONJG( (wk(1,1,izh2+1)-wk(1,1,izh2-1)) ) )/(4.*dz*dz)

       end do




       !Correction at top and bot for bz
       if(mype==0) then

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   wzsum_p(1) = wzsum_p(1) + real(  wk(ikx,iky,izbot2+1)   *CONJG(  wk(ikx,iky,izbot2+1) ))/(4.*dz*dz)
                end if

             enddo
          enddo

          wzsum_p(1) = wzsum_p(1) - 0.5*real( wk(1,1,izbot2+1)  *CONJG( wk(1,1,izbot2+1) ) )/(4.*dz*dz)


       elseif(mype==(npe-1)) then

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                    bsum_p(n3h0) =  bsum_p(n3h0) + real(  bk(ikx,iky,iztop2)    *CONJG(  bk(ikx,iky,iztop2) ))
                   wzsum_p(n3h0) = wzsum_p(n3h0) + real(wk(ikx,iky,iztop2-1)   *CONJG(wk(ikx,iky,iztop2-1) ))/(dz*dz)
                end if

             enddo
          enddo

           bsum_p(n3h0) = bsum_p(n3h0) - 0.5*real(  bk(1,1,iztop2)    *CONJG(  bk(1,1,iztop2) ))
          bzsum_p(n3h0) = 0.
           wsum_p(n3h0) = 0.
          wzsum_p(n3h0) = wzsum_p(n3h0)- 0.5*real(wk(1,1,iztop2-1)   *CONJG(wk(1,1,iztop2-1) ))/(dz*dz)


       end if





       !Compute the vertical scale at each vertical level
       do iz=1,n3h0
          hb_p(iz) = sqrt(bsum_p(iz)/bzsum_p(iz))
          hw_p(iz) = sqrt(wsum_p(iz)/wzsum_p(iz))
       end do


    !Share to other processors                                                                                                                             
    if(mype>0) call mpi_send(hb_p,n3h0,MPI_REAL,0,tag_hb2,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(hw_p,n3h0,MPI_REAL,0,tag_hw2,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part 
       do iz=1,n3h0
             hb(iz) = hb_p(iz)
             hw(iz) = hw_p(iz)
       end do


       !Receive other parts from other processors                                                                                                                             
       do nrec=1,npe-1

          call mpi_recv(hb_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hb2,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hb(izp) =  hb_r(iz)

          end do

          call mpi_recv(hw_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hw2,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                          
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hw(izp) =  hw_r(iz)

          end do

       end do

       !PLOT a z-profile!                                                                                                                                                       
     do iz=1,n3
        write(unit_hb2,fmt=*) za(iz),hb(iz),hw(iz)
     enddo
     write(unit_hb2,*) '           '
     call flush(unit_hb2)

     !Plot the tropopause value!
     write(unit=unit_hbt2,fmt=*) time,hb(n3/2),hw(n3/2)


   end if
       

     end subroutine buoy_vert_scale





  subroutine buoy_vert_scale_rot(wak,tk) !Version for unstaggered tk (or b_rot)  --- same as original version, but for izh1 field (b_rot)

       double complex, dimension(iktx,ikty,n3h1) :: tk
       double complex, dimension(iktx,ikty,n3h1) :: wak

       real, dimension(n3h0) :: bsum_p,bzsum_p,wsum_p,wzsum_p 
       real, dimension(n3h0) :: hb_p,hb_r,hw_p,hw_r
       real, dimension(n3)   :: hb,hw

       integer :: processor,izp,nrec


        bsum_p=0.
       bzsum_p=0.
        wsum_p=0.
       wzsum_p=0.

       do iz=1,n3h0
          izh1=iz+1

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                    bsum_p(iz) =  bsum_p(iz) + real( tk(ikx,iky,izh1)*CONJG( tk(ikx,iky,izh1) ) )
                   bzsum_p(iz) = bzsum_p(iz) + real( (tk(ikx,iky,izh1+1)-tk(ikx,iky,izh1-1))   *CONJG( (tk(ikx,iky,izh1+1)-tk(ikx,iky,izh1-1)) ) )/(4.*dz*dz)
                    wsum_p(iz) =  wsum_p(iz) + real( wak(ikx,iky,izh1)*CONJG( wak(ikx,iky,izh1) ) )
                   wzsum_p(iz) = wzsum_p(iz) + real( (wak(ikx,iky,izh1+1)-wak(ikx,iky,izh1-1))   *CONJG( (wak(ikx,iky,izh1+1)-wak(ikx,iky,izh1-1)) ) )/(4.*dz*dz)
                end if

             enddo
          enddo

          !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                                              
           bsum_p(iz) =  bsum_p(iz) - 0.5*real( tk(1,1,izh1)*CONJG( tk(1,1,izh1) ) )
          bzsum_p(iz) = bzsum_p(iz) - 0.5*real( (tk(1,1,izh1+1)-tk(1,1,izh1-1))   *CONJG( (tk(1,1,izh1+1)-tk(1,1,izh1-1)) ) )/(4.*dz*dz)
           wsum_p(iz) =  wsum_p(iz) - 0.5*real( wak(1,1,izh1)*CONJG( wak(1,1,izh1) ) )
          wzsum_p(iz) = wzsum_p(iz) - 0.5*real( (wak(1,1,izh1+1)-wak(1,1,izh1-1))   *CONJG( (wak(1,1,izh1+1)-wak(1,1,izh1-1)) ) )/(4.*dz*dz)

       end do




       !Correction at top and bot for bz
       if(mype==0) then
       
          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   bzsum_p(1) = bzsum_p(1) + real( tk(ikx,iky,izbot1+1)   *CONJG( tk(ikx,iky,izbot1+1) ))/(4.*dz*dz)
                   wzsum_p(1) = wzsum_p(1) + real(wak(ikx,iky,izbot1+1)   *CONJG(wak(ikx,iky,izbot1+1) ))/(4.*dz*dz)
                end if

             enddo
          enddo

          bzsum_p(1) = bzsum_p(1) - 0.5*real( tk(1,1,izbot1+1)   *CONJG( tk(1,1,izbot1+1) ))/(4.*dz*dz)
          wzsum_p(1) = wzsum_p(1) - 0.5*real(wak(1,1,izbot1+1)   *CONJG(wak(1,1,izbot1+1) ))/(4.*dz*dz)

       elseif(mype==(npe-1)) then

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   bzsum_p(n3h0) = bzsum_p(n3h0) + real( tk(ikx,iky,iztop1-1)   *CONJG( tk(ikx,iky,iztop1-1) ))/(dz*dz)
                   wzsum_p(n3h0) = wzsum_p(n3h0) + real(wak(ikx,iky,iztop1-1)   *CONJG(wak(ikx,iky,iztop1-1) ))/(dz*dz)
                end if

             enddo
          enddo

          bzsum_p(n3h0) = bzsum_p(n3h0) - 0.5*real( tk(1,1,iztop1-1)   *CONJG( tk(1,1,iztop1-1) ))/(dz*dz)
          wzsum_p(n3h0) = wzsum_p(n3h0) - 0.5*real(wak(1,1,iztop1-1)   *CONJG(wak(1,1,iztop1-1) ))/(dz*dz)

       end if


       !Compute the vertical scale at each vertical level
       do iz=1,n3h0
          hb_p(iz) = sqrt(bsum_p(iz)/bzsum_p(iz))
          hw_p(iz) = sqrt(wsum_p(iz)/wzsum_p(iz))
       end do


    !Share to other processors                                                                                                                             
    if(mype>0) call mpi_send(hb_p,n3h0,MPI_REAL,0,tag_hb,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(hw_p,n3h0,MPI_REAL,0,tag_hw,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part 
       do iz=1,n3h0
             hb(iz) = hb_p(iz)
             hw(iz) = hw_p(iz)
       end do


       !Receive other parts from other processors                                                                                                                             
       do nrec=1,npe-1

          call mpi_recv(hb_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hb,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hb(izp) =  hb_r(iz)

          end do

          call mpi_recv(hw_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hw,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                          
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hw(izp) =  hw_r(iz)

          end do

       end do

       !PLOT a z-profile!                                                                                                                                                       
     do iz=1,n3
        write(unit_hb,fmt=*) za(iz),hb(iz),hw(iz)
     enddo
     write(unit_hb,*) '           '
     call flush(unit_hb)

     !Plot the tropopause value!
     write(unit=unit_hbt,fmt=*) time,hb(n3/2),hw(n3/2)


   end if
       

 end subroutine buoy_vert_scale_rot



  subroutine buoy_vert_scale_original(wak,tk) !Version for unstaggered tk.

       double complex, dimension(iktx,ikty,n3h2) :: tk
       double complex, dimension(iktx,ikty,n3h1) :: wak

       real, dimension(n3h0) :: bsum_p,bzsum_p,wsum_p,wzsum_p 
       real, dimension(n3h0) :: hb_p,hb_r,hw_p,hw_r
       real, dimension(n3)   :: hb,hw

       integer :: processor,izp,nrec


        bsum_p=0.
       bzsum_p=0.
        wsum_p=0.
       wzsum_p=0.

       do iz=1,n3h0
          izh1=iz+1
          izh2=iz+2

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                    bsum_p(iz) =  bsum_p(iz) + real( tk(ikx,iky,izh2)*CONJG( tk(ikx,iky,izh2) ) )
                   bzsum_p(iz) = bzsum_p(iz) + real( (tk(ikx,iky,izh2+1)-tk(ikx,iky,izh2-1))   *CONJG( (tk(ikx,iky,izh2+1)-tk(ikx,iky,izh2-1)) ) )/(4.*dz*dz)
                    wsum_p(iz) =  wsum_p(iz) + real( wak(ikx,iky,izh1)*CONJG( wak(ikx,iky,izh1) ) )
                   wzsum_p(iz) = wzsum_p(iz) + real( (wak(ikx,iky,izh1+1)-wak(ikx,iky,izh1-1))   *CONJG( (wak(ikx,iky,izh1+1)-wak(ikx,iky,izh1-1)) ) )/(4.*dz*dz)
                end if

             enddo
          enddo

          !See diary Feb 10th 2014 or journal Oct 22nd 2013. Adjust sum to 2D dialiasing: substract half the kh=0 mode.                                                              
           bsum_p(iz) =  bsum_p(iz) - 0.5*real( tk(1,1,izh2)*CONJG( tk(1,1,izh2) ) )
          bzsum_p(iz) = bzsum_p(iz) - 0.5*real( (tk(1,1,izh2+1)-tk(1,1,izh2-1))   *CONJG( (tk(1,1,izh2+1)-tk(1,1,izh2-1)) ) )/(4.*dz*dz)
           wsum_p(iz) =  wsum_p(iz) - 0.5*real( wak(1,1,izh1)*CONJG( wak(1,1,izh1) ) )
          wzsum_p(iz) = wzsum_p(iz) - 0.5*real( (wak(1,1,izh1+1)-wak(1,1,izh1-1))   *CONJG( (wak(1,1,izh1+1)-wak(1,1,izh1-1)) ) )/(4.*dz*dz)

       end do




       !Correction at top and bot for bz
       if(mype==0) then
       
          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   bzsum_p(1) = bzsum_p(1) + real( tk(ikx,iky,izbot2+1)   *CONJG( tk(ikx,iky,izbot2+1) ))/(4.*dz*dz)
                   wzsum_p(1) = wzsum_p(1) + real(wak(ikx,iky,izbot1+1)   *CONJG(wak(ikx,iky,izbot1+1) ))/(4.*dz*dz)
                end if

             enddo
          enddo

          bzsum_p(1) = bzsum_p(1) - 0.5*real( tk(1,1,izbot2+1)   *CONJG( tk(1,1,izbot2+1) ))/(4.*dz*dz)
          wzsum_p(1) = wzsum_p(1) - 0.5*real(wak(1,1,izbot1+1)   *CONJG(wak(1,1,izbot1+1) ))/(4.*dz*dz)

       elseif(mype==(npe-1)) then

          do iky=1,ikty
             do ikx=1,iktx

                if(L(ikx,iky)==1) then
                   bzsum_p(n3h0) = bzsum_p(n3h0) + real( tk(ikx,iky,iztop2-1)   *CONJG( tk(ikx,iky,iztop2-1) ))/(dz*dz)
                   wzsum_p(n3h0) = wzsum_p(n3h0) + real(wak(ikx,iky,iztop1-1)   *CONJG(wak(ikx,iky,iztop1-1) ))/(dz*dz)
                end if

             enddo
          enddo

          bzsum_p(n3h0) = bzsum_p(n3h0) - 0.5*real( tk(1,1,iztop2-1)   *CONJG( tk(1,1,iztop2-1) ))/(dz*dz)
          wzsum_p(n3h0) = wzsum_p(n3h0) - 0.5*real(wak(1,1,iztop1-1)   *CONJG(wak(1,1,iztop1-1) ))/(dz*dz)

       end if





       !Compute the vertical scale at each vertical level
       do iz=1,n3h0
          hb_p(iz) = sqrt(bsum_p(iz)/bzsum_p(iz))
          hw_p(iz) = sqrt(wsum_p(iz)/wzsum_p(iz))
       end do


    !Share to other processors                                                                                                                             
    if(mype>0) call mpi_send(hb_p,n3h0,MPI_REAL,0,tag_hb,MPI_COMM_WORLD,ierror)
    if(mype>0) call mpi_send(hw_p,n3h0,MPI_REAL,0,tag_hw,MPI_COMM_WORLD,ierror)

    if(mype==0) then

       !Copy mype==0 part 
       do iz=1,n3h0
             hb(iz) = hb_p(iz)
             hw(iz) = hw_p(iz)
       end do


       !Receive other parts from other processors                                                                                                                             
       do nrec=1,npe-1

          call mpi_recv(hb_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hb,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                        
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hb(izp) =  hb_r(iz)

          end do

          call mpi_recv(hw_r,n3h0,MPI_REAL,MPI_ANY_SOURCE,tag_hw,MPI_COMM_WORLD,status,ierror)

          !Copy received share onto err                                                                                                                                          
          processor=status(MPI_SOURCE)

          do iz=1,n3h0
             izp=processor*n3h0+iz

             hw(izp) =  hw_r(iz)

          end do

       end do

       !PLOT a z-profile!                                                                                                                                                       
     do iz=1,n3
        write(unit_hb,fmt=*) za(iz),hb(iz),hw(iz)
     enddo
     write(unit_hb,*) '           '
     call flush(unit_hb)

     !Plot the tropopause value!
     write(unit=unit_hbt,fmt=*) time,hb(n3/2),hw(n3/2)


   end if
       

 end subroutine buoy_vert_scale_original



 

    SUBROUTINE  print_slab(uk,vk)

      double complex, dimension(iktx,ikty,n3h2) :: uk,vk
      
      double precision :: kh
      integer :: mode

      integer :: z_level = n3h0+2

      do ikx=1,iktx
         kx=kxa(ikx)
         do iky=1,ikty
            ky=kya(iky)
            kh2=kx*kx+ky*ky

            kh   = sqrt(1.D0*kh2)
            mode = ifix(real(kh*L1/twopi+0.5))

            if ((L(ikx,iky).eq.1) .and. (mode >= 10 ) .and. (mode <= 40 )) then
               write(unit_slab,fmt=*) real(DBLE(uk(ikx,iky,z_level))),real(DIMAG((uk(ikx,iky,z_level)))),real(DBLE(vk(ikx,iky,z_level))),real(DIMAG((vk(ikx,iky,z_level))))
            end if
         end do
      end do
      
      write(unit_slab,*) '           '
      call flush(unit_slab)
      
    end SUBROUTINE print_slab



    subroutine slab_klist !Print the list of (ikx,iky) couples as they appear in print_slab. Just for the first time step
      
      integer :: unit_klist = 12345
      character(len = 32) :: fname         

      double precision :: kh
      integer :: mode

      write (fname, "(A9)") "klist.dat"
      open (unit=unit_klist,file=fname,action="write",status="replace")

      do ikx=1,iktx
         kx=kxa(ikx)
         do iky=1,ikty
            ky=kya(iky)
            kh2=kx*kx+ky*ky

            kh   = sqrt(1.D0*kh2)
            mode = ifix(real(kh*L1/twopi+0.5))

            if ((L(ikx,iky).eq.1) .and. (mode >= 10 ) .and. (mode <= 40 )) then
               write(unit_klist,fmt=*) ikx,iky,kxa(ikx),kya(iky)
            end if
         end do
      end do


      close (unit=unit_klist)

    end subroutine slab_klist



    subroutine tropopause_meanders(uk,vk,wk,bk,ur,vr,wr,br)

      double complex,   dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk
      double precision, dimension(n1d,n2d,n3h2)   :: ur,vr,wr,br

      double complex, dimension(iktx,ikty,n3h2) :: umem,vmem,bmem

      integer, dimension(n1,n2) :: min_level_p,min_level_r

      real :: cost_function
      
      integer :: nrec
      character(len = 32) :: fname                !File name 

      !Tropopause-height spectrum
      integer :: z_eta !Vertical level corresponding to the tropopause

      double complex,   dimension(iktx,ikty,3) :: uk_tp
      double precision, dimension(n1d,n2d,3)   :: ur_tp, ur_tp_r

      real, dimension(0:ktx) :: spz        !spectrum for kinetic (1) and potential (2)                                                                                                                      
      real, dimension(0:ktx) :: num_modes  !number of modes for a given kh                                                                                                         
      integer :: mode                      !Mode no  
      double precision :: kh


      equivalence(uk_tp,ur_tp)


      !Save the inverse FFT by saving bk...
      bmem = bk

      !Move to r-space
      call fft_c2r(bk,br,n3h2)


      !Initialize with largest possible value
      min_level_p = n3 


      !For each processor (except the ones at the top and bottom to avoid errors with delta sheets...)
!      if(mype > 0 .and. mype < (npe-1) ) then  
      !For only the middle processors
      if(mype >= npe/4 .and. mype < 3*npe/4 ) then  

         do ix=1,n1
            do iy=1,n2
               do izh0=1,n3h0
                  izh2=izh0+2
                     
                  cost_function = r_1(izh2)*r_2(izh2) + (Fr*Fr/Ro)*( br(ix,iy,izh2+1) - br(ix,iy,izh2) )/dz  - 1.
                  
                  if(cost_function > 0) then
                     min_level_p(ix,iy) = mype*n3h0 + izh0
                     exit
                  end if
                  
               end do
            end do
         end do


         !MOVE THE TROPOPAUSE DOWN FOR TEST
!         do ix=1,n1
!            do iy=1,n2
!               min_level_p(ix,iy) = min_level_p(ix,iy) - 40
!            end do
!         end do



         !Send your minimum level to processor 0 for comparison
         call mpi_send(min_level_p,n1*n2,MPI_INTEGER,0,tag_eta,MPI_COMM_WORLD,ierror)

      end if 


      !Find the lowest grid point for which the cost_function > 0.
      if(mype==0) then
         
!         do nrec=1,npe-2! All processors sent except top and bottom (so npe-2 in total)
         do nrec=1,npe/2! Receive only from the middle processors
            
            call mpi_recv(min_level_r,n1*n2,MPI_INTEGER,MPI_ANY_SOURCE,tag_eta,MPI_COMM_WORLD,status,ierror)

          
            do ix=1,n1
               do iy=1,n2
                  if(min_level_r(ix,iy) < min_level_p(ix,iy)) min_level_p(ix,iy) = min_level_r(ix,iy)
               end do
            end do
            
         end do
         

         
         !Print eta
         
         write (fname, "(A3,I2,A4)") "eta",count_eta,".dat"
         open (unit=count_eta,file=fname,action="write",status="replace")
         
         do iy=1,n2
            write(unit=count_eta,fmt=334) (min_level_p(ix,iy),ix=1,n1)
            write(unit=count_eta,fmt=*) '           '
         enddo

334    format(1x,I5.2,1x)         
         
         close (unit=count_eta)
         
      end if


      count_eta=count_eta+1

      bk = bmem
      

      !If desired, plot the tropopause-height kinetic energy spectrum

      if(out_tspec == 1) then
      
         !Broadcast the tropopause height
         call mpi_bcast(min_level_p,n1*n2,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
         
         umem = uk
         vmem = vk
         bmem = wk
         
         call fft_c2r(uk,ur,n3h2)
         call fft_c2r(vk,vr,n3h2)
         call fft_c2r(wk,wr,n3h2)
         
         
         ur_tp = 0.D0
         
          if(mype > 0 .and. mype <(npe-1)) then 
         !For only the middle processors                                                                                                                                                                                               
!         if(mype >= npe/4 .and. mype < 3*npe/4 ) then
            
            do ix=1,n1
               do iy=1,n2
                  
                  !Find out which mype's got the tropopause for every (x,y)
                  if(min_level_p(ix,iy) > mype*n3h0 .and.  min_level_p(ix,iy) <= (mype+1)*n3h0) then
                     
                     z_eta = min_level_p(ix,iz) - mype*n3h0+2           !z-level for a 2-level halo field
                     
                     ur_tp(ix,iy,1) = ur(ix,iy,z_eta)
                     ur_tp(ix,iy,2) = vr(ix,iy,z_eta)
                     ur_tp(ix,iy,3) = wr(ix,iy,z_eta)
                     
                  end if
                  
               end do
            end do
            
            !Send your ur_tp to mype 0
            call mpi_send(ur_tp,n1*n2*3,MPI_DOUBLE,0,tag_eta2,MPI_COMM_WORLD,ierror)
            
         end if
         
         !Assemble the tropopause-height flow
         if(mype==0) then
            
            do nrec=1,npe-2! All processors sent except top and bottom (so npe-2 in total)                                                                                                                           
!            do nrec=1,npe/2 !middle processors only have sent stuff
               
               call mpi_recv(ur_tp_r,n1*n2*3,MPI_DOUBLE,MPI_ANY_SOURCE,tag_eta2,MPI_COMM_WORLD,status,ierror)
               
               do ix=1,n1
                  do iy=1,n2
                     if( ur_tp_r(ix,iy,1) /= 0.D0 )  then
                        if(ur_tp(ix,iy,1) /= 0.D0) then
                           write(*,*) "Problem in tropopause spec"
                           stop
                        end if
                        ur_tp(ix,iy,1) = ur_tp_r(ix,iy,1)
                        ur_tp(ix,iy,2) = ur_tp_r(ix,iy,2)
                        ur_tp(ix,iy,3) = ur_tp_r(ix,iy,3)
                     end if
                  end do
               end do
               
               
            end do
            
            !Now move to k-space:
            
            call fft_r2c(ur_tp,uk_tp,3)
            
            
            !Plot the tropopause-height kinetic energy spectrum
            spz=0.
            num_modes=0.
            

            do iky=1,ikty
               ky = kya(iky)
               do ikx=1,iktx
                  kx = kxa(ikx)
                  kh2  = kx*kx+ky*ky
                  kh   = sqrt(1.D0*kh2)

                  mode = ifix(real(kh*L1/twopi+0.5))

                  if (L(ikx,iky).eq.1) then

                     spz(mode)   = spz(mode) +     real( uk_tp(ikx,iky,1)*CONJG(uk_tp(ikx,iky,1)) )  
                     spz(mode)   = spz(mode) +     real( uk_tp(ikx,iky,2)*CONJG(uk_tp(ikx,iky,2)) )  
                     spz(mode)   = spz(mode) + Ar2*real( uk_tp(ikx,iky,3)*CONJG(uk_tp(ikx,iky,3)) )  
                     
                     num_modes(mode) = num_modes(mode) + 2

                  endif
               enddo
            enddo
          
            spz=0.5*spz

            do mode=0,ktx-1
               if (num_modes(mode).ne.0) then         !mode, kin energy, pot energy, num of modes                                                                                                                                           
                  write(unit_tspec,fmt=*) float(mode),spz(mode),num_modes(mode)    !There was a problem here.[r_1/r_2 multiplied a second time to spz...]                                                                           
               endif
            enddo
            write(unit_tspec,*) '           '
            call flush(unit_tspec)
            
            
         end if


         uk = umem
         vk = vmem
         wk = bmem
         
      end if

    end subroutine tropopause_meanders





  END MODULE diagnostics
