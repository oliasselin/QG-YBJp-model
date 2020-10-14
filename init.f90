MODULE init

USE parameters
USE mpi
USE diagnostics
USE special

IMPLICIT NONE

CONTAINS






    SUBROUTINE init_arrays

      integer :: dealiasing=1               !Dealias (1)  or not (0)                                                                                                        
      integer :: step1,step2                !for zath                                                                                                              

      double precision :: kh
     ! Initialize wavenumber arrays                                                                                                                              

      !Once created, the wavevector arrays will look like this (given that N is even and L1=L2=2pi)

      !kxa=[ 0 1 2 ... N/2]
      !kya=[ 0 1 2 ... N/2-1  0  -N/2+1 -N/2+2 ... -2 -1] 
  
      !One must notice the expected ky=-N/2 is put to 0 manually with the second if condition

      do  ikx = 1,iktx
         kxa(ikx) = ikx-1
      enddo
      do iky = 1,ikty
         jj = iky - 1
         if (iky.gt.kty)   jj = jj - 2*kty
         if (iky.eq.kty+1) jj = 0
         if (iky.gt.2*kty) jj = 0
         kya(iky) = jj
      enddo


      L  = 1

      if(dealiasing==1) then
     do ikx = 1,iktx
        kx = kxa(ikx)
        do iky = 1,ikty
           ky = kya(iky)
           kh = sqrt(real(kx*kx + ky*ky))


! Set L=0 where necessary:                                                                                                                                                                                                                  
           if (iky.eq.kty+1) then
               L(ikx,iky) = 0
!               write(*,*) "L", ikx, iky, "iky.eq.kty+1" !Basically when ky is manually set to 0 instead of -N/2 (repeated when n2=even) 

           else if (kx.lt.0) then
               L(ikx,iky) = 0
 !              write(*,*) "L", iky, iky,  "kx.lt.0"     !Deadbolt: iktx is such kx never goes under 0 (due to reality condition                                                                                                            

           else if (kx.eq.0 .and. ky.lt.0)   then
               L(ikx,iky) = 0
  !             write(*,*) "L", ikx, iky, "kx.eq.0 .and. ky.lt.0"  ! ???  WHY DON'T WE KEEP THESE NEGATIVE KY'S??????                                                                                                             

           else if (iky.gt.2*kty)  then
               L(ikx,iky) = 0
   !            write(*,*) "L", ikx, iky,  "iky.gt.2*kty" !Deadbolt again                                                                                                                                                                   


           else if ( (kh*L1/twopi).gt.ifix(float(n1)/3.+0.5)-0.5) then
               L(ikx,iky) = 0
    !           write(*,*) "L", ikx, iky,  "Dealiasing"                                                                                                                                                                                     
           else
     !          write(*,*) "L", ikx, iky,  "L=", L(ikx,iky)                                                                                                                                                                                 
           end if
        enddo
     enddo
     end if


! Initialize r-space dimensional positions x,y,z                                                                                                                                                                                
            
do ix=1,n1
  xa(ix)=L1*(ix-1)/n1
end do
do iy=1,n2
  ya(iy)=L2*(iy-1)/n2
end do
do iz=1,n3
  za(iz)=L3*iz/n3 !za(iz) = [dz 2dz 3dz... n3*dz=L]
end do




do izh2=1,n3h2

   if(mype==0) then
      iz=izh2-2
   elseif(mype==(npe-1)) then
      iz=izh2 + n3 - (n3h0+2)
   else
      iz=izh2 + mype*n3h0-2
   end if


   if(iz>=1 .and. iz<=n3) then
      zah2(izh2)=za(iz)
   else
      zah2(izh2)=-1.  !X-marked (won't be used)                                                                                                                                                                                          
   end if


end do
do izh1=1,n3h1

   if(mype==0) then
      iz=izh1-1
   elseif(mype==(npe-1)) then
      iz=izh1 + n3 - (n3h0+1)
   else
      iz=izh1 + mype*n3h0-1
   end if


   if(iz>=1 .and. iz<=n3) then
      zah1(izh1)=za(iz)
   else
      zah1(izh1)=-1.  !X-marked (won't be used)                                                                                                                                                                                             
   end if

end do

do izh0=1,n3h0
   iz=izh0+mype*n3h0
   zah0(izh0)=za(iz)
end do



!Staggered versions  zasX(iz)=zaX(iz) - dz/2

do iz=1,n3
  zas(iz)=za(iz) - dz/2      !zas=[dz/2   3dz/2   ...  N-dz/2]
end do

do izh2=1,n3h2
   if(zah2(izh2)==-1.) then
      zash2(izh2)=-1.
   else
      zash2(izh2)=zah2(izh2) - dz/2
   end if
end do

do izh1=1,n3h1
   if(zah1(izh1)==-1.) then
      zash1(izh1)=-1.
   else
      zash1(izh1)=zah1(izh1) - dz/2
   end if
end do

do izh0=1,n3h0
   zash0(izh0)=zah0(izh0) - dz/2
end do



END SUBROUTINE init_arrays





subroutine init_base_state !This subroutine defines the staggered and unstaggered versions of the base-state buoyancy profile N^2.

  double precision :: N2_nd,N2_nds,N2_ndst,N2_ndut   !nondimensional N^2, (un)staggered and (un)transposed

  do izh2=1,n3h2
   
     z=zah2(izh2)     !Unstaggered fields
     zs=zash2(izh2)   !Staggered   fields
     
     if(stratification==constant_N) then
        N2_nd   = 1.D0
        N2_nds  = 1.D0
     else if(stratification==skewed_gaussian) then
        N2_nd   = N12_sg*exp(-((z -z0_sg)**2)/(sigma_sg**2))*(1.+erf( alpha_sg*(z -z0_sg)/(sigma_sg*sqrt(2.))))+N02_sg
        N2_nds  = N12_sg*exp(-((zs-z0_sg)**2)/(sigma_sg**2))*(1.+erf( alpha_sg*(zs-z0_sg)/(sigma_sg*sqrt(2.))))+N02_sg
     else
        write(*,*) "Undefined stratification profile. Aborting."
        stop
     end if

     r_1(izh2)     =  1.D0
     r_1s(izh2)    =  1.D0
     r_2(izh2)     =  N2_nd
     r_2s(izh2)    =  N2_nds

     a_ell_u(izh2) = Bu/N2_nd
     a_ell(izh2)   = Bu/N2_nds
     
     rho_u(izh2)   = 1.D0
     rho_s(izh2)   = 1.D0
     
     r_3u(izh2)    = 0.D0
     r_3(izh2)     = 0.D0

  end do

   !Special case: I need r_1 at z=0.
  r_1(izbot2-1) = 1.D0
  
  
  !Transposed fields (r_3 rho_0 a_ell b_ell a_helm b_helm)                                                                                                                              

  !Print base-state!
  if (mype==0) open (unit = 154673, file = "rucoeff.dat")
  if (mype==0) open (unit = 154674, file = "rscoeff.dat")
  
                          
  do iz=1,n3
     
     z  = za(iz)
     zs = zas(iz)
     
     
     if(stratification==constant_N) then
        N2_ndut = 1.D0
        N2_ndst = 1.D0
     else if(stratification==skewed_gaussian) then
        N2_ndut  = N12_sg*exp(-((z -z0_sg)**2)/(sigma_sg**2))*(1.+erf( alpha_sg*(z -z0_sg)/(sigma_sg*sqrt(2.))))+N02_sg
        N2_ndst  = N12_sg*exp(-((zs-z0_sg)**2)/(sigma_sg**2))*(1.+erf( alpha_sg*(zs-z0_sg)/(sigma_sg*sqrt(2.))))+N02_sg
     else
        write(*,*) "Undefined stratification profile. Aborting."
        stop
     end if

     r_1ut(iz)    =  1.D0
     r_1st(iz)    =  1.D0
     
     r_2ut(iz)    =  N2_ndut
     r_2st(iz)    =  N2_ndst
     
     rho_ut(iz)= 1.D0
     rho_st(iz)= 1.D0
     
     r_3ut(iz)  = 0.D0
     r_3t(iz)  = 0.D0
     
     a_ell_ut(iz)= Bu/N2_ndut
     a_ell_t(iz) = Bu/N2_ndst

     write(154673,"(E12.5,E12.5,E12.5,E12.5,E12.5)") z ,r_1ut(iz),r_2ut(iz),r_3ut(iz),a_ell_ut(iz)
     write(154674,"(E12.5,E12.5,E12.5,E12.5,E12.5)") zs,r_1st(iz),r_2st(iz),r_3t(iz), a_ell_t(iz)

  end do
  
  a_helm = 1./Ar2
  b_helm = 0.

end subroutine init_base_state
























 SUBROUTINE init_psi_generic(uk,vk,wk,bk,psik,psir)

    double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk
    double complex, dimension(iktx,ikty,n3h1) :: psik     
    double precision, dimension(n1d,n2d,n3h1) :: psir

    !Unstaggered version of psi to get staggered b
    double complex, dimension(iktx,ikty,n3h1) :: psik_us     
    double precision, dimension(n1d,n2d,n3h1) :: psir_us

    equivalence(psir_us,psik_us)

    real, parameter :: amp_width=2.

    real :: phi(4*initial_k+1,4*initial_k+1)   


    double precision :: amplitude
    real :: phase,norm1,norm2
    double precision :: kh

    double precision :: kz,kk 
    
    
    !Set random amplitude to be broadcasted 
    if(mype==0) then

       CALL RANDOM_SEED (PUT=seed)

       do ikx=1,4*initial_k+1
          do iky=1,4*initial_k+1 

                call random_number(phase)
                phi(ikx,iky)=twopi*phase       !Random number between 0 and twopi

          enddo
       enddo

       

    end if 

    !Broadcast phi to everybody                                                                                                                                                                                                             
    call mpi_bcast(phi,(4*initial_k+1)*(4*initial_k+1),MPI_REAL,0,MPI_COMM_WORLD,ierror)

    

    !Now set psi!
    !-----------|
    
    psir   =0.D0
    psir_us=0.D0

    if(linear_vert_structure==1) then

       do ikx = -2*initial_k,2*initial_k
          do iky = -2*initial_k,2*initial_k
             
             kh2 = ikx*ikx + iky*iky
             
             if(kh2 >0) then 
                
                kk= sqrt(1.D0*kh2)                                                                                                                                
                amplitude = exp( - (kk - 1.D0*initial_k)*(kk - 1.D0*initial_k)/(2.D0*amp_width))         !Gaussian decaying away from k=initial_k  
                
                !R-space loop to set psi!                                                                                                                            
                do izh1=1,n3h1                                                                                                                                                   
                   do ix=1,n1d                                                                                                                                                   
                      do iy=1,n2d                                                                                                                                                
                         if(ix<=n1) then                                                                                                                                         
                            
                            psir(ix,iy,izh1) =  psir(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)  + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
                                                        
                         end if
                      end do
                   end do
                end do
             end if
             
          end do
       end do
       
       !Both psi have the same horizontal structure
       psir_us = psir 
       
       !Multiply by the vertical enveloppe!
       !----------------------------------!
       
          do izh1=1,n3h1
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      
                      psir(ix,iy,izh1) = (zash1(izh1)-z0 ) *    psir(ix,iy,izh1) 

                      !Correct for the bottom point (exactly z=0). Only exception in the code where I use X-marked regions.
                      if(mype==0 .and. izh1==1) then
                         psir_us(ix,iy,izh1) =  -z0  * psir_us(ix,iy,izh1) 
                      else
                         psir_us(ix,iy,izh1) = (zah1(izh1) -z0 ) * psir_us(ix,iy,izh1)                          
                      end if

                   end if
                end do
             end do
          end do






    elseif(linear_vert_structure==2) then 
       
       
       do ikx = -2*initial_k,2*initial_k
          do iky = -2*initial_k,2*initial_k
             
             kh2 = ikx*ikx + iky*iky
             kh  = sqrt(1.D0*kh2)
             
             kz  = 1.!kh/sqrt(Bu)     !In non-dim form, kz ~ NH/fL kh, Bu = (fL/NH)^2
             
             
             if(kh2 >0) then 
                
                kk= sqrt(1.D0*kh2)                                                                                                                                
                amplitude = exp( - (kk - 1.D0*initial_k)*(kk - 1.D0*initial_k)/(2.D0*amp_width))         !Gaussian decaying away from k=initial_k  
                
                !R-space loop to set psi!                                                                                                                            
                do izh1=1,n3h1                                                                                                                                                   
                   do ix=1,n1d                                                                                                                                                   
                      do iy=1,n2d                                                                                                                                                
                         if(ix<=n1) then                                                                                                                                         
                            
           psir(ix,iy,izh1) =  psir(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)  + 1.D0*kz*zash1(izh1) + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
                            
                      !Correct for the bottom point (exactly z=0). Only exception in the code where I use X-marked regions.                                                       
                      if(mype==0 .and. izh1==1) then
           psir_us(ix,iy,izh1) =  psir_us(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)                       + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
                      else
           psir_us(ix,iy,izh1) =  psir_us(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)  + 1.D0*kz*zah1(izh1) + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
                      end if    

                         end if
                      end do
                   end do
                end do
             end if
             
          end do
       end do
       
 

       !Multiply by the vertical enveloppe!
       !----------------------------------!
       
       if(enveloppe==1) then
          do izh1=1,n3h1
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      
                      psir(ix,iy,izh1) = 0.5*( tanh((zash1(izh1)-z_env)/sig_env) -  tanh((zash1(izh1)-(L3-z_env))/sig_env) ) * psir(ix,iy,izh1) 

                      if(mype==0 .and. izh1==1) then
                         psir_us(ix,iy,izh1) = 0.5*( tanh(-z_env/sig_env) -  tanh(-(L3-z_env)/sig_env) ) * psir_us(ix,iy,izh1) 
                      else
                         psir_us(ix,iy,izh1) = 0.5*( tanh((zah1(izh1)-z_env)/sig_env) -  tanh((zah1(izh1)-(L3-z_env))/sig_env) ) * psir_us(ix,iy,izh1) 
                      end if


                      
                   end if
                end do
             end do
          end do
       end if

    else !If not a linear psi
       
       
       do ikx = -2*initial_k,2*initial_k
          do iky = -2*initial_k,2*initial_k
             
             kh2 = ikx*ikx + iky*iky
             kh  = sqrt(1.D0*kh2)
             
             kz  = kh/sqrt(Bu)     !In non-dim form, kz ~ NH/fL kh, Bu = (fL/NH)^2
             
             
             if(kz>0. .and. kh2 >0) then 
                
                kk= sqrt(1.D0*kh2 + 1.D0*kz*kz)                                                                                                                                
                amplitude = exp( - (kk - 1.D0*initial_k)*(kk - 1.D0*initial_k)/(2.D0*amp_width))         !Gaussian decaying away from k=initial_k  
                
                !R-space loop to set psi!                                                                                                                            
                do izh1=1,n3h1                                                                                                                                                   
                   do ix=1,n1d                                                                                                                                                   
                      do iy=1,n2d                                                                                                                                                
                         if(ix<=n1) then                                                                                                                                         
                            
            psir(ix,iy,izh1) =  psir(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)  + 1.D0*kz*zash1(izh1) + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
 
            if(mype==0 .and. izh1==1) then
               psir_us(ix,iy,izh1) =  psir_us(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)                       + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
            else
               psir_us(ix,iy,izh1) =  psir_us(ix,iy,izh1) + amplitude*cos(1.D0*ikx*xa(ix)  + 1.D0*iky*ya(iy)  + 1.D0*kz*zah1(izh1) + phi(ikx+2*initial_k+1,iky+2*initial_k+1)) 
            end if
                           
                         end if
                      end do
                   end do
                end do
             end if
             
          end do
       end do
       
 

       !Multiply by the vertical enveloppe!
       !----------------------------------!
       

       if(enveloppe==1) then
          do izh1=1,n3h1
             do ix=1,n1d
                do iy=1,n2d
                   if(ix<=n1) then
                      
                      psir(ix,iy,izh1) = 0.5*( tanh((zash1(izh1)-z_env)/sig_env) -  tanh((zash1(izh1)-(L3-z_env))/sig_env) ) * psir(ix,iy,izh1) 
                      

                      if(mype==0 .and. izh1==1) then
                         psir_us(ix,iy,izh1) = 0.5*( tanh((-z_env)/sig_env) -  tanh(-(L3-z_env)/sig_env) ) * psir_us(ix,iy,izh1) 
                      else
                         psir_us(ix,iy,izh1) = 0.5*( tanh((zah1(izh1)-z_env)/sig_env) -  tanh((zah1(izh1)-(L3-z_env))/sig_env) ) * psir_us(ix,iy,izh1) 
                      end if


                   end if
                end do
             end do
          end do
       end if

    end if


    !Now move psi to k-space!
    !-----------------------!
    call fft_r2c(psir,psik,n3h1)
    call fft_r2c(psir_us,psik_us,n3h1)


    !Now compute u,v,w and t!
    !-----------------------!

     do izh1=1,n3h1
        izh2=izh1+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              if (L(ikx,iky).eq.1) then
                 uk(ikx,iky,izh2) =  - i*ky*psik(ikx,iky,izh1)
                 vk(ikx,iky,izh2) =    i*kx*psik(ikx,iky,izh1)
                 wk(ikx,iky,izh2) =    (0.D0,0.D0)
                 bk(ikx,iky,izh2) =    ( psik_us(ikx,iky,izh1) - psik_us(ikx,iky,izh1-1) )/(r_1s(izh2)*dz)
              else
                 uk(ikx,iky,izh2)   = (0.D0,0.D0)
                 vk(ikx,iky,izh2)   = (0.D0,0.D0)
                 wk(ikx,iky,izh2)   = (0.D0,0.D0)
                 bk(ikx,iky,izh2)   = (0.D0,0.D0)
                 psik(ikx,iky,izh1) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo




    !Normalize energy while keeping the ratio between u,v and psi: normalize by total energy rather than K and P separately

     if(normalize==1 .and. norm_energy==1) then
        call energy_linear(uk,vk,wk,bk,norm1,norm2)
        call mpi_bcast(norm1,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast(norm2,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        
        uk=uk*sqrt((k_init+p_init)/(norm1+norm2))
        vk=vk*sqrt((k_init+p_init)/(norm1+norm2))
    psik=psik*sqrt((k_init+p_init)/(norm1+norm2))
        wk=wk*sqrt((k_init+p_init)/(norm1+norm2))
        bk=bk*sqrt((k_init+p_init)/(norm1+norm2))
     end if

   END SUBROUTINE init_psi_generic


 SUBROUTINE init_q(qk,psik)

   !This subroutine simply computes the RHS of the elliptic equation for psi, knowing psi.                                                                                                                          
   !q = - kh2 psi + 1/rho d/dz (rho a_ell psi_z) with d psi./dz = 0 at bounds.                                                                                                                                                    
   !The output is q without its halo.                                                                                                                                                                                              

    double complex, dimension(iktx,ikty,n3h1) :: qk,psik

    qk = (0.D0,0.D0)
    do izh0=1,n3h0
       izh1=izh0+1
       izh2=izh0+2
       do iky=1,ikty
          ky = kya(iky)
          do ikx=1,iktx
             kx = kxa(ikx)
             kh2= kx*kx+ky*ky

             if (L(ikx,iky).eq.1) then
                qk(ikx,iky,izh1) =  (rho_u(izh2)*a_ell_u(izh2)/rho_s(izh2))*psik(ikx,iky,izh1+1)  - ( (rho_u(izh2)*a_ell_u(izh2) + rho_u(izh2-1)*a_ell_u(izh2-1) )/rho_s(izh2) + kh2*dz*dz)*psik(ikx,iky,izh1)  +  (rho_u(izh2-1)*a_ell_u(izh2-1)/rho_s(izh2))*psik(ikx,iky,izh1-1)

             endif

          end do
       end do
    end do

    !Boundary corrections                                                                                                                                                                                                    
    !At the top, d psi/dz = 0, so psik(top+1) = psi(top) and psi(bot-1) = psi(bot)                                                                                                                                              

    if(mype==0) then  !At the bottom                                                                                                                                                                                  

       do iky=1,ikty
          ky = kya(iky)
          do ikx=1,iktx
             kx = kxa(ikx)
             kh2= kx*kx+ky*ky

             if (L(ikx,iky).eq.1) then
                qk(ikx,iky,izbot1) =  (rho_u(izbot2)*a_ell_u(izbot2)/rho_s(izbot2))*psik(ikx,iky,izbot1+1)  - ( rho_u(izbot2)*a_ell_u(izbot2)/rho_s(izbot2) + kh2*dz*dz)*psik(ikx,iky,izbot1)
             endif

          end do
       end do
    else if(mype==(npe-1)) then !At the top   
       do iky=1,ikty
          ky = kya(iky)
          do ikx=1,iktx
             kx = kxa(ikx)
             kh2= kx*kx+ky*ky

             if (L(ikx,iky).eq.1) then
                qk(ikx,iky,iztop1) =   - ( rho_u(iztop2-1)*a_ell_u(iztop2-1)/rho_s(iztop2) + kh2*dz*dz)*psik(ikx,iky,iztop1)  +  (rho_u(iztop2-1)*a_ell_u(iztop2-1)/rho_s(iztop2))*psik(ikx,iky,iztop1-1)
             endif

          end do
       end do
    end if

    qk=qk/(dz*dz)

  end SUBROUTINE init_q







  !*************************************************************************!
  !*************************************************************************!
  !************   Test-type real-space initialization routine   ************!
  !*************************************************************************!
  !*************************************************************************!



    SUBROUTINE init_prog_vars(ur,vr,wr,br)  

      double precision, dimension(n1d,n2d,n3h2), intent(out) :: ur,vr,wr,br

      ! Initialize r-space velocity field  ur,vr,wr                                                                                                   
                                                                                                                   

      do ix=1,n1d
         if(ix<=n1) x=xa(ix)
         do iy=1,n2d
            y=ya(iy)
            do izh2=1,n3h2
               z=zah2(izh2)     !Use z for w,t
               zs=zash2(izh2)   !Use zs for u,v

               if(ix<=n1  .and. z>=0) then
                  ur(ix,iy,izh2)=                   cos(y)   !-sin(x)*cos(zs)        !USE zs
                  vr(ix,iy,izh2)=(2./3.)*cos(2*zs)* cos(3*y) !0.                     !USE zs  
                  wr(ix,iy,izh2)=        sin(2*z) * sin(3*y) ! cos(x)*sin(z)
                  br(ix,iy,izh2)=0.
                  
               else
                  ur(ix,iy,izh2)=0.
                  vr(ix,iy,izh2)=0.
                  wr(ix,iy,izh2)=0.
                  br(ix,iy,izh2)=0.
               end if
               
            end do
         end do
      end do
      
    END SUBROUTINE init_prog_vars







    SUBROUTINE generate_fields(f1,nz1,f2,nz2,f3,nz3)

      integer, intent(in) :: nz1,nz2,nz3                  !z-dimension size of the fields
      double precision, dimension(n1d,n2d,nz1), intent(out) :: f1
      double precision, dimension(n1d,n2d,nz2), intent(out) :: f2
      double precision, dimension(n1d,n2d,nz3), intent(out) :: f3
      
      double precision :: z1,z2,z3   !Corresponding value of z for each field (depends on halo)
      integer :: hlev1,hlev2,hlev3   !What's your halo?
      integer :: iz1,iz2,iz3


      hlev1=(nz1-n3h0)/2
      hlev2=(nz2-n3h0)/2
      hlev3=(nz3-n3h0)/2
      
do ix=1,n1d
   x=xa(ix)
   do iy=1,n2d
      y=ya(iy)
      do izh2=1,n3h2
         izh1=izh2-1
         izh0=izh2-2


         if(hlev1==2) then 
            z1=zah2(izh2)
            iz1=izh2
         elseif(hlev1==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z1=zah1(izh1)
            iz1=izh1
         elseif(hlev1==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z1=zah0(izh0)
            iz1=izh0
         else 
            z1=-1.   ! Don't enter the loop
         end if

         if(hlev2==2) then
            z2=zah2(izh2)
            iz2=izh2
         elseif(hlev2==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z2=zah1(izh1)
            iz2=izh1
         elseif(hlev2==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z2=zah0(izh0)
            iz2=izh0
         else
                z2=-1.   ! Don't enter the loop   
         end if


         if(hlev3==2) then 
            z3=zah2(izh2)
            iz3=izh2
         elseif(hlev3==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z3=zah1(izh1)
            iz3=izh1
         elseif(hlev3==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z3=zah0(izh0)
            iz3=izh0
         else
                z3=-1.   ! Don't enter the loop   
         end if

      if(ix<=n1) then
         if(z1>=0) f1(ix,iy,iz1)=0.!sin(x)*cos(x)! cos(x)*cos(y)
         if(z2>=0) f2(ix,iy,iz2)=0.!-sin(x)*sin(y)
         if(z3>=0) f3(ix,iy,iz3)=2*sin(2*z3)*cos(2*z3)!sin(z3)*cos(z3)
      else
         if(z1>=0) f1(ix,iy,iz1)=0.
         if(z2>=0) f2(ix,iy,iz2)=0.
         if(z3>=0) f3(ix,iy,iz3)=0.
      end if

      !In this case, the X-marked arrays (on the top and bottom mype will not be initialized at all. Shouldn't cause any problem since we never invoke them.

      end do
   end do
end do




  END SUBROUTINE generate_fields


    SUBROUTINE generate_fields_stag(f1s,nz1,f2s,nz2,f3s,nz3) !generate_fields for staggered fields (like u,v)

      integer, intent(in) :: nz1,nz2,nz3                  !z-dimension size of the fields
      double precision, dimension(n1d,n2d,nz1), intent(out) :: f1s
      double precision, dimension(n1d,n2d,nz2), intent(out) :: f2s
      double precision, dimension(n1d,n2d,nz3), intent(out) :: f3s
      
      double precision :: z1,z2,z3   !Corresponding value of z for each field (depends on halo)
      integer :: hlev1,hlev2,hlev3   !What's your halo?
      integer :: iz1,iz2,iz3


      hlev1=(nz1-n3h0)/2
      hlev2=(nz2-n3h0)/2
      hlev3=(nz3-n3h0)/2
      
do ix=1,n1d
   x=xa(ix)
   do iy=1,n2d
      y=ya(iy)
      do izh2=1,n3h2
         izh1=izh2-1
         izh0=izh2-2


         if(hlev1==2) then 
            z1=zash2(izh2)
            iz1=izh2
         elseif(hlev1==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z1=zash1(izh1)
            iz1=izh1
         elseif(hlev1==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z1=zash0(izh0)
            iz1=izh0
         else 
            z1=-1.   ! Don't enter the loop
         end if

         if(hlev2==2) then
            z2=zash2(izh2)
            iz2=izh2
         elseif(hlev2==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z2=zash1(izh1)
            iz2=izh1
         elseif(hlev2==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z2=zash0(izh0)
            iz2=izh0
         else
            z2=-1.   ! Don't enter the loop   
         end if


         if(hlev3==2) then 
            z3=zash2(izh2)
            iz3=izh2
         elseif(hlev3==1 .and. izh1>=1 .and. izh1<=n3h1) then
            z3=zash1(izh1)
            iz3=izh1
         elseif(hlev3==0 .and. izh0>=1 .and. izh0<=n3h0) then
            z3=zash0(izh0)
            iz3=izh0
         else
                z3=-1.   ! Don't enter the loop   
         end if

      if(ix<=n1) then
         
         !Set fields here. for fX, one must use the vertical index zX. ex. f2s(ix,iy,iz2)=cos(mmm*z2)  
         if(z1>=0) f1s(ix,iy,iz1)=sin(x)*sin(y)
         if(z2>=0) f2s(ix,iy,iz2)=0.
         if(z3>=0) f3s(ix,iy,iz3)=cos(8*z3/2.)
         
      else
         if(z1>=0) f1s(ix,iy,iz1)=0.
         if(z2>=0) f2s(ix,iy,iz2)=0.
         if(z3>=0) f3s(ix,iy,iz3)=0.
      end if

      !In this case, the X-marked arrays (on the top and bottom mype will not be initialized at all. Shouldn't cause any problem since we never invoke them.

      end do
   end do
end do




END SUBROUTINE generate_fields_stag









END MODULE init
