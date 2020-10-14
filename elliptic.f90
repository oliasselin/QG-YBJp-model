MODULE elliptic

USE parameters
USE fft
USE mpi

IMPLICIT NONE

CONTAINS
  

    SUBROUTINE helmholtzdouble(p,ft,b_bot,b_top)  

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: p    !Pressure in usual z-parallelization                                                                       
      double complex, dimension(iktx,n3d1, iktyp)             :: pt   !Transposed (ky-parallelization) pressure                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: ft   !Transposed (ky-parallelization) right-hand side               
      double complex, dimension(iktx,ikty), intent(in) :: b_bot,b_top   !Top and bottom buoyancy terms             

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value   
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                      

                                                                                                      
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3
                  
                  br(iz)=dz*dz*DBLE(  ft(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( ft(ikx,iz,ikyp) )  

               end do

               !--- For nonzero w* at the top/bot, we must add an extra to the RHS ---!
               !----------------------------------------------------------------------!

               iz=1
               
               br(iz) = br(iz) + (a_helm(iz) - 0.5*b_helm(iz)*dz)* DBLE(b_bot(ikx,iky))*dz
               bi(iz) = bi(iz) + (a_helm(iz) - 0.5*b_helm(iz)*dz)*DIMAG(b_bot(ikx,iky))*dz

               iz=n3
               
               br(iz) = br(iz) - (a_helm(iz) + 0.5*b_helm(iz)*dz)* DBLE(b_top(ikx,iky))*dz
               bi(iz) = bi(iz) - (a_helm(iz) + 0.5*b_helm(iz)*dz)*DIMAG(b_top(ikx,iky))*dz

               !----------------------------------------------------------------------!

 

               do iz=1,n3
                  if(iz==1) then
                     d(iz) = -a_helm(iz)-0.5*b_helm(iz)*dz-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  elseif(iz==n3) then
                     d(iz) = -a_helm(iz)+0.5*b_helm(iz)*dz-kh2*dz*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  else !1<iz<n3                                                                                                                                                                                                             
                     d(iz) = -2*a_helm(iz)-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  end if
               end do

               
               call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info
               
               do iz=1,n3
                  if(iz==1) then
                     d(iz) = -a_helm(iz)-0.5*b_helm(iz)*dz-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  elseif(iz==n3) then
                     d(iz) = -a_helm(iz)+0.5*b_helm(iz)*dz-kh2*dz*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  else !1<iz<n3                                                                                                                                                                                                             
                     d(iz) = -2*a_helm(iz)-kh2*dz*dz
                    du(iz) =  a_helm(iz)+0.5*b_helm(iz)*dz
                  dl(iz-1) =  a_helm(iz)-0.5*b_helm(iz)*dz
                  end if
               end do
               
               call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info


               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1) then
                     pt(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     pt(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
      
            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   pt(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO


      

      



      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(pt,iktx,n3d1,iktyp,p,ikty,n3h1)


    END SUBROUTINE HELMHOLTZDOUBLE






    SUBROUTINE psi_solver(psik,qt)

      ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi. 

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: psik   !Pressure in usual z-parallelization                                                                       
      double complex, dimension(iktx,n3d1, iktyp)             :: psit   !Transposed (ky-parallelization) pressure                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: qt       !Transposed (ky-parallelization) right-hand side               

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value   
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                      

                                                                                                                                                                          
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3
                  
                  br(iz)=dz*dz*DBLE(  qt(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( qt(ikx,iz,ikyp) )  

               end do

               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info
               
               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do

               call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info
               
               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1) then
                     psit(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     psit(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
               

            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   psit(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO

      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(psit,iktx,n3d1,iktyp,psik,ikty,n3h1)


    END SUBROUTINE PSI_SOLVER
    



    SUBROUTINE A_solver_ybj_plus_original(ak,bt)

      ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi. 

      double complex, dimension(iktx,ikty,n3h0), intent(out) :: ak   !A in usual z-parallelization                                                                       
      double complex, dimension(iktx,n3, iktyp)              :: at   !Transposed (ky-parallelization) A                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in)  :: bt   !Transposed (ky-parallelization) right-hand side               

      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value   
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                      

                                                                                                                                                                          
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3
                  
                  br(iz)=dz*dz*Bu*DBLE(  bt(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*Bu*DIMAG( bt(ikx,iz,ikyp) )  

               end do

               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz/4. )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info
               
               do iz=1,n3
                  if(iz==1) then                     
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz/4. )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do

               call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info
               
               !*** Put the solution in pt ***!    
               
               DO iz=1,n3
                  at(ikx,iz,ikyp)=br(iz)+i*bi(iz)
               END DO

            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO iz=1,n3
                   at(ikx,iz,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO

      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(at,iktx,n3,iktyp,ak,ikty,n3h0)


    END SUBROUTINE A_SOLVER_YBJ_PLUS_ORIGINAL
    


    SUBROUTINE A_solver_ybj_plus(ak,bt,ck)

       ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi.                                              

      double complex, dimension(iktx,ikty,n3h0), intent(out) :: ak   !A in usual z-parallelization                                                                                      
      double complex, dimension(iktx,n3, iktyp)              :: at   !Transposed (ky-parallelization) A                                                                                
      double complex, dimension(iktx,n3, iktyp), intent(in)  :: bt   !Transposed (ky-parallelization) right-hand side                                                                   
      double complex, dimension(iktx,ikty,n3h0), intent(out) :: ck   !C=A_z in usual z-parallelization                                                                                  
      double complex, dimension(iktx,n3, iktyp)              :: ct   !Transposed (ky-parallelization) C                                                                             


      integer :: info                                      ! Returns error for sgtsv                                                                                                    

      double precision :: d(n3),dl(n3-1),du(n3-1)          !diagonal value                                                                                                              
      double precision :: br(n3), bi(n3)                 !real and imaginary parts of rhs                                                                                              


       DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

             kh2=kx*kx + ky*ky

             if(kh2/=0 .and. L(ikx,iky)==1 ) then

                do iz=1,n3

                   br(iz)=dz*dz*Bu*DBLE(  bt(ikx,iz,ikyp) )
                  bi(iz)=dz*dz*Bu*DIMAG( bt(ikx,iz,ikyp) )

                end do

                do iz=1,n3
                  if(iz==1) then
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz/4. )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3                                                                                                                                                          
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do

                call DGTSV( n3, 1, dl, d, du, br, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble", info

                do iz=1,n3
                  if(iz==1) then
                     d(iz) = -( rho_ut(iz)*a_ell_ut(iz)/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  elseif(iz==n3) then
                     d(iz) = -( rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz) + kh2*dz*dz/4. )
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3                                                                                                                                                          
                     d(iz) = -( (rho_ut(iz)*a_ell_ut(iz) + rho_ut(iz-1)*a_ell_ut(iz-1))/rho_st(iz) + kh2*dz*dz/4. )
                    du(iz) = rho_ut(iz)*a_ell_ut(iz)/rho_st(iz)
                  dl(iz-1) = rho_ut(iz-1)*a_ell_ut(iz-1)/rho_st(iz)
                  end if
               end do

                call DGTSV( n3, 1, dl, d, du, bi, n3, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info

                !*** Put the solution in pt ***!                                                                                                                 

                DO iz=1,n3
                  at(ikx,iz,ikyp)=br(iz)+i*bi(iz)
               END DO

             else

                !set kh=0 modes to zero...                                                                                                                                  

                 DO iz=1,n3
                   at(ikx,iz,ikyp)=(0.,0.)
                END DO


              end if

          END DO
      END DO

      !Compute C=A_z while transposed!                                                                                                                                        
      !******************************!                                                                                                                                          

      DO ikx=1,iktx
         DO ikyp=1,iktyp

            do iz=1,n3-1

               ct(ikx,iz,ikyp) = ( at(ikx,iz+1,ikyp) - at(ikx,iz,ikyp) )/dz

            end do

            !Top: C = A_z = 0                                                                                                                                                   
            ct(ikx,n3,ikyp) = (0.D0,0.D0)

         end DO
      end DO

      !Transpose the result back to z-parallelized space!                                                                                                                        
      call mpitranspose(ct,iktx,n3,iktyp,ck,ikty,n3h0)


       !*********** Transposition to z-parallelized p ***************!                                                                                     
      call mpitranspose(at,iktx,n3,iktyp,ak,ikty,n3h0)


    END SUBROUTINE A_SOLVER_YBJ_PLUS


    SUBROUTINE omega_equation(wak,qt)

      ! This subroutines solves the elliptic equation a_ell(z) d^2 psi/dz^2 + b_ell(z) d psi/dz - kh2 psi = q for the streamfunction psi. 

      double complex, dimension(iktx,ikty,n3h1),  intent(out) :: wak   !vertical velocity in usual z-parallelization                                   
      double complex, dimension(iktx,n3d1, iktyp)             :: wat   !Transposed (ky-parallelization) wak                                                   
      double complex, dimension(iktx,n3, iktyp), intent(in) :: qt      !Transposed (ky-parallelization) right-hand side               

      integer :: step1,step2                               !Intermediate steps in calculation of indices
      integer :: info                                      ! Returns error for sgtsv                                                                                           

      double precision :: d(n3-1),dl(n3-2),du(n3-2)          !diagonal value   
      double precision :: br(n3-1), bi(n3-1)                 !real and imaginary parts of rhs                                                                      

                                                                                                                                                                          
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               
               do iz=1,n3-1
                  
                  br(iz)=dz*dz*DBLE(  qt(ikx,iz,ikyp) )  
                  bi(iz)=dz*dz*DIMAG( qt(ikx,iz,ikyp) )  

               end do

               do iz=1,n3-1
                  if(iz==1) then                     
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
!-(rho_ut(iz)/(rho_st(iz+1)+rho_st(iz)) + dz*dz*kh2/a_ell_ut(iz))      WRONG
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  elseif(iz==(n3-1)) then
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3-1
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3-1, 1, dl, d, du, br, n3-1, info )
               if(info/=0) write(*,*) "problem in helmdouble", info

               do iz=1,n3-1
                  if(iz==1) then                     
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  elseif(iz==(n3-1)) then
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  else !1<iz<n3-1
                     d(iz) = -(rho_ut(iz)*(rho_st(iz+1)+rho_st(iz))/(rho_st(iz+1)*rho_st(iz))  + dz*dz*kh2/a_ell_ut(iz))
                    du(iz) = rho_ut(iz+1)/rho_st(iz+1)
                  dl(iz-1) = rho_ut(iz-1)/rho_st(iz)
                  end if
               end do
               
               call DGTSV( n3-1, 1, dl, d, du, bi, n3-1, info )
               if(info/=0) write(*,*) "problem in helmdouble (imag)", info
               
               
               !*** Put the solution in pt ***!    
               
               do izth=1,n3d1
                  
                  step1=(izth-1)/n3h1   !mype                                                       
                  step2=izth-n3h1*step1 !position in mype                                           
                  iz=step1*n3h0+step2-1 !-1 factor since h1 fields                                  
                  
                  if(izth/=1 .and. izth/=n3d1 .and. iz/=n3) then    !Difference with psi_solver, if iz=n3, w=0. I didn't even solve the equation for that point...
                     wat(ikx,izth,ikyp)=br(iz)+i*bi(iz)
                  else
                     wat(ikx,izth,ikyp)=(0.,0.)
                  end if
                  
               end do
               

            else
            
               !set kh=0 modes to zero...                                                                                                                                    

                DO izth=1,n3d1
                   wat(ikx,izth,ikyp)=(0.,0.)
                END DO

             
             end if

         END DO
      END DO

      !*********** Transposition to z-parallelized p ***************!                                                                                                               

      call mpitranspose(wat,iktx,n3d1,iktyp,wak,ikty,n3h1)


    END SUBROUTINE OMEGA_EQUATION

END MODULE elliptic
