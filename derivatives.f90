MODULE derivatives

  USE parameters
  USE fft
  USE mpi

  IMPLICIT NONE

!  integer :: kx,ky,kh2
!  integer :: ikx,iky
!  integer :: izh0,izh1,izh2,izh3

  CONTAINS


    SUBROUTINE  convol(nuk,nvk,nwk,nbk,nur,nvr,nwr,nbr,uk,vk,wk,bk,ur,vr,wr,br)
      ! this subroutine computes 1/rho d/dxj( rho uj ui) in the divergence form on a staggered grid.
                                                                    
      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk   
      double complex, dimension(iktx,ikty,n3h1) :: nuk,nvk,nwk,nbk

      double precision,    dimension(n1d,n2d,n3h2) :: ur,vr,wr,br
      double precision,    dimension(n1d,n2d,n3h1) :: nur,nvr,nwr,nbr

      double complex, dimension(iktx,ikty,n3h2) :: umem,vmem,wmem,bmem


      !Suboptimal but simple.
      double complex,   dimension(iktx,ikty,n3h1) :: Auk,Buk,Cuk,Avk,Bvk,Cvk,Awk,Bwk,Cwk,Abk,Bbk,Cbk
      double precision, dimension(n1d,n2d,n3h1)   :: Aur,Bur,Cur,Avr,Bvr,Cvr,Awr,Bwr,Cwr,Abr,Bbr,Cbr

      equivalence(Aur,Auk)
      equivalence(Avr,Avk)
      equivalence(Awr,Awk)
      equivalence(Abr,Abk)

      equivalence(Bur,Buk)
      equivalence(Bvr,Bvk)
      equivalence(Bwr,Bwk)
      equivalence(Bbr,Bbk)

      equivalence(Cur,Cuk)
      equivalence(Cvr,Cvk)
      equivalence(Cwr,Cwk)
      equivalence(Cbr,Cbk)


      !DO A COPY OF VELOCITY, FASTER THAN TAKE A TRANSFORM BACK                                                                                     

      umem=uk
      vmem=vk
      wmem=wk
      bmem=bk

      call fft_c2r(uk,ur,n3h2)
      call fft_c2r(vk,vr,n3h2)
      call fft_c2r(wk,wr,n3h2)
      call fft_c2r(bk,br,n3h2)

      !Assign values to A, B and C for each of the three components of the nlt
      !nuk = ikx Auk + iky Buk + 1/2dz Cuk  (Same for nvk and nwk but u->v or ->w)
      !With A,B and C defined as:

      !Aur_n =  u_n^2
      !Bur_n =  v_n*u_n
      !Cur_n =  w_n*u_n+1 + w_n*u_n - u_n*w_n-1 - u_n-1*w_n-1

      !(same as for u)
      !Avr_n =  u_n*v_n
      !Bvr_n =  v_n^2
      !Cvr_n =  w_n*v_n+1 + w_n*v_n - w_n-1*v_n - w_n-1*v_n-1

      !Awr_n = 1/2 ( u_n+1*w_n + u_n*w_n)
      !Bwr_n = 1/2 ( v_n+1*w_n + v_n*w_n)
      !Cwr_n = 1/2 ( w_n+1^2 - w_n-1^2  ) + w_n * (w_n+1 - w_n-1)

      !Atr_n = 1/2 ( u_n+1*t_n + u_n*t_n)
      !Btr_n = 1/2 ( v_n+1*t_n + v_n*t_n)
      !Ctr_n = 1/2 [ w_n+1*t_n+1 + w_n+1*t_n + w_n*t_n+1 - w_n*t_n-1 - w_n-1*t_n - w_n-1*b_n-1 ]

      
      !DAVID'S CHEAT (NOT USED IN THE CURRENT VERSION): modify the point below the top and the one above the bottom for the unstaggered field w
      !David's cheat make the w-nonlinear term conservative in the momentum equation (triple integral of n_w = 0) but NOT in the energy equation 
      !(triple integral of w*n_w is NOT = 0). It is unlikely that this tweak makes energy conservation better, and it causes a O(dz) error at 
      !the modified grid points. 
      !See journal Feb 19th 2016 or code documentation. March 25h 2016: reverted to the original, non-David-cheated algorithm.


      do izh1=1,n3h1
         izh2=izh1+1
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   Aur(ix,iy,izh1) = ur(ix,iy,izh2)*ur(ix,iy,izh2)
                   Bur(ix,iy,izh1) = vr(ix,iy,izh2)*ur(ix,iy,izh2)
                   Cur(ix,iy,izh1) = rho_u(izh2)*wr(ix,iy,izh2)*ur(ix,iy,izh2+1) + rho_u(izh2)*wr(ix,iy,izh2)*ur(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*ur(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*ur(ix,iy,izh2-1) 

                   Avr(ix,iy,izh1) = ur(ix,iy,izh2)*vr(ix,iy,izh2)
                   Bvr(ix,iy,izh1) = vr(ix,iy,izh2)*vr(ix,iy,izh2)
                   Cvr(ix,iy,izh1) = rho_u(izh2)*wr(ix,iy,izh2)*vr(ix,iy,izh2+1) + rho_u(izh2)*wr(ix,iy,izh2)*vr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*vr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*vr(ix,iy,izh2-1) 

                   Awr(ix,iy,izh1) = 0.5 * ( ur(ix,iy,izh2+1)*wr(ix,iy,izh2)  +  ur(ix,iy,izh2)*wr(ix,iy,izh2) )    
                   Bwr(ix,iy,izh1) = 0.5 * ( vr(ix,iy,izh2+1)*wr(ix,iy,izh2)  +  vr(ix,iy,izh2)*wr(ix,iy,izh2) )    
                   Cwr(ix,iy,izh1) = 0.5 * ( rho_u(izh2+1)*wr(ix,iy,izh2+1)*wr(ix,iy,izh2+1) + rho_u(izh2+1)*wr(ix,iy,izh2+1)*wr(ix,iy,izh2) + rho_u(izh2)*wr(ix,iy,izh2)*wr(ix,iy,izh2+1) - rho_u(izh2)*wr(ix,iy,izh2)*wr(ix,iy,izh2-1) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*wr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*wr(ix,iy,izh2-1)  )


                   Abr(ix,iy,izh1) = ur(ix,iy,izh2)*br(ix,iy,izh2)
                   Bbr(ix,iy,izh1) = vr(ix,iy,izh2)*br(ix,iy,izh2)
                   Cbr(ix,iy,izh1) = rho_u(izh2)*wr(ix,iy,izh2)*br(ix,iy,izh2+1) + rho_u(izh2)*wr(ix,iy,izh2)*br(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*br(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*br(ix,iy,izh2-1) 

                   !Unstaggered buoyancy tk
!                   Atr(ix,iy,izh1) = 0.5 * ( ur(ix,iy,izh2+1)*tr(ix,iy,izh2)  +  ur(ix,iy,izh2)*tr(ix,iy,izh2) )    
!                   Btr(ix,iy,izh1) = 0.5 * ( vr(ix,iy,izh2+1)*tr(ix,iy,izh2)  +  vr(ix,iy,izh2)*tr(ix,iy,izh2) )    
!                   Ctr(ix,iy,izh1) = 0.5 * ( rho_u(izh2+1)*wr(ix,iy,izh2+1)*tr(ix,iy,izh2+1) + rho_u(izh2+1)*wr(ix,iy,izh2+1)*tr(ix,iy,izh2) + rho_u(izh2)*wr(ix,iy,izh2)*tr(ix,iy,izh2+1) - rho_u(izh2)*wr(ix,iy,izh2)*tr(ix,iy,izh2-1) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*tr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*tr(ix,iy,izh2-1)  )

                else

                   Aur(ix,iy,izh1) = 0.D0
                   Bur(ix,iy,izh1) = 0.D0
                   Cur(ix,iy,izh1) = 0.D0

                   Avr(ix,iy,izh1) = 0.D0
                   Bvr(ix,iy,izh1) = 0.D0
                   Cvr(ix,iy,izh1) = 0.D0

                   Awr(ix,iy,izh1) = 0.D0
                   Bwr(ix,iy,izh1) = 0.D0
                   Cwr(ix,iy,izh1) = 0.D0

                   Abr(ix,iy,izh1) = 0.D0
                   Bbr(ix,iy,izh1) = 0.D0
                   Cbr(ix,iy,izh1) = 0.D0

                end if

             end do
          end do
       end do

       !Boundary conditions! X-marked regions were assigned with insignificant values. We give the correct values here

       if(mype==0) then

         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   !Aur, Avr, Bur and Bur are fine already.

                   Cur(ix,iy,izbot1) = rho_u(izbot2)*wr(ix,iy,izbot2) * ( ur(ix,iy,izbot2+1) + ur(ix,iy,izbot2) ) 
                   Cvr(ix,iy,izbot1) = rho_u(izbot2)*wr(ix,iy,izbot2) * ( vr(ix,iy,izbot2+1) + vr(ix,iy,izbot2) ) 

                   Cwr(ix,iy,izbot1) = 0.5 * (rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2+1) + rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2) + rho_u(izbot2)*wr(ix,iy,izbot2)*wr(ix,iy,izbot2+1) )  !Original, non-Davidcheat
!                   Cwr(ix,iy,izbot1) = 0.5 * (rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2+1) + rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2) + rho_u(izbot2)*wr(ix,iy,izbot2)*wr(ix,iy,izbot2+1) +  rho_u(izbot2)*wr(ix,iy,izbot2)*wr(ix,iy,izbot2) )    !David's cheat

                   Cbr(ix,iy,izbot1) = rho_u(izbot2)*wr(ix,iy,izbot2) * ( br(ix,iy,izbot2+1) + br(ix,iy,izbot2) ) 
!                   Ctr(ix,iy,izbot1) = 0.5 * (wr(ix,iy,izbot2+1)*tr(ix,iy,izbot2+1) + wr(ix,iy,izbot2+1)*tr(ix,iy,izbot2) + wr(ix,iy,izbot2)*tr(ix,iy,izbot2+1) +  wr(ix,iy,izbot2)*tr(ix,iy,izbot2) )

                else

                   Cur(ix,iy,izbot1) = 0.D0
                   Cvr(ix,iy,izbot1) = 0.D0
                   Cwr(ix,iy,izbot1) = 0.D0
                   Cbr(ix,iy,izbot1) = 0.D0

                end if

             end do
          end do

       else if(mype==(npe-1)) then
          
          do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   !Aur, Avr, Awr, Bur, Bvr and Bwr  are fine already.                                                                                                                          
                   Cur(ix,iy,iztop1) = - rho_u(iztop2-1)*wr(ix,iy,iztop2-1) * ( ur(ix,iy,iztop2) + ur(ix,iy,iztop2-1) )
                   Cvr(ix,iy,iztop1) = - rho_u(iztop2-1)*wr(ix,iy,iztop2-1) * ( vr(ix,iy,iztop2) + vr(ix,iy,iztop2-1) )
                   Cbr(ix,iy,iztop1) = - rho_u(iztop2-1)*wr(ix,iy,iztop2-1) * ( br(ix,iy,iztop2) + br(ix,iy,iztop2-1) )

                   Awr(ix,iy,iztop1) = 0.D0
                   Bwr(ix,iy,iztop1) = 0.D0
                   Cwr(ix,iy,iztop1) = 0.D0

!                   Atr(ix,iy,iztop1) = 0.D0
!                   Btr(ix,iy,iztop1) = 0.D0
!                   Ctr(ix,iy,iztop1) = 0.D0

                   Cwr(ix,iy,iztop1-1) = - 0.5* (rho_u(iztop2-1)*wr(ix,iy,iztop2-1)*wr(ix,iy,iztop2-2) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-1) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-2) )  !Original
!                   Cwr(ix,iy,iztop1-1) = - 0.5* (rho_u(iztop2-1)*wr(ix,iy,iztop2-1)*wr(ix,iy,iztop2-1) + rho_u(iztop2-1)*wr(ix,iy,iztop2-1)*wr(ix,iy,iztop2-2) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-1) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-2) )    !David's cheat

!                   Ctr(ix,iy,iztop1-1) = - 0.5* (wr(ix,iy,iztop2-1)*tr(ix,iy,iztop2-1) + wr(ix,iy,iztop2-1)*tr(ix,iy,iztop2-2) + wr(ix,iy,iztop2-2)*tr(ix,iy,iztop2-1) + wr(ix,iy,iztop2-2)*tr(ix,iy,iztop2-2) ) 
                   
                else

                   Cur(ix,iy,iztop1) = 0.D0
                   Cvr(ix,iy,iztop1) = 0.D0
                   Cbr(ix,iy,iztop1) = 0.D0

                   Awr(ix,iy,iztop1) = 0.D0
                   Bwr(ix,iy,iztop1) = 0.D0
                   Cwr(ix,iy,iztop1) = 0.D0

!                   Atr(ix,iy,iztop1) = 0.D0
!                   Btr(ix,iy,iztop1) = 0.D0
!                   Ctr(ix,iy,iztop1) = 0.D0

                   Cwr(ix,iy,iztop1-1) = 0.D0
!                   Ctr(ix,iy,iztop1-1) = 0.D0


                end if

             end do
          end do


       end if



      !Move to h-space

      call fft_r2c(Aur,Auk,n3h1)
      call fft_r2c(Avr,Avk,n3h1)
      call fft_r2c(Awr,Awk,n3h1)
      call fft_r2c(Abr,Abk,n3h1)

      call fft_r2c(Bur,Buk,n3h1)
      call fft_r2c(Bvr,Bvk,n3h1)
      call fft_r2c(Bwr,Bwk,n3h1)
      call fft_r2c(Bbr,Bbk,n3h1)

      call fft_r2c(Cur,Cuk,n3h1)
      call fft_r2c(Cvr,Cvk,n3h1)
      call fft_r2c(Cwr,Cwk,n3h1)
      call fft_r2c(Cbr,Cbk,n3h1)


      !nuk = ikx Auk + iky Buk + 1/2dz Cuk  (Same for nvk, nwk and ntk  but u->v, ->w or -> t)

      do izh1=1,n3h1   
         izh2=izh1+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nuk(ikx,iky,izh1) =  i*kx*Auk(ikx,iky,izh1)  +  i*ky*Buk(ikx,iky,izh1)  +  Cuk(ikx,iky,izh1)/(2*dz*rho_s(izh2))
                  nvk(ikx,iky,izh1) =  i*kx*Avk(ikx,iky,izh1)  +  i*ky*Bvk(ikx,iky,izh1)  +  Cvk(ikx,iky,izh1)/(2*dz*rho_s(izh2))
                  nwk(ikx,iky,izh1) =  i*kx*Awk(ikx,iky,izh1)  +  i*ky*Bwk(ikx,iky,izh1)  +  Cwk(ikx,iky,izh1)/(2*dz*rho_u(izh2))
                  nbk(ikx,iky,izh1) =  i*kx*Abk(ikx,iky,izh1)  +  i*ky*Bbk(ikx,iky,izh1)  +  Cbk(ikx,iky,izh1)/(2*dz*rho_s(izh2))
!                  ntk(ikx,iky,izh1) =  i*kx*Atk(ikx,iky,izh1)  +  i*ky*Btk(ikx,iky,izh1)  +  Ctk(ikx,iky,izh1)/(2*dz*rho_u(izh2))
               else
                  nuk(ikx,iky,izh1) = (0.D0,0.D0)
                  nvk(ikx,iky,izh1) = (0.D0,0.D0)
                  nwk(ikx,iky,izh1) = (0.D0,0.D0)
                  nbk(ikx,iky,izh1) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      !Copy back velocity to avoid inverse FFTs!

      uk=umem
      vk=vmem
      wk=wmem
      bk=bmem


    END SUBROUTINE convol

    !The original version!
    SUBROUTINE  convol2(nuk,nvk,nwk,ntk,nur,nvr,nwr,ntr,uk,vk,wk,tk,ur,vr,wr,tr)
      ! this subroutine computes ((u.grad)u)i = d/dxj(uj ui) in the divergence form on a staggered grid.
                                                                    
      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,tk   
      double complex, dimension(iktx,ikty,n3h1) :: nuk,nvk,nwk,ntk

      double precision,    dimension(n1d,n2d,n3h2) :: ur,vr,wr,tr
      double precision,    dimension(n1d,n2d,n3h1) :: nur,nvr,nwr,ntr

      double complex, dimension(iktx,ikty,n3h2) :: umem,vmem,wmem,tmem


      !Suboptimal but simple.
      double complex,   dimension(iktx,ikty,n3h1) :: Auk,Buk,Cuk,Avk,Bvk,Cvk,Awk,Bwk,Cwk,Atk,Btk,Ctk
      double precision, dimension(n1d,n2d,n3h1)   :: Aur,Bur,Cur,Avr,Bvr,Cvr,Awr,Bwr,Cwr,Atr,Btr,Ctr

      equivalence(Aur,Auk)
      equivalence(Avr,Avk)
      equivalence(Awr,Awk)
      equivalence(Atr,Atk)

      equivalence(Bur,Buk)
      equivalence(Bvr,Bvk)
      equivalence(Bwr,Bwk)
      equivalence(Btr,Btk)

      equivalence(Cur,Cuk)
      equivalence(Cvr,Cvk)
      equivalence(Cwr,Cwk)
      equivalence(Ctr,Ctk)


      !DO A COPY OF VELOCITY, FASTER THAN TAKE A TRANSFORM BACK                                                                                     

      umem=uk
      vmem=vk
      wmem=wk
      tmem=tk

      call fft_c2r(uk,ur,n3h2)
      call fft_c2r(vk,vr,n3h2)
      call fft_c2r(wk,wr,n3h2)
      call fft_c2r(tk,tr,n3h2)

      !Assign values to A, B and C for each of the three components of the nlt
      !nuk = ikx Auk + iky Buk + 1/2dz Cuk  (Same for nvk and nwk but u->v or ->w)
      !With A,B and C defined as:

      !Aur_n =  u_n^2
      !Bur_n =  v_n*u_n
      !Cur_n =  w_n*u_n+1 + w_n*u_n - u_n*w_n-1 - u_n-1*w_n-1

      !(same as for u)
      !Avr_n =  u_n*v_n
      !Bvr_n =  v_n^2
      !Cvr_n =  w_n*v_n+1 + w_n*v_n - w_n-1*v_n - w_n-1*v_n-1

      !Awr_n = 1/2 ( u_n+1*w_n + u_n*w_n)
      !Bwr_n = 1/2 ( v_n+1*w_n + v_n*w_n)
      !Cwr_n = 1/2 ( w_n+1^2 - w_n-1^2  ) + w_n * (w_n+1 - w_n-1)

      !Atr_n = 1/2 ( u_n+1*t_n + u_n*t_n)
      !Btr_n = 1/2 ( v_n+1*t_n + v_n*t_n)
      !Ctr_n = 1/2 [ w_n+1*t_n+1 + w_n+1*t_n + w_n*t_n+1 - w_n*t_n-1 - w_n-1*t_n - w_n-1*b_n-1 ]

      !*****************************************************************************!
      !Attention: I am starting to doubt this is the best way to do it. Two comments:
      ! 1. Using the analytic rho_0 sometimes is not NECESSARLY better than using rho_0 at the points dictated by FDs
      ! 2. For the interpolation of fields (ex. wwz and wbz terms) there are other possibilities that would turn out simpler.

      !    For instance d\dz(rho w b)| idz could be just [rho w b|(i+1)dz - rho w b|(i-1)]/2dz
      !    and I am not sure if this is better or worse than what we use now with complicated cross terms...

      ! When code will be re-written, I want to make sure I got the best discretization possible.

      do izh1=1,n3h1
         izh2=izh1+1
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   Aur(ix,iy,izh1) = ur(ix,iy,izh2)*ur(ix,iy,izh2)
                   Bur(ix,iy,izh1) = vr(ix,iy,izh2)*ur(ix,iy,izh2)
                   Cur(ix,iy,izh1) = rho_u(izh2)*wr(ix,iy,izh2)*ur(ix,iy,izh2+1) + rho_u(izh2)*wr(ix,iy,izh2)*ur(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*ur(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*ur(ix,iy,izh2-1) 

                   Avr(ix,iy,izh1) = ur(ix,iy,izh2)*vr(ix,iy,izh2)
                   Bvr(ix,iy,izh1) = vr(ix,iy,izh2)*vr(ix,iy,izh2)
                   Cvr(ix,iy,izh1) = rho_u(izh2)*wr(ix,iy,izh2)*vr(ix,iy,izh2+1) + rho_u(izh2)*wr(ix,iy,izh2)*vr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*vr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*vr(ix,iy,izh2-1) 

                   Awr(ix,iy,izh1) = 0.5 * ( ur(ix,iy,izh2+1)*wr(ix,iy,izh2)  +  ur(ix,iy,izh2)*wr(ix,iy,izh2) )    
                   Bwr(ix,iy,izh1) = 0.5 * ( vr(ix,iy,izh2+1)*wr(ix,iy,izh2)  +  vr(ix,iy,izh2)*wr(ix,iy,izh2) )    
                   Cwr(ix,iy,izh1) = 0.5 * ( rho_u(izh2+1)*wr(ix,iy,izh2+1)*wr(ix,iy,izh2+1) + rho_u(izh2+1)*wr(ix,iy,izh2+1)*wr(ix,iy,izh2) + rho_u(izh2)*wr(ix,iy,izh2)*wr(ix,iy,izh2+1) - rho_u(izh2)*wr(ix,iy,izh2)*wr(ix,iy,izh2-1) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*wr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*wr(ix,iy,izh2-1)  )

                   Atr(ix,iy,izh1) = 0.5 * ( ur(ix,iy,izh2+1)*tr(ix,iy,izh2)  +  ur(ix,iy,izh2)*tr(ix,iy,izh2) )    
                   Btr(ix,iy,izh1) = 0.5 * ( vr(ix,iy,izh2+1)*tr(ix,iy,izh2)  +  vr(ix,iy,izh2)*tr(ix,iy,izh2) )    
                   Ctr(ix,iy,izh1) = 0.5 * ( rho_u(izh2+1)*wr(ix,iy,izh2+1)*tr(ix,iy,izh2+1) + rho_u(izh2+1)*wr(ix,iy,izh2+1)*tr(ix,iy,izh2) + rho_u(izh2)*wr(ix,iy,izh2)*tr(ix,iy,izh2+1) - rho_u(izh2)*wr(ix,iy,izh2)*tr(ix,iy,izh2-1) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*tr(ix,iy,izh2) - rho_u(izh2-1)*wr(ix,iy,izh2-1)*tr(ix,iy,izh2-1)  )

                else

                   Aur(ix,iy,izh1) = 0.D0
                   Bur(ix,iy,izh1) = 0.D0
                   Cur(ix,iy,izh1) = 0.D0

                   Avr(ix,iy,izh1) = 0.D0
                   Bvr(ix,iy,izh1) = 0.D0
                   Cvr(ix,iy,izh1) = 0.D0

                   Awr(ix,iy,izh1) = 0.D0
                   Bwr(ix,iy,izh1) = 0.D0
                   Cwr(ix,iy,izh1) = 0.D0

                   Atr(ix,iy,izh1) = 0.D0
                   Btr(ix,iy,izh1) = 0.D0
                   Ctr(ix,iy,izh1) = 0.D0

                end if

             end do
          end do
       end do

       !Boundary conditions! X-marked regions were assigned with insignificant values. We give the correct values here

       if(mype==0) then

         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   !Aur, Avr, Bur and Bur are fine already.

                   Cur(ix,iy,izbot1) = rho_u(izbot2)*wr(ix,iy,izbot2) * ( ur(ix,iy,izbot2+1) + ur(ix,iy,izbot2) ) 
                   Cvr(ix,iy,izbot1) = rho_u(izbot2)*wr(ix,iy,izbot2) * ( vr(ix,iy,izbot2+1) + vr(ix,iy,izbot2) ) 
                   Cwr(ix,iy,izbot1) = 0.5 * (rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2+1) + rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*wr(ix,iy,izbot2) + rho_u(izbot2)*wr(ix,iy,izbot2)*wr(ix,iy,izbot2+1) )  
                   Ctr(ix,iy,izbot1) = 0.5 * (rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*tr(ix,iy,izbot2+1) + rho_u(izbot2+1)*wr(ix,iy,izbot2+1)*tr(ix,iy,izbot2) + rho_u(izbot2)*wr(ix,iy,izbot2)*tr(ix,iy,izbot2+1) )   !1/4dz (W2*B2 + W2*B1 + W1*B2)

                else

                   Cur(ix,iy,izbot1) = 0.D0
                   Cvr(ix,iy,izbot1) = 0.D0
                   Cwr(ix,iy,izbot1) = 0.D0
                   Ctr(ix,iy,izbot1) = 0.D0

                end if

             end do
          end do

       else if(mype==(npe-1)) then
          
          do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   !Aur, Avr, Awr, Bur, Bvr and Bwr  are fine already.                                                                                                                          
                   Cur(ix,iy,iztop1) = - rho_u(iztop2-1)*wr(ix,iy,iztop2-1) * ( ur(ix,iy,iztop2) + ur(ix,iy,iztop2-1) )
                   Cvr(ix,iy,iztop1) = - rho_u(iztop2-1)*wr(ix,iy,iztop2-1) * ( vr(ix,iy,iztop2) + vr(ix,iy,iztop2-1) )

                   Awr(ix,iy,iztop1) = 0.D0
                   Bwr(ix,iy,iztop1) = 0.D0
                   Cwr(ix,iy,iztop1) = 0.D0

                   Atr(ix,iy,iztop1) = 0.D0
                   Btr(ix,iy,iztop1) = 0.D0
                   Ctr(ix,iy,iztop1) = 0.D0

                   Cwr(ix,iy,iztop1-1) = - 0.5* (rho_u(iztop2-1)*wr(ix,iy,iztop2-1)*wr(ix,iy,iztop2-2) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-1) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*wr(ix,iy,iztop2-2) )  !1/4dz (WN-1*WN-2 + WN-2*WN-1 + WN-2*WN-2)
                   Ctr(ix,iy,iztop1-1) = - 0.5* (rho_u(iztop2-1)*wr(ix,iy,iztop2-1)*tr(ix,iy,iztop2-2) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*tr(ix,iy,iztop2-1) + rho_u(iztop2-2)*wr(ix,iy,iztop2-2)*tr(ix,iy,iztop2-2) )  !1/4dz (WN-1*BN-2 + WN-2*BN-1 + WN-2*BN-2)
                   
                else

                   Cur(ix,iy,iztop1) = 0.D0
                   Cvr(ix,iy,iztop1) = 0.D0

                   Awr(ix,iy,iztop1) = 0.D0
                   Bwr(ix,iy,iztop1) = 0.D0
                   Cwr(ix,iy,iztop1) = 0.D0

                   Atr(ix,iy,iztop1) = 0.D0
                   Btr(ix,iy,iztop1) = 0.D0
                   Ctr(ix,iy,iztop1) = 0.D0

                   Cwr(ix,iy,iztop1-1) = 0.D0
                   Ctr(ix,iy,iztop1-1) = 0.D0


                end if

             end do
          end do


       end if



      !Move to h-space

      call fft_r2c(Aur,Auk,n3h1)
      call fft_r2c(Avr,Avk,n3h1)
      call fft_r2c(Awr,Awk,n3h1)
      call fft_r2c(Atr,Atk,n3h1)

      call fft_r2c(Bur,Buk,n3h1)
      call fft_r2c(Bvr,Bvk,n3h1)
      call fft_r2c(Bwr,Bwk,n3h1)
      call fft_r2c(Btr,Btk,n3h1)

      call fft_r2c(Cur,Cuk,n3h1)
      call fft_r2c(Cvr,Cvk,n3h1)
      call fft_r2c(Cwr,Cwk,n3h1)
      call fft_r2c(Ctr,Ctk,n3h1)


      !nuk = ikx Auk + iky Buk + 1/2dz Cuk  (Same for nvk, nwk and ntk  but u->v, ->w or -> t)

      do izh1=1,n3h1   
         izh2=izh1+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nuk(ikx,iky,izh1) =  i*kx*Auk(ikx,iky,izh1)  +  i*ky*Buk(ikx,iky,izh1)  +  Cuk(ikx,iky,izh1)/(2*dz*rho_s(izh2))
                  nvk(ikx,iky,izh1) =  i*kx*Avk(ikx,iky,izh1)  +  i*ky*Bvk(ikx,iky,izh1)  +  Cvk(ikx,iky,izh1)/(2*dz*rho_s(izh2))
                  nwk(ikx,iky,izh1) =  i*kx*Awk(ikx,iky,izh1)  +  i*ky*Bwk(ikx,iky,izh1)  +  Cwk(ikx,iky,izh1)/(2*dz*rho_u(izh2))
                  ntk(ikx,iky,izh1) =  i*kx*Atk(ikx,iky,izh1)  +  i*ky*Btk(ikx,iky,izh1)  +  Ctk(ikx,iky,izh1)/(2*dz*rho_u(izh2))
               else
                  nuk(ikx,iky,izh1) = (0.D0,0.D0)
                  nvk(ikx,iky,izh1) = (0.D0,0.D0)
                  nwk(ikx,iky,izh1) = (0.D0,0.D0)
                  ntk(ikx,iky,izh1) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      !Copy back velocity to avoid inverse FFTs!

      uk=umem
      vk=vmem
      wk=wmem
      tk=tmem


    END SUBROUTINE convol2



    SUBROUTINE  convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)
      ! this subroutine computes ((u.grad)q) = d/dxj(uj q) in the divergence form on a staggered grid.
      ! Notice that this routine outputs the r-space velocity fields.
                                                                    
      double complex, dimension(iktx,ikty,n3h2) :: uk,vk   
      double complex, dimension(iktx,ikty,n3h1) :: qk   
      double complex, dimension(iktx,ikty,n3h0) :: nqk   

      double complex, dimension(iktx,ikty,n3h1) :: qmem   

      double precision,    dimension(n1d,n2d,n3h2) :: ur,vr
      double precision,    dimension(n1d,n2d,n3h1) :: qr
      double precision,    dimension(n1d,n2d,n3h0) :: nqr

      !Suboptimal but simple.
      double complex,   dimension(iktx,ikty,n3h0) :: Aqk,Bqk
      double precision, dimension(n1d,n2d,n3h0)   :: Aqr,Bqr

      equivalence(Aqr,Aqk)
      equivalence(Bqr,Bqk)

      !I think we don't need to keep a copy of u in qg, right?
      call fft_c2r(uk,ur,n3h2)
      call fft_c2r(vk,vr,n3h2)
     
      qmem = qk

      call fft_c2r(qk,qr,n3h1)

      !Assign values to A and B
      !nqk = ikx Aqk + iky Bqk 
      !With A,B and C defined as:

      !Aqr_n =  u_n*q_n
      !Bqr_n =  v_n*q_n

      Aqr=0.
      Bqr=0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   Aqr(ix,iy,izh0) = ur(ix,iy,izh2)*qr(ix,iy,izh1)
                   Bqr(ix,iy,izh0) = vr(ix,iy,izh2)*qr(ix,iy,izh1)

                end if
             end do
          end do
       end do

      !Move to h-space

      call fft_r2c(Aqr,Aqk,n3h0)
      call fft_r2c(Bqr,Bqk,n3h0)

      !nqk = ikx Aqk + iky Bqk

      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nqk(ikx,iky,izh0) =  i*kx*Aqk(ikx,iky,izh0)  +  i*ky*Bqk(ikx,iky,izh0)
               else
                  nqk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      qk = qmem

    END SUBROUTINE convol_q



    SUBROUTINE  convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
      ! this subroutine computes ((u.grad)q) = d/dxj(uj q) in the divergence form on the staggered grid.
      ! Notice that this routine outputs the r-space velocity fields.
                                                                    
      double complex, dimension(iktx,ikty,n3h2) :: uk,vk   
      double complex, dimension(iktx,ikty,n3h1) :: qk   
      double complex, dimension(iktx,ikty,n3h0) :: nqk   
      double complex, dimension(iktx,ikty,n3h0) :: BRk, BIk
      double complex, dimension(iktx,ikty,n3h0) :: nBRk, nBIk

      double complex, dimension(iktx,ikty,n3h1) :: qmem   
      double complex, dimension(iktx,ikty,n3h0) :: BRmem   
      double complex, dimension(iktx,ikty,n3h0) :: BImem   

      double precision, dimension(n1d,n2d,n3h2) :: ur,vr
      double precision, dimension(n1d,n2d,n3h1) :: qr
      double precision, dimension(n1d,n2d,n3h0) :: nqr
      double precision, dimension(n1d,n2d,n3h0) :: BRr, BIr
      double precision, dimension(n1d,n2d,n3h0) :: nBRr, nBIr


      !Terms to be differentiated
      double complex,   dimension(iktx,ikty,n3h0) :: utermk,vtermk
      double precision, dimension(n1d,n2d,n3h0)   :: utermr,vtermr

      equivalence(utermr,utermk)
      equivalence(vtermr,vtermk)

      !I think we don't need to keep a copy of u in qg, right?
      call fft_c2r(uk,ur,n3h2)
      call fft_c2r(vk,vr,n3h2)
     
      qmem  = qk
      BRmem = BRk
      BImem = BIk

      call fft_c2r(qk,qr,n3h1)
      call fft_c2r(BRk,BRr,n3h0)
      call fft_c2r(BIk,BIr,n3h0)
      


      ! ---- J(psi,q) ---- !

      !Assign values to uterm and vterm
      !nqk = ikx utermk + iky vtermk 
      !With A,B and C defined as:

      !utermr_n =  u_n*q_n
      !utermr_n =  v_n*q_n

      utermr=0.
      vtermr=0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   utermr(ix,iy,izh0) = ur(ix,iy,izh2)*qr(ix,iy,izh1)
                   vtermr(ix,iy,izh0) = vr(ix,iy,izh2)*qr(ix,iy,izh1)

                end if
             end do
          end do
       end do

      !Move to k-space

      call fft_r2c(utermr,utermk,n3h0)
      call fft_r2c(vtermr,vtermk,n3h0)

      !nqk = ikx utermk + iky vtermk

      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nqk(ikx,iky,izh0) =  i*kx*utermk(ikx,iky,izh0)  +  i*ky*vtermk(ikx,iky,izh0)
               else
                  nqk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      qk = qmem






      ! ---- J(psi,BR) ---- !

      utermr=0.
      vtermr=0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   utermr(ix,iy,izh0) = ur(ix,iy,izh2)*BRr(ix,iy,izh0)
                   vtermr(ix,iy,izh0) = vr(ix,iy,izh2)*BRr(ix,iy,izh0)

                end if
             end do
          end do
       end do

      !Move to k-space

      call fft_r2c(utermr,utermk,n3h0)
      call fft_r2c(vtermr,vtermk,n3h0)

      !nqk = ikx utermk + iky vtermk

      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nBRk(ikx,iky,izh0) =  i*kx*utermk(ikx,iky,izh0)  +  i*ky*vtermk(ikx,iky,izh0)
               else
                  nBRk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      BRk = BRmem





      ! ---- J(psi,BI) ---- !

      utermr=0.
      vtermr=0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   utermr(ix,iy,izh0) = ur(ix,iy,izh2)*BIr(ix,iy,izh0)
                   vtermr(ix,iy,izh0) = vr(ix,iy,izh2)*BIr(ix,iy,izh0)

                end if
             end do
          end do
       end do

      !Move to k-space

      call fft_r2c(utermr,utermk,n3h0)
      call fft_r2c(vtermr,vtermk,n3h0)

      !nqk = ikx utermk + iky vtermk

      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then                                                               
                  nBIk(ikx,iky,izh0) =  i*kx*utermk(ikx,iky,izh0)  +  i*ky*vtermk(ikx,iky,izh0)
               else
                  nBIk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      BIk = BImem










    END SUBROUTINE convol_waqg






    SUBROUTINE  refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)
      ! this subroutine computes the refractive term: ~B*zeta
      ! Notice that this routine outputs the r-space velocity fields.

      double complex, dimension(iktx,ikty,n3h1) :: psik   
      double complex, dimension(iktx,ikty,n3h0) :: zetak 
      double complex, dimension(iktx,ikty,n3h0) :: BRk, BIk
      double complex, dimension(iktx,ikty,n3h0) :: rBRk, rBIk


      double precision, dimension(n1d,n2d,n3h1) :: psir
      double precision, dimension(n1d,n2d,n3h0) :: zetar
      double precision, dimension(n1d,n2d,n3h0) :: BRr, BIr
      double precision, dimension(n1d,n2d,n3h0) :: rBRr, rBIr

      double complex, dimension(iktx,ikty,n3h0) :: BRmem   
      double complex, dimension(iktx,ikty,n3h0) :: BImem   

      equivalence(zetar , zetak)

      !Compute vertical vorticity from the streamfunction: zetak = - kh2 * psik
      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if(L(ikx,iky)==1) then
                  zetak(ikx,iky,izh0) = -kh2*psik(ikx,iky,izh1)
               else
                  zetak(ikx,iky,izh0) = (0.D0,0.D0)
               end if
            enddo
         enddo
      end do

      call fft_c2r(zetak,zetar,n3h0)

      BRmem = BRk
      BImem = BIk

      call fft_c2r(BRk,BRr,n3h0)
      call fft_c2r(BIk,BIr,n3h0)
      


      ! ---- B*zeta ---- !

      rBRr=0.
      rBIr=0.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   rBRr(ix,iy,izh0) = zetar(ix,iy,izh0)*BRr(ix,iy,izh0)
                   rBIr(ix,iy,izh0) = zetar(ix,iy,izh0)*BIr(ix,iy,izh0)

                end if
             end do
          end do
       end do

      !Move to k-space

      call fft_r2c(rBRr,rBRk,n3h0)
      call fft_r2c(rBIr,rBIk,n3h0)

      !Dealias

      do izh0=1,n3h0 
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.0) then                                                               
                  rBRk(ikx,iky,izh0) = (0.D0,0.D0)
                  rBIk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      BRk = BRmem
      BIk = BImem



    END SUBROUTINE refraction_waqg





    SUBROUTINE  compute_qw(qwk,BRk,BIk,qwr,BRr,BIr)
      ! this subroutine computes the feedback of near-inertial waves onto QGPV

      double complex,   dimension(iktx,ikty,n3h0) :: qwk
      double precision, dimension(n1d,n2d,n3h0)   :: qwr

      double complex, dimension(iktx,ikty,n3h0) :: BRk, BIk
      double precision, dimension(n1d,n2d,n3h0) :: BRr, BIr

      double complex, dimension(iktx,ikty,n3h0) :: BRmem   
      double complex, dimension(iktx,ikty,n3h0) :: BImem   

      !Can probably be optimized
      double complex,   dimension(iktx,ikty,n3h0) :: BRxk,BRyk,BIxk,BIyk
      double precision, dimension(n1d,n2d,n3h0)   :: BRxr,BRyr,BIxr,BIyr

      !This could be one of the above: to optimize
      double complex,   dimension(iktx,ikty,n3h0) :: tempk
      double precision, dimension(n1d,n2d,n3h0)   :: tempr

      equivalence(BRxr,BRxk)
      equivalence(BRyr,BRyk)
      equivalence(BIxr,BIxk)
      equivalence(BIyr,BIyk)

      equivalence(tempr,tempk)

      BRxk = (0.D0,0.D0)
      BRyk = (0.D0,0.D0)
      BIxk = (0.D0,0.D0)
      BRyk = (0.D0,0.D0)

      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  BRxk(ikx,iky,izh0) =  i*kx*BRk(ikx,iky,izh0)
                  BRyk(ikx,iky,izh0) =  i*ky*BRk(ikx,iky,izh0)
                  BIxk(ikx,iky,izh0) =  i*kx*BIk(ikx,iky,izh0)
                  BIyk(ikx,iky,izh0) =  i*ky*BIk(ikx,iky,izh0)
               endif
            enddo
         enddo
      enddo

      call fft_c2r(BRxk,BRxr,n3h0)
      call fft_c2r(BRyk,BRyr,n3h0)
      call fft_c2r(BIxk,BIxr,n3h0)
      call fft_c2r(BIyk,BIyr,n3h0)


      ! ---- 1st term: (i/2) J(B*,B) ---- !

      ! Note that (i/2) J(B*,B) = (i/2)*(2i) J(Br,Bi) = - J(Br,Bi) = (Bry Bix - Brx Biy)


      qwr=0.
      

      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   qwr(ix,iy,izh0) = BRyr(ix,iy,izh0)*BIxr(ix,iy,izh0) - BRxr(ix,iy,izh0)*BIyr(ix,iy,izh0)

                end if
             end do
          end do
       end do


      ! ---- 2nd term: |B|^2 ---- !       

      BRmem = BRk
      BImem = BIk

      call fft_c2r(BRk,BRr,n3h0)
      call fft_c2r(BIk,BIr,n3h0)

      tempr = 0.D0

      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   tempr(ix,iy,izh0) = BRr(ix,iy,izh0)*BRr(ix,iy,izh0) + BIr(ix,iy,izh0)*BIr(ix,iy,izh0)

                end if
             end do
          end do
       end do      

      call fft_r2c(tempr,tempk,n3h0)


      ! --- Add the 1st and 2nd term to get qw --- !

      call fft_r2c(qwr,qwk,n3h0)
            
      
      do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2 = kx*kx + ky*ky
               if (L(ikx,iky).eq.1) then
                  qwk(ikx,iky,izh0) = qwk(ikx,iky,izh0) - 0.25*kh2*tempk(ikx,iky,izh0)
               else
                  qwk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo


      ! --- Get the right dimensions --- !
      qwk = qwk*Ro*W2F  

      BRk = BRmem
      BIk = BImem

    END SUBROUTINE compute_qw




    SUBROUTINE  compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)
      ! this subroutine computes the vertical integral of A at every wavenumber.

      double complex, dimension(iktx,ikty,2)               :: sigma_to_reduce     !This is the sum local to each processor. Last dimension: 1=real 2=imag
      double complex, dimension(iktx,ikty,2), intent(out)  :: sigma               !This is the global sum after all processors shared theirs

      !**** n = nonlinear advection term J(psi,B) **** r = refractive term ~ B*vort                                                                                            
      double complex,   dimension(iktx,ikty,n3h0), intent(in) :: nBRk, nBIk, rBRk, rBIk

      sigma_to_reduce = (0.D0,0.D0)
      sigma           = (0.D0,0.D0)

      !There is the Coriolis parameter missing: figure dimensions out
      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2 = kx*kx + ky*ky
               if ((L(ikx,iky).eq.1) .and. kh2 > 0) then
                  sigma_to_reduce(ikx,iky,1) = sigma_to_reduce(ikx,iky,1) + (rBRk(ikx,iky,izh0) + 2*nBIk(ikx,iky,izh0))/(1.D0*kh2) 
                  sigma_to_reduce(ikx,iky,2) = sigma_to_reduce(ikx,iky,2) + (rBIk(ikx,iky,izh0) - 2*nBRk(ikx,iky,izh0))/(1.D0*kh2) 
               endif
            enddo
         enddo
      enddo

      call mpi_allreduce(sigma_to_reduce,sigma,iktx*ikty*2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierror)

      ! Get the right dimensions !
      sigma = sigma*Bu*Ro

    END SUBROUTINE compute_sigma



    SUBROUTINE  sumB(BRk, BIk)
      ! this subroutine computes the vertical integral of B at every wavenumber.

      double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk

      double complex, dimension(iktx,ikty,2) :: sum_to_reduce     !This is the sum local to each processor. Last dimension: 1=real 2=imag
      double complex, dimension(iktx,ikty,2) :: aveB              !This is the global sum after all processors shared theirs. To be divided by n3 for the average B at each (kx,ky)

      sum_to_reduce = (0.D0,0.D0)
      aveB          = (0.D0,0.D0)

      !There is the Coriolis parameter missing: figure dimensions out
      do izh0=1,n3h0 
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2 = kx*kx + ky*ky
               if ((L(ikx,iky).eq.1) .and. kh2 > 0) then
                  sum_to_reduce(ikx,iky,1) = sum_to_reduce(ikx,iky,1) + BRk(ikx,iky,izh0)
                  sum_to_reduce(ikx,iky,2) = sum_to_reduce(ikx,iky,2) + BIk(ikx,iky,izh0)
               endif
            enddo
         enddo
      enddo

      call mpi_allreduce(sum_to_reduce,aveB,iktx*ikty*2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierror)

      aveB = aveB/n3

      do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2 = kx*kx + ky*ky
               if ((L(ikx,iky).eq.1) .and. kh2 > 0) then
                  BRk(ikx,iky,izh0) = BRk(ikx,iky,izh0) - aveB(ikx,iky,1)
                  BIk(ikx,iky,izh0) = BIk(ikx,iky,izh0) - aveB(ikx,iky,2)
               endif
            enddo
         enddo
      enddo
   
         
    END SUBROUTINE sumB




















    SUBROUTINE compute_A(ARk,AIK,BRkt,BIkt,CRk,CIk,sigma)
 
      double complex, dimension(iktx,ikty,2), intent(in)  :: sigma               !This is the global sum after all processors shared theirs                                    
      double complex,   dimension(iktx,ikty,n3h0), intent(out) :: ARk, AIk
      double complex,   dimension(iktx,ikty,n3h0), intent(out) :: CRk, CIk       !C=A_z on the unstaggered grid
 
       double complex, dimension(iktx,n3, iktyp), intent(in) :: BRkt          !Transposed (ky-parallelization) BRk
       double complex, dimension(iktx,n3, iktyp), intent(in) :: BIkt          !Transposed (ky-parallelization) BIk 
 
      !This transposed field exists only inside this subroutine. (Temporary solution for non-haloed A: compute while parallelized in ky)
      double complex, dimension(iktx,n3, iktyp) :: CRkt          !Transposed (ky-parallelization) BRk
      double complex, dimension(iktx,n3, iktyp) :: CIkt          !Transposed (ky-parallelization) BIk 

       !Just for internal usage: A to-be-transposed version.
      double complex, dimension(iktx,n3, iktyp) :: ARkt          !Transposed (ky-parallelization) BIk                                                            
      double complex, dimension(iktx,n3, iktyp) :: AIkt          !Transposed (ky-parallelization) BIk                                                               
 
      double complex, dimension(iktx, iktyp) :: sumAR, sumAI, sumBR, sumBI
      
      !-Initialize to 0-!
      ARkt  = (0.D0,0.D0)
      AIkt  = (0.D0,0.D0)
      sumAR = (0.D0,0.D0)
      sumAI = (0.D0,0.D0)
      sumBR = (0.D0,0.D0)
      sumBI = (0.D0,0.D0)

      CRkt  = (0.D0,0.D0)
      CIkt  = (0.D0,0.D0)



      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)
            kh2=kx*kx + ky*ky
            if(kh2/=0 .and. L(ikx,iky)==1 ) then
               sumBR(ikx,ikyp) = BRkt(ikx,1,ikyp)
               sumBI(ikx,ikyp) = BIkt(ikx,1,ikyp)
            end if
         end DO
      end DO

      !Compute \tilde{A}, which is \hat{A} up to an arbitrary constant (\hat{A} at z = dz/2 set to 0)
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then

               do iz=2,n3

                  ARkt(ikx,iz,ikyp) = ARkt(ikx,iz-1,ikyp) + sumBR(ikx,ikyp)*r_2ut(iz-1)*dz*dz
                  AIkt(ikx,iz,ikyp) = AIkt(ikx,iz-1,ikyp) + sumBI(ikx,ikyp)*r_2ut(iz-1)*dz*dz

                  sumAR(ikx,ikyp) = sumAR(ikx,ikyp) + ARkt(ikx,iz,ikyp)
                  sumAI(ikx,ikyp) = sumAI(ikx,ikyp) + AIkt(ikx,iz,ikyp)

                  sumBR(ikx,ikyp) = sumBR(ikx,ikyp) + BRkt(ikx,iz,ikyp)
                  sumBI(ikx,ikyp) = sumBI(ikx,ikyp) + BIkt(ikx,iz,ikyp)

               end do
            end if

         end DO
      end DO


      !Compute \hat{A}, the actual solution we are looking for
      DO ikx=1,iktx
         kx=kxa(ikx)
         DO ikyp=1,iktyp
            iky=ikyp+iktyp*mype
            ky=kya(iky)

            kh2=kx*kx + ky*ky

            if(kh2/=0 .and. L(ikx,iky)==1 ) then

               do iz=1,n3

                  ARkt(ikx,iz,ikyp) = ARkt(ikx,iz,ikyp) + ( sigma(ikx,iky,1) - sumAR(ikx,ikyp) )/n3
                  AIkt(ikx,iz,ikyp) = AIkt(ikx,iz,ikyp) + ( sigma(ikx,iky,2) - sumAI(ikx,ikyp) )/n3

               end do
            end if

         end DO
      end DO


      !Temporary solution for potential energy: compute C=A_z while in ky-parallelized space.
      DO ikx=1,iktx
         DO ikyp=1,iktyp

               do iz=1,n3-1

                  CRkt(ikx,iz,ikyp) = ( ARkt(ikx,iz+1,ikyp) - ARkt(ikx,iz,ikyp) )/dz 
                  CIkt(ikx,iz,ikyp) = ( AIkt(ikx,iz+1,ikyp) - AIkt(ikx,iz,ikyp) )/dz 

               end do

               !Top: C = A_z = 0

               CRkt(ikx,n3,ikyp) = (0.D0,0.D0)
               CIkt(ikx,n3,ikyp) = (0.D0,0.D0)

         end DO
      end DO


      !Transpose A back to the regular z-parallelized world
      call mpitranspose(ARkt,iktx,n3,iktyp,ARk,ikty,n3h0)
      call mpitranspose(AIkt,iktx,n3,iktyp,AIk,ikty,n3h0)

      !Transpose C to the regular z-parallelized world
      call mpitranspose(CRkt,iktx,n3,iktyp,CRk,ikty,n3h0)
      call mpitranspose(CIkt,iktx,n3,iktyp,CIk,ikty,n3h0)

    END SUBROUTINE compute_A
      


    SUBROUTINE compute_streamfunction(uk,vk,psik)  !Computes the QG streamfunction (rotational part of the horizontal flow)

      double complex,   dimension(iktx,ikty,n3h2), intent(in)  :: uk,vk
      double complex,   dimension(iktx,ikty,n3h1), intent(out) :: psik 

      !Compute psi                                                                                                                                                                                                                                
      do izh1=1,n3h1
         izh2=izh1+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2  = kx*kx+ky*ky
               
               psik(ikx,iky,izh1) = i * ( ky*uk(ikx,iky,izh2) - kx*vk(ikx,iky,izh2) )/kh2
               
               if(kh2==0) psik(ikx,iky,izh1) = (0.D0,0.D0)
               
            end do
         end do
      end do
         

    end SUBROUTINE compute_streamfunction



    SUBROUTINE compute_velo(uk,vk,wk,bk,psik)  !Computes the QG velocity and buoyancy fields from the streamfunction, assuming psi_z=0 at the top/bottom. bk is staggered. 
      
      double complex,   dimension(iktx,ikty,n3h2), intent(out) :: uk,vk,wk,bk  
      double complex,   dimension(iktx,ikty,n3h1), intent(in)  :: psik         


      !Compute velocity/buoyancy fields                                                                                                                                                                        
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
                  bk(ikx,iky,izh2) =    ( psik(ikx,iky,izh1+1) - psik(ikx,iky,izh1-1) )/(r_1s(izh2)*2.*dz)    !1/r_1 d psi/dz                                                                                   
               else
                  uk(ikx,iky,izh2) = (0.D0,0.D0)
                  vk(ikx,iky,izh2) = (0.D0,0.D0)
                  wk(ikx,iky,izh2) = (0.D0,0.D0)
                  bk(ikx,iky,izh2) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo
      
      !Explicit treatment of boundary conditions                                                                                                                                                                        
      if(mype==0) then
         do iky=1,ikty
            do ikx=1,iktx
               bk(ikx,iky,izbot2) = ( psik(ikx,iky,izbot1+1) - psik(ikx,iky,izbot1) )/(r_1s(izbot2)*2.*dz)    !1/r_1 d psi/dz                                                                                             
            enddo
         enddo
      elseif(mype==(npe-1)) then
         do iky=1,ikty
            do ikx=1,iktx
               bk(ikx,iky,iztop2) = ( psik(ikx,iky,iztop1) - psik(ikx,iky,iztop1-1) )/(r_1s(iztop2)*2.*dz)    !1/r_1 d psi/dz                                                                                           
            enddo
         enddo
      end if

    end SUBROUTINE compute_velo






    SUBROUTINE omega_verification2(wak,lhs)

      double complex,   dimension(iktx,ikty,n3h1) :: wak                                       
      double complex,   dimension(iktx,ikty,n3h0) :: lhs                                       

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,izh0) = (wak(ikx,iky,izh1+1) -2.*wak(ikx,iky,izh1) +wak(ikx,iky,izh1-1))/(dz*dz) + r_3u(izh2)*(wak(ikx,iky,izh1+1) -wak(ikx,iky,izh1-1))/(2.*dz)  +     (r_5u(izh2)-kh2/a_ell_u(izh2))*wak(ikx,iky,izh1)
               else
                  lhs(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      if(mype==0) then
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,1) = (wak(ikx,iky,izbot1+1) -2.*wak(ikx,iky,izbot1) )/(dz*dz) + r_3u(izbot2)*wak(ikx,iky,izbot1+1)/(2.*dz)  +     (r_5u(izbot2)-kh2/a_ell_u(izbot2))*wak(ikx,iky,izbot1)
               else
                  lhs(ikx,iky,1) = (0.D0,0.D0)
               endif
            enddo
         enddo
      end if


      if(mype==(npe-1)) then
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,n3h0-1) = (wak(ikx,iky,iztop1-2) -2.*wak(ikx,iky,iztop1-1))/(dz*dz) - r_3u(iztop2-1)*wak(ikx,iky,iztop1-2)/(2.*dz)  +     (r_5u(iztop2-1)-kh2/a_ell_u(iztop2-1))*wak(ikx,iky,iztop1-1)
                  lhs(ikx,iky,n3h0) = (0.D0,0.D0)
               else
                  lhs(ikx,iky,n3h0-1) = (0.D0,0.D0)
                  lhs(ikx,iky,n3h0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      end if

     

    end SUBROUTINE omega_verification2




    SUBROUTINE omega_verification(wak,lhs)

      double complex,   dimension(iktx,ikty,n3h1) :: wak                                       
      double complex,   dimension(iktx,ikty,n3h0) :: lhs                                       

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,izh0) = (rho_u(izh2+1)/rho_s(izh2+1))*wak(ikx,iky,izh1+1) - (rho_u(izh2)*(rho_s(izh2+1)+rho_s(izh2))/(rho_s(izh2+1)*rho_s(izh2)) + r_1(izh2)*r_2(izh2)*dz*dz*kh2/Bu)*wak(ikx,iky,izh1) +  (rho_u(izh2-1)/rho_s(izh2))*wak(ikx,iky,izh1-1)
               else
                  lhs(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      if(mype==0) then
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,1) = (rho_u(izbot2+1)/rho_s(izbot2+1))*wak(ikx,iky,izbot1+1) - (rho_u(izbot2)*(rho_s(izbot2+1)+rho_s(izbot2))/(rho_s(izbot2+1)*rho_s(izbot2)) + r_1(izbot2)*r_2(izbot2)*dz*dz*kh2/Bu)*wak(ikx,iky,izbot1) 
               else
                  lhs(ikx,iky,1) = (0.D0,0.D0)
               endif
            enddo
         enddo
      end if


      if(mype==(npe-1)) then
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               kh2=kx*kx+ky*ky
               if (L(ikx,iky).eq.1) then
                  lhs(ikx,iky,n3h0-1) = - (rho_u(iztop2-1)*(rho_s(iztop2)+rho_s(iztop2-1))/(rho_s(iztop2)*rho_s(iztop2-1)) + r_1(iztop2-1)*r_2(iztop2-1)*dz*dz*kh2/Bu)*wak(ikx,iky,iztop1-1) +  (rho_u(iztop2-2)/rho_s(iztop2-1))*wak(ikx,iky,iztop1-2)
                  lhs(ikx,iky,n3h0) = (0.D0,0.D0)
               else
                  lhs(ikx,iky,n3h0-1) = (0.D0,0.D0)
                  lhs(ikx,iky,n3h0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      end if

      lhs=lhs/(dz*dz)

    end SUBROUTINE omega_verification




    SUBROUTINE  omega_eqn_rhs(rhs,rhsr,psik)
      ! this subroutine computes the right hand side of the omega equation = 2 J(psi_z,nabla^2 psi) 
      ! Notice that this routine outputs the r-space velocity fields.

      double complex,   dimension(iktx,ikty,n3h1) :: psik                                       
      double complex,   dimension(iktx,ikty,n3h0) :: rhs                             
      double precision, dimension(n1d,n2d,n3h0)   :: rhsr

      !Suboptimal but simple.
      double complex,   dimension(iktx,ikty,n3h0) :: bxk,byk,xxk,xyk 
      double precision, dimension(n1d,n2d,n3h0)   :: bxr,byr,xxr,xyr

      equivalence(bxr,bxk)
      equivalence(byr,byk)
      equivalence(xxr,xxk)
      equivalence(xyr,xyk)

      !First compute psik_z and store (temporarily) in rhs

      rhs = (0.D0,0.D0)

      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  rhs(ikx,iky,izh0) =  ( psik(ikx,iky,izh1+1) - psik(ikx,iky,izh1) )/dz
               else
                  rhs(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo
      
     if(mype==(npe-1)) then
         do iky=1,ikty
            do ikx=1,iktx
               rhs(ikx,iky,n3h0) = (0.D0,0.D0)
            enddo
        enddo
     end if



     !Now compute psi_z and lap(psi)'s horizontal derivatives...

     bxk = (0.D0,0.D0)
     byk = (0.D0,0.D0)
     xxk = (0.D0,0.D0)
     xyk = (0.D0,0.D0)

      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  bxk(ikx,iky,izh0) = i*kx*rhs(ikx,iky,izh0)
                  byk(ikx,iky,izh0) = i*ky*rhs(ikx,iky,izh0)
                  xxk(ikx,iky,izh0) = -i*kx*(kx*kx+ky*ky)*0.5*( psik(ikx,iky,izh1+1) + psik(ikx,iky,izh1) )   !Need to interpolate so everything is unstaggered.
                  xyk(ikx,iky,izh0) = -i*ky*(kx*kx+ky*ky)*0.5*( psik(ikx,iky,izh1+1) + psik(ikx,iky,izh1) )
               else
                  bxk(ikx,iky,izh0) = (0.D0,0.D0)
                  byk(ikx,iky,izh0) = (0.D0,0.D0)
                  xxk(ikx,iky,izh0) = (0.D0,0.D0)
                  xyk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      if(mype==(npe-1)) then
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  bxk(ikx,iky,n3h0) = (0.D0,0.D0)
                  byk(ikx,iky,n3h0) = (0.D0,0.D0)
                  xxk(ikx,iky,n3h0) = -i*kx*(kx*kx+ky*ky)*psik(ikx,iky,iztop1)
                  xyk(ikx,iky,n3h0) = -i*ky*(kx*kx+ky*ky)*psik(ikx,iky,iztop1)
               else
                  bxk(ikx,iky,n3h0) = (0.D0,0.D0)
                  byk(ikx,iky,n3h0) = (0.D0,0.D0)
                  xxk(ikx,iky,n3h0) = (0.D0,0.D0)
                  xyk(ikx,iky,n3h0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      end if

      call fft_c2r(bxk,bxr,n3h0)
      call fft_c2r(byk,byr,n3h0)
      call fft_c2r(xxk,xxr,n3h0)
      call fft_c2r(xyk,xyr,n3h0)

      rhsr = 0.D0

      !Now compute the r-space RHS.

      do izh0=1,n3h0
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then
                   rhsr(ix,iy,izh0) = 2.D0*( bxr(ix,iy,izh0)*xyr(ix,iy,izh0) - byr(ix,iy,izh0)*xxr(ix,iy,izh0) )
                else
                   rhsr(ix,iy,izh0) = 0.D0
                end if

             end do
          end do
       end do

      call fft_r2c(rhsr,rhs,n3h0)


    END SUBROUTINE omega_eqn_rhs















   SUBROUTINE  compute_w(psik,psi_old,uk,vk,wak,ur,vr)

      !This subroutine computes the first order vertical velocity w_1 = -1/(r1r2) (Fr/Ro)^2 D/Dt psi_z

      ! this subroutine computes ((u.grad)psiz) = d/dxj(uj psiz) in the divergence form on a unstaggered grid.                                                                                                      

      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wak
      double complex, dimension(iktx,ikty,n3h0) :: psizk
      double complex,   dimension(iktx,ikty,n3h1) :: psik   

      double complex, dimension(iktx,ikty,n3h1) :: psi_old
      double complex, dimension(iktx,ikty,n3h0) :: psiz_old

      double precision,    dimension(n1d,n2d,n3h2) :: ur,vr
      double precision,    dimension(n1d,n2d,n3h0) :: psizr


      double complex, dimension(iktx,ikty,n3h2) :: umem,vmem

      !Suboptimal but simple.                                                                                                                                                                                                                
      double complex,   dimension(iktx,ikty,n3h0) :: Aqk,Bqk
      double precision, dimension(n1d,n2d,n3h0)   :: Aqr,Bqr

      equivalence(Aqr,Aqk)
      equivalence(Bqr,Bqk)
      equivalence(psizr,psizk)


      !First compute psi_z at current and past time steps:

      
      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  psizk(ikx,iky,izh0)    = ( psik(ikx,iky,izh1+1)    - psik(ikx,iky,izh1)    )/dz 
                  psiz_old(ikx,iky,izh0) = ( psi_old(ikx,iky,izh1+1) - psi_old(ikx,iky,izh1) )/dz 
              else
                  psizk(ikx,iky,izh0)    = (0.,0.)
                  psiz_old(ikx,iky,izh0) = (0.,0.)
              endif
           enddo
        enddo
     enddo


     if(mype==(npe-1)) then !Correction at the top... psi_z = 0 right at the top and bot...
         do iky=1,ikty
            do ikx=1,iktx
                  psizk(ikx,iky,n3h0)    = (0.,0.)
                  psiz_old(ikx,iky,n3h0) = (0.,0.)
           enddo
        enddo
     end if

     !Compute the Eularian time derivative of psi_z using Euler \partial_t psi_z = [ psi_z(n+1) - psi_z(n) ] / dt
     !Overwrite the psiz_old

     psiz_old = ( psizk - psiz_old )/delt
     
     !Now compute the horizontal advection of psi_z (at new time step n+1)

      !I think we don't need to keep a copy of u in qg, right? We sure do here.
      umem = uk
      vmem = vk

      call fft_c2r(uk,ur,n3h2)
      call fft_c2r(vk,vr,n3h2)
      call fft_c2r(psizk,psizr,n3h0)

      !Assign values to A and B                                                                                                                                                                                                
      !nqk = ikx Aqk + iky Bqk                                                                                                                                                                                                      
      !With A and B defined as:                                                                                                                                                                                               

      !Aqr_n =  0.5*(u_n+u_n+1)*psiz_n                                                                                                                                                                                          
      !Bqr_n =  0.5*(v_n+v_n+1)*psiz_n                                                                                                                                                                                                  

      !Indeed, psiz is unstag while u,v are stag

      Aqr=0.
      Bqr=0.

      do izh0=1,n3h0
         izh2=izh0+2
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

!                   Aqr(ix,iy,izh0) = ur(ix,iy,izh2)*psizr(ix,iy,izh0)
!                   Bqr(ix,iy,izh0) = vr(ix,iy,izh2)*psizr(ix,iy,izh0)

                   Aqr(ix,iy,izh0) = 0.5*(ur(ix,iy,izh2)+ ur(ix,iy,izh2+1))*psizr(ix,iy,izh0)
                   Bqr(ix,iy,izh0) = 0.5*(vr(ix,iy,izh2)+ vr(ix,iy,izh2+1))*psizr(ix,iy,izh0)

                end if
             end do
          end do
       end do

       if(mype==(npe-1)) then   !Correction at the top (u(N)=u(N+1)), but psiz = 0 anyway...
         do ix=1,n1d
             do iy=1,n2d
                if(ix<=n1) then

                   Aqr(ix,iy,n3h0) = (0.,0.)
                   Bqr(ix,iy,n3h0) = (0.,0.)

                end if
             end do
          end do
       end if


      !Move to h-space                                                                                                                                                                                                               

      call fft_r2c(Aqr,Aqk,n3h0)
      call fft_r2c(Bqr,Bqk,n3h0)

      !nqk = ikx Aqk + iky Bqk    
      !Overwrite Aqk instead of creating new field... Aqk = u_g cdot grad(psi_z) after the loop
      do izh0=1,n3h0
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  Aqk(ikx,iky,izh0) =  i*kx*Aqk(ikx,iky,izh0)  +  i*ky*Bqk(ikx,iky,izh0)
               else
                  Aqk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      uk = umem
      vk = vmem


      !Now bring everything together to compute w.

      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  wak(ikx,iky,izh2)    = -(1./(r_1(izh2)*r_2(izh2)))*(Fr/Ro)*(Fr/Ro)*(psiz_old(ikx,iky,izh0) +  Aqk(ikx,iky,izh0) )
              else
                  wak(ikx,iky,izh2)    = (0.,0.)
              endif
           enddo
        enddo
     enddo


   END SUBROUTINE compute_w







    SUBROUTINE divergence(usk,vsk,wsk,rhs)

      !Computes the divergence of n3h1 u* (not general for now) and put result in the n3h0 field RHS

      double complex,   dimension(iktx,ikty,n3h1) :: usk,vsk,wsk   !u*!                                                                                                     
      double complex,   dimension(iktx,ikty,n3h0) :: rhs       !pressure, and rhs of pressure equation! 


      do izh0=1,n3h0
         izh1=izh0+1
         izh2=izh0+2
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  rhs(ikx,iky,izh0) =  i*kx*usk(ikx,iky,izh1) + i*ky*vsk(ikx,iky,izh1) + (0.5*r_3(izh2) + 1./dz)*wsk(ikx,iky,izh1) + (0.5*r_3(izh2) - 1./dz)* wsk(ikx,iky,izh1-1)  
               else
                  rhs(ikx,iky,izh0) = (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      !Since wsk is defined properly at the top and bottom (z=0), no need to adjust manually


      !If b=0 at bnds:
      !Boundary conditions: must manually set that of the bottom, top not necessary
!      if(mype==0) then

!         do iky=1,ikty
!            ky = kya(iky)
!            do ikx=1,iktx
!               kx = kxa(ikx)
!               if (L(ikx,iky).eq.1) then
!                  rhs(ikx,iky,1) =  i*kx*usk(ikx,iky,izbot1) + i*ky*vsk(ikx,iky,izbot1)  + (0.5*r_3(izbot2) + 1./dz)*wsk(ikx,iky,izbot1)    
!               else
!                  rhs(ikx,iky,1) = (0.D0,0.D0)
!               endif
!            enddo
!         enddo

!      elseif(mype==(npe-1)) then

!         do iky=1,ikty
!            ky = kya(iky)
!            do ikx=1,iktx
!               kx = kxa(ikx)
!               if (L(ikx,iky).eq.1) then
!                  rhs(ikx,iky,n3h0) =  i*kx*usk(ikx,iky,iztop1) + i*ky*vsk(ikx,iky,iztop1) + (0.5*r_3(iztop2) - 1./dz)*wsk(ikx,iky,iztop1-1)  
!               else
!                  rhs(ikx,iky,n3h0) = (0.D0,0.D0)
!               endif
!            enddo
!         enddo

!      end if

    END SUBROUTINE divergence


    SUBROUTINE gradient(phi,phi_x,phi_y,phi_z)

      !Computes the gradient of phi (n3h1) and stores it in phi_x,y,z (n3h0).
      !Not essential to encapsulate as a subroutine, but cleaner this way.

      double complex,   dimension(iktx,ikty,n3h1) :: phi                     !pressure, and rhs of pressure equation!                                  
      double complex,   dimension(iktx,ikty,n3h0) :: phi_x,phi_y,phi_z       !grad(phi)

      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               if (L(ikx,iky).eq.1) then
                  phi_x(ikx,iky,izh0) =  i*kx*phi(ikx,iky,izh1) 
                  phi_y(ikx,iky,izh0) =  i*ky*phi(ikx,iky,izh1) 
                  phi_z(ikx,iky,izh0) =  (phi(ikx,iky,izh1+1) - phi(ikx,iky,izh1))/dz 
               else
                  phi_x(ikx,iky,izh0) =  (0.D0,0.D0) 
                  phi_y(ikx,iky,izh0) =  (0.D0,0.D0)
                  phi_z(ikx,iky,izh0) =  (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo

      !phi_z( top ) is crap but it is not necessary to compute anyway...


      !Special treatment of boundaries: must set top boundary dphi/dz to 0
!      if(mype==(npe-1)) then

!         do iky=1,ikty
!            do ikx=1,iktx
               
!               phi_z(ikx,iky,n3h0) = (0.D0,0.D0)
               
!            enddo
!         enddo

!      end if

    END SUBROUTINE gradient





    SUBROUTINE  dissipation_q_nv(dqk,qok)

      ! This subroutine computes the VERTICAL dissipation of q at ts n-1 (lagged), assuming that just like u,v,psi, q is symmetric about the top and bot (dq/dz=0) so Q_N+1 = Q_N and Q_0 = Q_1.

      double complex, dimension(iktx,ikty,n3h1) :: qok
      double complex, dimension(iktx,ikty,n3h0) :: dqk


      do izh0=1,n3h0
         izh1=izh0+1
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  dqk(ikx,iky,izh0) =  nuz *(qok(ikx,iky,izh1+1) -2.*qok(ikx,iky,izh1) + qok(ikx,iky,izh1-1))/(dz*dz)  
               else
                  dqk(ikx,iky,izh0) = (0.D0,0.D0)
               endif
           enddo
        enddo
     enddo

     if(mype==0) then

         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  dqk(ikx,iky,1) =  nuz *( qok(ikx,iky,izbot1+1) - qok(ikx,iky,izbot1) )/(dz*dz)  
               else
                  dqk(ikx,iky,1) = (0.D0,0.D0)
               endif
           enddo
        enddo


     else if(mype==(npe-1)) then

         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  dqk(ikx,iky,n3h0) =  nuz*( qok(ikx,iky,iztop1-1) - qok(ikx,iky,iztop1) )/(dz*dz)  
               else
                  dqk(ikx,iky,n3h0) = (0.D0,0.D0)
               endif
           enddo
        enddo

     end if


   END SUBROUTINE dissipation_q_nv



    subroutine vort(uk,vk,wk,zxk,zyk,zzk)

      double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk
      double complex, dimension(iktx,ikty,n3h1) :: zxk,zyk,zzk

      !Nondim vort: Ar2 are added where w's appear
      !u,v stag and w unstag                                                                                                                                                                                                    
      !so that zx,zy unstag and zz stag                                                                                                                                                                                           

      do izh1=1,n3h1
         izh2=izh1+1
         
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
               
               if(L(ikx,iky)==1) then
                  
                  zxk(ikx,iky,izh1) = i*ky*Ar2*wk(ikx,iky,izh2) - (vk(ikx,iky,izh2+1)-vk(ikx,iky,izh2))/dz
                  zyk(ikx,iky,izh1) = (uk(ikx,iky,izh2+1)-uk(ikx,iky,izh2))/dz - i*kx*Ar2*wk(ikx,iky,izh2)
                  zzk(ikx,iky,izh1) = i*kx*vk(ikx,iky,izh2) - i*ky*uk(ikx,iky,izh2)
                  
               else
                  
                  zxk(ikx,iky,izh1) = (0.D0,0.D0)
                  zyk(ikx,iky,izh1) = (0.D0,0.D0)
                  zzk(ikx,iky,izh1) = (0.D0,0.D0)
                  
               end if
               
            enddo
         enddo
      end do

      !Boundary conditions. Might want to double-check, written after bike trip.

      if(mype==(npe-1)) then
         
         do iky=1,ikty
            do ikx=1,iktx
               
                  zxk(ikx,iky,iztop1) = (0.D0,0.D0)
                  zyk(ikx,iky,iztop1) = (0.D0,0.D0)
        
            enddo
         enddo

      end if


    end subroutine vort

    SUBROUTINE gradient2(tk,tk_x,tk_y,tk_z)

      !Would be nice to generalize gradient subroutine. This is used for PV calculation only for now.                                                                                                                                
      !Computes the gradient of tk (n3h2) and stores it in tk_x,y,z (n3h1).                                                                                                                                                          
      !Not essential to encapsulate as a subroutine, but cleaner this way.                                                                                                                                                         
      !tk is defined on UNstag grid points i*delta z with i=1,2,...,N and is =0 at the bounds.                                                                                                                                     
      !Its 3rd component (ONLY) is defined on STAG grid (i-1/2)delta z with again i=1,...N.  

      double complex,   dimension(iktx,ikty,n3h2) :: tk                    !pressure, and rhs of pressure equation!                                                                                                                    
      double complex,   dimension(iktx,ikty,n3h1) :: tk_x,tk_y,tk_z       !grad(phi)                                                                                                                                                     

      do izh1=1,n3h1
          izh2=izh1+1
         do iky=1,ikty
            ky = kya(iky)
            do ikx=1,iktx
               kx = kxa(ikx)
                      if (L(ikx,iky).eq.1) then
                  tk_x(ikx,iky,izh1) =  i*kx*tk(ikx,iky,izh2)
                  tk_y(ikx,iky,izh1) =  i*ky*tk(ikx,iky,izh2)
                  tk_z(ikx,iky,izh1) =  (tk(ikx,iky,izh2) - tk(ikx,iky,izh2-1))/dz        !Different than gradient (other subroutine)                                                                                                  
               else
                  tk_x(ikx,iky,izh1) =  (0.D0,0.D0)
                  tk_y(ikx,iky,izh1) =  (0.D0,0.D0)
                  tk_z(ikx,iky,izh1) =  (0.D0,0.D0)
                         endif
            enddo
         enddo
      enddo

      !Special treatment of boundaries                                                                                                                                                                                                   

      if(mype==0) then

         do iky=1,ikty
            do ikx=1,iktx

               tk_z(ikx,iky,izbot1) = tk(ikx,iky,izbot2) / dz        !Different than gradient (other subroutine)                                                                                                                             

            enddo
         enddo

      end if

    END SUBROUTINE gradient2



    SUBROUTINE compute_rot(uk,vk,wk,bk,wak,psik,u_a,v_a,w_a,b_a)

       double complex, dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk    
       double complex, dimension(iktx,ikty,n3h1) :: u_a,v_a,w_a,b_a
       double complex, dimension(iktx,ikty,n3h1) :: psik
       double complex, dimension(iktx,ikty,n3h1) :: u_rot,v_rot,b_rot,wak

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

       call generate_halo_q(b_rot)


       !Now, compute the ageostrophic part:  

       ! u = u_0 + Ro u_1 + ...
       ! w = 0   + Ro w_1 + ...
       ! b = b_0 + Ro b_1 + ...

       ! w_1 is wak, the vertical velocity retrieved from the omega equation. We could define wg=Ro*wak if wanted.
       ! u_0 and b_0 are the geostrophic components of uk and bk. Define the leftover as ageotrophic.

      do izh1=1,n3h1
         izh2=izh1+1
         do iky=1,ikty
            do ikx=1,iktx
               if (L(ikx,iky).eq.1) then
                  u_a(ikx,iky,izh1) = uk(ikx,iky,izh2) - u_rot(ikx,iky,izh1) 
                  v_a(ikx,iky,izh1) = vk(ikx,iky,izh2) - v_rot(ikx,iky,izh1) 
                  w_a(ikx,iky,izh1) = wk(ikx,iky,izh2) - Ro*wak(ikx,iky,izh1) 
                  b_a(ikx,iky,izh1) = bk(ikx,iky,izh2) - b_rot(ikx,iky,izh1)
               else 
                  u_a(ikx,iky,izh1) =  (0.D0,0.D0)
                  v_a(ikx,iky,izh1) =  (0.D0,0.D0)
                  w_a(ikx,iky,izh1) =  (0.D0,0.D0)
                  b_a(ikx,iky,izh1) =  (0.D0,0.D0)
               endif
            enddo
         enddo
      enddo


    end SUBROUTINE compute_rot
    END MODULE derivatives
