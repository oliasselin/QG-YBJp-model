PROGRAM main

  USE parameters
  USE mpi
  USE fft
  USE init
  USE derivatives
  USE elliptic
  USE diagnostics
  USE files
  USE IO_ncf

  !********************** Declaring variables *****************************!

  double precision, dimension(n1d,n2d,n3h2)   :: ur,vr,wr,br                     !Velocity and potential temperature fields (r-space)
  double complex,   dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk                     !Velocity and potential temperature fields (k-space)

  double precision, dimension(n1d,n2d,n3h1)   :: war       !1st order vertical velocity as computed in QG (w^1)
  double complex,   dimension(iktx,ikty,n3h1) :: wak       

  double precision, dimension(n1d,n2d,n3h1)   :: qr         !QGPV          
  double complex,   dimension(iktx,ikty,n3h1) :: qk        
  double complex,   dimension(iktx,ikty,n3h1) :: qok        !QGPV at previous time step 
  double complex,   dimension(iktx,ikty,n3h1) :: qtempk     !Temporary storage for filtering   

  !**** B = L+A, and both A and B are decomposed into their real and imag parts (ex.: A = AR + iAI)
  double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk, ARk, AIk
  double precision, dimension(n1d,n2d,n3h0)   :: BRr, BIr, ARr, AIr
  double complex,   dimension(iktx,ikty,n3h0) :: BRok, BIok            !B at the old time step
  double complex,   dimension(iktx,ikty,n3h0) :: BRtempk, BItempk      !B before filering

  !**** C = Az and is decomposed into real and imag parts (ex.: C = CR + iCI) even though in Fourier-space both CRk and CIk are complex. // Only used for diagnostics.
  double complex,   dimension(iktx,ikty,n3h0) :: CRk, CIk          
  double precision, dimension(n1d,n2d,n3h0)   :: CRr, CIr

  !**** n = nonlinear advection term J(psi,B) **** r = refractive term ~ B*vort
  double complex,   dimension(iktx,ikty,n3h0) :: nBRk, nBIk, rBRk, rBIk
  double precision, dimension(n1d,n2d,n3h0)   :: nBRr, nBIr, rBRr, rBIr

  !**** qw, the wave-averaged feedback onto QGPV ****!
  double complex,   dimension(iktx,ikty,n3h0) :: qwk
  double precision, dimension(n1d,n2d,n3h0)   :: qwr

  double complex,   dimension(iktx,ikty,n3h0) :: dqk         !dissipation
  double complex,   dimension(iktx,ikty,n3h1) :: psik        !QG streamfunction!
  double precision, dimension(n1d,n2d,n3h1)   :: psir
  double complex,   dimension(iktx,ikty,n3h1) :: psi_old     !For computing w...

  double complex,   dimension(iktx,ikty,n3h0) :: rhs         !RHS of elliptic equation (n3h0 version of q at n+1)
  double precision, dimension(n1d,n2d,n3h0)   :: rhsr

  double complex, dimension(iktx,n3, iktyp) :: qt            !Transposed (ky-parallelization) right-hand side   
  double complex, dimension(iktx,n3, iktyp) :: BRkt          !Transposed (ky-parallelization) BRk (this array can most likely be recycled)
  double complex, dimension(iktx,n3, iktyp) :: BIkt          !Transposed (ky-parallelization) BIk (this array can most likely be recycled)

  double precision, dimension(n1d,n2d,n3h0)   :: nqr         !Nonlinear q advection, J(psi,q)         
  double complex,   dimension(iktx,ikty,n3h0) :: nqk        

  double complex, dimension(iktx,ikty,2) :: sigma            !Vertial integral of A(kx,ky), 1=real part, 2=imag part. Only needed for regular YBJ (ybj_plus==0)

  !Integrating factor to account for horizontal hyperdiffusion. One for the flow, one for the waves (w)                                                                                                                                                   
  double precision :: int_factor,int_factor_w

  !For in-place Fourier transforms: Xr and Xk occupy the same storage space.
  equivalence(ur,uk)
  equivalence(vr,vk)
  equivalence(wr,wk)
  equivalence(br,bk)
  equivalence(war,wak)
  equivalence(rhsr,rhs)
  equivalence(psir,psik)
  equivalence(qr,qk)
  equivalence(qwr,qwk)
  equivalence(nqr,nqk)
  equivalence(BRr,BRk)
  equivalence(BIr,BIk)
  equivalence(ARr,ARk)
  equivalence(AIr,AIk)
  equivalence(CRr,CRk)
  equivalence(CIr,CIk)
  equivalence(nBRr,nBRk)
  equivalence(nBIr,nBIk)
  equivalence(rBRr,rBRk)
  equivalence(rBIr,rBIk)

  !FFT initialization dummy arrays
  double precision, dimension(n1d,n2d) :: array2dr
  double complex,   dimension(iktx,ikty) :: array2di
  double precision, dimension(n3)   :: fr_even,fk_even
  double precision, dimension(n3-1) :: fr_odd ,fk_odd
  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )
  equivalence(array2dr,array2di)

  !Additional arrays for diagnostics only!
  double complex, dimension(iktx,ikty,n3h1) :: u_rot     !Rotational part of u for slice...                                                                                                                                                                             
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr
  equivalence(u_rotr,u_rot)

  !For the Eulerian frequency diagnostic only
  double complex,   dimension(iktx,ikty,n3h0) :: dBRk, dBIk
  double precision, dimension(n1d,n2d,n3h0)   :: dBRr, dBIr
  equivalence(dBRr,dBRk)
  equivalence(dBIr,dBIk)

  !For YBJ terms magnitude only                                                                                                                            
  double complex,   dimension(iktx,ikty,n3h0) :: dARk, dAIk
  double precision, dimension(n1d,n2d,n3h0)   :: dARr, dAIr
  equivalence(dARr,dARk)
  equivalence(dAIr,dAIk)

  
  !********************!
  !*** Initializing ***!
  !********************!

  call initialize_mpi                                                   !Initialize parallel computation stuff
  call init_files                                                       !Define necessary (text) files for various diagnostics
  call initialize_fftw(array2dr,array2di,fr_even,fk_even,fr_odd,fk_odd) !Set up Fourier transform
  call init_arrays                                                      !Initialize all r-space, k-space arrays and the dealiasing matrix L
  call init_base_state                                                  !Set stratification and associated coefficients
  if(mype==0)  call validate_run                                        !Make sure all is good to run

  !Initialize physical fields: by default manually-set fields are used.
  call generate_fields_stag(psir,n3h1,ARr,n3h0,BRr,n3h0)      !Manually initialize 3 real-space fields with analytical formulae. Make modifications in 'init.f90'
  call fft_r2c(psir,psik,n3h1)                                !Transform to k-space all these manually-set r-space fields
  call fft_r2c(ARr,ARk,n3h0)                                  !Transform to k-space all these manually-set r-space fields
  call fft_r2c(BRr,BRk,n3h0)                                  !Transform to k-space all these manually-set r-space fields
  call sumB(BRk, BIk)                                         !Call only if the wave initial condition is horizontally uniform: this routine forces int(L+A)=0. Shouldn't be performed if A has horizontal structure

  !Initialize other fields to zero.   
  AIk = (0.D0,0.D0)                                             
  BIk = (0.D0,0.D0)
  CRk = (0.D0,0.D0)
  CIk = (0.D0,0.D0)

  !If desired, overwrite the manual initialization of psi and B by reading netCDF inputs.
  if(init_ncf_psi==1) call ncread_psi(psik,psir) 
  if(init_ncf_la ==1) call ncread_la(BRk,BRr,BIk,BIr)

  !Set q from psi, or set to 0 in the case of fixed flow (such that q = 0 for all times)
  if(fixed_flow == 0) then
     call init_q(qk,psik)
  else
     qk = (0.D0,0.D0)
  end if

  !Compute velocities and set necessary halos
  call compute_velo(uk,vk,wk,bk,psik) 
  call generate_halo(uk,vk,wk,bk)
  call generate_halo_q(qk) 

  !Dump initial condition with NetCDF
  if(out_psi ==1) call ncdump_psi(psik,psir,dump_count_psi)
  if(out_la  ==1) call ncdump_la(BRk,BRr,BIk,BIr,dump_count_la)
  if(out_n2  ==1) call ncdump_n2


 !***************************!
 !*** Initial diagnostics ***!
 !***************************!

  !For Eulerian Freq onlt
  dBRk = (0.D0,0.D0)
  dBIk = (0.D0,0.D0)

 !Compute war/wak if desired                                                                                                                                                     
 if(out_omega==1)  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
 end if

 if(out_etot ==1) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)

 do id_field=1,nfields                                            
    if(out_slice ==1)  call slices(ARk,AIK,ARr,AIr,BRk,BIk,BRr,BIr,CRk,CIk,CRr,CIr,dBRk,dBIk,dBRr,dBIr,id_field)
 end do

 do id_field=1,nfields2                                            
    if(out_slice2==1)  call slices2(uk,vk,wak,bk,psik,ur,vr,war,br,psir,id_field)
 end do
 
 do iz=1,num_spec
    if(out_hspecw ==1) call hspec_waves(BRk,BIk,CRk,CIk,iz)
 end do

 if(out_we   ==1) call wave_energy(ARk,AIk,BRk,BIk,CRk,CIk)
 if(out_wvave==1) call we_vave(BRk,BIk,BRr,BIr)
 if(out_conv ==1) call we_conversion(ARk, AIk, nBRk, nBIk, rBRk, rBIk, nBRr, nBIr, rBRr, rBIr)


 !************************************************************************!
 !*** 1st time timestep using the projection method with Forward Euler ***!
 !************************************************************************!
 
 time=delt
 if(itermax>0) then
 iter=1

 !Compute refractive and wave advective terms 
 call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
 call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)

 !Compute dissipation 
 call dissipation_q_nv(dqk,qok)

 !Set fields to zero in special cases
 if(inviscid==1) then  
    dqk=(0.D0,0.D0)
 end if

 if(linear==1) then
    nqk=(0.D0,0.D0)
   nBRk=(0.D0,0.D0)
   nBIk=(0.D0,0.D0)
 end if

 if(no_dispersion==1) then
    ARk=(0.D0,0.D0)
    AIk=(0.D0,0.D0)
 end if

 if(passive_scalar==1) then
    ARk = (0.D0,0.D0)
    AIk = (0.D0,0.D0)
   rBRk = (0.D0,0.D0)
   rBIk = (0.D0,0.D0)
end if

 !Store initial time step
  qok = qk
 BRok = BRk
 BIok = BIk

 !Compute q^1 and B^1 with Forward Euler  
 do izh0=1,n3h0
    izh1=izh0+1
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          kh2=kx*kx+ky*ky

          !Integrating factor for horizontal diffusion                                                                                                                                                                                                                                                                       
          int_factor   = delt* ( nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 )) )
          int_factor_w = delt* ( nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w)) )

          if (L(ikx,iky).eq.1) then

             qk(ikx,iky,izh1) = (  qok(ikx,iky,izh1) - delt* nqk(ikx,iky,izh0)  + delt*dqk(ikx,iky,izh0) )*exp(-int_factor)
            BRk(ikx,iky,izh0) = ( BRok(ikx,iky,izh0) - delt*nBRk(ikx,iky,izh0)  - delt*(0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) + delt*0.5*rBIk(ikx,iky,izh0) )*exp(-int_factor_w)
            BIk(ikx,iky,izh0) = ( BIok(ikx,iky,izh0) - delt*nBIk(ikx,iky,izh0)  + delt*(0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) - delt*0.5*rBRk(ikx,iky,izh0) )*exp(-int_factor_w)

          else

             qk(ikx,iky,izh1) = (0.D0,0.D0)
            BRk(ikx,iky,izh0) = (0.D0,0.D0)
            BIk(ikx,iky,izh0) = (0.D0,0.D0)

          endif
       enddo
    enddo
 enddo


 !Generate halo for q
 call generate_halo_q(qk)

!Recover the new streamfunction if the flow is evolving!
if(fixed_flow==0) then

   if(no_feedback == 1) then
      qwk = (0.D0,0.D0)
   else
      call compute_qw(qwk,BRk,BIk,qwr,BRr,BIr)           ! Compute qw                                                                                                                     
   end if

   do izh0=1,n3h0                                     ! Compute q* = q - qw
      izh1=izh0+1
      do iky=1,ikty
         do ikx=1,iktx
            if (L(ikx,iky).eq.1) then
               qwk(ikx,iky,izh0)=  qk(ikx,iky,izh1) - qwk(ikx,iky,izh0)
            endif
         enddo
      enddo
   enddo

   call mpitranspose(qwk,iktx,ikty,n3h0,qt,n3,iktyp)  !Transpose q*                                                                                      
   call psi_solver(psik,qt)                           !Solve the QGPV equation L(phi)=q*, assuming psi_z = 0 at top/bot (homogeneous problem)                    

end if

!Recover A from B unless in passive scalar mode (for which dispersion = 0)
if(passive_scalar==0) then

   if(ybj_plus==0) call sumB(BRk,BIk)                            !Resets the vertical sum of B to zero

   call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space 
   call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space 

   if(ybj_plus==1) then  !YBJ+ case
      call A_solver_ybj_plus(ARk,BRkt,CRk)
      call A_solver_ybj_plus(AIk,BIkt,CIk)
   else   !Normal YBJ solver 
      call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A
      call compute_A(ARk,AIK,BRkt,BIkt,CRk,CIK,sigma)               !Compute A!
   end if

end if

!Compute the corresponding u,v,w and t (u and v to be used in convol)                                    
call compute_velo(uk,vk,wk,bk,psik)
call generate_halo(uk,vk,wk,bk)

end if



 !********************************************************!
 !*** Subsequent timesteps using leapfrog timestepping ***!
 !********************************************************!


do iter=2,itermax
   
   time=iter*delt

   !Compute refractive and wave advective terms
   call convol_waqg(nqk,nBRk,nBIk,nqr,nBRr,nBIr,uk,vk,qk,BRk,BIk,ur,vr,qr,BRr,BIr)
   call refraction_waqg(rBRk,rBIk,rBRr,rBIr,BRk,BIk,psik,BRr,BIr,psir)

   !Compute dissipation from qok
   call dissipation_q_nv(dqk,qok)

   !Set fields to zero in special cases
   if(inviscid==1) then
      dqk=(0.D0,0.D0)
   end if
   
   if(linear==1) then
      nqk=(0.D0,0.D0)
      nBRk=(0.D0,0.D0)
      nBIk=(0.D0,0.D0)
   end if
   
   if(no_dispersion==1) then
      ARk=(0.D0,0.D0)
      AIk=(0.D0,0.D0)
   end if
   
   if(passive_scalar==1) then
      ARk = (0.D0,0.D0)
      AIk = (0.D0,0.D0)
      rBRk = (0.D0,0.D0)
      rBIk = (0.D0,0.D0)
   end if
   

     !Compute q^n+1 and B^n+1 using leap-frog
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky

              !Integrating factor for horizontal diffusion                                                                                                                                                                                                                                                             
              int_factor   = delt* ( nuh1 *((1.*kx)**(2.*ilap1 ) + (1.*ky)**(2.*ilap1 )) + nuh2 *((1.*kx)**(2.*ilap2 ) + (1.*ky)**(2.*ilap2 )) )
              int_factor_w = delt* ( nuh1w*((1.*kx)**(2.*ilap1w) + (1.*ky)**(2.*ilap1w)) + nuh2w*((1.*kx)**(2.*ilap2w) + (1.*ky)**(2.*ilap2w)) )

              if (L(ikx,iky).eq.1) then
                 qtempk(ikx,iky,izh1) =  qok(ikx,iky,izh1)*exp(-2*int_factor) - 2*delt*nqk(ikx,iky,izh0)*exp(-int_factor)  + 2*delt*dqk(ikx,iky,izh0)*exp(-2*int_factor)
                BRtempk(ikx,iky,izh0) = BRok(ikx,iky,izh0)*exp(-2*int_factor_w) - 2*delt*(nBRk(ikx,iky,izh0) + (0.5/(Bu*Ro))*kh2*AIk(ikx,iky,izh0) - 0.5*rBIk(ikx,iky,izh0) )*exp(-int_factor_w)
                BItempk(ikx,iky,izh0) = BIok(ikx,iky,izh0)*exp(-2*int_factor_w) - 2*delt*(nBIk(ikx,iky,izh0) - (0.5/(Bu*Ro))*kh2*ARk(ikx,iky,izh0) + 0.5*rBRk(ikx,iky,izh0) )*exp(-int_factor_w)
              else
                 qtempk(ikx,iky,izh1) = (0.D0,0.D0)
                BRtempk(ikx,iky,izh0) = (0.D0,0.D0)
                BItempk(ikx,iky,izh0) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo

     !For diagnostic of Eulerian Frequency only!
     !d(LA)/dt = [LA^(n+1)-LA^(n-1)]/2dt
     dBRk=(BRtempk-BRok)/(2.*delt)
     dBIk=(BItempk-BIok)/(2.*delt)
     !*****************************************!


     !Apply Robert-Asselin filter to damp the leap-frog computational mode
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           do ikx=1,iktx
              if (L(ikx,iky).eq.1) then
                 qok(ikx,iky,izh1) =  qk(ikx,iky,izh1) + gamma * (  qok(ikx,iky,izh1) - 2 *  qk(ikx,iky,izh1) +  qtempk(ikx,iky,izh1) )
                BRok(ikx,iky,izh0) = BRk(ikx,iky,izh0) + gamma * ( BRok(ikx,iky,izh0) - 2 * BRk(ikx,iky,izh0) + BRtempk(ikx,iky,izh0) )
                BIok(ikx,iky,izh0) = BIk(ikx,iky,izh0) + gamma * ( BIok(ikx,iky,izh0) - 2 * BIk(ikx,iky,izh0) + BItempk(ikx,iky,izh0) )
              else
                 qok(ikx,iky,izh1) = (0.D0,0.D0)
                BRok(ikx,iky,izh0) = (0.D0,0.D0)
                BIok(ikx,iky,izh0) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo

!Overwrite the new field uk with u^{n+1} 
 qk =  qtempk
BRk = BRtempk
BIk = BItempk

 !Generate halo for q
 call generate_halo_q(qk)
 call generate_halo_q(qok)
 
!Recover the new streamfunction if the flow is evolving!  
if(fixed_flow==0) then

   if(no_feedback == 1) then
      qwk = (0.D0,0.D0)
   else
      call compute_qw(qwk,BRk,BIk,qwr,BRr,BIr)           ! Compute qw                                                                                                                     
   end if
   
   do izh0=1,n3h0                                     ! Compute q* = q - qw                                                                                 
      izh1=izh0+1
      do iky=1,ikty
         do ikx=1,iktx
            if (L(ikx,iky).eq.1) then
               qwk(ikx,iky,izh0)=  qk(ikx,iky,izh1) - qwk(ikx,iky,izh0)
            endif
         enddo
      enddo
   enddo
   
   call mpitranspose(qwk,iktx,ikty,n3h0,qt,n3,iktyp)  !Transpose rhs -> ft                                                                            
   call psi_solver(psik,qt)                           !Solve the pressure equation laplacian(phi)=f                                                              

end if

!Recover A from B unless in passive scalar mode (for which dispersion = 0)
if(passive_scalar==0) then

   if(ybj_plus==0) call sumB(BRk,BIk)                           !Resets the vertical sum of B to zero

   call mpitranspose(BRk,iktx,ikty,n3h0,BRkt,n3,iktyp)           !Transpose BR to iky-parallelized space 
   call mpitranspose(BIk,iktx,ikty,n3h0,BIkt,n3,iktyp)           !Transpose BK to iky-parallelized space 

   if(ybj_plus==1) then  !YBJ+ case
      call A_solver_ybj_plus(ARk,BRkt,CRk)
      call A_solver_ybj_plus(AIk,BIkt,CIk)
   else   !Normal YBJ solver 
      call compute_sigma(sigma,nBRk, nBIk, rBRk, rBIk)              !Compute the sum of A
      call compute_A(ARk,AIK,BRkt,BIkt,CRk,CIK,sigma)               !Compute A!
   end if
   
end if


 !Compute the corresponding u,v,w and t 
 call compute_velo(uk,vk,wk,bk,psik)
 call generate_halo(uk,vk,wk,bk) 

!Dump r-space output!
if(out_psi ==1 .and. mod(iter,freq_psi)==0) call ncdump_psi(psik,psir,dump_count_psi)
if(out_la  ==1 .and. mod(iter,freq_la )==0) call ncdump_la(BRk,BRr,BIk,BIr,dump_count_la)


 !*******************!
 !*** Diagnostics ***!
 !*******************!

 !Compute w if desired
 if(out_omega==1 .and. (mod(iter,freq_omega) ==0))  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
    call generate_halo_q(wak)
 end if
 
if(out_etot ==1 .and. mod(iter,freq_etot )==0) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)

do id_field=1,nfields
   if(out_slice ==1 .and. mod(iter,freq_slice)==0 .and. count_slice(id_field)<max_slices)  call slices(ARk,AIK,ARr,AIr,BRk,BIk,BRr,BIr,CRk,CIk,CRr,CIr,dBRk,dBIk,dBRr,dBIr,id_field)
end do

do id_field=1,nfields3
   if(out_slice3==1 .and. mod(iter,freq_slice3)==0 .and. count_slice3(id_field)<max_slices)  call slices3(ARk,AIK,ARr,AIr,dBRk,dBIk,dBRr,dBIr,nBRk,nBIk,nBRr,nBIr,rBRk,rBIk,rBRr,rBIr,id_field)
end do

do id_field=1,nfields2
   if(out_slice2==1 .and. mod(iter,freq_slice2)==0 .and. count_slice2(id_field)<max_slices)  call slices2(uk,vk,wak,bk,psik,ur,vr,war,br,psir,id_field)
end do
                                                                          
 do iz=1,num_spec
    if(out_hspecw ==1  .and. mod(iter,freq_hspecw)==0 ) call hspec_waves(BRk,BIk,CRk,CIk,iz)
 end do

 if(out_we ==1   .and. mod(iter,freq_we   )==0)  call wave_energy(ARk,AIk,BRk,BIk,CRk,CIk)
 if(out_wvave==1 .and. mod(iter,freq_wvave)==0)  call we_vave(BRk,BIk,BRr,BIr)
 if(out_conv ==1 .and. mod(iter,freq_conv )==0)  call we_conversion(ARk, AIk, nBRk, nBIk, rBRk, rBIk, nBRr, nBIr, rBRr, rBIr)
 if(out_gamma==1 .and. mod(iter,freq_gamma )==0) call gamma_conversion(ARk, AIk, BRk, BIk, nBRk, nBIk, rBRk, rBIk, nBRr, nBIr, rBRr, rBIr)
     

 !**************************************************************************!

 if(time>maxtime) EXIT
end do !End loop

!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
