module IO_ncf
  use parameters
  use mpi
  use fft
  use netcdf
  implicit none
  
  integer :: idin,idout !Files ID
  integer :: idx,idy,idz,idt !Dimension IDs
  integer :: idpsi,idlar,idlai,idtime    !Variable IDs
  integer :: idn2

CONTAINS

subroutine ncdump_psi(psik,psir,dump_count_psi)

  !Create netcdf file and write for the geostrophic streamfunction psi (output)
  implicit none

  double complex, dimension(iktx,ikty,n3h1) ::   psik
  double precision, dimension(n1d,n2d,n3h1)   :: psir

  double precision,    dimension(iktx,ikty,n3h1) :: psi_save  !Save psi to avoid inverse FFT
  real, dimension(n1,n2,n3h0) :: psi_clean                    !suboptimal: de-haloed, de-extra x-dim values version of psi, real (not double). Dimensional units (m^2/s)

  integer, dimension(3) :: ncdims
  integer, dimension(3) :: nccount,ncstart

  character *13 :: file_name                                  !Name of the netCDF file 
  real :: time_seconds                                        !Print time in seconds => nondimensional time * L_scale/U_scale
  
  integer :: dump_count_psi                                   !Count the number of psi files printed for numbering netCDF file

  !Calculate time in seconds!
  time_seconds = time*L_scale/U_scale

  !Move to r-space!
  psi_save = psik !Save to avoid inverse fft
  call fft_c2r(psik,psir,n3h1)

  !---------------------------------------------------!
  !Suboptimal: get rid of halos and extra x-dim values, and make dimensional!

  do izh0=1,n3h0
     izh1=izh0+1
     do iy=1,n2
        do ix=1,n1

           psi_clean(ix,iy,izh0) = real(psir(ix,iy,izh1)*U_scale*L_scale)

        end do
     end do
  end do
  !Recover the k-space psi without fft back
  psik=psi_save
  !---------------------------------------------------!

!---------------------------------------------------------------------
!     DEFINING THE OUTPUT NETCDF FILE 

!!! prep restart file (output) (define the variables and dims)
  !Set file name
  write(file_name,"(A3,I0.3,A7)") "psi",dump_count_psi,".out.nc"

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF,
  call check( nf90_create(file_name,ior(NF90_NETCDF4,NF90_MPIIO),idout, comm = MPI_COMM_WORLD,info = MPI_INFO_NULL) )

  ! Define the dimensions: x, y, z. time. 
  call check( nf90_def_dim(idout, "x",int(n1,4),idx) )
  call check( nf90_def_dim(idout, "y",int(n2,4),idy) )
  call check( nf90_def_dim(idout, "z",int(n3,4),idz) )
  call check( nf90_def_dim(idout, "t",        1,idt) )
     
  ! ncdimsk array used to pass the dimension IDs to the variables
  ncdims = (/idx,idy,idz/)
  
  ! Define the variables which are the fields that going to be stored
  call check( nf90_def_var(idout,"psi" ,NF90_FLOAT,ncdims,idpsi ) )
  call check( nf90_def_var(idout,"time",NF90_FLOAT,idt   ,idtime) )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(idout) )

!---------------------------------------------------------------------
  !     WRITING VARIABLES IN NETCDF FILE
  
  ! how many to count in each dimension when writing files
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0

  ! where to start on the output file
  ncstart(1) = 1
  ncstart(2) = 1
  
  ! write time (time is written only one time so just root process is used
  if (mype.eq.0)     call check( nf90_put_var(idout, idtime, time_seconds) )

  ! in the z-direction mype=0 is in the first place
  ncstart(3) = mype*n3h0+1
  call check( nf90_put_var(idout, idpsi, psi_clean,ncstart,nccount))    

  
!---------------------------------------------------------------------
!     CLOSING the  NETCDF FILE 
  call check (nf90_close(idout))

  dump_count_psi=  dump_count_psi + 1

end subroutine ncdump_psi





subroutine ncdump_la(BRk,BRr,BIk,BIr,dump_count_la)

  !Create netcdf file and write for the geostrophic streamfunction psi (output)

  implicit none

  double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk
  double precision, dimension(n1d,n2d,n3h0)   :: BRr, BIr

  double precision,    dimension(iktx,ikty,n3h0) :: BR_save,BI_save  !Save psi to avoid inverse FFT  
  real, dimension(n1,n2,n3h0) :: BR_clean, BI_clean                  !suboptimal: de-haloed, de-extra x-dim values version of psi, real (not double)

  integer, dimension(3) :: ncdims
  integer, dimension(3) :: nccount,ncstart

  character *12 :: file_name
  real :: time_seconds        !Print time in seconds => nondimensional time * L_scale/U_scale
  
  integer :: dump_count_la

  !Calculate time in seconds!
  time_seconds = time*L_scale/U_scale

  !Move to r-space!
  BR_save = BRk !Save to avoid inverse fft
  BI_save = BIk !Save to avoid inverse fft
  call fft_c2r(BRk,BRr,n3h0)
  call fft_c2r(BIk,BIr,n3h0)

  !---------------------------------------------------!
  !Suboptimal: get rid of halos and extra x-dim values and make dimensional!

  do izh0=1,n3h0
     do iy=1,n2
        do ix=1,n1

           BR_clean(ix,iy,izh0) = real(BRr(ix,iy,izh0)*Uw_scale)
           BI_clean(ix,iy,izh0) = real(BIr(ix,iy,izh0)*Uw_scale)

        end do
     end do
  end do
  !Recover the k-space psi without fft back
  BRk=BR_save
  BIk=BI_save
  !---------------------------------------------------!

!---------------------------------------------------------------------
!     DEFINING THE OUTPUT NETCDF FILE 

!!! prep restart file (output) (define the variables and dims)
  !Set file name
  write(file_name,"(A2,I0.3,A7)") "la",dump_count_la,".out.nc"

  ! Create the netCDF file. 
  call check( nf90_create(file_name,ior(NF90_NETCDF4,NF90_MPIIO),idout, comm = MPI_COMM_WORLD,info = MPI_INFO_NULL) )

  ! Define the dimension: kx, ky, kz. time. RI used as another dim to distinguish between real and imag parts
  call check( nf90_def_dim(idout, "x",int(n1,4),idx) )
  call check( nf90_def_dim(idout, "y",int(n2,4),idy) )
  call check( nf90_def_dim(idout, "z",int(n3,4),idz) )
  call check( nf90_def_dim(idout, "t",        1,idt) )
     
  ! ncdimsk array used to pass the dimension IDs to the variables
  ncdims = (/idx,idy,idz/)
  
  ! Define the variables which are the fields that going to be stored
  call check( nf90_def_var(idout,"LAr" ,NF90_FLOAT,ncdims,idlar ) )
  call check( nf90_def_var(idout,"LAi" ,NF90_FLOAT,ncdims,idlai ) )
  call check( nf90_def_var(idout,"time",NF90_FLOAT,idt   ,idtime) )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(idout) )

!---------------------------------------------------------------------
  !     WRITING VARIABLES IN NETCDF FILE
  
  ! how many to count in each dimension when writing files
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0

  ! where to start on the output file
  ncstart(1) = 1
  ncstart(2) = 1
  
  ! write time (time is written only one time so just root process is used
  if (mype.eq.0)     call check( nf90_put_var(idout, idtime, time_seconds) )

  ! in the z-direction mype=0 is in the first place
  ncstart(3) = mype*n3h0+1
  call check( nf90_put_var(idout, idlar, BR_clean,ncstart,nccount))    
  call check( nf90_put_var(idout, idlai, BI_clean,ncstart,nccount))    

  
!---------------------------------------------------------------------
!     CLOSING the  NETCDF FILE 
  call check (nf90_close(idout))

  dump_count_la= dump_count_la+1

end subroutine ncdump_la





subroutine ncread_psi(psik,psir)
! Read the netcdf data from an input file (psi.in.ncf)                                                                
  implicit none

  double complex, dimension(iktx,ikty,n3h1) ::   psik   !k-space, model-usable streamfunciton with halos
  double precision, dimension(n1d,n2d,n3h1)   :: psir   !r-space version

  real, dimension(n1,n2,n3h0) :: psi_clean              !Dimensional psi directly from the read file

  integer :: n1in,n2in,n3in                             !Dimensions read in the file. Must match n1,n2,n3 for now
  integer, dimension(3) :: nccount,ncstart              !Number of values to read for each dimension and where to start

  real :: time_seconds                                  !Time in seconds. For now we ignore this and start at time=0

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access                                                            
  print*,'Yo! reading from netcdf restart file for psi, BTW im mype',mype
  call check( nf90_open(init_ncf_psi_filename, NF90_NOWRITE,idin) )        

  ! Get the dimensions IDs based on their name                                                                                   
  call check( nf90_inq_dimid(idin, "x", idx) )
  call check( nf90_inq_dimid(idin, "y", idy) )
  call check( nf90_inq_dimid(idin, "z", idz) )
  call check( nf90_inq_dimid(idin, "t", idt) )

  ! Get the dimension length and check if the grid resolution matches                                                   
  call check( nf90_inquire_dimension(idin,idx,len=n1in))
  call check( nf90_inquire_dimension(idin,idy,len=n2in))
  call check( nf90_inquire_dimension(idin,idz,len=n3in))

  if (n1in.ne.n1 .or. n2in.ne.n2 .or. n3in.ne.n3) then
    print*,'Sorry, do not know how to change resolution.'
    stop
  endif

  ! Get the variables IDs                                                                                                                                                      
  call check( nf90_inq_varid(idin,"time",idtime))
  call check( nf90_inq_varid(idin, "psi",idpsi) )

  ! prepare for reading variables                                                                                                                                              
  ncstart    = 1
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0
  ncstart(3)= int(mype*n3h0+1)

  !Get variables
  call check(  nf90_get_var(idin,idtime,time_seconds))
  call check( nf90_get_var(idin,idpsi,psi_clean,ncstart,nccount))

  !Print time
  if (mype.eq.0) print*,"File read starts at t = ",time_seconds,"s"

  !Back to the streamfunction: get its real-space, non dimensional version first!
  psir = 0.D0
  do izh0=1,n3h0
     izh1=izh0+1
     do iy=1,n2
        do ix=1,n1

          psir(ix,iy,izh1) =  psi_clean(ix,iy,izh0)/(U_scale*L_scale) 

        end do
     end do
  end do
  call fft_r2c(psir,psik,n3h1)   !Move to k-space
  call generate_halo_q(psik)     !Generate 1-lvl halo
  !------------!

  call check (nf90_close(idin))

end subroutine ncread_psi




subroutine ncread_la(BRk,BRr,BIk,BIr)
! Read the netcdf data from an input file (psi.in.ncf)                                                                
  implicit none

  double complex,   dimension(iktx,ikty,n3h0) :: BRk, BIk            !k-space, model-usable L+A with halos
  double precision, dimension(n1d,n2d,n3h0)   :: BRr, BIr            !r-space version
  real, dimension(n1,n2,n3h0) :: BR_clean, BI_clean                  !suboptimal: de-haloed, de-extra x-dim values version of psi, real (not double)    

  integer :: n1in,n2in,n3in                             !Dimensions read in the file. Must match n1,n2,n3 for now
  integer, dimension(3) :: nccount,ncstart              !Number of values to read for each dimension and where to start

  real :: time_seconds                                  !Time in seconds. For now we ignore this and start at time=0

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access                                                            
  print*,'Yo! reading from netcdf restart file for LA, BTW im mype',mype
  call check( nf90_open(init_ncf_la_filename, NF90_NOWRITE,idin) )        

  ! Get the dimensions IDs based on their name                                                                                   
  call check( nf90_inq_dimid(idin, "x", idx) )
  call check( nf90_inq_dimid(idin, "y", idy) )
  call check( nf90_inq_dimid(idin, "z", idz) )
  call check( nf90_inq_dimid(idin, "t", idt) )

  ! Get the dimension length and check if the grid resolution matches                                                   
  call check( nf90_inquire_dimension(idin,idx,len=n1in))
  call check( nf90_inquire_dimension(idin,idy,len=n2in))
  call check( nf90_inquire_dimension(idin,idz,len=n3in))

  if (n1in.ne.n1 .or. n2in.ne.n2 .or. n3in.ne.n3) then
    print*,'Sorry, do not know how to change resolution.'
    stop
  endif

  ! Get the variables IDs                                                                                                                                                      
  call check( nf90_inq_varid(idin,"time",idtime))
  call check( nf90_inq_varid(idin, "LAr",idlar) )
  call check( nf90_inq_varid(idin, "LAi",idlai) )


  ! prepare for reading variables                                                                                                                                              
  ncstart    = 1
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0
  ncstart(3)= int(mype*n3h0+1)

  !Get variables
  call check( nf90_get_var(idin,idtime,time_seconds))
  call check( nf90_get_var(idin,idlar,BR_clean,ncstart,nccount))
  call check( nf90_get_var(idin,idlai,BI_clean,ncstart,nccount))

  !Print time
  if (mype.eq.0) print*,"File read starts at t = ",time_seconds,"s"

  !Back to the streamfunction: get its real-space, non dimensional version first!
  BRr = 0.D0
  BIr = 0.D0
  do izh0=1,n3h0
     do iy=1,n2
        do ix=1,n1

          BRr(ix,iy,izh0) =  BR_clean(ix,iy,izh0)/(Uw_scale) 
          BIr(ix,iy,izh0) =  BI_clean(ix,iy,izh0)/(Uw_scale) 

        end do
     end do
  end do
  call fft_r2c(BRr,BRk,n3h0)   !Move to k-space
  call fft_r2c(BIr,BIk,n3h0)   !Move to k-space
  !------------!

  call check (nf90_close(idin))

end subroutine ncread_la




subroutine ncdump_N2

  !Create netcdf file and write for the geostrophic streamfunction psi (output)
  implicit none

  integer, dimension(1) :: ncdims
  integer, dimension(1) :: nccount,ncstart

  character *9 :: file_name                                  !Name of the netCDF file 

!---------------------------------------------------------------------
!     DEFINING THE OUTPUT NETCDF FILE 

!!! prep restart file (output) (define the variables and dims)
  !Set file name
  write(file_name,"(A9)") "N2.out.nc"

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF,
  call check( nf90_create(file_name,ior(NF90_NETCDF4,NF90_MPIIO),idout, comm = MPI_COMM_WORLD,info = MPI_INFO_NULL) )

  ! Define the dimensions: x, y, z. time. 
  call check( nf90_def_dim(idout, "z",int(n3,4),idz) )
     
  ! ncdimsk array used to pass the dimension IDs to the variables
  ncdims = (/idz/)
  
  ! Define the variables which are the fields that going to be stored
  call check( nf90_def_var(idout,"N2" ,NF90_FLOAT,ncdims,idn2 ) )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(idout) )

!---------------------------------------------------------------------
  !     WRITING VARIABLES IN NETCDF FILE
  
  ! how many to count in each dimension when writing files
  nccount(1) = n3

  ! where to start on the output file
  ncstart(1) = 1
  
  ! write time (time is written only one time so just root process is used
!  if (mype.eq.0)     call check( nf90_put_var(idout, idtime, time_seconds) )

  ! in the z-direction mype=0 is in the first place
!  ncstart(3) = mype*n3h0+1
!  call check( nf90_put_var(idout, idpsi, psi_clean,ncstart,nccount))    

  if(mype .eq. 0)  call check( nf90_put_var(idout, idn2, r_2st*N0*N0,ncstart,nccount))    
  
!---------------------------------------------------------------------
!     CLOSING the  NETCDF FILE 
  call check (nf90_close(idout))



end subroutine ncdump_N2



subroutine check(status)
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop "Stopped!!!"
  end if
end subroutine check



end module IO_ncf
