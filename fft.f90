MODULE fft

  USE parameters

  IMPLICIT NONE

  
  
  integer, parameter :: fftw_forward         = -1, fftw_backward                  = 1
  integer, parameter :: fftw_real_to_complex = -1, fftw_complex_to_real           = 1
  integer, parameter :: fftw_estimate        =  0, fftw_measure=1,fftw_use_wisdom = 16
  integer, parameter :: fftw_in_place        =  8

  !FFT plans                                                                                                                                                                                                                      
  integer*8          :: plan_r2c
  integer*8          :: plan_c2r

  integer, parameter :: norm_rc = n1*n2


  !For z-spectra.
  integer*8 :: plan_even_z
  integer*8 :: plan_odd_z

  integer, parameter :: fftw_redft10         =  5, fftw_rodft00                   = 7

 ! https://code.google.com/p/fdtd3d/source/browse/trunk/fortran/fftw3.f


  CONTAINS

    SUBROUTINE initialize_fftw(array2dr,array2di,fr_even,fk_even,fr_odd,fk_odd)
      
      !2D horizontal transforms                                                                                                                                                                                        
      double precision, dimension(n1d,n2d) :: array2dr
      double complex,   dimension(iktx,ikty) :: array2di

      !For z-spectra
      double precision, dimension(n3)   :: fr_even,fk_even
      double precision, dimension(n3-1) :: fr_odd ,fk_odd


      !2D horizontal transforms                                                                                                                                                                                              
      call dfftw_plan_dft_r2c_2d(plan_r2c,n1,n2,array2dr,array2di,fftw_estimate+fftw_in_place)
      call dfftw_plan_dft_c2r_2d(plan_c2r,n1,n2,array2di,array2dr,fftw_estimate+fftw_in_place)

      !For z-spectra
      call dfftw_plan_r2r_1d(plan_even_z,n3  ,fr_even,fk_even,  FFTW_REDFT10,     fftw_estimate+fftw_in_place)
      call dfftw_plan_r2r_1d(plan_odd_z ,n3-1,fr_odd ,fk_odd ,  FFTW_RODFT00,     fftw_estimate+fftw_in_place)


    END SUBROUTINE initialize_fftw


    SUBROUTINE fft_r2c(rvar,kvar,nz)

        double precision :: rvar(n1d,n2d,nz)
        double complex   :: kvar(iktx,ikty,nz)
        integer :: nz

        do iz=1,nz
           call dfftw_execute_dft_r2c(plan_r2c,rvar(1,1,iz),kvar(1,1,iz))
        end do
        kvar=kvar/norm_rc

      END SUBROUTINE fft_r2c




    SUBROUTINE fft_c2r(kvar,rvar,nz)

        double precision :: rvar(n1d,n2d,nz)
        double complex   :: kvar(iktx,ikty,nz)
        integer :: nz

        if(dealiasing==1) call realit(kvar,nz)
        
        do iz=1,nz
           call dfftw_execute_dft_c2r(plan_c2r,kvar(1,1,iz),rvar(1,1,iz))
        end do

      END SUBROUTINE fft_c2r



    
    SUBROUTINE realit(kvar,nz) !This routine adds the right modes (satisfying the reality condition) in (0,ky-) so that the r-field is the right one.

      double complex, dimension(iktx,ikty,nz) :: kvar
      integer :: nz
      integer :: iz !dummy index
      integer :: inky !index for negative ky's

      !Modes (0,-ky) are just the conjugate of (0,ky+) modes
      !Recall kya=[ 0 1 2 ... N/2-1  0  -N/2+1 -N/2+2 ... -2 -1] 

      do iz=1,nz
         do iky=2,ikty/2 !n2 must be even, otherwise it won't work.
            inky=ikty-iky+2

            kvar(1,inky,iz)=dble(kvar(1,iky,iz)) - i*dimag(kvar(1,iky,iz))  ! u(0,ky-)=u*(0,ky+)

         end do
      end do

    END SUBROUTINE realit




    SUBROUTINE kill_fftw

        !*** r-space <=> h-space ***!                                                                                                                                                                                                        

        call dfftw_destroy_plan(plan_r2c)
        call dfftw_destroy_plan(plan_c2r)

    END SUBROUTINE kill_fftw

END MODULE fft
