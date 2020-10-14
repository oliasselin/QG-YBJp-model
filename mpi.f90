MODULE mpi

  USE parameters

  IMPLICIT NONE

  include 'mpif.h'

  integer :: mype,size,ierror
  integer :: status(MPI_STATUS_SIZE) !for mpitranspose...

  integer :: tag_nlt = 123
  integer :: tag_kzs = 124
  integer :: tag_kzu = 125
  integer :: tag_pzs = 126
  integer :: tag_dz  = 127

  integer :: tag_kzs2 = 128
  integer :: tag_kzu2 = 129
  integer :: tag_pzs2 = 130

  integer :: tag_r1  = 131
  integer :: tag_r2  = 132
  integer :: tag_p0  = 133

  integer :: tag_slice_xz(nfields) = [141,142,143,144,145,146,147,148]  !Up to 14X where X is the number of fields to slice
  integer :: tag_slice_xz2(nfields2) = [1141,1142,1143,1144,1145,1146,1147,1148,1149]  !Up to 114X where X is the number of fields to slice
  integer :: tag_slice_xz3(nfields3) = [11141,11142,11143,11144,11145,11146,11147,11148]  !Up to 114X where X is the number of fields to slice

  integer :: tag_nzs = 151
  integer :: tag_nzu = 152
  integer :: tag_lzs = 153

  integer :: tag_ezs = 154
  integer :: tag_ezu = 155

  integer :: tag_rzs = 156
  integer :: tag_rzu = 157
  integer :: tag_hb  = 158
  integer :: tag_hw  = 159

  integer :: tag_hb2 = 1582
  integer :: tag_hw2 = 1592

  integer :: tag_ns = 160
  integer :: tag_us = 161
  integer :: tag_vs = 162
  integer :: tag_ws = 163
  integer :: tag_ts = 164


  integer :: tag_cond =165
  integer :: tag_condwz=166

  integer :: tag_psir=167
  integer :: tag_psik=168

  integer :: tag_cont = 169

  integer :: tag_eta  = 170
  integer :: tag_eta2 = 171

  integer :: tag_grow = 172

  CONTAINS

    SUBROUTINE initialize_mpi

         call mpi_init(ierror)
         call mpi_comm_size(mpi_comm_world,size,ierror)
         call mpi_comm_rank(mpi_comm_world,mype,ierror)
!         if (mype.eq.0) print*,'size = ',size
         if (size.ne.npe) then
         print*,size,npe,'Wrong number of processors!'
         stop
         end if

    END SUBROUTINE initialize_mpi


    SUBROUTINE kill_mpi

        call mpi_finalize(ierror)

    END SUBROUTINE kill_mpi

    SUBROUTINE mpitranspose(x,d1,d2,d3p,xt,d3,d2p)

      integer :: d1,d2,d3,d2p,d3p  !I change the names to avoid confusion...                                                                                                                                                               
      integer :: incr,j,k,ip,ito,ifrom,jblock,kblock,koff,n123p

      double complex, dimension(d1,d2,d3p), intent(in)  :: x
      double complex, dimension(d1,d3,d2p), intent(out) :: xt
      double complex, dimension(d1,d2p,d3p,npe) :: xsend,xrecv


      n123p=d1*d2p*d3p

      do ip=1,npe
         do k=1,d3p
            do j=1,d2p
               jblock=j+(ip-1)*d2p
                      do incr=1,d1
                  xsend(incr,j,k,ip) = x(incr,jblock,k)
                         enddo
            enddo
         enddo
      enddo

  call mpi_alltoall(xsend,n123p,MPI_DOUBLE_COMPLEX,xrecv,n123p,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierror)

  do ip=1,npe
     do j=1,d2p
        do k=1,d3p
           do incr=1,d1
              xt(incr,k+(ip-1)*d3p,j) = xrecv(incr,j,k,ip)
           enddo
        enddo
     enddo
  enddo


    END SUBROUTINE  mpitranspose



    SUBROUTINE generate_halo(uk,vk,wk,tk) 

      !This subroutine distributes the data to create the proper halo for velocity fields. On the input, fields are also n3h2, but the values in the halo (first and last 2 levels) are insignificant. 
      !The way it proceeds is that it first save the 2 first and last significant levels for each processor on a chunk (chun_send_..._...). This chunk will become the previous and next processor's new top and bottom halo.
      !So first, the even taged (mype 0, 2, 4 ...) processor send both of their chunk (execpt for mype 0 and npe-1 that respectively just send their top and bottom) to the right (previous or next) processor.
      !Then, we do the inverse: odd processors (mype = 1, 3, ...) send their chunk to the even ones.
      !This subroutine is based on the fact that the number of processors is even and that n3h0 is at least 2, and must be modified otherwise.
      

      integer :: prev,next,tag
      integer :: nrec,maxrec

      integer :: ikx,iky,incrz

      integer :: chunksize
      integer, parameter :: tag_u=1,tag_v=2,tag_w=3,tag_t=4

      double complex, dimension(iktx,ikty,2) :: chunk_send_bot_u,chunk_send_top_u,chunk_send_bot_v,chunk_send_top_v,chunk_send_bot_w,chunk_send_top_w,chunk_send_bot_t,chunk_send_top_t
      double complex, dimension(iktx,ikty,2) :: chunk_recv_u,chunk_recv_v,chunk_recv_w,chunk_recv_t
      double complex, dimension(iktx,ikty,n3h2), intent(inout) :: uk,vk,wk,tk

      !maxsend is the number of processors from which mype receives velocity fields

      if(mype==0) then
         maxrec=1
      elseif(mype==(npe-1)) then
         maxrec=1
      else
         maxrec=2
      end if

      !Send to which processor                                                                                                                                                                                                               
      prev=mype-1
      next=mype+1

      !Create the chunks to send                                                                                                                                                                                                                                  

      do ikx=1,iktx
         do iky=1,ikty
            do incrz=1,2
 
               if(mype/=0) then
               chunk_send_bot_u(ikx,iky,3-incrz)=uk(ikx,iky,2+incrz)
               chunk_send_bot_v(ikx,iky,3-incrz)=vk(ikx,iky,2+incrz)
               chunk_send_bot_w(ikx,iky,3-incrz)=wk(ikx,iky,2+incrz)
               chunk_send_bot_t(ikx,iky,3-incrz)=tk(ikx,iky,2+incrz)
               end if

               if(mype/=(npe-1)) then
               chunk_send_top_u(ikx,iky,3-incrz)=uk(ikx,iky,n3h0+3-incrz)
               chunk_send_top_v(ikx,iky,3-incrz)=vk(ikx,iky,n3h0+3-incrz)
               chunk_send_top_w(ikx,iky,3-incrz)=wk(ikx,iky,n3h0+3-incrz)
               chunk_send_top_t(ikx,iky,3-incrz)=tk(ikx,iky,n3h0+3-incrz)
               end if

            end do
         end do
      end do

      chunksize=iktx*ikty*2

      !Even processors send their chunks

      if(mod( (mype+2),2)==0) then

         if(mype/=0) then
         call mpi_send(chunk_send_bot_u,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_u,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_v,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_v,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_w,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_w,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_t,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_t,MPI_COMM_WORLD,ierror)
         end if

         call mpi_send(chunk_send_top_u,chunksize,MPI_DOUBLE_COMPLEX,next,tag_u,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_v,chunksize,MPI_DOUBLE_COMPLEX,next,tag_v,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_w,chunksize,MPI_DOUBLE_COMPLEX,next,tag_w,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_t,chunksize,MPI_DOUBLE_COMPLEX,next,tag_t,MPI_COMM_WORLD,ierror)

      else

      !Odd processors receive chunks from the even ones

      do nrec=1,maxrec

         call mpi_recv(chunk_recv_u,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_u,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_v,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_v,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_w,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_w,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_t,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_t,MPI_COMM_WORLD,status,ierror)

         if(status(MPI_SOURCE)==next) then

            do ikx=1,iktx
               do iky=1,ikty
                  do incrz=1,2

                     uk(ikx,iky,n3h2+1-incrz)=chunk_recv_u(ikx,iky,incrz)
                     vk(ikx,iky,n3h2+1-incrz)=chunk_recv_v(ikx,iky,incrz)
                     wk(ikx,iky,n3h2+1-incrz)=chunk_recv_w(ikx,iky,incrz)
                     tk(ikx,iky,n3h2+1-incrz)=chunk_recv_t(ikx,iky,incrz)

                  end do
               end do
            end do

         elseif(status(MPI_SOURCE)==prev) then

            do ikx=1,iktx
               do iky=1,ikty
                  do incrz=1,2

                     uk(ikx,iky,incrz)=chunk_recv_u(ikx,iky,incrz)
                     vk(ikx,iky,incrz)=chunk_recv_v(ikx,iky,incrz)
                     wk(ikx,iky,incrz)=chunk_recv_w(ikx,iky,incrz)
                     tk(ikx,iky,incrz)=chunk_recv_t(ikx,iky,incrz)

                  end do
               end do
            end do

         else
            write(*,*) "There is a tag/source mismatch in mype=",mype
         end if

      end do

   end if


     !Now, the inverse, odd processors send their chunks

      if(mod( (mype+2),2)/=0) then

         call mpi_send(chunk_send_bot_u,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_u,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_v,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_v,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_w,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_w,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_bot_t,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_t,MPI_COMM_WORLD,ierror)

         if(mype/=(npe-1)) then
         call mpi_send(chunk_send_top_u,chunksize,MPI_DOUBLE_COMPLEX,next,tag_u,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_v,chunksize,MPI_DOUBLE_COMPLEX,next,tag_v,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_w,chunksize,MPI_DOUBLE_COMPLEX,next,tag_w,MPI_COMM_WORLD,ierror)
         call mpi_send(chunk_send_top_t,chunksize,MPI_DOUBLE_COMPLEX,next,tag_t,MPI_COMM_WORLD,ierror)
         end if

      else

         !And even ones receive chunks from odds

         do nrec=1,maxrec

         call mpi_recv(chunk_recv_u,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_u,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_v,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_v,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_w,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_w,MPI_COMM_WORLD,status,ierror)
         call mpi_recv(chunk_recv_t,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_t,MPI_COMM_WORLD,status,ierror)


          if(status(MPI_SOURCE)==next) then

            do ikx=1,iktx
               do iky=1,ikty
                  do incrz=1,2

                     uk(ikx,iky,n3h2+1-incrz)=chunk_recv_u(ikx,iky,incrz)
                     vk(ikx,iky,n3h2+1-incrz)=chunk_recv_v(ikx,iky,incrz)
                     wk(ikx,iky,n3h2+1-incrz)=chunk_recv_w(ikx,iky,incrz)
                     tk(ikx,iky,n3h2+1-incrz)=chunk_recv_t(ikx,iky,incrz)

                  end do
               end do
            end do

         elseif(status(MPI_SOURCE)==prev) then

            do ikx=1,iktx
               do iky=1,ikty
                  do incrz=1,2

                     uk(ikx,iky,incrz)=chunk_recv_u(ikx,iky,incrz)
                     vk(ikx,iky,incrz)=chunk_recv_v(ikx,iky,incrz)
                     wk(ikx,iky,incrz)=chunk_recv_w(ikx,iky,incrz)
                     tk(ikx,iky,incrz)=chunk_recv_t(ikx,iky,incrz)

                  end do
               end do
            end do

         else
            write(*,*) "There is a tag/source mismatch in mype=",mype
         end if




      end do

   end if

 END SUBROUTINE generate_halo





    SUBROUTINE generate_halo_q(qk) 

      !This subroutine distributes the data to create the proper halo for qk. On the input, q is also n3h1, but the values in the halo (first and last levels) are insignificant. 
      !The way it proceeds is that it first save the first and last significant levels for each processor on a chunk (chun_send_..._...). This chunk will become the previous and next processor's new top and bottom halo.
      !So first, the even taged (mype 0, 2, 4 ...) processor send both of their chunk (execpt for mype 0 and npe-1 that respectively just send their top and bottom) to the right (previous or next) processor.
      !Then, we do the inverse: odd processors (mype = 1, 3, ...) send their chunk to the even ones.

      !This subroutine is based on the fact that the number of processors is even and that n3h0 is at least 2, and must be modified otherwise.
      

      integer :: prev,next,tag
      integer :: nrec,maxrec

      integer :: ikx,iky,incrz

      integer :: chunksize
      integer, parameter :: tag_q=5

      double complex, dimension(iktx,ikty) :: chunk_send_bot_q,chunk_send_top_q
      double complex, dimension(iktx,ikty) :: chunk_recv_q
      double complex, dimension(iktx,ikty,n3h1), intent(inout) :: qk

      !maxsend is the number of processors from which mype receives velocity fields

      if(mype==0) then
         maxrec=1
      elseif(mype==(npe-1)) then
         maxrec=1
      else
         maxrec=2
      end if

      !Send to which processor                                                                                                                                                                                                               
      prev=mype-1
      next=mype+1

      !Create the chunks to send                                                                                                                                                                                                                                  

      do ikx=1,iktx
         do iky=1,ikty
 
               if(mype/=0) then
               chunk_send_bot_q(ikx,iky)=qk(ikx,iky,2)
               end if

               if(mype/=(npe-1)) then
               chunk_send_top_q(ikx,iky)=qk(ikx,iky,n3h0+1)
               end if

         end do
      end do

      chunksize=iktx*ikty

      !Even processors send their chunks

      if(mod( (mype+2),2)==0) then

         if(mype/=0) then
         call mpi_send(chunk_send_bot_q,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_q,MPI_COMM_WORLD,ierror)
         end if

         call mpi_send(chunk_send_top_q,chunksize,MPI_DOUBLE_COMPLEX,next,tag_q,MPI_COMM_WORLD,ierror)
      else

      !Odd processors receive chunks from the even ones

      do nrec=1,maxrec

         call mpi_recv(chunk_recv_q,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_q,MPI_COMM_WORLD,status,ierror)

         if(status(MPI_SOURCE)==next) then

            do ikx=1,iktx
               do iky=1,ikty

                     qk(ikx,iky,n3h1)=chunk_recv_q(ikx,iky)

               end do
            end do

         elseif(status(MPI_SOURCE)==prev) then

            do ikx=1,iktx
               do iky=1,ikty

                     qk(ikx,iky,1)=chunk_recv_q(ikx,iky)

               end do
            end do

         else
            write(*,*) "There is a tag/source mismatch in mype=",mype
         end if

      end do

   end if


     !Now, the inverse, odd processors send their chunks

      if(mod( (mype+2),2)/=0) then

         call mpi_send(chunk_send_bot_q,chunksize,MPI_DOUBLE_COMPLEX,prev,tag_q,MPI_COMM_WORLD,ierror)

         if(mype/=(npe-1)) then
         call mpi_send(chunk_send_top_q,chunksize,MPI_DOUBLE_COMPLEX,next,tag_q,MPI_COMM_WORLD,ierror)
         end if

      else

         !And even ones receive chunks from odds

         do nrec=1,maxrec

         call mpi_recv(chunk_recv_q,chunksize,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,tag_q,MPI_COMM_WORLD,status,ierror)

          if(status(MPI_SOURCE)==next) then

            do ikx=1,iktx
               do iky=1,ikty

                     qk(ikx,iky,n3h1)=chunk_recv_q(ikx,iky)

               end do
            end do

         elseif(status(MPI_SOURCE)==prev) then

            do ikx=1,iktx
               do iky=1,ikty

                     qk(ikx,iky,1)=chunk_recv_q(ikx,iky)

               end do
            end do

         else
            write(*,*) "There is a tag/source mismatch in mype=",mype
         end if




      end do

   end if

 END SUBROUTINE generate_halo_q




END MODULE mpi
