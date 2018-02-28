

program FFT_MPI_3D
use, intrinsic :: iso_c_binding
implicit none
include 'mpif.h'
include 'fftw3-mpi.f03'
integer(C_INTPTR_T), parameter :: L =1
integer(C_INTPTR_T), parameter :: M =1
integer(C_INTPTR_T), parameter :: N =8
type(C_PTR) :: planfft,planifft, cdata,data2
complex(C_DOUBLE_COMPLEX), pointer :: fdata(:,:,:)
integer(C_INTPTR_T) :: alloc_local, local_N, local_N_offset
integer(C_INTPTR_T) :: i, j,k
complex(C_DOUBLE_COMPLEX) :: fout
integer :: ierr, myid, nproc
real,parameter :: pi2=2*3.141592654
real rph
! Initialize
call mpi_init(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call fftw_mpi_init()
!local_N=8
!local_N_offset=N/local_N
! get local data size and allocate (note dimension reversal)
alloc_local = fftw_mpi_local_size_3d(N, M, L, MPI_COMM_WORLD, local_N, local_N_offset)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
print*,local_N, local_N_offset
cdata = fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata, fdata, [L, M, local_N])
! create MPI plan for in-place forward DFT (note dimension reversal)
planfft = fftw_mpi_plan_dft_3d(N, M, L, fdata, fdata, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE)
planifft = fftw_mpi_plan_dft_3d( N,M,l , Fdata, Fdata, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE)
! initialize data to some function my_function(i,j)
do k = 1, local_N
do j = 1, M
do i = 1, L
        !CALL RANDOM_NUMBER(rph)
        fdata(i, j,k) = cmplx(k+local_N_offset,0)  !fout
enddo
enddo
enddo

! compute transform (as many times as desired)
if(myid.eq.0) write(*,100) fdata(1,3,3)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
call fftw_mpi_execute_dft(planfft, fdata, fdata)!
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
if(myid.eq.0) write(*,100) fdata
call fftw_mpi_execute_dft(planifft,  fdata, fdata)!
if(myid.eq.0) write(*,100) fdata(1,3,3)/L/M/N
! deallocate and destroy plans
call fftw_destroy_plan(planfft)
call fftw_destroy_plan(planifft)
call fftw_mpi_cleanup()
call fftw_free(cdata)
call fftw_free(data2)
call mpi_finalize(ierr)
100 format(6f15.6)
end



