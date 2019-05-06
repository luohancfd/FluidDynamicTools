!
!*****************************************************************************
!
SUBROUTINE SET_ZIGGURAT
!
!--author: isebasti
!--calculations required by Ziggurate Method
!--call this subroutine with a different jsrl value for a different set of seeds
!
USE CALC
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: nth,idtl,ishr3  !local variables
!
nth=1                      !for serial compilation
!$omp parallel
!$omp master
!$    nth=omp_get_num_threads( )
!$    write(*,*) 'Available threads     ', omp_get_max_threads( )
!$    write(*,*) 'Threads in use        ', nth
!$omp end master
!$omp end parallel
!
allocate(iseed(0:nth-1))
!
do idtl=0,nth-1
call ishr(ishr3)
iseed(idtl)=ishr3         !set a seed for each thread
!$    write(*,*) 'Thread',idtl,'Seed',iseed(idtl)
end do
!
RETURN
END SUBROUTINE SET_ZIGGURAT
