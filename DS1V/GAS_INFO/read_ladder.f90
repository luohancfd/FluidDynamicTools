program main
implicit none
type :: aho_dt
  logical :: iaho=.false.
  integer :: vmax, mnphase
  real(8),pointer :: vphase(:,:,:)
end type aho_dt
type(aho_dt) :: aho_dat
integer :: i

open(10,file='on_vphase.bin',form='unformatted', action='read')
read(10) aho_dat%vmax
write(*,*) aho_dat%vmax
read(10) aho_dat%mnphase
write(*,*) aho_dat%mnphase
allocate(aho_dat%vphase(aho_dat%mnphase,2,aho_dat%vmax+1))
read(10) aho_dat%vphase(:,1,:)
read(10) aho_dat%vphase(:,2,:)
do i=0,aho_dat%vmax,1
    write(*,*) aho_dat%vphase(5,1,i+1)
end do
! write(*,*) aho_dat(i)%vphase(aho_dat(i)%mnphase,2,aho_dat(i)%vmax+1)
close(10)

end program main
