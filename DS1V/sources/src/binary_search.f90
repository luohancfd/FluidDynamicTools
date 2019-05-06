
! subroutine to search index ii and ij in x meets
!       x(ii) < x0 < x(ij)
! if x0 equals to some value in x, ii=ij=index of the element
      subroutine binary_search(ii,ij,x0,x,lx,rx,ierror)
        implicit none
        integer :: ii, ij
        integer,intent(in) ::  rx,lx
        double precision,intent(in) :: x0,x(lx:rx)
        logical :: ierror
        integer :: L,R,M

        ierror = .false.
        L = lx
        R = rx
        M = (L+R)/2

        if ( (x0 .lt. x(L)) .or. (x0 .gt. x(R)))then
          write(6,*) 'No correct phase is found'
          ierror = .true.
          return
        elseif (x0 .eq. x(L)) then
          ii = L
          ij = L
          return
        elseif (x0 .eq. x(R)) then
          ii = R
          ij = R
          return
        endif

        do while ( (x(L) .le. x0) .and. (x(R) .ge. x0) .and. (R-L)>1 )
          if (x0 .eq. x(M)) then
            ii = M
            ij = M
            R = M
            L =M
            exit
          elseif (x0 .gt. x(M))then
            L = M
            M = (L+R)/2
          elseif (x0 .le. x(M)) then
            R = M
            M = (L+R)/2
          endif
        enddo
        if ( x(R) .eq. x0) then
          ii = R
          ij = R
        elseif ( x(L) .eq. x0) then
          ii = L
          ij = L
        else
          ii = L
          ij = R
        end if

      end subroutine binary_search
