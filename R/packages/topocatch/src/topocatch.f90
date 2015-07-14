
!USE FOR PRINTING: See section 'Printing from FORTRAN' in R-exts.pdf
!  subroutine dblepr(label, nchar, data, ndata )
!  subroutine realpr(label, nchar, data, ndata )
!  subroutine intpr (label, nchar, data, ndata )

!###############################################################################
!### BEGIN OF MODULE ###########################################################
!###############################################################################

module topocatchInternals
implicit none

!Generic interfaces
interface adjacentCells
  module procedure adjacentCellsDbl, adjacentCellsInt
end interface adjacentCells

contains

!-------------------------------------------------------------------------------
!Get the values of the neighboring cell of a cell. If the cell is located
!at the margins of the matrix, the neighbors outside are assumed to be nodata.
!Position 1 refers to the upper left (north-west) corner and counts clockwise.
!  1  2  3
!  8  x  4
!  7  6  5
!There is an implementation for matrices of type double and integer
function adjacentCellsDbl(matrix,col,row,nodata) result (values)
implicit none
  !input
  real(8),dimension(:,:),intent(in):: matrix
  integer,intent(in):: col,row
  real(8),intent(in):: nodata
  !output
  real(8),dimension(8):: values
  !local
  integer:: ncols, nrows
  !code
  ncols= size(matrix, dim=1)
  nrows= size(matrix, dim=2)
  values(:)= nodata
  if ((col .lt. 1) .or. (col .gt. ncols)) return
  if ((row .lt. 1) .or. (row .gt. nrows)) return
  if (col .gt. 1) then
    values(8)= matrix(col-1,row)
    if (row .gt. 1) values(1)= matrix(col-1,row-1)
    if (row .lt. nrows) values(7)= matrix(col-1,row+1)
  end if
  if (col .lt. ncols) then
    values(4)= matrix(col+1,row)
    if (row .gt. 1) values(3)= matrix(col+1,row-1)
    if (row .lt. nrows) values(5)= matrix(col+1,row+1)
  end if  
  if (row .gt. 1) values(2)= matrix(col,row-1)
  if (row .lt. nrows) values(6)= matrix(col,row+1)
end function adjacentCellsDbl
function adjacentCellsInt(matrix,col,row,nodata) result (values)
implicit none
  !input
  integer,dimension(:,:),intent(in):: matrix
  integer,intent(in):: col,row
  integer,intent(in):: nodata
  !output
  integer,dimension(8):: values
  !local
  integer:: ncols, nrows
  !code
  ncols= size(matrix, dim=1)
  nrows= size(matrix, dim=2)
  values(:)= nodata
  if ((col .lt. 1) .or. (col .gt. ncols)) return
  if ((row .lt. 1) .or. (row .gt. nrows)) return
  if (col .gt. 1) then
    values(8)= matrix(col-1,row)
    if (row .gt. 1) values(1)= matrix(col-1,row-1)
    if (row .lt. nrows) values(7)= matrix(col-1,row+1)
  end if
  if (col .lt. ncols) then
    values(4)= matrix(col+1,row)
    if (row .gt. 1) values(3)= matrix(col+1,row-1)
    if (row .lt. nrows) values(5)= matrix(col+1,row+1)
  end if  
  if (row .gt. 1) values(2)= matrix(col,row-1)
  if (row .lt. nrows) values(6)= matrix(col,row+1)
end function adjacentCellsInt

!-------------------------------------------------------------------------------
!Returns the direction of the steepest downward gradient when looking from a
!cell into its 8-cell environment.
!Return value:
! - The 'nodata_out' value, if any of the cells in the 8-cell environment is nodata_in
! - The 'undefined' value, if the area is flat or the cell in the center is sink
! - A value in 1:8, encoding a direction 1==nw,2==nn,ne,ee,se,ss,sw,8==ww.
!   In the case that the max. gradient orrcurs in more than one direction,
!   only the code of the FIRST direction is returned.
function maxDownDir(center,nw,nn,ne,ee,se,ss,sw,ww,nodata_in,nodata_out,undefined) result(dir)
implicit none
  !input
  real(8),intent(in):: center,nw,nn,ne,ee,se,ss,sw,ww
  real(8),intent(in):: nodata_in
  integer,intent(in):: nodata_out, undefined
  !output
  integer:: dir
  !local
  real(8),parameter:: R2= sqrt(2.0d0), UN= 1.0d0
  real(8),dimension(8):: grad,marg,dist
  integer:: i
  !code
  marg= (/nw,nn,ne,ee,se,ss,sw,ww/)
  dist= (/R2,UN,R2,UN,R2,UN,R2,UN/)
  if (any(marg .eq. nodata_in)) then
    dir= nodata_out
  else
    grad= ((/(center,i=1,8)/) - marg(:)) / dist(:)
    if (minval(marg) .eq. center) then
      dir= undefined
    else
      dir= maxloc(grad,dim=1)
    end if
  end if
end function maxDownDir

!-------------------------------------------------------------------------------
! Distance of two points
function dist(x1,y1,x2,y2) result(d)
implicit none
  !input
  real(8),intent(in):: x1,y1,x2,y2
  !output
  real(8):: d
  !local
  real(8),parameter:: TWO= 2.0d0
  !code
  d= sqrt((x1-x2)**TWO + (y1-y2)**TWO)
end function dist

!-------------------------------------------------------------------------------
! Position of first match in an array. Returns zero if no match was found.
integer function match(x, array)
  integer,intent(in):: x
  integer,dimension(:),intent(in):: array
  integer:: i
  match= 0
  do i=1,size(array)
    if (array(i) .eq. x) then
      match= i
      exit      
    end if
  end do
end function match

!-------------------------------------------------------------------------------
!Functions to return the coordinates of a grid cell's margins
!Cell(1,1) is the cell in the upper left corner of the grid.
!Cell(maxcol,maxrow) is the cell in the lower right corner of the grid.
real(8) function cellXmin(xllcorner,cellsize,col)
implicit none
real(8),intent(in):: xllcorner, cellsize
integer,intent(in):: col
  cellXmin= xllcorner + dble(col-1)*cellsize
end function cellXmin
real(8) function cellXmax(xllcorner,cellsize,col)
implicit none
real(8),intent(in):: xllcorner, cellsize
integer,intent(in):: col
  cellXmax= xllcorner + dble(col)*cellsize
end function cellXmax
real(8) function cellYmin(yllcorner,cellsize,nrow,row)
implicit none
real(8),intent(in):: yllcorner, cellsize
integer,intent(in):: nrow,row
  cellYmin= yllcorner + dble(nrow-row)*cellsize
end function cellYmin
real(8) function cellYmax(yllcorner,cellsize,nrow,row)
implicit none
real(8),intent(in):: yllcorner, cellsize
integer,intent(in):: nrow,row
  cellYmax= yllcorner + dble(nrow-(row-1))*cellsize
end function cellYmax
!-------------------------------------------------------------------------------
!Functions to return the coordinates of a grid cell's center
!Cell(1,1) is the cell in the upper left corner of the grid.
!Cell(maxcol,maxrow) is the cell in the lower right corner of the grid.
real(8) function cellXavg (col,xllcorner,cellsize)
  integer,intent(in):: col
  real(8),intent(in):: xllcorner, cellsize
  real(8),parameter:: HALF= 0.5d0
  cellXavg= xllcorner + (col-1) * cellsize + HALF * cellsize
end function cellXavg
real(8) function cellYavg (row,nrows,yllcorner,cellsize)
  integer,intent(in):: row, nrows
  real(8),intent(in):: yllcorner, cellsize
  real(8),parameter:: HALF= 0.5d0
  cellYavg= yllcorner + (nrows-row) * cellsize + HALF * cellsize
end function cellYavg

!-------------------------------------------------------------------------------
!Returns the column and row index of a grid cell which contains a point defined
!by x/y coordinates. If the point is located outside the grid the indices are
!set zo zero.
subroutine whichCell(ncol, nrow, xll, yll, cellsize, x, y, colindex, rowindex)
implicit none
  !input
  integer,intent(in):: ncol, nrow      !number of cols and rows of the grid
  real(8),intent(in):: xll, yll        !coordinates of lower-left grid corner
  real(8),intent(in):: cellsize        !cellsize of the grid
  real(8),intent(in):: x, y            !coordinates of the point
  !output
  integer,intent(out):: colindex, rowindex  !indices of matching cell (or zero)
  !local
  integer:: col,row
  !CODE
  colindex= 0
  rowindex= 0
  !Find the target row
  do row=1,nrow
    if (((y .ge. cellYmin(yll,cellsize,nrow,row)) .and. &
         (y .lt. cellYmax(yll,cellsize,nrow,row))) .or. &
        ((row .eq. 1) .and. (y .eq. cellYmax(yll,cellsize,nrow,row)))) then
      rowindex= row
      exit
    end if
  end do
  !Find the target column
  do col=1,ncol
    if (((x .ge. cellXmin(xll,cellsize,col)) .and. &
         (x .lt. cellXmax(xll,cellsize,col))) .or. &
        ((col .eq. ncol) .and. (x .eq. cellXmax(xll,cellsize,col)))) then
      colindex= col
      exit
    end if
  end do
end subroutine whichCell

!-------------------------------------------------------------------------------
!Find col/row indices of the nearest neighbor of a cell.
! - Parameter 'maxdist' controls the (maximum) search DISTANCE.
! - The parameters 'range_lbnd', 'range_ubnd', and 'in_range' define the
!   VALUE which a cell must have to be accepted as a neighbor.
subroutine findNN(ncols,nrows,matrix,cellsize,nodata,col,row,maxdist, &
  range_lbnd,range_ubnd,in_range,col_neighbor,row_neighbor)
implicit none
  type t_index
    integer:: col,row
  end type t_index
  !input
  integer,intent(in):: ncols,nrows
  integer,dimension(1:ncols,1:nrows),intent(in):: matrix
  integer,intent(in):: nodata
  real(8),intent(in):: cellsize
  integer,intent(in):: col,row
  real(8),intent(in):: maxdist
  integer,intent(in):: range_lbnd,range_ubnd
  logical,intent(in):: in_range
  !output
  integer,intent(out):: col_neighbor,row_neighbor
  !local
  integer:: i,pass,lcol,rcol,trow,brow,ncol_pass,nrow_pass
  type(t_index),dimension(:),allocatable:: xlist
  real(8):: dist_neighbor,temp
  !code
  pass=0
  col_neighbor= 0
  row_neighbor= 0
  dist_neighbor= huge(dist_neighbor)
  do
    !Set parameters of search ring
    pass= pass + 1
    trow= max(1,row-pass)
    brow= min(nrows,row+pass)
    lcol= max(1,col-pass)
    rcol= min(ncols,col+pass)  
    !Stop if maximum search distance is exceeded
    if ((dble(max(col-lcol,rcol-col,row-trow,brow-row)) * cellsize) .gt. maxdist) exit
    ncol_pass= rcol-lcol+1
    nrow_pass= brow-trow+1
    allocate(xlist(2*ncol_pass+2*nrow_pass))
    xlist(1:nrow_pass)%col=                         (/(lcol,i=trow,brow)/)
    xlist(nrow_pass+1:2*nrow_pass)%col=                 (/(rcol,i=trow,brow)/)
    xlist(2*nrow_pass+1:2*nrow_pass+ncol_pass)%col=         (/(i,i=lcol,rcol)/)
    xlist(2*nrow_pass+ncol_pass+1:2*nrow_pass+2*ncol_pass)%col= (/(i,i=lcol,rcol)/)
    xlist(1:nrow_pass)%row=                         (/(i,i=trow,brow)/)
    xlist(nrow_pass+1:2*nrow_pass)%row=                 (/(i,i=trow,brow)/)
    xlist(2*nrow_pass+1:2*nrow_pass+ncol_pass)%row=         (/(trow,i=lcol,rcol)/)
    xlist(2*nrow_pass+ncol_pass+1:2*nrow_pass+2*ncol_pass)%row= (/(brow,i=lcol,rcol)/)
    !Search and remember value of nearest cell
    do i=1,size(xlist)
      !Skip nodata cells
      if (matrix(xlist(i)%col,xlist(i)%row) .eq. nodata) goto 1000
      if (in_range .and. ((matrix(xlist(i)%col,xlist(i)%row) .lt. range_lbnd) .or. &
        (matrix(xlist(i)%col,xlist(i)%row) .gt. range_ubnd))) goto 1000
      if (.not. in_range .and. ((matrix(xlist(i)%col,xlist(i)%row) .ge. range_lbnd) .and. &
        (matrix(xlist(i)%col,xlist(i)%row) .le. range_ubnd))) goto 1000
        !Remember position if this is the nearest neighbor
        temp= sqrt(real(col-xlist(i)%col)**2.0+real(row-xlist(i)%row)**2.0)
        if (temp .lt. dist_neighbor) then
          col_neighbor= xlist(i)%col
          row_neighbor= xlist(i)%row
          dist_neighbor= temp
        end if
      1000 continue
    end do
    deallocate(xlist)
    !Stop search if successful
    if ((col_neighbor*row_neighbor) .gt. 0) exit  
    !Stop if no neighbors were found
    if ((trow.eq.1).and.(brow.eq.nrows).and.(lcol.eq.1).and.(rcol.eq.ncols)) exit
  end do
end subroutine findNN

! Print progress info
! NOTE: C-call used for printing!
subroutine progress (first, last, current, delta)
  implicit none
  interface
    subroutine printstate(x)
      real(8):: x
    end subroutine printstate
  end interface
  real(8), intent(in):: first, last, current, delta
  real(8):: percDone_now, percDone_old
  integer:: i
  integer,dimension(7),parameter:: p=(/1,2,5,10,20,50,100/)
  percDone_now= (100.*(current-first))/(last-first)
  percDone_old= (100.*(current-first-delta))/(last-first)
  do i=1, size(p)
    if ((percDone_now .ge. dble(p(i))) .and. (percDone_old .lt. dble(p(i)))) then
      call printstate(dble(p(i)))
      exit
    end if
  end do
end subroutine progress


end module topocatchInternals

!###############################################################################
!### END OF MODULE #############################################################
!###############################################################################

!-------------------------------------------------------------------------------
! Reclassification of values in a vector using a reclass table
subroutine vectreclass(vectlen, vect_in, nclasses, class_min, class_max, &
  class_val, include_min, include_max, vect_out)
implicit none
  !input
  integer,intent(in):: vectlen
  real(8),dimension(vectlen),intent(in):: vect_in
  integer,intent(in):: nclasses
  real(8),dimension(nclasses),intent(in):: class_min
  real(8),dimension(nclasses),intent(in):: class_max
  real(8),dimension(nclasses),intent(in):: class_val
  integer,intent(in):: include_min
  integer,intent(in):: include_max
  !output
  real(8),dimension(vectlen),intent(out):: vect_out
  !local
  integer:: i,k
  !code
  vect_out= vect_in
  if ((include_min .ne. 0) .and. (include_max .ne. 0)) then
    do i=1,vectlen
      do k=1,nclasses
        if ((vect_in(i) .ge. class_min(k)) .and. (vect_in(i) .le. class_max(k))) then
          vect_out(i)= class_val(k)
          exit
        end if
      end do
    end do
  else if ((include_min .ne. 0) .and. (include_max .eq. 0)) then
    do i=1,vectlen
      do k=1,nclasses
        if ((vect_in(i) .ge. class_min(k)) .and. (vect_in(i) .lt. class_max(k))) then
          vect_out(i)= class_val(k)
          exit
        end if
      end do
    end do
  else if ((include_min .eq. 0) .and. (include_max .ne. 0)) then
    do i=1,vectlen
      do k=1,nclasses
        if ((vect_in(i) .gt. class_min(k)) .and. (vect_in(i) .le. class_max(k))) then
          vect_out(i)= class_val(k)
          exit
        end if
      end do
    end do
  else
    do i=1,vectlen
      do k=1,nclasses
        if ((vect_in(i) .gt. class_min(k)) .and. (vect_in(i) .lt. class_max(k))) then
          vect_out(i)= class_val(k)
          exit
        end if
      end do
    end do
  end if
end subroutine vectreclass

!-------------------------------------------------------------------------------
! Filling of sinks in an elevation model using a non-iterative approach.
! Notes:
!  - The elevation matrix must be passed as a vector. This vector must contain
!    the concatenated columns of the matrix. Thus, if the original matrix is
!       1 4 7
!       2 5 8
!       3 6 9
!    the input vector must be: 1,2,3,4,5,6,7,8,9
!  - The elevation model is processes after scaling and rounding to reduce the
!    effort. This is because the filling process is
!    done for every UNIQUE value being present in the input.
!  - The processed elevation model is returned as a vector. This vector contains
!    the concatenated columns of the corresponding matrix, thus, if a matrix
!    is to be constructed from that vector this has to be done 'by column'.
subroutine sinkfill(ncol, nrow, vect_in, nodata, ndigits, silent, vect_out)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  real(8),dimension(ncol*nrow),intent(in):: vect_in
  real(8),intent(in):: nodata
  integer,intent(in):: ndigits
  integer,intent(in):: silent
  !output
  real(8),dimension(ncol*nrow),intent(out):: vect_out
  !local
  real(8),parameter:: ZERO= 0.0_8, ONE= 1.0_8, TEN=10.0_8
  integer:: z
  integer:: mini, maxi
  integer:: newsum
  integer:: icol,irow
  integer,dimension(1:ncol,1:nrow):: matrix_in
  integer:: nodata_tmp
  !Temporary grid and matrix with a two-cell nodata margin
  integer,dimension(1:(ncol+4),1:(nrow+4)):: tmp
  logical,dimension(1:(ncol+4),1:(nrow+4)):: wetmatrix
  !code
  !Transform input vector (concatenated columns) into matrix,
  !and scale and round data
  do icol=1,ncol
    matrix_in(icol,:)= nint(vect_in( (icol-1)*nrow+1 : icol*nrow ) * (TEN**ndigits))
  end do
  !Determine value limits
  mini= nint(minval(vect_in, mask=(vect_in .ne. nodata)))
  maxi= nint(maxval(vect_in, mask=(vect_in .ne. nodata)))
  !Set additional margin in tmp to nodata and copy central values from input
  nodata_tmp= mini - TEN
  tmp(:,:)= nodata_tmp
  tmp(3:(size(tmp,dim=1)-2),3:(size(tmp,dim=2)-2))= matrix_in(:,:)
  !Initialize wetmatrix=true where dem=nodata
  wetmatrix(:,:)= .false.
  where (tmp .eq. nodata_tmp)
    wetmatrix= .true.
  end where
  !Loop over elevation values (in increasing order)
  do z=mini, maxi

    !Progress info
    !if (silent .eq. 0) call intpr("Filling at scaled level (current / min / max):", -1, (/z,mini,maxi/), 3)
    if (silent .eq. 0) call progress(dble(mini), dble(maxi), dble(z), 1.0_8)

    !Find connected wet cells at current level
    newsum= 1
    do while (newsum .gt. 0)
      newsum= 0
      !Loop through grid
      do icol=2,(size(tmp,dim=1)-1)
        do irow=2,(size(tmp,dim=2)-1)
          if ((.not. wetmatrix(icol,irow)) .and. (tmp(icol,irow) .lt. z)) then
            if (wetmatrix(icol-1,irow-1)) goto 1
            if (wetmatrix(icol+1,irow+1)) goto 1
            if (wetmatrix(icol-1,irow+1)) goto 1
            if (wetmatrix(icol+1,irow-1)) goto 1
            if (wetmatrix(icol  ,irow-1)) goto 1
            if (wetmatrix(icol  ,irow+1)) goto 1
            if (wetmatrix(icol-1,irow  )) goto 1
            if (wetmatrix(icol+1,irow  )) goto 1
            goto 2
            1 wetmatrix(icol,irow)= .true.
            newsum= newsum + 1
            2 continue
          end if
        end do
      end do
    end do
    !Determine unconnected wet cells at current level
    where ((.not. wetmatrix) .and. (tmp .lt. z))
      tmp= z
    end where
  end do
  !Transform result matrix into output vector (of concatenated colums)
  !If the contents of tmp is
  !   0 0 0 0 0 0 0
  !   0 0 0 0 0 0 0
  !   0 0 1 3 5 0 0
  !   0 0 2 4 6 0 0
  !   0 0 0 0 0 0 0
  !   0 0 0 0 0 0 0
  !the resulting vect_out is: 1,2,3,4,5,6
  !Also undo the scaling
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= tmp((icol+2), 3:(size(tmp,dim=2)-2)) / dble(TEN**ndigits)
  end do
  !Set nodata
  where (vect_in .eq. nodata)
    vect_out= nodata
  end where 
end subroutine sinkfill

!-------------------------------------------------------------------------------
! Analyze flow directions in a digital elevation model.
! Sinks in the DEM should have been filled before.
! The direction codes in the output are defined as follows (1 == Northwest):
!   1  2  3
!   8  x  4
!   7  6  5
! The special code of 0 is returned for cells where no direction could be
! computed. These may be cells which
! - represent sinks (not removed in earlier step of DEM processing)
! - represent flat areas at the margin of the DEM
subroutine flowdir(ncol, nrow, vect_in, nodata, silent, vect_out)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  real(8),dimension(ncol*nrow),intent(in):: vect_in
  real(8),intent(in):: nodata
  integer,intent(in):: silent
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_out
  !local
  real(8),parameter:: ZERO= 0.0d8
  integer,parameter:: LINK_MISSING= 0
  integer,parameter:: NODATA_OUT=-9999
  real(8),dimension(1:ncol,1:nrow):: matrix_in
  integer,dimension(1:ncol,1:nrow):: matrix_out
  integer,dimension(1:ncol,1:nrow):: matrix_tmp
  integer,dimension(8),parameter:: colshift= (/-1,0,1,1,1,0,-1,-1/)
  integer,dimension(8),parameter:: rowshift= (/-1,-1,-1,0,1,1,1,0/)
  integer,dimension(8),parameter:: directions= (/2,4,6,8,1,3,5,7/) ! straight links before diagonal first
  integer:: icol,irow,c,r,i
  logical:: change
  !code
  !Transform input vector (concatenated columns) into matrix
  do icol=1,ncol
    matrix_in(icol,:)= vect_in( (icol-1)*nrow+1 : icol*nrow ) 
  end do
  !Initialize output matrix to nodata
  matrix_out= NODATA_OUT
  !1st pass: Process central part of the grid (margin remains nodata)
  if (silent .eq. 0) call intpr("Processing central part...",-1, 0, 0)      
  do icol=2,(ncol-1)
    if (silent .eq. 0) call progress(dble(2), dble(ncol-1), dble(icol), dble(1))
    do irow=2,(nrow-1)
      matrix_out(icol,irow)= maxDownDir( &
        matrix_in(icol,irow), &     
        matrix_in(icol-1,irow-1),matrix_in(icol  ,irow-1), &   !nw,n
        matrix_in(icol+1,irow-1),matrix_in(icol+1,irow  ), &   !ne,e
        matrix_in(icol+1,irow+1),matrix_in(icol  ,irow+1), &   !se,s
        matrix_in(icol-1,irow+1),matrix_in(icol-1,irow  ), &   !sw,w
        nodata,NODATA_OUT,LINK_MISSING)
    end do
  end do
  !2nd pass: Handling of flat areas
  if (silent .eq. 0) call intpr("Processing flat areas...",-1, 0, 0)      
  change= .true.
  do while (change)
    change= .false.
    matrix_tmp(:,:)= NODATA_OUT
    do icol=2,(ncol-1)
      do irow=2,(nrow-1)
        !Find a neighboring cell into which unlinked cells could drain
        if (matrix_out(icol,irow) .eq. LINK_MISSING) then
          !Check if neighboring cell is linked and its elevation is not higher.
          !We accept the FIRST cell (of the 8 possible cells) that fits.
          do i=1,8
            c= icol+colshift(directions(i))
            r= irow+rowshift(directions(i))
            !Flat cells next to the margin - allow drainage to the margin
            if ((c .eq. 1) .or. (c .eq. ncol) .or. (r .eq. 1) .or. (r .eq. nrow)) then
              if (matrix_in(c,r) .le. matrix_in(icol,irow)) then
                matrix_tmp(icol,irow)= directions(i)
                change= .true.
                exit
              end if
            !In central parts
            else
              if ((matrix_out(c,r) .ne. LINK_MISSING) .and. &
                  (matrix_out(c,r) .ne. NODATA_OUT) .and. &
                  (matrix_in(c,r) .le. matrix_in(icol,irow))) then
                matrix_tmp(icol,irow)= directions(i)
                change= .true.
                exit
              end if
            end if
          end do
        end if
      end do
    end do
    !Save newly linked cells
    where (matrix_tmp .ne. NODATA_OUT)
      matrix_out= matrix_tmp
    end where
  end do
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix_out(icol, 1:nrow)
  end do
end subroutine flowdir

!-------------------------------------------------------------------------------
! Compute flow accumulation based on flow direction codes
! The result is the number of upstream cells (integer), NOT the area!
subroutine flowacc(ncol, nrow, vect_in, nodata, silent, vect_out, errorlevel)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  integer,dimension(ncol*nrow),intent(in):: vect_in
  integer,intent(in):: nodata
  integer,intent(in):: silent
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_out
  integer,intent(out):: errorlevel
  !local
  integer,dimension(1:ncol,1:nrow):: matrix_in, matrix_flowacc
  integer,dimension(8),parameter:: colshift= (/-1,0,1,1,1,0,-1,-1/)
  integer,dimension(8),parameter:: rowshift= (/-1,-1,-1,0,1,1,1,0/)
  integer:: icol, irow, icol_start, irow_start, tmp1, tmp2
  !Transform input vector (concatenated columns) into matrix
  do icol=1,ncol
    matrix_in(icol,:)= vect_in( (icol-1)*nrow+1 : icol*nrow ) 
  end do
  ! (1) Compute flow accumulation grid from a grid of flow direction codes.
  ! Valid flow direction codes are in range 1...8.
  !Initialize output matrix to a value of one (or nodata)
  where (matrix_in .eq. nodata)
    matrix_flowacc= 0
  elsewhere
    matrix_flowacc= 1
  end where
  !Start a flow path at each valid cell
  if (silent .eq. 0) call intpr("Computing flow accumulation...",-1, 0, 0)      
  do icol_start=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol_start), dble(1))
    do irow_start=1,nrow
      if ((matrix_in(icol_start,irow_start) .ge. 1) .and. (matrix_in(icol_start,irow_start) .le. 8)) then
        !Walk along the flow paths
        icol= icol_start
        irow= irow_start
        do
          tmp1= icol + colshift(matrix_in(icol,irow))
          tmp2= irow + rowshift(matrix_in(icol,irow))
          icol= tmp1
          irow= tmp2
          if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
            exit   ! Stop at grid border
          else if ((matrix_in(icol,irow) .lt. 1) .or. (matrix_in(icol,irow) .gt. 8)) then
            exit   ! Stop at nodata or invalid codes (e.g. code 0)
          else
            matrix_flowacc(icol,irow)= matrix_flowacc(icol,irow) + 1   ! Add current flow path
            if (matrix_flowacc(icol,irow) .gt. (ncol*nrow)) then
              goto 1   !We are lost in a sink!
            end if
          end if
        end do
      end if
    end do
  end do
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix_flowacc(icol, 1:nrow)
  end do
  !Set error level
  errorlevel= 0; return
  1 errorlevel= 1; return
end subroutine flowacc

!-------------------------------------------------------------------------------
! Create flow paths (vectors) from the flow direction grid
subroutine flowpaths(ncol, nrow, vect_flowdir, vect_flowacc, xllcorner, yllcorner, cellsize, &
  area_crit, x_incatch, y_incatch, silent, ncoords, path_id, path_x, path_y, id_outlet, errorlevel)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  integer,dimension(ncol*nrow),intent(in):: vect_flowdir
  integer,dimension(ncol*nrow),intent(in):: vect_flowacc
  real(8),intent(in):: xllcorner, yllcorner, cellsize
  real(8),intent(in):: area_crit
  real(8),intent(in):: x_incatch, y_incatch
  integer,intent(in):: silent
  !output
  integer,intent(out):: ncoords
  integer,dimension(ncol*nrow*2),intent(out):: path_id        !dimension = max. possible size
  real(8),dimension(ncol*nrow*2),intent(out):: path_x, path_y  !dimension = max. possible size
  integer,intent(out):: id_outlet
  integer,intent(out):: errorlevel
  !local
  real(8),parameter:: ZERO= 0.0_8, TWO= 2.0_8
  integer,dimension(8),parameter:: colshift= (/-1,0,1,1,1,0,-1,-1/)
  integer,dimension(8),parameter:: rowshift= (/-1,-1,-1,0,1,1,1,0/)
  integer:: icol,irow, icol_start, irow_start, tmp1, tmp2
  integer,dimension(1:ncol,1:nrow):: matrix_flowdir, matrix_flowacc, matrix_npaths
  logical,dimension(1:ncol,1:nrow):: matrix_start_current, matrix_start_new, matrix_start_used
  integer:: id_path
  integer:: npaths_start
  integer:: ncells_crit
  integer:: icol_end, irow_end
  integer:: i,k
  real(8):: dx1,dx2,dy1,dy2
  logical,dimension(1:ncol,1:nrow):: matrix_rivnet, matrix_occupied
  logical,dimension(ncol*nrow*2):: path_keep
  integer:: ioresult
  !code
  !Transform input vectors (concatenated columns) into matrices
  do icol=1,ncol
    matrix_flowdir(icol,:)= vect_flowdir( (icol-1)*nrow+1 : icol*nrow )
    matrix_flowacc(icol,:)= vect_flowacc( (icol-1)*nrow+1 : icol*nrow )
  end do
  ! (1) Compute flow accumulation grid from a grid of flow direction codes.
  ! --> Has become a subroutine input: Range of values is 0....Inf
  ! (2) Set start cells at upstream ends of head reaches
  if (silent .eq. 0) call intpr("Setting source cells...",-1, 0, 0)      
  ! 2a. Select river cells
  ncells_crit= nint(area_crit / (cellsize**TWO))
  if (ncells_crit .lt. 1) then
    goto 1
  end if
  where (matrix_flowacc .ge. ncells_crit)
    matrix_rivnet= .true.
  else where
    matrix_rivnet= .false.
  end where
  ! 2b. Set start cells
  matrix_start_current= .false.
  matrix_occupied= .false.
  do icol_start=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol_start), dble(1))
    do irow_start=1,nrow
      if (matrix_rivnet(icol_start,irow_start)) then
        if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. (matrix_flowdir(icol_start,irow_start) .le. 8)) then
          !Walk along the flow paths
          icol= icol_start
          irow= irow_start
          if (.not. matrix_occupied(icol,irow)) then
            matrix_start_current(icol,irow)= .true.
          end if
          matrix_occupied(icol,irow)= .true. 
          do
            tmp1= icol + colshift(matrix_flowdir(icol,irow))
            tmp2= irow + rowshift(matrix_flowdir(icol,irow))
            icol= tmp1
            irow= tmp2
            if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
              exit   ! Stop at grid border
            else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
              exit   ! Stop at nodata or invalid codes (e.g. code 0)
            else
                matrix_start_current(icol,irow)= .false.
                matrix_occupied(icol,irow)= .true. 
            end if
          end do
        end if
      end if
    end do
  end do
  ! 2c. Find system outlet corresponding to the user-supplied in-catchment coordinates
  if (silent .eq. 0) call intpr("Identifying the system's outlet...",-1, 0, 0)      
  call whichCell(ncol, nrow, xllcorner, yllcorner, cellsize, &
    x_incatch, y_incatch, icol_start, irow_start)
  !Handle bad coordinates 
  if ((icol_start * irow_start) .eq. 0) then
    goto 2
  end if
  if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. (matrix_flowdir(icol_start,irow_start) .le. 8)) then
    !Walk along the flow paths
    icol= icol_start
    irow= irow_start
    do
      tmp1= icol + colshift(matrix_flowdir(icol,irow))
      tmp2= irow + rowshift(matrix_flowdir(icol,irow))
      icol= tmp1
      irow= tmp2
      if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
        exit   ! Stop at grid border
      else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
        exit   ! Stop at nodata or invalid codes (e.g. code 0)
      else
        continue ! Just walk on
      end if
    end do
    icol_end= icol
    irow_end= irow
  else
    goto 3
  end if
  ! 2d. Remove start cells for rivers outside the area of interest
  if (silent .eq. 0) call intpr("Identifying area of interest...",-1, 0, 0)      
  do icol_start=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol_start), dble(1))
    do irow_start=1,nrow
      if (matrix_start_current(icol_start,irow_start)) then
        if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. &
            (matrix_flowdir(icol_start,irow_start) .le. 8)) then
          !Walk along the flow paths
          icol= icol_start
          irow= irow_start        
          do
            tmp1= icol + colshift(matrix_flowdir(icol,irow))
            tmp2= irow + rowshift(matrix_flowdir(icol,irow))
            icol= tmp1
            irow= tmp2
            if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
              exit   ! Stop at grid border
            else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
              exit   ! Stop at nodata or invalid codes (e.g. code 0)
            else
              continue !Just move on
            end if
          end do          
        end if
      end if
      if ((icol .ne. icol_end) .or. (irow .ne. irow_end)) then
        matrix_start_current(icol_start,irow_start)= .false.   !Remove
      end if
    end do
  end do
  ! (3) Compute the number of flow paths / cell
  if (silent .eq. 0) call intpr("Computing path/cell...",-1, 0, 0)      
  ! Walk downhill from the initial start points and register the number of paths / cell
  matrix_npaths= 0
  do icol_start=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol_start), dble(1))
    do irow_start=1,nrow
      if (matrix_start_current(icol_start,irow_start)) then
        if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. &
            (matrix_flowdir(icol_start,irow_start) .le. 8)) then
          !Walk along the flow paths
          icol= icol_start
          irow= irow_start
          matrix_npaths(icol,irow)= matrix_npaths(icol,irow) + 1   ! Add current path
          do
            tmp1= icol + colshift(matrix_flowdir(icol,irow))
            tmp2= irow + rowshift(matrix_flowdir(icol,irow))
            icol= tmp1
            irow= tmp2
            if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
              exit   ! Stop at grid border
            else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
              exit   ! Stop at nodata or invalid codes (e.g. code 0)
            else
              matrix_npaths(icol,irow)= matrix_npaths(icol,irow) + 1   ! Add current path
            end if
          end do          
        end if
      end if
    end do
  end do
!!  open(unit=10, file="/home/dkneis/npaths.asc", status="replace",action="write")
!!  write(10,'(a,i0)')"ncols ",ncol; write(10,'(a,i0)')"nrows ",nrow
!!  write(10,'(a,f0.1)')"xllcorner ",xllcorner; write(10,'(a,f0.1)')"yllcorner ",yllcorner
!!  write(10,'(a,f0.1)')"cellsize ",cellsize; write(10,'(a,i0)')"NODATA_value ",nodata
!!  do irow=1,nrow
!!    write(10,*)(matrix_npaths(icol,irow), icol=1, ncol)
!!  end do
!!  close(10)
  ! (4) Create reach vectors:
  !     Walk downhill from the start points again and register the coordinates.
  !     Where the number of paths / cell changes (i.e. at junctions):
  !     - let the current path end
  !     - generate a start point for another path
  !     After this, update the start points and repeat until no more start points exist.
  !  In addition, identify the reach being the system's outlet.
  if (silent .eq. 0) call intpr("Assembling reach vectors...",-1, 0, 0)      
  !Initialize output vectors
  ncoords= 0
  path_id= 0
  path_x= ZERO
  path_y= ZERO
  !Create vectors
  id_path= 0
  id_outlet= 0  
  matrix_start_used= .false.
  where (matrix_start_current) matrix_start_used= .true.
  do while (count(matrix_start_current) > 0)
    matrix_start_new= .false.
    do icol_start=1,ncol
      do irow_start=1,nrow
        if (matrix_start_current(icol_start,irow_start)) then
          if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. &
              (matrix_flowdir(icol_start,irow_start) .le. 8)) then
            !Walk along the flow paths
            icol= icol_start
            irow= irow_start
            id_path= id_path + 1
            npaths_start= matrix_npaths(icol,irow)
            !Register start coordinates
            ncoords= ncoords + 1
            path_id(ncoords)= id_path
            path_x(ncoords)= cellXavg(icol,xllcorner,cellsize)
            path_y(ncoords)= cellYavg(irow,nrow,yllcorner,cellsize)
            do
              tmp1= icol + colshift(matrix_flowdir(icol,irow))
              tmp2= irow + rowshift(matrix_flowdir(icol,irow))
              icol= tmp1
              irow= tmp2
              !Register coordinates
              ncoords= ncoords + 1
              path_id(ncoords)= id_path
              path_x(ncoords)= cellXavg(icol,xllcorner,cellsize)
              path_y(ncoords)= cellYavg(irow,nrow,yllcorner,cellsize)
              if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
                exit   ! Stop at grid border
              else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
                exit   ! Stop at nodata or invalid codes (e.g. code 0)
              else
                !At a junction: Register new start point (if not used) and let path end here
                if (matrix_npaths(icol,irow) .gt. npaths_start) then  ! Are we at a junction?
                  if (.not. matrix_start_used(icol,irow)) then        ! Does the startpoint already exist?
                    matrix_start_new(icol,irow)= .true.
                    matrix_start_used(icol,irow)= .true.
                  end if
                  exit
                end if
              end if
            end do
            !Register outlet ID
            if ((icol .eq. icol_end) .and. (irow .eq. irow_end)) id_outlet= id_path
          end if
        end if
      end do
    end do
    !Update start points
    matrix_start_current= matrix_start_new
  end do
  !Delete unnecessary interior points (i.e. points in the middle of straight lines) from the flow paths
  if (silent .eq. 0) call intpr("Simplifying reach vectors...",-1, 0, 0)      
  !Step 1: Find points to be deleted
  path_keep= .false.
  path_keep(1:ncoords)= .true.
  if (ncoords .gt. 2) then
    do i=2, (ncoords-1)
      !If points belong to the same path
      if ((path_id(i-1) .eq. path_id(i)) .and. (path_id(i+1) .eq. path_id(i))) then
        dx1= path_x(i) - path_x(i-1)
        dy1= path_y(i) - path_y(i-1)
        dx2= path_x(i+1) - path_x(i)
        dy2= path_y(i+1) - path_y(i)         
        !If points are on a line parallel to x- or y-axis or points are on a diagonal
        if (((dx1 .eq. ZERO) .and. (dx2 .eq. ZERO)) .or. &
            ((dy1 .eq. ZERO) .and. (dy2 .eq. ZERO)) .or. &
            ((dx1 .eq. dx2) .and. (dy1 .eq. dy2))) then
          path_keep(i)= .false.
        end if
      end if
    end do    
  end if
  !Step 2: Remove points from result arrays
  path_id= pack(path_id, mask=path_keep, vector=(/(0,k=1,size(path_id))/))
  path_x= pack(path_x, mask=path_keep, vector=(/(ZERO,k=1,size(path_id))/))
  path_y= pack(path_y, mask=path_keep, vector=(/(ZERO,k=1,size(path_id))/))
  ncoords= count(path_keep)
  !Set error level
  errorlevel= 0; return
  1 errorlevel= 1; return
  2 errorlevel= 2; return
  3 errorlevel= 3; return
end subroutine flowpaths

!-------------------------------------------------------------------------------
! Compute a concentration time index for each cell based on Manning's Eqn.
! 
!   Manning's eqn.:     u = 1/n * sqrt(I) * R^(2/3)
!
!             with      u: Velocity
!                       n: Roughness
!                       I: Energy gradient
!                       R: Hydraulic radius (== width if flow height << width)
!
!   Definition of velocity:   u = ds / dt
!
!                     with    ds: Length of flow path
!                             dt: Travel time
!
!   Concentration time index:  cti = ds / srqt(I)
!
!   For each cell, the cti is the sum of the ctis of the DOWNstream cells
!   in the respective flow path.
!
!   --> Value of cti is high for (1) long flow paths and (2) low gradients
!
subroutine conctimeindex(ncol, nrow, vect_flowdir, nodata_flowdir, &
  vect_dem, nodata_dem, vect_flowacc, cellsize, area_crit, dz_min, silent, &
  vect_out, errorlevel)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  integer,dimension(ncol*nrow),intent(in):: vect_flowdir
  integer,intent(in):: nodata_flowdir
  real(8),dimension(ncol*nrow),intent(in):: vect_dem
  real(8),intent(in):: nodata_dem
  integer,dimension(ncol*nrow),intent(in):: vect_flowacc
  real(8),intent(in):: cellsize
  real(8),intent(in):: area_crit
  real(8),intent(in):: dz_min
  integer,intent(in):: silent
  !output
  real(8),dimension(ncol*nrow),intent(out):: vect_out
  integer,intent(out):: errorlevel
  !local
  real(8),parameter:: ZERO= 0.0d0, TWO= 2.0d0, R2= sqrt(2.0d0), UN= 1.0d0, NODATA_OUT=-9999.d0
  integer,dimension(8),parameter:: colshift= (/-1,0,1,1,1,0,-1,-1/)
  integer,dimension(8),parameter:: rowshift= (/-1,-1,-1,0,1,1,1,0/)
  real(8),dimension(8):: distances
  integer:: icol,irow, icol_start, irow_start, tmp1, tmp2
  integer,dimension(1:ncol,1:nrow):: matrix_flowdir, matrix_flowacc
  real(8),dimension(1:ncol,1:nrow):: matrix_dem, matrix_out
  logical,dimension(1:ncol,1:nrow):: matrix_rivnet
  real(8):: dst
  real(8):: z1, z2, dz
  integer:: ncells_crit
  real(8):: cti
  !code
  distances= (/R2,UN,R2,UN,R2,UN,R2,UN/)
  !Check args
  ! (a) Critical number of cells must be >= 1
  ncells_crit= max(1, nint(area_crit / (cellsize**TWO)))
  if (ncells_crit .lt. 1) then
    goto 1
  end if
  ! (b) Lower limit of elevation difference must be positive
  if (dz_min .le. ZERO) then
    goto 2
  end if
  ! (c) Cells with valid flow direction must also have a valid elevation
  if (any((matrix_flowdir .ne. nodata_flowdir) .and. (matrix_dem .eq. nodata_dem))) then
    goto 3
  end if
  !Transform input vectors (concatenated columns) into matrices
  do icol=1,ncol
    matrix_flowdir(icol,:)= vect_flowdir( (icol-1)*nrow+1 : icol*nrow ) 
    matrix_dem(icol,:)= vect_dem( (icol-1)*nrow+1 : icol*nrow ) 
    matrix_flowacc(icol,:)= vect_flowacc( (icol-1)*nrow+1 : icol*nrow ) 
  end do
  ! (1) Compute flow accumulation grid from a grid of flow direction codes.
  ! --> Has become a subroutine input: Range of values is 0....Inf
  ! (2) Select river cells
  where (matrix_flowacc .ge. ncells_crit)
    matrix_rivnet= .true.
  else where
    matrix_rivnet= .false.
  end where
  ! (3) Compute the concentration time index (cti) for each flow path and
  !     assign the value to the start cell. Flow paths end at river cells.
  if (silent .eq. 0) call intpr("Computing CT index...",-1, 0, 0)      
  matrix_out= NODATA_OUT
  do icol_start=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol_start), dble(1))
    do irow_start=1,nrow
      if ((matrix_flowdir(icol_start,irow_start) .ge. 1) .and. (matrix_flowdir(icol_start,irow_start) .le. 8)) then
        !Initialize (Assume at least a flow distance of half a cell. This is
        !better than assuming ZERO. Using ZERO might lead to numerical problems
        !when used in a model.)
        cti= 0.5*cellsize / sqrt(dz_min / (0.5*cellsize))
        if (.not. matrix_rivnet(icol_start,irow_start)) then
          !Walk along the flow paths
          icol= icol_start
          irow= irow_start
          do
            dst= distances(matrix_flowdir(icol,irow))
            z1= matrix_dem(icol,irow)
            tmp1= icol + colshift(matrix_flowdir(icol,irow))
            tmp2= irow + rowshift(matrix_flowdir(icol,irow))
            icol= tmp1
            irow= tmp2
            z2= matrix_dem(icol,irow)
            if ((icol .lt. 1) .or. (icol .gt. ncol) .or. (irow .lt. 1) .or. (irow .gt. nrow)) then
              exit   ! Stop at grid border
            else if ((matrix_flowdir(icol,irow) .lt. 1) .or. (matrix_flowdir(icol,irow) .gt. 8)) then
              exit   ! Stop at nodata or invalid codes (e.g. code 0)
            else
              dz= max(dz_min, (z1 - z2))
              cti= cti + dst*cellsize / sqrt(dz / (dst*cellsize))
              if (matrix_rivnet(icol,irow)) then
                exit     !End of flow path encountered
              end if
            end if
          end do
        end if
        matrix_out(icol_start,irow_start)= cti
      end if
    end do
  end do
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix_out(icol, 1:nrow)
  end do
  !Set error level
  errorlevel= 0; return
  1 errorlevel= 1; return
  2 errorlevel= 2; return
  3 errorlevel= 3; return
end subroutine conctimeindex

!-------------------------------------------------------------------------------
!Convert lines to a grid.
!The lines are defined by the 3 vectors: x, y, ID. The ID values are integers.
!Grid cells not touched by any line are set to 'code_undefined'.
!Grid cells touched by multiple lines are set to 'code_conflict'.
subroutine lines2raster(ncol, nrow, xllcorner, yllcorner, cellsize, &
  code_undefined, code_conflict, nbuffer, silent, &
  length, x, y, id, vect_out)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol, nrow
  real(8),intent(in):: xllcorner, yllcorner
  real(8),intent(in):: cellsize
  integer,intent(in):: code_undefined
  integer,intent(in):: code_conflict
  integer,intent(in):: nbuffer
  integer,intent(in):: silent
  integer,intent(in):: length
  real(8),dimension(length),intent(in):: x
  real(8),dimension(length),intent(in):: y
  integer,dimension(length),intent(in):: id 
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_out
  !local
  real(8),parameter:: FACT= 0.5d0
  integer,dimension(8):: near
  integer:: conflict_marker
  integer,dimension(1:ncol,1:nrow):: matrix
  integer,dimension(1:ncol,1:nrow):: tmp
  real(8):: dx, dy, xsub, ysub
  integer:: nsub
  integer:: i, k
  integer:: icol,irow
  !code
  conflict_marker= minval(id) - 1
  !Init matrix to nodata
  matrix= code_undefined
  !Loop through the coordinates of the lines
  if (silent .eq. 0) call intpr("Gridding...",-1, 0, 0)      
  do i=1, (length-1)
    if (silent .eq. 0) call progress(dble(1), dble(length-1), dble(i), dble(1))
    !Process a segment of a particular line (id changes where new line starts)
    if (id(i+1) .eq. id(i)) then
      !Define sequence of points along the segment with sufficient density
      nsub= int(dist(x(i),y(i),x(i+1),y(i+1)) / (FACT * cellsize)) + 1
      dx= (x(i+1)-x(i)) / dble(nsub)
      dy= (y(i+1)-y(i)) / dble(nsub)
      !Identify cells corresonding to the points
      do k=0,nsub
        xsub= x(i) + dble(k) * dx
        ysub= y(i) + dble(k) * dy
        !Find col/row indices (zero if outside grid extent)
        call whichCell(ncol, nrow, xllcorner, yllcorner, cellsize, &
          xsub, ysub, icol, irow)
        if ((icol * irow) .ne. 0) then
          !If the cell is virgin, assign value
          if (matrix(icol,irow) .eq. code_undefined) then
            matrix(icol,irow)= id(i)
          !If the cell already has this id, do nothing
          else if (matrix(icol,irow) .eq. id(i)) then
            continue
          !If the cell already has another ID, remove this and mark the conflict
          else
            matrix(icol,irow)= conflict_marker
          end if
        end if 
      end do     
    end if
  end do
  !Create buffer if requested
  if (nbuffer .gt. 0) then
    if (silent .eq. 0) call intpr("Creating buffers...",-1, 0, 0)      
    do i=1, nbuffer
      tmp= code_undefined
      do icol=1,ncol
        do irow=1,nrow
          if (matrix(icol,irow) .eq. code_undefined) then
            near= adjacentCells(matrix,icol,irow,code_undefined)
            do k=1,8
              if (near(k) .eq. conflict_marker) then
                tmp(icol,irow)= conflict_marker
                exit
              end if
              if (near(k) .ne. code_undefined) then
                if (tmp(icol,irow) .eq. code_undefined) then
                  tmp(icol,irow)= near(k)
                else if (tmp(icol,irow) .eq. near(k)) then
                  continue
                else
                  tmp(icol,irow)= conflict_marker
                  exit
                end if
              end if
            end do
          end if
        end do
      end do
      where (tmp .ne. code_undefined) matrix= tmp
    end do
  end if
  !Replace conflict marker by desired code
  where (matrix .eq. conflict_marker) matrix= code_conflict
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix(icol, 1:nrow)
  end do
end subroutine lines2raster

!-------------------------------------------------------------------------------
! Creates catchments based on
! - an initial matrix
! - a matrix of flow directions
! These two matrices must be of the same dimensions.
subroutine buildcatchments(ncol, nrow, nodata_init, conflict_init, vect_init, &
  nodata_fldir, vect_fldir, undef_catch, vect_catch)
implicit none
  !input
  integer,intent(in):: ncol
  integer,intent(in):: nrow 
  integer,dimension(ncol*nrow),intent(in):: vect_init, vect_fldir
  integer,intent(in):: conflict_init
  integer,intent(in):: nodata_init, nodata_fldir
  integer,intent(in):: undef_catch
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_catch
  !local
  integer,dimension(1:ncol,1:nrow):: matrix_fldir
  integer,dimension(1:ncol,1:nrow):: matrix_catch
  integer,dimension(8),parameter:: colshift= (/-1,0,1,1,1,0,-1,-1/)
  integer,dimension(8),parameter:: rowshift= (/-1,-1,-1,0,1,1,1,0/)
  integer:: UNDEF_FLDIR
  integer:: icol,irow
  integer:: dir, catchID
  logical:: change
  !code
  !Convert flow direction vector to matrix
  do icol=1,ncol
    matrix_fldir(icol,:)= vect_fldir( (icol-1)*nrow+1 : icol*nrow ) 
  end do
  !Identify cells with unknown flow direction (neither in 1...8, nor nodata)
  UNDEF_FLDIR= nodata_fldir - 1
  where ((matrix_fldir .lt. 1) .and. (matrix_fldir .gt. 8) .and. &
    (matrix_fldir .ne. nodata_fldir))
    matrix_fldir= UNDEF_FLDIR
  end where
  !Initialize catchment matrix
  do icol=1,ncol
    matrix_catch(icol,:)= vect_init( (icol-1)*nrow+1 : icol*nrow ) 
  end do
  !Recursively find new linked cells that belong to a catchment based on the flow direction
  change= .true.
  do while (change)
    change= .false.
    do icol=1,ncol
      do irow=1,nrow
        !Set only cells that haven't been set in an earlier pass
        if (matrix_catch(icol,irow) .eq. nodata_init) then
          !Set only cells that are linked to another one
          dir= matrix_fldir(icol,irow)
          if ((dir .ne. nodata_fldir) .and. (dir .ne. UNDEF_FLDIR)) then
            !Assign new cell only if linked cell was already set
            catchID= matrix_catch(icol+colshift(dir),irow+rowshift(dir))
            if (catchID .ne. nodata_init) then
              matrix_catch(icol,irow)= catchID
              change= .true.
            end if
          end if
        end if
      end do
    end do
  end do
  !Set undefined cells 1: Cells with undefined low direction
  where (matrix_fldir .eq. UNDEF_FLDIR)
    matrix_catch= undef_catch
  end where
  !Set undefined cells 2: Catchment of cells marked with the conflict code in the init grid
  where (matrix_catch .eq. conflict_init)
    matrix_catch= undef_catch
  end where
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_catch( (icol-1)*nrow+1 : icol*nrow )= matrix_catch(icol, 1:nrow)
  end do
end subroutine buildcatchments

!-------------------------------------------------------------------------------
!Fill grid cells having a certain value by a nearest neighbor approach
subroutine nnfill(ncol,nrow,cellsize,nodata,vect_in,removevalue,maxdist,silent,vect_out)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol,nrow
  real(8),intent(in):: cellsize
  integer,intent(in):: nodata
  integer,intent(in):: removevalue
  integer,dimension(ncol*nrow),intent(in):: vect_in
  real(8),intent(in):: maxdist
  integer,intent(in):: silent
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_out
  !local
  integer,dimension(1:ncol,1:nrow):: matrix_in
  integer,dimension(1:ncol,1:nrow):: matrix_out
  integer:: icol,irow,col_neighbor,row_neighbor
  !code
  !Convert input vector to matrix
  do icol=1,ncol
    matrix_in(icol,:)= vect_in( (icol-1)*nrow+1 : icol*nrow )
  end do
  !Copy input to output
  matrix_out= matrix_in
  !Fill
  do icol=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol), dble(1))
    do irow=1,nrow
      if (matrix_in(icol,irow) .eq. removevalue) then
        call findNN(ncol,nrow,matrix_in,cellsize,nodata,icol,irow,maxdist, &
          range_lbnd=removevalue,range_ubnd=removevalue,in_range=.false., &
          col_neighbor=col_neighbor,row_neighbor=row_neighbor)
        if ((row_neighbor * col_neighbor) .ne. 0) then
          matrix_out(icol,irow)= matrix_in(col_neighbor,row_neighbor)
        end if
      end if    
    end do
  end do
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix_out(icol, 1:nrow)
  end do
end subroutine nnfill

!-------------------------------------------------------------------------------
!Returns the cell values for a set of points
subroutine gridvaluesatpoints(ncol, nrow, xllcorner, yllcorner, cellsize, &
  nodata, vect_grid, npoints, x, y, silent, values)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol, nrow
  real(8),intent(in):: xllcorner, yllcorner, cellsize, nodata
  real(8),dimension(ncol*nrow),intent(in):: vect_grid
  integer,intent(in):: npoints
  real(8),dimension(npoints),intent(in):: x, y
  integer,intent(in):: silent
  !output
  real(8),dimension(npoints),intent(out):: values
  !local
  integer:: icol
  integer:: colind, rowind
  integer:: k
  real(8),dimension(1:ncol,1:nrow):: matrix_grid
  !code
  !Transform input vector (concatenated columns) into matrix
  do icol=1,ncol
    matrix_grid(icol,:)= vect_grid( (icol-1)*nrow+1 : icol*nrow )
  end do
  !Determine grid values for all point cordinates
  do k=1,npoints
    if (silent .eq. 0) call progress(dble(1), dble(npoints), dble(k), dble(1))
    call whichCell(ncol,nrow,xllcorner,yllcorner,cellsize,x(k),y(k),colind,rowind)
    if ((colind*rowind) .eq. 0) then
      values(k)= nodata
    else
      values(k)= matrix_grid(colind,rowind)
    end if
  end do
end subroutine gridvaluesatpoints

!-------------------------------------------------------------------------------
!Returns a grid where the cell values represent the index of the nearest point in a set of points
subroutine gridnearestpoints(ncol, nrow, xllcorner, yllcorner, cellsize, &
  npoints, x, y, silent, vect_out)
use topocatchInternals
implicit none
  !input
  integer,intent(in):: ncol, nrow
  real(8),intent(in):: xllcorner, yllcorner, cellsize
  integer,intent(in):: npoints
  real(8),dimension(npoints),intent(in):: x, y
  integer,intent(in):: silent
  !output
  integer,dimension(ncol*nrow),intent(out):: vect_out
  !local
  integer:: icol, irow
  integer:: k
  integer,dimension(1:ncol,1:nrow):: matrix_out
  real(8):: d, dmin
  !code
  matrix_out= 0
  ! Find the nearest point for each cell
  do icol=1,ncol
    if (silent .eq. 0) call progress(dble(1), dble(ncol), dble(icol), dble(1))
    do irow=1,nrow
      dmin= huge(dmin)
      do k=1,npoints
	d= dist(x(k),y(k),cellXavg(icol,xllcorner,cellsize), &
	  cellYavg(irow,nrow,yllcorner,cellsize))
	if (d .lt. dmin) then
          dmin= d
          matrix_out(icol,irow)= k
	end if
      end do
    end do
  end do
  !Transform result matrix into output vector (of concatenated colums)
  do icol=1,ncol
    vect_out( (icol-1)*nrow+1 : icol*nrow )= matrix_out(icol, 1:nrow)
  end do
end subroutine gridnearestpoints
