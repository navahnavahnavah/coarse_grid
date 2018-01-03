MODULE globals
  IMPLICIT NONE
  SAVE



  ! JDF PARAMS WITH FRACTURE AUGUST 9
  INTEGER, PARAMETER :: testInt = 31, xn =1801, yn = 52, altnum = 136,  cell = 1 !50000
  INTEGER, PARAMETER :: cellx = 90, celly = 1
  INTEGER, PARAMETER ::  g_pri = 5, g_sec = 80, g_sol = 15, g_med = 7, g_iso = 2, cstep = 125, ar = 800 ! cstep = 1000
  INTEGER, PARAMETER :: tn = 1000000, mstep = 25, wscale = 1, ison = 10, inertn = 10! ison = 10000, inertn = 100000
  INTEGER, PARAMETER :: write_factor = 25, res_factor = 300 ! res factor is how many flushes happen within t_max
  INTEGER :: active_cells
  INTEGER, PARAMETER :: particle_sat = 1, inert_sat = 10
  REAL(4) :: tb_res = 1.0e10
  REAL(4) :: cstep_int
  INTEGER :: cstep_num
  REAL(4) :: u_1d
  REAL(4) :: x_min = 0.0, x_max = 90000.0, y_min = -1250.0, y_max = 25.0 !y_min = -1650.0, y_max = 25.0
  REAL(4) :: t_min = 0.0, t_max = 1.02e14!23.55e13 !9.42e13 !
  REAL(4) :: ki = .76, ra = 100.0, viscosity = .001, cp = 1173.0, alpha =4.0e-4, k, calc0, psi_round, psi_round2 !alpha =4.0e-5 !cp = 1175.0
  INTEGER :: thresh=0, theta0
  REAL(4) :: in_left1, in_left2, in_right1, in_right2, in_left3, in_right3, max_left, max_right, fac_left, fac_right
  REAL(4) :: scope! = -5.0e-10
  REAL(4) :: scope1, bh5, bh10, ah5, ah10, psi50(xn), f50a(xn), f50b(xn), f50c(xn), hflat(xn), maskflat(xn)
  REAL(4) :: first = 0.56, factor = -.000206!-10.0
  REAL(4) :: psi5(yn), f5a(yn), f5b(yn), f5c(yn)
  REAL(4) :: psi10(yn), f10a(yn), f10b(yn), f10c(yn)
  REAL(4) :: top5, top10, bottom5, bottom10, sumplace
  REAL(4) :: ksed = 1e-16
  REAL(4) :: fix_b, dpd, dt_bit
  REAL(4) :: psi_bl, psi_br
  REAL(4) :: permf, phi_calc_denom
  REAL(4) :: h_base, y_base, h_top, y_top, h_adjacent
  INTEGER :: jj_base, jj_top, leng, i_count
  LOGICAL, DIMENSION(cellx,celly) :: i_mask
  INTEGER :: th_bool(.TRUE.:.FALSE.)
  LOGICAL :: sig_bool_a, sig_bool_b
  INTEGER :: tag

  REAL(4) :: dt, dx, dy, dt0 = 0.001
  REAL(4) :: dPsi, psiLast(xn,yn)
  INTEGER :: loop
  REAL(4) :: lambda = 1.14
  REAL(4) :: grav = 9.8
  REAL(4) :: rho_fluid = 1000.0
  CHARACTER(len=300) :: s_i
  REAL(4) :: rho_total=2530.0
  CHARACTER(len=20) :: kinetics!=" precipitate_only"
  CHARACTER(len=25) :: s_pressure = "250.0"
  CHARACTER(len=25) :: si_hematite
  INTEGER :: i_unit, code
  REAL(4) :: something, refL, refR
  INTEGER :: vfe1=24, vfe2 = 16
  INTEGER :: diss_toggle
  REAL(4) :: precip_th=0.0
  !CHARACTER(len=25) :: ph_fix = "10.0^-7.5"
  REAL(4) :: ph_count
  REAL(4) :: ph_sum
  INTEGER :: m_count
  INTEGER :: par_rounds, par_round, par_current, last_active
  INTEGER :: active_coarse, end_loop


  INTEGER :: slave_vector((3*(xn-1)/cellx)*(yn/(2*celly)))
  INTEGER :: slave_count, index_count
  INTEGER :: dabs_loop

  INTEGER :: se_toggle
  INTEGER :: se_factor = 100
  INTEGER :: j_root
  INTEGER :: se_spinup = 100
  INTEGER :: age_toggle
  INTEGER :: age_root
  CHARACTER(len=25) :: my_id_s
  REAL(4) :: vol_th
  CHARACTER(len=25) :: s_precip_nope= "100.0"
  CHARACTER(len=25) :: fixed_volume
  CHARACTER(len=25) :: grow_bit
  CHARACTER(len=25) :: surf_scale_1_string, surf_scale_2_string, surf_scale_gen_string
  REAL(4) :: surf_scale_1, surf_scale_2
  REAL(4) :: glob_t_vol_s = 0.008
  REAL(4) :: glob_t_vol_a = 0.0075
  REAL(4) :: glob_t_vol_b = 0.0005

CONTAINS







  ! ----------------------------------------------------------------------------------%%
  !
  ! CHECK SOMETHING
  !
  ! ----------------------------------------------------------------------------------%%

  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    !if(status /= nf90_noerr) then
    !print *, trim(nf90_strerror(status))
    !stop "Stopped"
    !end if

  END SUBROUTINE check
  ! ----------------------------------------------------------------------------------%%
  !
  ! BAND
  !
  ! SUMMARY: Gaussian elimination for the banded case, when you know the structure of
  ! 		   the band, which you usually do
  !
  ! INPUTS: a(n,m) : original matrix, but just the band; n rows deep, m columns wide
  !		  		   relation between matrix A and band a is A(i,j) = a(i,j-i+(m+1)/2)
  !		  m : width of the band (MUST BE ODD! no 4-wide bands. makes sense)
  !		  n : order of the matrix A or # of rows in a
  !
  !
  ! RETURNS: band(n,m) : which contains upper triangular and lower triangular in the
  !       			   columns that are not the central column
  !
  ! ----------------------------------------------------------------------------------%%


  FUNCTION band(a,m,n)

    IMPLICIT NONE
    INTEGER g,h,i,j,k,m,n,r
    REAL(4) :: a(n,m), band(n,m)
    REAL eps

    WRITE(*,*) m !57 = 2*(yn/2 - 2) + 1 = 2*28 + 1 = 56 + 1
    !
    WRITE(*,*) n !3304 = 28 * 118

    band = a
    r = (m+1)/2
    eps = 1.0e-6

    WRITE(*,*) "made it?"



    DO 20 k = 1,n
       ! 	write(*,*) "row", k
       ! 	write(*,*) sum(band(k,:))
       IF (ABS(band(k,r)) .LE. eps) GOTO 99
       band(k,r) = 1.0/band(k,r)
       h = r-1
       i = k+1
10     IF (h .LT. 1 .OR. i .GT. n) GOTO 20
       band(i,h) = band(i,h) * band(k,r)
       j = h+1
       g = r+1
30     IF (g .GT. m .OR. j .GT. (r+n-i)) GOTO 40
       band(i,j) = band(i,j) - band(i,h)*band(k,g)
       j = j+1
       g = g+1
       GOTO 30

40     CONTINUE
       i = i+1
       h = h-1
       GOTO 10
20     CONTINUE
       RETURN


       WRITE(*,*) "made it??"

99     m=0
       RETURN
     END FUNCTION band

     ! subroutine band(a,m,n)
     !
     ! implicit none
     ! integer :: g,h,i,j,k,m,n,r
     ! real(4) :: a(n,m)!, band(n,m)
     ! real(4) :: eps
     !
     !
     !
     !
     ! ! do 20 k = 1,n
     ! ! 	10 if (abs(band(k,r)) .le. eps) then
     ! ! 			go to 99
     ! ! 		else
     ! ! 			band(k,r) = 1.0/band(k,r)
     ! ! 				h = r-1
     ! ! 				i = k+1
     ! ! 		end if
     ! !
     ! !
     ! ! 	if (h .lt. 1 .or. i .gt. n) then
     ! ! 		go to 20
     ! ! 	else
     ! ! 		band(i,h) = band(i,h) * band(k,r)
     ! ! 		j = h+1
     ! ! 		g = r+1
     ! ! 	end if
     ! !
     ! !
     ! !
     ! ! do while (g .le. m .and. j .le. (r+n-i))
     ! ! 	band(i,j) = band(i,j) - band(i,h)*band(k,g)
     ! ! 	j = j+1
     ! ! 	g = g+1
     ! ! end do
     ! !
     ! ! do while (g .gt. m .and. j .gt. (r+n-i))
     ! ! 	continue
     ! ! end do
     ! !
     ! !
     ! ! !if (g .gt. m .or. j .gt. (r+n-i)) then
     ! ! 	i = i+1
     ! ! 	h = h-1
     ! ! 	go to 10
     ! ! 20 continue
     ! ! return
     ! !
     ! ! 99 m=0
     ! ! return
     !
     !
     ! !band = a
     ! r = (m+1)/2
     ! eps = 1.0e-6
     !
     !
     ! do 20 k = 1,n
     ! if (abs(a(k,r)) .le. eps) goto 99
     ! a(k,r) = 1.0/a(k,r)
     ! 	h = r-1
     ! 	i = k+1
     !
     ! 10 	if (h .lt. 1 .or. i .gt. n) goto 20
     ! 	a(i,h) = a(i,h) * a(k,r)
     ! 		j = h+1
     ! 		g = r+1
     !
     ! 30		if (g .gt. m .or. j .gt. (r+n-i)) goto 40
     ! 		a(i,j) = a(i,j) - a(i,h)*a(k,g)
     ! 			j = j+1
     ! 			g = g+1
     ! 			goto 30
     !
     ! 40 continue
     ! 	i = i+1
     ! 	h = h-1
     ! 	goto 10
     !
     ! 20 continue
     ! !return
     !
     ! 99 return
     ! !m = 0
     ! !write(*,*) "got here?"
     ! !return
     !
     ! end subroutine band




     FUNCTION solve(a,b,m,n)

       IMPLICIT NONE
       INTEGER i,j,k,m,n,r
       REAL(4) :: a(n,m), b(n), solve(n)

       solve = b
       r = (m+1)/2

       DO 100 k=1,n-1
          i = k+1
          j = r-1
110       IF (j .LT. 1 .OR. i .GT. n) GOTO 100
          solve(i) = solve(i) - a(i,j)*solve(k)
          i = i+1
          j = j-1
          GOTO 110
100       CONTINUE

          DO 120 k=n,1,-1
             i = k+1
             j = r+1
130          IF (j .GT. m .OR. i .GT. n) GOTO 140
             solve(k) = solve(k) - a(k,j)*solve(i)
             i = i+1
             j=j+1
             GOTO 130

140          CONTINUE
             solve(k) = solve(k)*a(k,r)

120          CONTINUE
             RETURN
           END FUNCTION

           ! ----------------------------------------------------------------------------------%%
           !
           ! LINSPACE
           !
           ! SUMMARY: Makes a vector of dimension n, starting with a_first, ending with a_last
           !          with n-2 evenly spaced items along the way.
           !
           ! INPUTS: n : number of elements
           !         a_first : value of first element
           !         a_last : value of last element
           !
           ! RETURNS: linspace(a)
           !
           ! ----------------------------------------------------------------------------------%%

           FUNCTION linspace ( n, a_first, a_last )

             IMPLICIT NONE
             INTEGER, INTENT(in) :: n
             REAL (4) :: linspace(n)
             INTEGER :: i
             REAL (4) :: a(n)
             REAL (4) , INTENT(in) ::  a_first, a_last

             IF ( n == 1 ) THEN
                a(1) = ( a_first + a_last ) / 2.0
             ELSE
                DO i = 1, n
                   a(i) = ( REAL ( n - i,     kind = 8 ) * a_first &
                        + REAL (     i - 1, kind = 8 ) * a_last ) &
                        / REAL ( n     - 1, kind = 8 )
                END DO
             END IF
             linspace = a
             RETURN

           END FUNCTION linspace







           ! ----------------------------------------------------------------------------------%%
           !
           ! FOUND GAUSSIAN ELIMINATION FUNCTION
           !
           ! SUMMARY: Solve Ax = b style linear system using gaussian elimination
           !
           ! INPUTS: a : A concatenated to b
           !         row : n-by-n dimension
           !
           ! RETURNS: Gsselm(row)
           !
           ! ----------------------------------------------------------------------------------%%

           FUNCTION Gsselm(a,row)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: row
             REAL(4) , INTENT(IN OUT)  ::  a(:,:)   	!Assume shape (:)
             REAL(4) , DIMENSION(row) :: Gsselm

             INTEGER i,j,k
             INTEGER, DIMENSION(2) :: shap
             REAL(4) , ALLOCATABLE :: swap_ik(:)
             REAL(4)  :: tmp

             ALLOCATE (swap_ik(row+1))

             ! 	Initialise
             swap_ik(:) = 0.0d0            ! Whole vector initialized to zero
             tmp = 0.0d0

             ! Check dimensions of input matrix
             shap = SHAPE(a)
             IF ( (shap(1) .NE. row) .OR.  (shap(2) .NE. row+1) ) THEN
                !  	call break()
             END IF


             !/*   Gaussian Elimination - Row Reduction of matrix  */
             DO k=1, row-1                             ! total of row-1 operations

                !/*  Pivotal strategy - SWAP rows to make pivotal element a[k][k] have the
                !    greatest magnitude in its column. This prevents unnecessary division by
                !    a small number.                           */
                !   	do i = k+1, row
                !      	if ( (dabs(a(i,k))-dabs(a(k,k))).gt.eps  ) then
                !         	do j = k, row+1                     !/* If pivotal element is not */
                !            	swap_ik(j) = a(k,j)              !/* the highest then  */
                !          	   a(k,j) = a(i,j)                  !/* swap i'th and k'th rows */
                !            	a(i,j) = swap_ik(j)
                !         	end do 		!j-loop
                !         end if
                !   	end do 				!i-loop


                !/*   If the Matrix is SINGULAR then EXIT program      */
                !  		IF ( dabs(a(k,k)) < EPS ) then
                !   		print *,'After swapping rows to make the pivotal element become the'
                !      	print *,'highest magnitude element in its column, its magnitude is'
                !	      print *,'still extremely small.'
                !	      print *,'Hence this is a SINGULAR MATRIX - no unique solution or '
                !	      print *,'check that the input dimensions are correct.'
                !   	   call break()
                !	   END if



                !/*      Perform row-reduction with pivotal element a[k][k]     */
		DO i = k+1, row
                   DO j = row+1, k, -1			!/* starting from end of column */
                      a(i,j) = a(i,j) - a(k,j) / a(k,k) * a(i,k)
                   END DO 							!/* end of j loop     */
		END DO 	 							!/* end of 2nd i loop */

             END DO 									!/* end of k loop     */
             !  At this point, the bottom triangle is Zero


             !/*   Back Substitution - Solutions of equations   */
             Gsselm(row) = a(row,row+1) / a(row,row)
             DO k = row-1, 1, -1
                tmp = 0.0d0
                DO j = k+1, row
                   tmp = tmp + a(k,j)*Gsselm(j)
                END DO 							!j-loop
                Gsselm(k) = ( a(k,row+1) - tmp ) / a(k,k)
             END DO 								!k-loop

             DEALLOCATE (swap_ik)
             RETURN
           END FUNCTION Gsselm



           ! ----------------------------------------------------------------------------------%%
           !
           ! TRIDIAG
           !
           ! SUMMARY: Solves Ax = b style linear system where A is tridiagonal matrix
           !
           ! INPUTS: aa : tridiagonal matrix
           !         nn : nn-by-nn dimension
           !
           ! RETURNS: tridiag(nn)
           !
           ! ----------------------------------------------------------------------------------%%





           FUNCTION tridiag(a,b,c,d,nn)
             INTEGER :: nn, k, km1, i
             REAL(4) :: aa(nn,nn), a(nn), b(nn), c(nn), d(nn), tridiag(nn), xm
             !      dimension a(nn),b(nn),c(nn),d(nn)
             !tridiag = aa(:,nn+1)
             tridiag = d
             a(1) = 0.0
             !b(1) = aa(1,1)
             !c(1) = aa(1,2)
             !a(nn) = aa(nn,nn-1)
             !b(nn) = aa(nn,nn)
             c(nn) = 0.0
             !do k = 2,nn-1
             !	a(k) = aa(k,k-1)
             !	b(k) = aa(k,k)
             !	c(k) = aa(k,k+1)
             !end do


             !      if(nn .eq. 1) then
             !      tridiag(1)=d(1)/b(1)
             !      write(*,*) "shit"
             !      end if

             DO k = 2,nn
                km1 = k - 1

                IF(b(k-1) .EQ. 0.0) THEN
                   WRITE(*,*) "what"
                END IF

                xm  = a(k)/b(km1)
                b(k)  = b(k) - xm*c(km1)
                tridiag(k)  = tridiag(k) - xm*tridiag(km1)
             END DO

             tridiag(nn)   = tridiag(nn)/b(nn)

             k = nn
             DO i = 2,nn
                k = nn + 1 - i
                tridiag(k) = (tridiag(k) - c(k)*tridiag(k+1))/b(k)
             END DO

             !      format(/3x,'diagonal element .eq. 0 in tridag at k = ',i2/)
             RETURN
           END FUNCTION tridiag











           ! ----------------------------------------------------------------------------------%%
           !
           ! TIMESTAMP
           !
           ! SUMMARY: Write a timestamp
           !
           ! RETURNS: timestamp
           !
           ! ----------------------------------------------------------------------------------%%

           SUBROUTINE timestamp ( )

             IMPLICIT NONE
             CHARACTER ( len = 8 ) ampm
             INTEGER ( kind = 4 ) d
             INTEGER ( kind = 4 ) h
             INTEGER ( kind = 4 ) m
             INTEGER ( kind = 4 ) mm
             CHARACTER ( len = 9 ), PARAMETER, DIMENSION(12) :: month = (/ &
                  'January  ', 'February ', 'March    ', 'April    ', &
                  'May      ', 'June     ', 'July     ', 'August   ', &
                  'September', 'October  ', 'November ', 'December ' /)
             INTEGER ( kind = 4 ) n
             INTEGER ( kind = 4 ) s
             INTEGER ( kind = 4 ) values(8)
             INTEGER ( kind = 4 ) y
             CALL DATE_AND_TIME ( values = values )

             y = values(1)
             m = values(2)
             d = values(3)
             h = values(5)
             n = values(6)
             s = values(7)
             mm = values(8)

             IF ( h < 12 ) THEN
                ampm = 'AM'
             ELSE IF ( h == 12 ) THEN
                IF ( n == 0 .AND. s == 0 ) THEN
                   ampm = 'Noon'
                ELSE
                   ampm = 'PM'
                END IF
             ELSE
                h = h - 12
                IF ( h < 12 ) THEN
                   ampm = 'PM'
                ELSE IF ( h == 12 ) THEN
                   IF ( n == 0 .AND. s == 0 ) THEN
                      ampm = 'Midnight'
                   ELSE
                      ampm = 'AM'
                   END IF
                END IF
             END IF

             WRITE ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
                  d, TRIM ( month(m) ), y, h, ':', n, ':', s, '.', mm, TRIM ( ampm )

             RETURN
           END SUBROUTINE timestamp


           ! ----------------------------------------------------------------------------------%%
           !
           ! GIVE ME A NUMBER I AM NOT USING SO I CAN USE IT TO MAKE MYSELF SOME FILES
           !
           ! ----------------------------------------------------------------------------------%%

           FUNCTION get_unit ( )

             IMPLICIT NONE
             INTEGER :: i, ios, get_unit
             LOGICAL lopen
             get_unit = 0
             DO i = 1, 99
                IF ( i /= 5 .AND. i /= 6 .AND. i /= 9 ) THEN
                   INQUIRE ( unit = i, opened = lopen, iostat = ios )
                   IF ( ios == 0 ) THEN
                      IF ( .NOT. lopen ) THEN
                         get_unit = i
                         RETURN
                      END IF
                   END IF
                END IF
             END DO
             RETURN
           END FUNCTION get_unit



!!!!! NEW LINSPACE !!!!!!!!

           FUNCTION f_linspace(x_1, x_n, n)
             INTEGER :: n, i
             REAL(4) :: x_1, x_n, f_linspace(n)

             DO i = 1,n
                f_linspace(i) = x_1 + (x_n - x_1)*(i-1) / (n-1)
             END DO
           END FUNCTION f_linspace



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine for smoothing of a time-series data .......            !!!
!!! Input to the subroutine is 1-d array and output is also the same  !!!
!!! Input:                                                            !!!
!!!             x       = one dimensional array; x(n)                 !!!
!!!             n       = no. of elements in input array              !!!
!!!             m       = no. of points smoothing (should be odd no.) !!!
!!! Output:                                                           !!!
!!!             y       = one dimensional array; y(n)                 !!!
!!! Version: March 24, 2010 (ATJ/DBS)                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           FUNCTION moving_average (x1, n, m)

             IMPLICIT NONE
             REAL(4) :: x1(n), moving_average(n), sumx
             INTEGER :: m, i, j, i1, i2, k1, k2, iflag, nd, n

             !        Dimension x1(n), moving_average(n)

             IF (MOD(m,2) == 0) m = m + 1

             DO i = 1, n

                k2 = (m/2)

151             CONTINUE

                iflag = 0

                i1 = i - k2
                i2 = i + k2

                IF ((i1 <= 0).OR.(i2 > n)) iflag = 1
                IF ((i1 <= 0).OR.(i2 > n)) k2 = k2 -1

                IF (iflag == 1) Go to 151

                sumx = 0.0

                DO j = i1, i2
                   sumx = sumx + x1(j)
                ENDDO

                nd = (i2 - i1) + 1

                moving_average(i) = sumx/float(nd)

             ENDDO

             RETURN

           END FUNCTION moving_average

         END MODULE globals
