! mesh-blunt-flaredcone-withShock generator V1.0, blunted flared cone with small medium-radius
! outline is a parabolic line conjucted with a linear line

!  Version: Beta V1.0
!  Main Function: generate the whole mesh block  of the blunthead-straightcone-flaredcone model
!  Exceptional Mode: increase the mesh density on both sides of the input local shock curve

! New feature of V2.0:
! destance between wall to the first mesh is constant and are given in control file
! Continues mesh step between conjuction point (between cone and body).
! New feature of V3.0:
! Step of the first mesh to the wall can be adjust by user. (linear function for hy1 (i=1) to hy2 (i=Ny) )
! IF IFLAG_MESHY_TYPE == 1, V1.0 TYPE mesh, 2: V2.0 type mesh
!--------------------------------------------------------------

! Developed by Zhenghao Feng on October 2022

    ! ! Parameter declaration
    ! implicit doubleprecision (a - h, o - z)
    ! real*8, parameter:: PI = 3.1415926535897932
    ! real*8, allocatable:: sx(:), sy(:, :), xa(:), xb(:), ya(:), yb(:), SL(:), delty(:)
    ! real*8, allocatable:: ss(:), ss_new(:), xb_new(:), yb_new(:)
    ! real*8, allocatable:: xx(:,:), yy(:,:)
    ! common/para/ b2, d2, seta2, ys

    ! open(66, file = 'grid2d-withleading.in')
    ! read(66, *)
    ! read(66, *) 
    ! read(66, *) nx, ny
    ! read(66, *)
    ! read(66, *) r1, seta1, b2, d2, seta2
    ! read(66, *)
    ! ! read(66, *) S0, alfay
    ! read(66, *) x_end, alfax, nx_buff, alfax_buff, nxconjunction_buffer
    ! read(66, *)
    ! read(66, *) IFLAG_MESHY_TYPE, deltyw_i_begin, deltyw_i_end, alfay1, alfay2 
    ! read(66,*)
    ! close(66)

    ! Parameter declaration
    implicit doubleprecision (a - h, o - z)
    ! implicit none
    real*8, parameter:: PI = 3.1415926535897932
    real*8, parameter:: MAX_INF = 99999.d0
    real*8, parameter:: MIN_INF = -99999.d0
    integer, parameter:: USER_PARA = 100, USER_LEN = 2000
    !!!!!!!!!!!!!!!!!!!!! 2022/10/07
    !real*8, allocatable:: sx(:), sy(:, :), xa(:), xb(:), ya(:), yb(:), SL(:), delty(:)
    !real*8, allocatable:: ss(:), ss_new(:), xb_new(:), yb_new(:)
    !real*8, allocatable:: xx(:,:), yy(:,:)
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6
    real*8:: xa(USER_LEN), xb(USER_LEN), xc(USER_LEN), ya(USER_LEN), yb(USER_LEN), yc(USER_LEN), thetaa(USER_LEN), thetac(USER_LEN)
    real*8:: ss_b(USER_LEN), ss_c(USER_LEN), ss_b_new(USER_LEN), ss_c_new(USER_LEN)
    real*8:: xb_new(USER_LEN), yb_new(USER_LEN), xc_new(USER_LEN), yc_new(USER_LEN)
    ! real*8:: rn, theta1, R, x3c, x4c
    ! real*8:: a, b, c, xm, xt, ym, yt
    real*8:: x1, x2, x3, x4
    real*8:: dev_X, dev_Y, alfax_buff, delta_yc_trans
    real*8:: nx_ratio_array(USER_PARA), ny_ratio_array(USER_PARA), &
             parax_array(USER_PARA), eta_X(USER_PARA), As_X(USER_PARA), &
             paray_array(USER_PARA), eta_Y(USER_PARA), As_Y(USER_PARA)
    ! real*8:: parayw_array(USER_PARA), parays_array(USER_PARA)
    integer nx_tot, ny_tot, num_x_part, num_y_part, Mesh_X_trans, Mesh_Y_trans, nx_buff, nxconjuction
    integer Mesh_X_TYPE(USER_PARA), Mesh_X_dense(USER_PARA), Mesh_Y_TYPE(USER_PARA), Mesh_Y_dense(USER_PARA)
    ! integer Mesh_Y_normgrowth_TYPE(USER_PARA), Mesh_Y_firstlayer_TYPE(USER_PARA), Mesh_Y_denseshock_TYPE(USER_PARA)
    integer ID_X(USER_LEN), ID_Y(USER_LEN, USER_LEN), ID_x_flag(USER_PARA + 2), ID_y_flag(USER_LEN, USER_PARA + 2), ID_y_flag_tmp(USER_PARA + 2)
    integer delta_y_first_min_id, delta_y_first_max_id, delta_y_final_min_id, delta_y_final_max_id
    
    real*8:: SLX_part(USER_PARA), SLX_prop(USER_LEN), sx_wall_tot(USER_LEN), SLY_part_tmp(USER_PARA), sy_tot_tmp(USER_LEN)
    real*8:: SLX_total
    real*8:: SLY_part(USER_PARA, USER_LEN), sy_tot(USER_LEN, USER_LEN)
    real*8:: SLY_total(USER_LEN)
    real*8:: xx(USER_LEN, USER_LEN), yy(USER_LEN, USER_LEN), xx_tmp(USER_LEN), yy_tmp(USER_LEN)
    real*8:: delta_first_array(USER_PARA), delta_final_array(USER_PARA)
    real*8:: delta_y_first_min(USER_PARA), delta_y_first_max(USER_PARA), delta_y_final_min(USER_PARA), delta_y_final_max(USER_PARA)
    
    ! read the input parameters
    open(66, file = 'grid2d-blunt-flaredcone-withShock-V2.0.in')
    read(66, *)
    read(66, *)
    read(66, *)
    read(66, *)
    read(66, *)
    read(66, *)
    read(66, *)
    read(66, *)
    ! [Model geometry]
    ! [Global parameters]
    read(66, *) rn, theta1, R, x3c, x4c
    read(66, *)
    read(66, *)
    ! [Boundary Definition - the shock wave fitting curve]
    read(66, *) a, b, c, xm, xt, ym, yt, setac
    read(66, *)
    read(66, *)
    ! [Boundary Definition - farfield curve for computation]
    read(66, *) b2, d2, setab
    read(66, *)
    read(66, *)
    read(66, *)
    ! [Grid parameters]
    ! [Global grid number & Partition strategy]
    read(66, *) nx_tot, ny_tot
    read(66, *)
    read(66, *) num_x_part, num_y_part, (nx_ratio_array(i), i = 1, num_x_part), (ny_ratio_array(j), j = 1, num_y_part)
    read(66, *)
    read(66, *)
    ! [Grid Type - x - streamwise direction]
    read(66, *) (Mesh_X_TYPE(i), i = 1, num_x_part), Mesh_X_trans, (parax_array(i), i = 1, num_x_part), dev_X
    read(66, *)
    read(66, *) (Mesh_X_dense(i), i = 1, num_x_part), (eta_X(i), i = 1, num_x_part), (As_X(i), i = 1, num_x_part)
    read(66, *)
    read(66, *)
    ! [Buffer section]
    read(66, *) nx_buff, alfax_buff, nxconjuction
    read(66, *)
    read(66, *)
    ! [Grid Type - y - wall-normal direction] interpolation involved
    !read(66, *) (Mesh_Y_normgrowth_TYPE(j), j = 1, num_y_part), Mesh_Y_trans, dev_Y
    !read(66, *)
    !read(66, *) (Mesh_Y_firstlayer_TYPE(j), j = 1, num_y_part), (parayw_array(i), j = 1, num_y_part)
    !read(66, *)
    !read(66, *) (Mesh_Y_denseshock_TYPE(j), j = 1, num_y_part), (parays_array(i), j = 1, num_y_part)
    read(66, *) (Mesh_Y_TYPE(i), i = 1, num_y_part), Mesh_Y_trans, (paray_array(i), i = 1, num_y_part), dev_Y
    read(66, *)
    read(66, *) (Mesh_Y_dense(i), i = 1, num_y_part), (eta_Y(i), i = 1, num_y_part), (As_Y(i), i = 1, num_y_part)
    close(66)

    ! allocate(SLX_part(:), SLX_prop(:), sx_wall_tot(:))
    
    !allocate(sx(nx), sy(nx, ny), xa(nx), xb(nx),  &
    !         ya(nx), yb(nx), SL(nx), delty(nx), xx(nx, ny), yy(nx, ny),  &
    !         ss(nx), ss_new(nx), xb_new(nx), yb_new(nx))
 
    ! Module 1: construct the basic model
    ! x1 - the coor. pf the blunt head stagnation point
    ! x2 - the coor. of the transition point from the head to the straight frustum
    ! x3 - the beginning coor. of the flared section
    !      i.e. the ending coor. of the straight section
    ! x4 - the ending coor. of the flared section
    !      i.e. the length of the whole model

    theta1 = theta1 * (PI / 180.d0)
    x1 = - 1.d0 * rn
    x2 = - 1.d0 * rn * sin(theta1)
    x3 = x3c - (rn / sin(theta1))
    x4 = x4c - (rn / sin(theta1))
    
    setab = setab * (PI / 180.d0)
    setac = setac * (PI / 180.d0)

    ! Section 1: the blunt head [x1, x2]
    SLX_part(1) = (PI / 2. - theta1) * rn
    ! Section 2: the straight cone frustum (x2, x3]
    SLX_part(2) = (x3 + rn * sin(theta1)) / cos(theta1)
    ! Section 3: the flared cone section (x3, x4]
    h1 = x3c * tan(theta1)
    Ax4c = (x4c - x3c) * (2 * R * sin(theta1) + (x4c - x3c))
    h2 = h1 + R * cos(theta1) - sqrt(R * R * cos(theta1) * cos(theta1) - Ax4c)
    ! the change of theta along the flared arc
    theta_arc = atan((R * sin(theta1) + (x4c - x3c)) / (R * cos(theta1) - (h2 - h1))) - theta1
    SLX_part(3) = R * theta_arc

    ! SLX_total = (PI / 2. - theta1) * rn + (x3 + rn * sin(theta1)) / cos(theta1)
    ! sum up to get the total length of the wall (from the blunt head to the flared rear)
    ! besides, the length proportions of each section "SLX_prop" are also calculated.
    SLX_total = 0
    do k = 1, num_x_part
        SLX_total = SLX_total + SLX_part(k)
    enddo
    do k = 1, num_x_part
        SLX_prop(k) = SLX_part(k) / SLX_total
    enddo
       
    ! generate the grid in streamwise direction (x direction)
    ! sx_wall is always between [0, 1]
    open(99, file = "mesh_generation_info.log")
    call getsx_wall(nx_tot, num_x_part, nx_ratio_array, nx_buff, alfax_buff, &
                    Mesh_X_TYPE, Mesh_X_trans, parax_array, dev_X, &
                    Mesh_X_dense, eta_X, As_X, SLX_part, SLX_total, 1, 1, &
                    sx_wall_tot, ID_x_flag, delta_first_array, delta_final_array)

    ! x grid space = sx(k + 1) - sx(k), k = 1, 2, 3, ..., nx - 1
    ! x_delta_1 = 
    write(*, *)
    write(*, *) "Writing the sx [0, 1] on the model wall ......"
    write(99, *)
    write(99, *) "Writing the sx [0, 1] on the model wall ......"
    open(55, file = "sx_wall.dat")
    do k = 1, nx_tot - 1
        write(55, *) k, sx_wall_tot(k), sx_wall_tot(k + 1) - sx_wall_tot(k)
    enddo
    write(*, *) "Finish writing!"
    write(99, *) "Finish writing!"
    
    ! cal. the corresponding points on the three characteristic boundaries
    ! Boundary 1: (xa, ya) - wall surface
    
    ! Boundary 2: (xb, yb) - far field - parabolic-line curve
    ! parabolic line equ. x = b2 * y * y + d2
    ! (xb, yb) is the coordination of the outline of the computational regime 
    ! (xb, yb) is located on the parabolic line
    
    ! (xbs, ybs) - the transition point from the parabolic curve to the line
    ybs = 1.d0 / (2.d0 * b2 * tan(setab))
	xbs = b2 * ybs * ybs + d2
    
    ! Boundary 3: (xc, yc) - shock wave - fitting curve
    !#  the formula of the fitting shock curve: y = f(x) = a * ((x + 1) ** b - 1.d0) ** (1.d0 / c)
    !#  using the 'lsqcurvefit' algorithm in MATLAB
    !#  a - the 1st fitting coeff.;
    !   b - the 2nd fitting coeff.;
    !   c - the 3rd fitting coeff.;
    !   xm - the MAX x for the scaled fitting range;
    !   xt - the transformed distance of x;
    !   ym - the MAX y for the scaled fitting range;
    !   yt - the transformed distance of y;
    !#  the formula with the real coor.: y_real = (a * ((((((x_real / rn) + xt) / xm)) ** b - 1.d0) ** (1.d0 / c)) * ym - yt) * rn    
    
    ! coeff. list
    c1 = a * rn * ym
    c2 = 1.d0 / (rn * xm)
    c3 = ((xt * 1.d0) / (xm * 1.d0)) + 1.d0
    c4 = b * 1.d0
    c5 = 1.d0 / (c * 1.d0)
    c6 = (- 1.d0) * rn * yt
    
    deltax0 = (- 1.d0) * xt * rn
    deltay0 = c1 * (((c2 * deltax0 + c3) ** c4 - 1.d0) ** c5) + c6
    
    c6 = c6 - deltay0
    
    ! search (xcs, ycs) - the transition point from the fitting curve to the line
    ! ycs = 1.d0 / (2.d0 * b2 * tan(setab))
	! xcs = b2 * ycs * ycs + d2
    find_xcs_gap = 5.d-2
    find_xcs_s = (- 1.d0) * xt * rn
    find_xcs_t = x4c
    find_xcs_err_min = MAX_INF
    N_xcs = ((find_xcs_t - find_xcs_s) / (find_xcs_gap * 1.d0)) + 1
    
    ! do xcs_tmp = find_xcs_s, find_xcs_gap, find_xcs_t
    do k = 1, N_xcs
        xcs_tmp = find_xcs_s + (k - 1) * find_xcs_gap
        find_xcs_err_tmp = abs((c1 * c2 * c4 * c5) * (((c2 * xcs_tmp + c3) ** c4 - 1.d0) ** (c5 - 1.d0)) * ((c2 * xcs_tmp + c3) ** (c4 - 1.d0)) - tan(setac))
        if (find_xcs_err_min > find_xcs_err_tmp) then
            find_xcs_err_min = find_xcs_err_tmp
            xcs = xcs_tmp
        endif
    enddo
    ycs = c1 * (((c2 * xcs + c3) ** c4 - 1.d0) ** c5) + c6
    
    ! Note: Make sure the grid lines are normal to the previous cuvre (shock wave to wall surface, farfield to shock wave)
    ! scan the curve procedure
    do i = 1, nx_tot

        ! natural coor. s (ncs), i.e. the arc length from the stagnation to the present pos. 
        s = sx_wall_tot(i) * SLX_total

        ! Section 1: blunt head section
	    if (s .le. SLX_part(1)) then
            
            ID_X(i) = 1

            ! (xa, ya) - wall surface
            ! alfa - circular center angle
            alfa = s / rn
            xa(i) = - rn * cos(alfa)
            ya(i) = rn * sin(alfa)
            
            ! (xc, yc) - shock wave
            thetaa(i) = (PI / 2.d0) - alfa
            ! common module
            call getxc_shock_wall_normal(xa(i), ya(i), thetaa(i), xc_tmp, yc_tmp, thetac_tmp)
            xc(i) = xc_tmp
            yc(i) = yc_tmp
            thetac(i) = thetac_tmp
            ! yc(1) = 0.d0
            
            !if (xc(i) .le. xcs) then
            !    thetac(i) = atan( (c1 * c2 * c4 * c5) * (((c2 * xc(i) + c3) ** c4 - 1.d0) ** (c5 - 1.d0)) * ((c2 * xc(i) + c3) ** (c4 - 1.d0)) )
            !else
            !    thetac(i) = setac
            !endif
            
            !! (xb, yb) - parabolic-line farfield
            !xb(i) = (1.d0 - sqrt(1.d0 - 4.d0 * b2 * d2 * tan(alfa) ** 2)) / (2.d0 * b2 * tan(alfa) ** 2)
            !yb(i) = - xb(i) * tan(alfa)
            !! equal to...
            !! yb(i) = (- 1.d0 + sqrt(1.d0 - 4.d0 * b2 * d2 * tan(alfa) ** 2)) / (2.d0 * b2 * tan(alfa))
            !! xb(i) = b2 * yb(i) * yb(i) + d2
            !if (xb(i) .gt. xbs) then
            !    ! (xb, yb) is located on the linear line
	           ! xb(i) = - (ybs - tan(setab) * xbs) / (tan(alfa) + tan(setab))
	           ! yb(i) = - tan(alfa) * xb(i)
            !endif
            
            ! (xb, yb) - parabolic-line farfield
            call getxb_farfield_shock_normal(xc(i), yc(i), thetac(i), xb_tmp, yb_tmp)
            xb(i) = xb_tmp
            yb(i) = yb_tmp
	    
        ! Section 2: straight cone section
        else if ((s .gt. SLX_part(1)) .and. (s .le. (SLX_part(1) + SLX_part(2)))) then
            
            ID_X(i) = 2

            ! (xa, ya) - wall surface
	        sa = s - SLX_part(1)
	        xa(i) = - rn * sin(theta1) + sa * cos(theta1)
	        ya(i) = rn * cos(theta1) + sa * sin(theta1)
            
            ! (xc, yc) - shock wave
            thetaa(i) = theta1
            ! common module
            call getxc_shock_wall_normal(xa(i), ya(i), thetaa(i), xc_tmp, yc_tmp, thetac_tmp)
            xc(i) = xc_tmp
            yc(i) = yc_tmp
            thetac(i) = thetac_tmp

         !   ! (xb, yb) - parabolic-line farfield
	        !at = b2
	        !bt = tan(theta1)
	        !ct = d2 - xa(i) - ya(i) * tan(theta1)
	        !yb(i) = (- bt + sqrt(bt * bt - 4. * at * ct)) / (2. * at)
	        !xb(i) = b2 * yb(i) * yb(i) + d2
         !   if (xb(i) .gt. xs) then
         !       xb(i) = (ya(i) - ybs + tan(setab) * xbs - tan(PI / 2. + theta1) * xa(i)) /  &
         !               (tan(setab) - tan(PI / 2. + theta1))
	        !    yb(i) = ybs + tan(setab) * (xb(i) - xbs)
         !   endif
            
            ! (xb, yb) - parabolic-line farfield
            call getxb_farfield_shock_normal(xc(i), yc(i), thetac(i), xb_tmp, yb_tmp)
            xb(i) = xb_tmp
            yb(i) = yb_tmp
            
        ! Section 3: flared cone section
        else if ((s .gt. (SLX_part(1) + SLX_part(2))) .and. (s .le. SLX_total)) then
            
            ID_X(i) = 3
            
            ! (xa, ya) - wall surface
            sa = s - (SLX_part(1) + SLX_part(2))
            thetaa(i) = (sa * 1.d0) / (R * 1.d0) + theta1
            xa(i) = R * (sin(thetaa(i)) - sin(theta1)) - rn * sin(theta1) + SLX_part(2) * cos(theta1)
            ya(i) = R * (cos(theta1) - cos(thetaa(i))) + rn * cos(theta1) + SLX_part(2) * sin(theta1)

            ! (xc, yc) - shock wave
            ! common module
            call getxc_shock_wall_normal(xa(i), ya(i), thetaa(i), xc_tmp, yc_tmp, thetac_tmp)
            xc(i) = xc_tmp
            yc(i) = yc_tmp
            thetac(i) = thetac_tmp

            ! (xb, yb) - parabolic-line farfield
            call getxb_farfield_shock_normal(xc(i), yc(i), thetac(i), xb_tmp, yb_tmp)
            xb(i) = xb_tmp
            yb(i) = yb_tmp
            
        ! Section 4: buffer section (straight line)
        else if (s .gt. SLX_total) then
            
            ID_X(i) = 4
            
            ! (xa, ya) - wall surface
            sa = s - SLX_total
            thetaa(i) = thetaa(i - 1)
            xa(i) = xa(i - 1) + sa * cos(thetaa(i))
            ya(i) = ya(i - 1) + sa * sin(thetaa(i))
            
            ! (xc, yc) - shock wave
            ! common module
            call getxc_shock_wall_normal(xa(i), ya(i), thetaa(i), xc_tmp, yc_tmp, thetac_tmp)
            xc(i) = xc_tmp
            yc(i) = yc_tmp
            thetac(i) = thetac_tmp
            
            ! (xb, yb) - parabolic-line farfield
            call getxb_farfield_shock_normal(xc(i), yc(i), thetac(i), xb_tmp, yb_tmp)
            xb(i) = xb_tmp
            yb(i) = yb_tmp
            
        endif

        ! cal. the arc length of outter-line (integal along the y direction)
	    call simpson_curve(ss_b(i), 0.d0, yb(i), 1)
        call simpson_curve(ss_c(i), 0.d0, yc(i), 2)
    
    enddo
    
    !!!!!!!!!!!!!!!!!!!!
    ! Test module
    ! Adjust mesh points near the conjunction points to smooth the connection

    ! i_conjunction - find the point nearest to the conjection point
    do i = 1, nx_tot
        s = sx_wall_tot(i) * SLX_total
	    if (s .gt. SLX_part(1)) then
           i_conjunction = i
           goto 101
		endif
	enddo

101 continue
    write(*, *)
    write(*, *) "The point found nearest to the conjection point: ID = ", i_conjunction
    write(99, *)
    write(99, *) "The point found nearest to the conjection point: ID = ", i_conjunction
    
    ! important module
    ! Adjust mesh points near the conjunction points to smooth the connection
	call adjust_outter(nx_tot, ss_b, ss_b_new, i_conjunction, nxconjuction)
    write(*, *)
    write(*, *) "Finish to adjust the outer parabolic-line boundary!"
    write(99, *)
    write(99, *) "Finish to adjust the outer parabolic-line boundary!"
    call adjust_outter(nx_tot, ss_c, ss_c_new, i_conjunction, nxconjuction)
    write(*, *)
    write(*, *) "Finish to adjust the shock wave boundary!"
    write(99, *)
    write(99, *) "Finish to adjust the shock wave boundary!"
     
	do i = 1, nx_tot
	    call get_y_form_ss(ss_b_new(i), yb_new(i), 1.d0, 1)
        call get_y_form_ss(ss_c_new(i), yc_new(i), 1.d0, 2)
    enddo

    open(33, file = 'yb_comp.dat')
    do i = 1, nx_tot
	    write(33, '(6f16.8)') i * 1., yb(i), yb_new(i)
    enddo
    open(33, file = 'yc_comp.dat')
    do i = 1, nx_tot
	    write(33, '(6f16.8)') i * 1., yc(i), yc_new(i)
    enddo
    
    do i = 1, nx_tot

        ! located on the parabolic line
        if (yb_new(i) .lt. ybs) then
            xb_new(i) = b2 * yb_new(i) * yb_new(i) + d2
		else
            xb_new(i) = xbs + (yb_new(i) - ybs) / tan(setab)
        endif
        
        ! located on the shock wave
        if (yc_new(i) .lt. ycs) then
            xc_new(i) = ((1.d0 + ((yc_new(i) - c6) / c1) ** (1.d0 / c5)) ** (1.d0 / c4) - c3) / (c2 * 1.d0)
		else
            xc_new(i) = xcs + (yc_new(i) - ycs) / tan(setac)
		endif
	
    enddo

    open(55, file = "grid1d_comp.dat")
    do i = 1, nx_tot
        write(55, "(10f15.6)") xa(i), ya(i), xb(i), yb(i), xb_new(i), yb_new(i), xc(i), yc(i), xc_new(i), yc_new(i)
    enddo
    close(55)
    
    ! Update the mesh (xb, yb) --> new mesh (xb_new, yb_new)
    !                 (xc, yc) --> new mesh (xc_new, yc_new)
    xb = xb_new
	yb = yb_new
    xc = xc_new
	yc = yc_new
    
    !!!!!!!!!!!!!!!!!!!!
    
    write(*, *)
    write(*, *) "Writing the grids on the three boudaries: wall, shock wave, farfield ......"
    write(99, *)
    write(99, *) "Writing the grids on the three boudaries: wall, shock wave, farfield ......"
    !open(77, file = "x_curves_coor.dat")
    !do i = 1, nx_tot
    !    write(77, *) i, xa(i), ya(i), xc(i), yc(i), xb(i), yb(i)
    !enddo
    write(*, *) "Writing the coor. xa ..."
    write(99, *) "Writing the coor. xa ..."
    open(77, file = "xa.dat")
    do i = 1, nx_tot
        write(77, *) i, ID_X(i), xa(i), ya(i)
    enddo
    write(*, *) "Writing the coor. xb ..."
    write(99, *) "Writing the coor. xb ..."
    open(77, file = "xb.dat")
    do i = 1, nx_tot
        write(77, *) i, ID_X(i), xb(i), yb(i)
    enddo
    write(*, *) "Writing the coor. xc ..."
    write(99, *) "Writing the coor. xc ..."
    open(77, file = "xc.dat")
    do i = 1, nx_tot
        write(77, *) i, ID_X(i), xc(i), yc(i) !, thetac(i)
    enddo
    write(*, *) "Writing the angle theta_xa ..."
    write(99, *) "Writing the angle theta_xa ..."
    open(77, file = "theta_xa.dat")
    do i = 1, nx_tot
        write(77, *) i, ID_X(i), thetaa(i)
    enddo
    write(*, *) "Writing the angle theta_xc ..."
    write(99, *) "Writing the angle theta_xc ..."
    open(77, file = "theta_xc.dat")
    do i = 1, nx_tot
        write(77, *) i, ID_X(i), thetac(i)
    enddo
    write(*, *) "Finish writing!"
    write(*, *)
    write(99, *) "Finish writing!"
    write(99, *)
    
    
    ! handle the shock wave curve
    ! ensure the grid lines are normal to the previous cuvre (wall surface)
    ! the shock wave expression
    ! the formula with the real coor.: y_real = (a * ((((((x_real / rn) + xt) / xm)) ** b - 1.d0) ** (1.d0 / c)) * ym - yt) * rn
    
    ! cal. the first grid height along the wall surface
    ! cal. the distribution of the control parameter
    !do k = 1, num_x_part
    !    
    !    if (k == 1) then
    !        ID_x_s = ID_x_flag(k)
    !    else
    !        ID_x_s = ID_x_flag(k) + 1
    !    endif
    !    ID_x_t = ID_x_flag(k + 1)
    !    
    !    ! vary para. linear distribution
    !    if (Mesh_Y_firstlayer_vary(k) == 0) then
    !        
    !        if (Mesh_Y_firstlayer_TYPE(k) == 1) then
    !        parayw_array
    !        else if (Mesh_Y_firstlayer_TYPE(k) == 2) then
    !        
    !        else if (Mesh_Y_firstlayer_TYPE(k) == 3) then
    !    
    !    ! fixed para. linear distribution
    !    else if (Mesh_Y_firstlayer_vary(k) == 1)
    !        
    !    endif
    !    
    !enddo
    
    delta_y_first_min = MAX_INF
    delta_y_first_max = MIN_INF
    delta_y_final_min = MAX_INF
    delta_y_final_max = MIN_INF
        
    do i = 1, nx_tot
        
        ! cal. the length of the curve-normal grid lines
        !SLX_total = 0
        !do k = 1, num_x_part
        !    SLX_total = SLX_total + SLX_part(k)
        !enddo
        !do k = 1, num_x_part
        !    SLX_prop(k) = SLX_part(k) / SLX_total
        !enddo
        
        SLY_part(1, i) = sqrt((xc(i) - xa(i)) ** 2 + (yc(i) - ya(i)) ** 2)
        SLY_part(2, i) = sqrt((xb(i) - xc(i)) ** 2 + (yb(i) - yc(i)) ** 2)
        SLY_part_tmp(1) = SLY_part(1, i)
        SLY_part_tmp(2) = SLY_part(2, i)
        SLY_total(i) = SLY_part_tmp(1) + SLY_part_tmp(2)
        
        !do k = 1, num_y_part
        !    SLY_prop(k) = SLY_part_tmp(k) / SLY_total(i)
        !enddo
        
        call getsx_wall(ny_tot, num_y_part, ny_ratio_array, 0, 0.d0, &
                        Mesh_Y_TYPE, Mesh_Y_trans, paray_array, dev_Y, &
                        Mesh_Y_dense, eta_Y, As_Y, SLY_part_tmp, SLY_total(i), 2, i, &
                        sy_tot_tmp, id_y_flag_tmp, delta_first_array, delta_final_array)
        
        do k = 1, num_y_part
            
            if (delta_y_first_min(k) > delta_first_array(k)) then
                delta_y_first_min(k) = delta_first_array(k)
                delta_y_first_min_id = i
            endif
            if (delta_y_first_max(k) < delta_first_array(k)) then
                delta_y_first_max(k) = delta_first_array(k)
                delta_y_first_max_id = i
            endif
    
            if (Mesh_Y_dense(k) == 1) then
                if (delta_y_final_min(k) > delta_final_array(k)) then
                    delta_y_final_min(k) = delta_final_array(k)
                    delta_y_final_min_id = i
                endif
                if (delta_y_final_max(k) < delta_final_array(k)) then
                    delta_y_final_max(k) = delta_final_array(k)
                    delta_y_final_max_id = i
                endif
            endif
    
        enddo
        
        do j = 1, ny_tot
            
            sy_tot(j, i) = sy_tot_tmp(j)
            ! id_y_flag(j, i) = id_y_flag_tmp(j)
            
            ! natural wall-normal coor. s (ncs)
            s = sy_tot_tmp(j) * SLY_total(i)

            ! Section 1: blunt head section
	        if (s .le. SLY_part_tmp(1)) then
                xx_tmp(j) = xa(i) + (xc(i) - xa(i)) * ((s * 1.d0) / (SLY_part_tmp(1) * 1.d0))
                yy_tmp(j) = ya(i) + (yc(i) - ya(i)) * ((s * 1.d0) / (SLY_part_tmp(1) * 1.d0))
                xx(i, j) = xx_tmp(j)
                yy(i, j) = yy_tmp(j)
            !! judge condition too strict!!!
            else if ((s .gt. SLY_part_tmp(1)) .and. (s .le. SLY_total(i) + 1.d-6))  then
                xx_tmp(j) = xc(i) + (xb(i) - xc(i)) * ((s - SLY_part_tmp(1) * 1.d0) / (SLY_part_tmp(2) * 1.d0))
                yy_tmp(j) = yc(i) + (yb(i) - yc(i)) * ((s - SLY_part_tmp(1) * 1.d0) / (SLY_part_tmp(2) * 1.d0))
                xx(i, j) = xx_tmp(j)
                yy(i, j) = yy_tmp(j)
            endif
            
        enddo
        
    enddo
    
    ! ! generate the mesh along the wall-normal direction
    ! do i = 1, nx
    !    do j = 1, ny
	!        xx(i, j) = xb(i) + (xa(i) - xb(i)) * sy(i, j)
	!        yy(i, j) = yb(i) + (ya(i) - yb(i)) * sy(i, j)
	!    enddo
	!enddo
    
    ! output deltay info.
    write(*, *)
    write(*, *) "Write Y Grid info ..."
    write(99, *)
    write(99, *) "Write Y Grid info ..."
    do k = 1, num_y_part
        write(*, *) "Y Section", k, "MIN. first deltay = ", delta_y_first_min(k)
        write(*, *) " at ID = ", delta_y_first_min_id
        write(*, *) "Y Section", k, "MAX. first deltay = ", delta_y_first_max(k)
        write(*, *) " at ID = ", delta_y_first_max_id
        write(99, *) "Y Section", k, "MIN. first deltay = ", delta_y_first_min(k)
        write(99, *) " at ID = ", delta_y_first_min_id
        write(99, *) "Y Section", k, "MAX. first deltay = ", delta_y_first_max(k)
        write(99, *) " at ID = ", delta_y_first_max_id
        if (Mesh_Y_dense(k) == 1) then
            write(*, *) "Y Section", k, "MIN. final deltay = ", delta_y_final_min(k)
            write(*, *) " at ID = ", delta_y_final_min_id
            write(*, *) "Y Section", k, "MAX. final deltay = ", delta_y_final_max(k)
            write(*, *) " at ID = ", delta_y_final_max_id
            write(99, *) "Y Section", k, "MIN. final deltay = ", delta_y_final_min(k)
            write(99, *) " at ID = ", delta_y_final_min_id
            write(99, *) "Y Section", k, "MAX. final deltay = ", delta_y_final_max(k)
            write(99, *) " at ID = ", delta_y_final_max_id
        endif
    enddo
    
    ! Write the final mesh
    write(*, *)
    write(*, *) "Writing the final mesh ..."
    write(99, *)
    write(99, *) "Writing the final mesh ..."
    open(33, file = 'grid_xy.dat')
    write(33, *) 'variables = x, y'
    write(33, *) 'zone i = ', nx_tot, 'j = ', ny_tot
    do j = 1, ny_tot
        do i = 1, nx_tot
            write(33, '(7f15.6)') xx(i, j), yy(i, j)
        enddo
    enddo
    close(33)
    open(33, file = 'grid_xy_plot.dat')
    write(33, *) nx_tot, ny_tot
    do j = 1, ny_tot
        do i = 1, nx_tot
            write(33, '(7f15.6)') xx(i, j), yy(i, j)
        enddo
    enddo
    close(33)
    write(*, *) "Finish writing mesh! "
    write(*, *)
    write(99, *) "Finish writing mesh! "
    write(99, *)
    
    write(*, *)
    write(*, *) "Cal. the Jacobian martix from mesh ..."
    write(99, *)
    write(99, *) "Cal. the Jacobian martix from mesh ..."
    ! call get_Jacobian(nx_tot, ny_tot, yy_new, xx_new)
    write(*, *) "Finish writing the Jacobian martix!"
    write(99, *) "Finish writing the Jacobian martix!"
    
    write(*, *)
    write(*, *) "This is the END of the PROGRAM! Please check the quality of the mesh!"
    write(*, *)
    write(99, *)
    write(99, *) "This is the END of the PROGRAM! Please check the quality of the mesh!"
    write(99, *)

    end
    
    ! cal. the grid points on the shock wave curve which satisfies the wall-normal requirements
    subroutine getxc_shock_wall_normal(xa, ya, thetaa, xc, yc, thetac)
    implicit doubleprecision (a - h, o - z)
    real*8, parameter:: PI = 3.1415926535897932
    real*8, parameter:: MAX_INF = 99999.d0
    ! common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6
    
        eps = 1.d-6
        
        if (xa .le. 0.d0) then
            
         !   ! search (xc, yc) using the slope relation
         !   ! ycs = 1.d0 / (2.d0 * b2 * tan(setab))
	        !! xcs = b2 * ycs * ycs + d2
         !   find_xc_gap = 1.d-5
         !   find_xc_s = (- 1.d0) * xt * rn
         !   find_xc_t = 0.d0
         !   find_xc_err_min = MAX_INF
         !   N_xc = ((find_xc_t - find_xc_s) / (find_xc_gap * 1.d0)) + 1
         !
         !   ! do xcs_tmp = find_xcs_s, find_xcs_gap, find_xcs_t
         !   do k = 1, N_xc
         !       xc_tmp = find_xc_s + (k - 1) * find_xc_gap
         !       yc_tmp = c1 * (((c2 * xc_tmp + c3) ** c4 - 1.d0) ** c5) + c6
         !       find_xc_err_tmp = abs((yc_tmp - ya) - (tan(thetaa + (PI / 2.d0)) * (xc_tmp - xa)))
         !       if (find_xc_err_min > find_xc_err_tmp) then
         !           find_xc_err_min = find_xc_err_tmp
         !           xc_ans = xc_tmp
         !           yc_ans = yc_tmp
         !       endif
         !   enddo
         !   ! yc = c1 * (((c2 * xc + c3) ** c4 - 1.d0) ** c5) + c6
            
            ! search (xc, yc) using the slope relation
            find_yc_gap = 1.d-5
            find_yc_s = 0.d0
            find_yc_t = c1 * (((c2 * xa + c3) ** c4 - 1.d0) ** c5) + c6
            find_yc_err_min = MAX_INF
            N_yc = ((find_yc_t - find_yc_s) / (find_yc_gap * 1.d0)) + 1
    
            ! do xcs_tmp = find_xcs_s, find_xcs_gap, find_xcs_t
            do k = 1, N_yc
                yc_tmp = find_yc_s + (k - 1) * find_yc_gap
                xc_tmp = ((1.d0 + ((yc_tmp - c6) / c1) ** (1.d0 / c5)) ** (1.d0 / c4) - c3) / (c2 * 1.d0)
                find_yc_err_tmp = abs((yc_tmp - ya) - (tan(thetaa + (PI / 2.d0)) * (xc_tmp - xa)))
                if (find_yc_err_min > find_yc_err_tmp) then
                    find_yc_err_min = find_yc_err_tmp
                    yc_ans = yc_tmp
                    xc_ans = xc_tmp
                endif
            enddo
            ! yc = c1 * (((c2 * xc + c3) ** c4 - 1.d0) ** c5) + c6
            
            if (xc_ans .le. xcs) then
                xc = xc_ans
                yc = yc_ans
            else
                dev_upper = (ya - ycs) + (xcs * tan(setac)) - (xa * tan(thetaa + (PI / 2.d0)))
                dev_lower = tan(setac) - tan(thetaa + (PI / 2.d0)) + eps
                xc = (1.d0 * dev_upper) / (1.d0 * dev_lower)
                yc = (tan(setac) * (xc - xcs)) + ycs
            endif
    
        else
            
            ! initial guess coefficient
            ! if (xa .le. 0.d0) then
            !    xc = xa - 1.d0
            !else
            xc = xa
    
            ! using Newton iterative method get the real xc
            flag_find_sol = 0
            do while (flag_find_sol == 0)

                pcs = c1 * (((c2 * xc + c3) ** c4 - 1.d0) ** c5) + c6
                pcsx = (c1 * c2 * c4 * c5) * (((c2 * xc + c3) ** c4 - 1.d0) ** (c5 - 1.d0)) * ((c2 * xc + c3) ** (c4 - 1.d0))
        
                fcs = pcs - (tan(thetaa + (PI / 2.d0)) * xc) + ((tan(thetaa + (PI / 2.d0)) * xa) - ya)
                fcsx = pcsx - tan(thetaa + (PI / 2.d0))

                xc_new = xc - fcs / fcsx
                xc_err = abs(xc - xc_new)
            
                if (xc_err .gt. 1.d-6) then
                    xc = xc_new
                else
                    xc_final = xc_new
                    flag_find_sol = 1
                endif
        
            enddo
        
            xc_tmp = xc_final
            yc_tmp = c1 * (((c2 * xc_tmp + c3) ** c4 - 1.d0) ** c5) + c6
            
            if (xc_tmp .le. xcs) then
                xc = xc_tmp
                yc = yc_tmp
            else
                dev_upper = (ya - ycs) + (xcs * tan(setac)) - (xa * tan(thetaa + (PI / 2.d0)))
                dev_lower = tan(setac) - tan(thetaa + (PI / 2.d0)) + eps
                xc = (1.d0 * dev_upper) / (1.d0 * dev_lower)
                yc = (tan(setac) * (xc - xcs)) + ycs
            endif
            
        endif
        
        if (xc .le. xcs) then
            thetac = atan( (c1 * c2 * c4 * c5) * (((c2 * xc + c3) ** c4 - 1.d0) ** (c5 - 1.d0)) * ((c2 * xc + c3) ** (c4 - 1.d0)) )
        else
            thetac = setac
        endif
    
    end subroutine
    
    ! cal. the grid points on the farfield curve which satisfies the shock wave-normal requirements
    subroutine getxb_farfield_shock_normal(xc, yc, thetac, xb, yb)
    implicit doubleprecision (a - h, o - z)
    real*8, parameter:: PI = 3.1415926535897932
    real*8, parameter:: MAX_INF = 99999.d0
    ! common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6
        
        eps = 1.d-6
    
        if (xc .le. 0.) then
            
         !   ! search (xb, yb) using the slope relation
         !   ! ybs = 1.d0 / (2.d0 * b2 * tan(setab))
	        !! xbs = b2 * ybs * ybs + d2
         !   find_xb_gap = 1.d-5
         !   find_xb_s = d2
         !   find_xb_t = 0.
         !   find_xb_err_min = MAX_INF
         !   N_xb = ((find_xb_t - find_xb_s) / (find_xb_gap * 1.d0)) + 1
         !
         !   ! do xcs_tmp = find_xcs_s, find_xcs_gap, find_xcs_t
         !   do k = 1, N_xb
         !       xb_tmp = find_xb_s + (k - 1) * find_xb_gap
         !       yb_tmp = sqrt((xb_tmp - d2) / (b2 * 1.d0))
         !       find_xb_err_tmp = abs((yb_tmp - yc) - (tan(thetac + (PI / 2.d0)) * (xb_tmp - xc)))
         !       if (find_xb_err_min > find_xb_err_tmp) then
         !           find_xb_err_min = find_xb_err_tmp
         !           xb_ans = xb_tmp
         !           yb_ans = yb_tmp
         !       endif
         !   enddo
         !   ! yb = sqrt((xb - d2) / (b2 * 1.d0))
            
            ! search (xb, yb) using the slope relation
            ! ybs = 1.d0 / (2.d0 * b2 * tan(setab))
	        ! xbs = b2 * ybs * ybs + d2
            find_yb_gap = 1.d-5
            find_yb_s = 0.d0
            find_yb_t = sqrt((xc - d2) / (b2 * 1.d0))
            find_yb_err_min = MAX_INF
            N_yb = ((find_yb_t - find_yb_s) / (find_yb_gap * 1.d0)) + 1
    
            ! do xcs_tmp = find_xcs_s, find_xcs_gap, find_xcs_t
            do k = 1, N_yb
                yb_tmp = find_yb_s + (k - 1) * find_yb_gap
                xb_tmp = b2 * yb_tmp * yb_tmp + d2
                find_yb_err_tmp = abs((yb_tmp - yc) - (tan(thetac + (PI / 2.d0)) * (xb_tmp - xc)))
                if (find_yb_err_min > find_yb_err_tmp) then
                    find_yb_err_min = find_yb_err_tmp
                    yb_ans = yb_tmp
                    xb_ans = xb_tmp
                endif
            enddo
            ! yb = sqrt((xb - d2) / (b2 * 1.d0))
            
            if (xb_ans .le. xbs) then
                xb = xb_ans
                yb = yb_ans
            else
                dev_upper = (yc - ybs) + (xbs * tan(setab)) - (xc * tan(thetac + (PI / 2.d0)))
                dev_lower = tan(setab) - tan(thetac + (PI / 2.d0)) + eps
                xb = (1.d0 * dev_upper) / (1.d0 * dev_lower)
                yb = (tan(setab) * (xb - xbs)) + ybs
            endif
            
        else
            
            p = b2 * tan(thetac + (PI / 2.d0))
            q = yc + (tan(thetac + (PI / 2.d0)) * (d2 - xc))
            yb_tmp = (1.d0 - sqrt(1.d0 - (4.d0 * p * q))) / (2.d0 * p + eps)
            xb_tmp = b2 * yb_tmp * yb_tmp + d2
            
            if (xb_tmp .le. xbs) then
                xb = xb_tmp
                yb = yb_tmp
            else
                dev_upper = (yc - ybs) + (xbs * tan(setab)) - (xc * tan(thetac + (PI / 2.d0)))
                dev_lower = tan(setab) - tan(thetac + (PI / 2.d0)) + eps
                xb = (1.d0 * dev_upper) / (1.d0 * dev_lower)
                yb = (tan(setab) * (xb - xbs)) + ybs
            endif
            
        endif
    
    end subroutine


    ! cal. s distance along the wall from the stagnation point of each grid
    subroutine getsx_wall(nx, num_x_part, nx_ratio_array, nx_buff, alfax_buff, &
                          Mesh_X_TYPE, Mesh_X_trans, parax_array, dev_X, &
                          Mesh_X_dense, eta_X, As_X, SLX_part, SLX_total, IFLAG_X, IFLAG_CNT, &
                          sx_wall_tot, id_x, delta_first_array, delta_final_array)
    implicit doubleprecision (a - h, o - z)
    ! implicit none
    !real*8, parameter:: PI = 3.1415926535897932
    !real*8, parameter:: MAX_INF = 99999.d0
    !integer, parameter:: USER_PARA = 100, USER_LEN = 2000
    !! define the allocable arrays
    !integer nx, num_x_part, nx_buff, Mesh_X_trans
    !integer Mesh_X_TYPE(USER_PARA), Mesh_X_dense(USER_PARA)
    !! real*8:: alfax_buff, dev_X, SLX_total, alfax_get, alfat_get, betax_get, betat_get
    !real*8:: nx_ratio_array(USER_PARA), parax_array(USER_PARA)
    !real*8:: eta_X(USER_PARA), As_X(USER_PARA), SLX_part(USER_LEN)
    
    ! real*8, parameter:: PI = 3.1415926535897932
    real*8, parameter:: MAX_INF = 99999.d0
    real*8, parameter:: MIN_INF = -99999.d0
    ! integer, parameter:: USER_PARA = 100, USER_LEN = 2000
    ! define the allocable arrays
    ! integer nx, num_x_part, nx_buff, Mesh_X_trans
    integer Mesh_X_TYPE(num_x_part), Mesh_X_dense(num_x_part)
    ! real*8:: alfax_buff, dev_X, SLX_total, alfax_get, alfat_get, betax_get, betat_get
    real*8:: nx_ratio_array(num_x_part), parax_array(num_x_part)
    real*8:: eta_X(num_x_part), As_X(num_x_part), SLX_part(nx)
    
    real*8:: sx_wall_part(nx, num_x_part)
    real*8:: sx_wall_tmp(nx), sx_wall_tot(nx)
    real*8:: delta_first_array(num_x_part), delta_final_array(num_x_part)
    integer id_x(num_x_part + 2), num_dx(num_x_part)
    
    real*8, parameter:: PI = 3.1415926535897932
    
	! real*8:: sx_wall_part(nx, num_x_part)
    ! real*8:: sx_wall_tmp(nx), sx_wall_tot(nx)
    ! integer id_x(num_x_part + 1), num_dx(num_x_part)

        ! cal. the number of grids of each section
        nx1 = nx - nx_buff
        ! id_x: the boundary id. of the grid in x direction
        id_x(1) = 1
        id_x(num_x_part + 1) = nx1
        id_x(num_x_part + 2) = nx
        do i = 1, (num_x_part - 1)
            num_dx(i) = floor(nx1 * nx_ratio_array(i)) - 1
            id_x(i + 1) = id_x(i) + num_dx(i)
        enddo

        ! first, generate grids except the buffer region id_x = [nx1 + 1, nx]
        do k = 1, num_x_part

            ! get the starting and ending id. of the grid in the present section k
            id_x_s = id_x(k)
            id_x_t = id_x(k + 1)
            ! the number of nodes in the present section k
            nx_tmp = (id_x_t - id_x_s) + 1

            ! Mode 1: Equal grid space
            if (Mesh_X_TYPE(k) == 1) then
                
                ! X direction
                if (IFLAG_X == 1) then
                    write(*, *) "Enter Mesh_X global generation mode: Equal grid space"
                    write(99, *) "Enter Mesh_X global generation mode: Equal grid space"
                else if (IFLAG_X == 2) then
                    write(*, *) "Enter Mesh_Y", IFLAG_CNT, "global generation mode: Equal grid space"
                    write(99, *) "Enter Mesh_Y", IFLAG_CNT, "global generation mode: Equal grid space"
                endif

                ! the equally divided grid interval in the present section k
                dx_tmp = 1.d0 / (nx_tmp - 1.d0)

                ! cal. the "sx_wall_part" of the present section
                do j = 1, nx_tmp

                    ! sx_wall_part is always between [0, 1]
                    sx_wall_tmp(j) = (j - 1) * dx_tmp

                enddo
                delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))

            ! Mode 2: Expoential grid stretch - "alfax" in "parax_array" - specify the coeff. of Exp. relation
            else if (Mesh_X_TYPE(k) == 2) then
                
                if (IFLAG_X == 1) then
                    write(*, *) "Enter Mesh_X global generation mode: Expoential grid stretch"
                    write(99, *) "Enter Mesh_X global generation mode: Expoential grid stretch"
                else if (IFLAG_X == 2) then
                    write(*, *) "Enter Mesh_Y", IFLAG_CNT, "global generation mode: Expoential grid stretch"
                    write(99, *) "Enter Mesh_Y", IFLAG_CNT, "global generation mode: Expoential grid stretch"
                endif

                ! choose mode according to the grid transition flag
                ! 0 - respectively generate grids
                !     use "alfax" in "parax_array", the coeff. of the exp. growth, to cal. grid distribution of each section;
                if (Mesh_X_trans == 0) then

                    call gets_exp_free(nx_tmp, parax_array(k), dev_X, sx_wall_tmp)
                    delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        write(99, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                    else if (IFLAG_X == 2) then
                        write(*, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        write(99, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                    endif
                    ! sx_wall_part(j, k) = sx_wall_tmp
                    ! cal. and save the last grid space

                ! 1 - ensure the smoothness of the intersectional transition of grids
                else if (Mesh_X_trans == 1) then
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "Detailed X mode: smooth transition"
                        write(99, *) "Detailed X mode: smooth transition"
                    else if (IFLAG_X == 2) then
                        write(*, *) "Detailed Y mode: smooth transition"
                        write(99, *) "Detailed Y mode: smooth transition"
                    endif

                    ! use the given "alfax" in "parax_array" to cal. the first grid step by default and cal. distribution,
                    ! i.e. the first "deltax" =  ((exp(alfax / (n - 1)) - 1) / (exp(alfax) - 1));
                    if (k == 1) then

                        call gets_exp_free(nx_tmp, parax_array(k), dev_X, sx_wall_tmp)
                        delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))
                        
                        if (IFLAG_X == 1) then
                            write(*, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                            write(99, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        else if (IFLAG_X == 2) then
                            write(*, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                            write(99, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        endif
                        
                    ! in the case that the present part is NOT the first blunt head (i.e. Part_ID > 1), 
                    ! then cal. new "alfax" to ensure the first grid step is equal to the final grid step of the previous part;
                    else if (k > 1) then

                        ! the first grid that meets the requirement is equal to the last grid in the previous section
                        nx_pre = (id_x(k) - id_x(k - 1)) + 1
                        delta_target = (SLX_part(k - 1) * (sx_wall_part(nx_pre, k - 1) - sx_wall_part(nx_pre - 1, k - 1))) / SLX_part(k)
                        delta_first = SLX_part(k) * delta_target
                        
                        call gets_exp_fixdelta(nx_tmp, delta_target, dev_X, sx_wall_tmp, alfax_get)
                        
                        if (IFLAG_X == 1) then
                            write(*, *) "X Section", k, "Mode: Exp. + fixed delta alfat = ", alfax_get
                            write(99, *) "X Section", k, "Mode: Exp. + fixed delta alfat = ", alfax_get
                        else if (IFLAG_X == 2) then
                            write(*, *) "Y Section", k, "Mode: Exp. + fixed delta alfat = ", alfax_get
                            write(99, *) "Y Section", k, "Mode: Exp. + fixed delta alfat = ", alfax_get
                        endif

                    endif

                endif

                ! employ grid dense strategy
                if (Mesh_X_dense(k) == 1) then
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "Detailed X mode: densed grids on both sides"
                        write(99, *) "Detailed X mode: densed grids on both sides"
                    else if (IFLAG_X == 2) then
                        write(*, *) "Detailed Y mode: densed grids on both sides"
                        write(99, *) "Detailed Y mode: densed grids on both sides"
                    endif
                    
                    call gets_exp_densegrid(nx_tmp, parax_array(k), dev_X, eta_X(k), As_X(k), sx_wall_tmp, alfat_get)
                    delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))
                    delta_final = SLX_part(k) * (sx_wall_tmp(nx_tmp) - sx_wall_tmp(nx_tmp - 1))
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "X Section", k, "Mode: Exp. + densed grid alfat = ", alfat_get
                        write(99, *) "X Section", k, "Mode: Exp. + densed grid alfat = ", alfat_get
                    else if (IFLAG_X == 2) then
                        write(*, *) "Y Section", k, "Mode: Exp. + densed grid alfat = ", alfat_get
                        write(99, *) "Y Section", k, "Mode: Exp. + densed grid alfat = ", alfat_get
                    endif

                endif

            ! Mode 3: Linear grid stretch
            else if (Mesh_X_TYPE(k) == 3) then
                
                if (IFLAG_X == 1) then
                    write(*, *) "Enter Mesh_X generation mode: Linear grid stretch"
                    write(99, *) "Enter Mesh_X generation mode: Linear grid stretch"
                else if (IFLAG_X == 2) then
                    write(*, *) "Enter Mesh_Y", IFLAG_CNT, "generation mode: Linear grid stretch"
                    write(99, *) "Enter Mesh_Y", IFLAG_CNT, "generation mode: Linear grid stretch"
                endif

                ! choose mode according to the grid transition flag
                ! 0 - respectively generate grids
                !     use "betax" in "parax_array", the coeff. of the exp. growth, to cal. grid distribution of each section;
                if (Mesh_X_trans == 0) then

                    call gets_lin_free(nx_tmp, parax_array(k), sx_wall_tmp)
                    delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "X Section", k, "Mode: Lin. + free given alfat = ", parax_array(k)
                        write(99, *) "X Section", k, "Mode: Lin. + free given alfat = ", parax_array(k)
                    else if (IFLAG_X == 2) then
                        write(*, *) "Y Section", k, "Mode: Lin. + free given alfat = ", parax_array(k)
                        write(99, *) "Y Section", k, "Mode: Lin. + free given alfat = ", parax_array(k)
                    endif
                    ! sx_wall_part(k) = sx_wall_tmp
                    ! cal. and save the last grid space

                ! 1 - ensure the smoothness of the intersectional transition of grids
                else if (Mesh_X_trans == 1) then
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "Detailed X mode: smooth transition"
                        write(99, *) "Detailed X mode: smooth transition"
                    else if (IFLAG_X == 2) then
                        write(*, *) "Detailed Y mode: smooth transition"
                        write(99, *) "Detailed Y mode: smooth transition"
                    endif

                    ! use the given "alfax" in "parax_array" to cal. the first grid step by default and cal. distribution,
                    ! i.e. the first "deltax" =  ((exp(alfax / (n - 1)) - 1) / (exp(alfax) - 1));
                    if (k == 1) then

                        call gets_lin_free(nx_tmp, parax_array(k), sx_wall_tmp)
                        delta_target = (SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))) / SLX_part(k)
                        
                        if (IFLAG_X == 1) then
                            write(*, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                            write(99, *) "X Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        else if (IFLAG_X == 2) then
                            write(*, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                            write(99, *) "Y Section", k, "Mode: Exp. + free given alfat = ", parax_array(k)
                        endif
                        
                    ! in the case that the present part is NOT the first blunt head (i.e. Part_ID > 1), 
                    ! then cal. new "betax" to ensure the first grid step is equal to the final grid step of the previous part;
                    else if (k > 1) then

                        ! the first grid that meets the requirement is equal to the last grid in the previous section
                        nx_pre = (id_x(k) - id_x(k - 1)) + 1
                        delta_target = (SLX_part(k - 1) * (sx_wall_part(nx_pre, k - 1) - sx_wall_part(nx_pre - 1, k - 1))) / SLX_part(k)
                        delta_first = SLX_part(k) * delta_target
                        
                        call gets_lin_fixdelta(nx_tmp, delta_target, sx_wall_tmp, betax_get)
                        
                        if (IFLAG_X == 1) then
                            write(*, *) "X Section ", k, "Mode: Lin. + fixed delta alfat = ", betax_get
                            write(99, *) "X Section ", k, "Mode: Lin. + fixed delta alfat = ", betax_get
                        else if (IFLAG_X == 2) then
                            write(*, *) "Y Section ", k, "Mode: Lin. + fixed delta alfat = ", betax_get
                            write(99, *) "Y Section ", k, "Mode: Lin. + fixed delta alfat = ", betax_get
                        endif

                    endif
                    
                endif

                if (Mesh_X_dense(k) == 1) then
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "Detailed X mode: densed grids on both sides"
                        write(99, *) "Detailed X mode: densed grids on both sides"
                    else if (IFLAG_X == 2) then
                        write(*, *) "Detailed Y mode: densed grids on both sides"
                        write(99, *) "Detailed Y mode: densed grids on both sides"
                    endif

                    call gets_lin_densegrid(nx_tmp, parax_array(k), eta_X(k), As_X(k), sx_wall_tmp, betat_get)
                    delta_first = SLX_part(k) * (sx_wall_tmp(2) - sx_wall_tmp(1))
                    delta_final = SLX_part(k) * (sx_wall_tmp(nx_tmp) - sx_wall_tmp(nx_tmp - 1))
                    
                    if (IFLAG_X == 1) then
                        write(*, *) "X Section", k, "Mode: Lin. + densed grid alfat = ", betat_get
                        write(99, *) "X Section", k, "Mode: Lin. + densed grid alfat = ", betat_get
                    else if (IFLAG_X == 2) then
                        write(*, *) "Y Section", k, "Mode: Lin. + densed grid alfat = ", betat_get
                        write(99, *) "Y Section", k, "Mode: Lin. + densed grid alfat = ", betat_get
                    endif

                endif

            endif

            do j = 1, nx_tmp
                sx_wall_part(j, k) = sx_wall_tmp(j)
            enddo
            
            ! delta_first_min = MAX_INF
            ! delta_first_max = MIN_INF
            
            ! output info.
            if (IFLAG_X == 1) then
                
                write(*, *) "X Section", k, "First grid height = ", delta_first
                write(99, *) "X Section", k, "First grid height = ", delta_first
                delta_first_array(k) = delta_first
                
                if (Mesh_X_dense(k) == 1) then
                    write(*, *) "X Section", k, "Last grid height = ", delta_final
                    write(99, *) "X Section", k, "Last grid height = ", delta_final
                    delta_final_array(k) = delta_final
                endif
                
            else if (IFLAG_X == 2) then
                
                write(*, *) "Y Section", k, "First grid height = ", delta_first
                write(99, *) "Y Section", k, "First grid height = ", delta_first
                delta_first_array(k) = delta_first
                
                if (Mesh_X_dense(k) == 1) then
                    write(*, *) "Y Section", k, "Last grid height = ", delta_final
                    write(99, *) "Y Section", k, "Last grid height = ", delta_final
                    delta_final_array(k) = delta_final
                endif
                
            endif

        enddo

        ! combined all the fractions into one part
        ! sx_wall_part: [0, 1]
        nx_cnt = 0
        do k = 1, num_x_part
            do j = 1, (id_x(k + 1) - id_x(k)) + 1
                
                if ((k > 1) .and. (j == 1)) then
                    continue
                else
                    nx_cnt = nx_cnt + 1
                endif

                if (k == 1) then
                    sx_wall_tot(nx_cnt) = (sx_wall_part(j, k) * SLX_part(k)) / SLX_total
                else
                    sx_wall_tot(nx_cnt) = ((sx_wall_part(j, k) * SLX_part(k)) / SLX_total) + sx_wall_tot(id_x(k))
                endif
                
            enddo
        enddo

        do i = nx1 + 1, nx
            sx_wall_tot(i) = sx_wall_tot(i - 1) + alfax_buff * (sx_wall_tot(i - 1) - sx_wall_tot(i - 2))
        enddo

        ! x grid space = sx(k + 1) - sx(k), k = 1, 2, 3, ..., nx - 1
        ! x_delta_1 = 
        ! open(55, file = "sx.dat")
        ! do k = 1, nx - 1
        !     write(55, *) k, sx(k), sx(k + 1) - sx(k)
        ! enddo

    end subroutine

    subroutine gets_exp_free(nx, alfax, dev, sx)
    implicit doubleprecision(a - h, o - z)
    dimension sx(nx)

        dx = 1.d0 / (nx - 1.d0)
        A = 1.d0 / (exp(alfax) - 1.d0)
        do i = 1, nx
            s = (i - dev) * dx
            sx(i) = A * (exp(s * alfax) - 1.d0)
        enddo

    end subroutine

    subroutine gets_exp_fixdelta(nx, deltax, dev, sx, alfax_get)
    implicit doubleprecision (a - h, o - z)
    dimension sx(nx)

        dx = 1.d0 / (nx - 1.d0)
        ! initial guess coefficient
        ! b = 3.5d0
        b = 1.005d0
        ! deltax = delta / SL(i)
    
        ! using Newton method get coefficient b
        flag_find_sol = 0
        do while (flag_find_sol == 0)

            !  fb = exp(b / (ny - 1.d0)) - 1.d0 - delta * (exp(b) - 1.d0)
            !  fbx = exp(b / (ny - 1.d0)) / (ny - 1.d0) - delta * exp(b)
            fb = (exp(b / (nx - 1.d0)) - dev) / (exp(b) - 1.d0) - deltax
            fbx = (exp(b / (nx - 1.d0)) / (nx - 1.d0) * (exp(b) - 1.d0) -  &
                  (exp(b / (nx - 1.d0)) - dev) * exp(b)) / ((exp(b) - 1.d0)) ** 2

            bnew = b - fb / fbx
            
            if (abs(b - bnew) .gt. 1.d-6) then
                b = bnew
            else
                b_final = bnew
                flag_find_sol = 1
            endif
        
        enddo
        
        ! itertive result coefficient b
        A = 1.d0 / (exp(b_final) - 1.d0)
        do i = 1, nx
            s = (i - dev) * dx
            sx(i) = A * (exp(s * b_final) - 1.d0)
        enddo
        
        alfax_get = b_final

        ! !-------------- revise ----------
        !  reverse the array 
        ! do j = 1, ny
        !     j1 = ny + 1 - j
        !     sy(j) = sy1(ny) - sy1(j1)
        ! enddo

    end subroutine

    subroutine gets_exp_densegrid(nx, alfas, dev, eta, As, sx, alfat_get)
    implicit doubleprecision (a - h, o - z)
    real*8, parameter:: PI = 3.1415926535897932
    real*8, parameter:: MAX_INF = 99999.d0
    integer, parameter:: USER_PARA = 100, USER_LEN = 2000
    
    dimension sx(nx), sx_err(nx)
    ! real*8:: sx_err(USER_LEN)

        sx_err_min = MAX_INF
        do i = 1, nx
            ! sx_err(i) = abs(As * ((exp(((i - 1) / (n - 1)) * alfax) - 1) / (exp(alfax) - 1)) - eta);
            sx_err(i) = abs(As * ((exp(((i - 1.d0) / (nx - 1.d0)) * alfas) - dev) / (exp(alfas) - 1.d0)) - eta);
            if (sx_err(i) < sx_err_min) then
                nd = i
                sx_err_min = sx_err(i)
            endif

        enddo

        an = ((nx - nd) * 1.d0) / (nx - 1.d0)

        ! get the desired alfat using the Newton iterative method
        ! set the initial value 
        alfab = alfas

        flag_find_sol = 0
        do while (flag_find_sol == 0)

            fb = (exp(an * alfab) - dev) / (exp(alfab) - 1.d0) - ((1.d0 - eta) / As)
            fbx = ((an * exp(an * alfab) * (exp(alfab) - 1.d0)) -  &
                   (exp(alfab) * (exp(an * alfab) - dev))) / ((exp(alfab) - 1) ** 2)

            alfab_new = alfab - (fb / fbx)

            if (abs(alfab - alfab_new) > 1.d-6) then
                alfab = alfab_new
            else
                alfat = alfab_new
                flag_find_sol = 1
            endif

        enddo

        do i = 1, nx
            if (i <= nd) then
                sx(i) = As * ((exp(((i - 1.d0) / (nx - 1.d0)) * alfas) - dev) / (exp(alfas) - 1.d0))
            else
                sx(i) = (1.d0 - As * ((exp((((nx - i) * 1.d0) / (nx - 1.d0)) * alfat) - dev) / (exp(alfat) - 1.d0)))
            endif
        enddo
        
        alfat_get = alfat

    end subroutine

    subroutine gets_lin_free(nx, betax, sx)
    implicit doubleprecision (a - h, o - z)
    dimension sx(nx)

        deltax = (1.d0 - betax) / (1.d0 - betax ** (nx - 1.d0))
        do i = 1, nx
            sx(i) = deltax * ((1.d0 - betax ** (i - 1.d0)) / (1.d0 - betax))
        enddo

    end subroutine

    subroutine gets_lin_fixdelta(nx, deltax, sx, betax_get)
    implicit doubleprecision (a - h, o - z)
    dimension sx(nx)

        ! initial guess coefficient
        b = 1.005d0

        ! deltax = delta / SL(i)
    
        ! using Newton method get coefficient b
        flag_find_sol = 0
        do while (flag_find_sol == 0)

            !  fb = exp(b / (ny - 1.d0)) - 1.d0 - delta * (exp(b) - 1.d0)
            !  fbx = exp(b / (ny - 1.d0)) / (ny - 1.d0) - delta * exp(b)
            fb = deltax * ((1.d0 - b ** (nx - 1.d0)) / (1.d0 - b)) - 1.d0
            fbx = (deltax * (1.d0 - ((nx - 1.d0) * (b ** (nx - 2.d0))) + ((nx - 2.d0) * (b ** (nx - 1.d0))))) / ((1.d0 - b) ** 2)

            bnew = b - fb / fbx
            
            if (abs(b - bnew) .gt. 1.d-6) then
                b = bnew
            else
                beta_final = bnew
                flag_find_sol = 1
            endif
        
        enddo
        
        ! itertive result coefficient b
        do i = 1, nx
            sx(i) = deltax * ((1.d0 - beta_final ** (i - 1.d0)) / (1.d0 - beta_final))
        enddo
        
        betax_get = beta_final

    end subroutine

    subroutine gets_lin_densegrid(nx, betas, eta, As, sx, betat_get)
    implicit doubleprecision (a - h, o - z)
    dimension sx(nx), sx_err(nx)

        sx_err_min = MAX_INF
        do i = 1, nx

            sx_err(i) = abs((As * (((1.d0 - betas) ** (i - 1.d0)) / (1.d0 - betas))) - eta);
            if (sx_err(i) < sx_err_min) then
                nd = i
                sx_err_min = sx_err(i)
            endif

        enddo

        an = nx - nd - 1.d0

        ! get the desired betat using the Newton iterative method
        ! set the initial value
        betab = betas

        flag_find_sol = 0
        do while (flag_find_sol == 0)

            fb =  ((1.d0 - betab ** (i - 1.d0)) / (1.d0 - betab)) - ((1.d0 - eta) / As)
            fbx = (1.d0 - (an * (betab ** (an - 1.d0))) + ((an - 1.d0) * (betab ** an))) / ((1.d0 - betab) ** 2)

            betab_new = betab - (fb / fbx)

            if (abs(betab - betab_new) > 1.d-6) then
                betab = betab_new
            else
                betat = betab_new
                flag_find_sol = 1
            endif

        enddo

        do i = 1, nx
            if (i <= nd) then
                sx(i) = As * ((1.d0 - betat ** (i - 1.d0)) / (1.d0 - betat))
            else
                sx(i) = (1.d0 - As * ((1.d0 - betat ** (nx - i)) / (1.d0 - betat)))
            endif
        enddo
        
        betat_get = betat

    end subroutine
    
    !-------------------------------------------
    ! complex simpson integral formula
    subroutine simpson_curve(s, ya, yb, IFLAG_curve)
    implicit doubleprecision (a - h, o - z)
    integer, parameter:: ny = 100
    ! common/para/ b2, d2, setab, xbs, ybs
    ! common/para/ c1, c2, c3, c4, c5, c6, xcs, ycs, setac
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6

        dy = (yb - ya) / (ny - 1.d0)
        s1 = 0.d0
        s2 = 0.d0
        do k = 2, ny
            y = ya + (k - 1.d0) * dy - (dy / 2.d0)
            if (IFLAG_curve == 1) then
                s1 = s1 + funb(y)
            else if (IFLAG_curve == 2) then
                s1 = s1 + func(y)
            endif
        enddo
        do k = 2, ny - 1
            y = ya + (k - 1.d0) * dy
            if (IFLAG_curve == 1) then
                s2 = s2 + funb(y)
            else if (IFLAG_curve == 2) then
                s2 = s2 + func(y)
            endif
        enddo
        if (IFLAG_curve == 1) then
            s = (funb(ya) + funb(yb) + 4.d0 * s1 + 2.d0 * s2) * dy / 6.d0
        else if (IFLAG_curve == 2) then
            s = (func(ya) + func(yb) + 4.d0 * s1 + 2.d0 * s2) * dy / 6.d0
        endif

    end

    ! Farfield: parabolic-line curve dsy = sqrt(1 + x' ** 2)
    ! cal. the arc length differential
    ! Note: integal by x direction
	function funb(y)
	implicit doubleprecision (a - h, o - z)
    ! common/para/ b2, d2, setab, xbs, ybs
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6

        ! funs = cos(x)
        ! parabolic: x = b2 * y * y + d2
        if (y .lt. ybs) then
            ! xd = 2.d0 * b2 * y
            funb = sqrt(1.d0 + ((2.d0 * b2 * y) ** 2))  ! x is r ; located in the parabolic line
        else
            funb = 1.d0 / sin(setab)                 ! linear  
        endif

    end
    
    ! Shock wave: dsy = sqrt(1 + x' ** 2)
    ! cal. the arc length differential
    ! Note: integal by x direction
	function func(y)
	implicit doubleprecision (a - h, o - z)
    ! common/para/ c1, c2, c3, c4, c5, c6, xcs, ycs, setac
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6

        ! funs = cos(x)
        ! parabolic: y = b * x * x + d
        if (y .lt. ycs) then
            ! y = c1 * (((c2 * x + c3) ** c4 - 1.d0) ** c5) + c6
            xd = (((1 + ((y - c6) / c1) ** (1.d0 / c5)) ** ((1.d0 / c4) - 1)) * (((y - c6) / c1) ** ((1.d0 / c5) - 1))) / (c1 * c2 * c4 * c5 * 1.d0)
            func = sqrt(1.d0 + (xd ** 2))            ! shock wave fitting curve
        else
            func = 1.d0 / sin(setac)                 ! linear section
        endif

	end

    !---------------------------------------------
    subroutine get_y_form_ss(sa, y, y0, IFLAG_curve)
    implicit doubleprecision (a - h, o - z)
    ! real*8, parameter:: dx = 0.001d0
    ! common/para/ d2, xt, rn
    common/para/ rn, theta1, R, x3c, x4c, a, b, c, xm, xt, ym, yt, b2, d2, setab, setac
    common/para/ xbs, ybs, xcs, ycs, c1, c2, c3, c4, c5, c6
	
        y = y0

100     continue
		
        if (IFLAG_curve == 1) then
            call simpson_curve(s, 0.d0, y, IFLAG_curve)
            sy = funb(y)
        else if (IFLAG_curve == 2) then
            call simpson_curve(s, 0.d0, y, IFLAG_curve)
            sy = func(y)
        endif
        ynew = y - (s - sa) / sy
        
        if (abs(y - ynew) .gt. 1.d-6) then
            y = ynew
		    goto 100
        endif
        
        y = ynew

    end
    
    !--------------------------------------------
    ! 调整外边界的网格，使得连接点附近弧长间距光滑过渡
    ! 采用三次函数来过渡
    subroutine adjust_outter(nx, ss, ss_new, i_conjunction, nxconjunction_buffer)
    implicit doubleprecision (a - h, o - z)
    real*8 ss(nx), ss_new(nx)

        ss_new = ss
        
        n2 = nxconjunction_buffer / 2
        ia = i_conjunction - n2
        ib = i_conjunction + n2

        det1 = ss(ia) - ss(ia - 1)
        det2 = ss(ib + 1) - ss(ib)
        ds = 1.d0 / (2.d0 * n2)
        
        do i = 1, 2 * n2 + 1
        
            s = (i - 1.d0) * ds
            xh0 = (1.d0 + 2.d0 * s) * (s - 1.d0) ** 2
            xh1 = (1.d0 + 2.d0 * (1.d0 - s)) * s * s
            xhf0 = s * (s - 1.d0) ** 2
            xhf1 = (s - 1.d0) * s * s

            ss_new(i + ia - 1) = ss(ia) * xh0 + ss(ib) * xh1 + det1 / ds * xhf0 + det2 / ds * xhf1
        
        enddo

     !   open(55, file = "ss.dat")
     !   do k = 2, nx
     !       write(55, "(I3,4f16.8)") k, ss(k), ss_new(k), & 
     !                                ss(k) - ss(k - 1),  &
     !                                ss_new(k) - ss_new(k - 1)
	    !enddo

	end
    
    subroutine get_Jacobian(nx, ny, xx, yy)
    implicit doubleprecision(a - h, o - z)
    parameter (LAP = 3)
    dimension xx(nx, ny), yy(nx, ny)
	dimension Akx(nx, ny), Aky(nx, ny), Aix(nx, ny), Aiy(nx, ny), Ajac(nx, ny)
	dimension xk(nx, ny), xi(nx, ny), yk(nx, ny), yi(nx, ny)
	dimension xx1(1 - LAP : nx, ny), yy1(1 - LAP : nx, ny)
    dimension d(nx, ny), u(nx, ny), v(nx, ny), T(nx, ny)

        hx = 1.d0 / (nx - 1.d0)
        hy = 1.d0 / (ny - 1.d0)

        do j = 1, ny
            do i = 1, nx
                xx1(i, j) = xx(i, j)
                yy1(i, j) = yy(i, j)
            enddo
        enddo
        do j = 1, ny
            do i = 1, LAP
                i1 = 1 - i
                xx1(i1, j) = - xx(i, j)
                yy1(i1, j) = yy(i, j)
            enddo
        enddo

        call dx0(xx1, xk, nx, ny, hx, LAP)
        call dx0(yy1, yk, nx, ny, hx, LAP)
        call dy0(xx, xi, nx, ny, hy)
        call dy0(yy, yi, nx, ny, hy)
        
        Ajac = xk * yi - xi * yk
        Ajac = 1. / Ajac
        Akx = Ajac * yi
        Aky = - Ajac * xi
        Aix = - Ajac * yk
        Aiy = Ajac * xk

        open(55, file = 'OCFD2d-Jacobi.dat', form = 'unformatted')
        write(55) xx
        write(55) yy
        write(55) Akx
        write(55) Aky
        write(55) Aix
        write(55) Aiy
        write(55) Ajac
        close(55)

        do j = 1, ny
            do i = 1, nx
                if (Ajac(i, j) .lt. 1.e-5) print*, i, j, Ajac(i, j)
            enddo
        enddo

        ! Write the final mesh
        open(33, file = 'grid.dat')
        write(33, *) 'variables = x, y, Akx, Aky, Aix, Aiy, Ajac'
        write(33, *) 'zone i = ', nx , 'j = ', ny
        do j = 1, ny
            do i = 1, nx
                write(33, '(7f15.6)') yy(i, j), xx(i, j), Akx(i, j), Aky(i, j),  &
                                      Aix(i, j), Aiy(i, j), Ajac(i, j)
            enddo
        enddo
        close(33)

        ! Generate initial data ... 
        Istep = 0
        tt = 0.
        d = 1.
        u = 0.
        v = 1.
        T = 1.

        open(44, file = 'cone2d0.dat', form = 'unformatted')
        write(44) Istep, tt
        write(44) d
        write(44) u
        write(44) v
        write(44) T
        close(44)

    end
    
    subroutine dx0(f, fx, nx, ny, hx, LAP)
    implicit doubleprecision (a - h, o - z)
    dimension f(1 - LAP : nx, ny), fx(nx, ny)

        b1 = 8.d0 / (12.d0 * hx)
        b2 = 1.d0 / (12.d0 * hx)
        a1 = 1.d0 / (60.d0 * hx)
        a2 = - 3.d0 / (20.d0 * hx)
        a3 = 3.d0 / (4.d0 * hx)

        do j = 1, ny
            do i = 1, nx - 3
                fx(i, j) = a1 * (f(i + 3, j) - f(i - 3, j)) + a2 * (f(i + 2, j) - f(i - 2, j)) + a3 * (f(i + 1, j) - f(i - 1, j))
            enddo
        enddo

        do j = 1, ny
            fx(nx - 2, j) = b1 * (f(nx - 1, j) - f(nx - 3, j)) - b2 * (f(nx, j) - f(nx - 4, j))
            fx(nx - 1, j) = (f(nx - 3, j) - 6.d0 * f(nx - 2, j) + 3.d0 * f(nx - 1, j) + 2.d0 * f(nx, j)) / (6.d0 * hx)
            fx(nx, j) = (f(nx - 2, j) - 4.d0 * f(nx - 1, j) + 3.d0 * f(nx, j)) / (2.d0 * hx)
        enddo

    end

    !----------------------------------------------
    subroutine dy0(f, fy, nx, ny, hy)
    implicit doubleprecision (a - h, o - z)
    dimension f(nx, ny), fy(nx, ny)
    
    b1 = 8.d0 / (12.d0 * hy)
    b2 = 1.d0 / (12.d0 * hy)
    a1 = 1.d0 / (60.d0 * hy)
    a2 = - 3.d0 / (20.d0 * hy)
    a3 = 3.d0 / (4.d0 * hy)

    do j = 4, ny - 3
        do i = 1, nx
            fy(i, j) = a1 * (f(i, j + 3) - f(i, j - 3)) + a2 * (f(i, j + 2) - f(i, j - 2)) + a3 * (f(i, j + 1) - f(i, j - 1))
        enddo
    enddo

    do i = 1, nx

        fy(i, 1) = (- 3.d0 * f(i, 1) + 4.d0 * f(i, 2) - f(i, 3)) / (2.d0 * hy)            
        fy(i, 2) = (- 2.d0 * f(i, 1) - 3.d0 * f(i, 2) + 6.d0 * f(i, 3) - f(i, 4)) / (6.d0 * hy)
        fy(i, 3) = b1 * (f(i, 4) - f(i, 2)) - b2*(f(i,5)-f(i,1))

        fy(i, ny - 2) = b1 * (f(i, ny - 1) - f(i, ny - 3)) - b2 * (f(i, ny) - f(i, ny - 4))
        fy(i, ny - 1) = (f(i, ny - 3) - 6.d0 * f(i, ny - 2) + 3.d0 * f(i, ny - 1) + 2.d0 * f(i, ny)) / (6.d0 * hy)
        fy(i, ny) = (f(i, ny - 2) - 4.d0 * f(i, ny - 1) + 3.d0 * f(i, ny)) / (2.d0 * hy)

    enddo

    end
    
    ! subroutine find_n_id(n_array)

    ! alfax = 2.5;
    ! n = 500;
    ! eta = 0.5;
    ! err = 1e-6;
    ! As = 2.0;

    ! sx = zeros(n, 1);
    ! sx1 = zeros(n - 1, 1);
    ! for i = 1 : n
    !     sx(i) = As * ((exp(((i - 1) / (n - 1)) * alfax) - 1) / (exp(alfax) - 1));
    ! end

    ! nd = find(abs(sx - eta) == min(abs(sx - eta)))
    ! an = (n - nd) / (n - 1)

    ! % dx = 1 / (n - 1);
    ! alfab = alfax;

    ! flag_find_nd = 0;
    ! while (flag_find_nd == 0)

    !     fb = (exp(an * alfab) - 1) / (exp(alfab) - 1) - ((1 - eta) / As);
    !     fbx = ((an * exp(an * alfab) * (exp(alfab) - 1)) - (exp(alfab) * (exp(an * alfab) - 1))) / ((exp(alfab) - 1) ^ 2);

    !     alfab_new = alfab - (fb / fbx);

    !     if (abs(alfab - alfab_new) > err)
    !         alfab = alfab_new;
    !     else
    !         alfab_final = alfab_new
    !         flag_find_nd = 1;
    !     end

    ! end

    ! for i = 1 : n
    !     if (i <= nd)
    !         sx(i) = As * ((exp(((i - 1) / (n - 1)) * alfax) - 1) / (exp(alfax) - 1));
    !     else
    !         sx(i) = (1 - As * ((exp(((n - i) / (n - 1)) * alfab_final) - 1) / (exp(alfab_final) - 1)));
    !     end
    ! end

    ! for i = 1 : n - 1
    !     sx1(i) = sx(i + 1) - sx(i);
    ! end

    ! sx(nd) - sx(nd - 1)
    ! sx1(end)






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!