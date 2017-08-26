MODULE initialize

  USE globals
  IMPLICIT NONE

  SAVE

  INTEGER :: g, gg, f, ff, long, longP, maskBit
  REAL(4) :: x(xn), y(yn), t(tn)
  REAL(4) :: rho0(xn,yn)
  REAL(4) :: bcx0(xn,2), bcy0(2,yn), bcxPsi(xn,2), bcyPsi(2,yn), ic0(xn,yn)
  REAL(4) :: kMat(xn,yn), lambdaMat(xn,yn), lambdaMatW(xn,yn), porosity(xn,yn), reactive(xn,yn), reactiveCoarse((xn-1)/cellx,yn/celly)
  REAL(4) :: reactiveCoarse1((xn-1)/cellx,yn/celly), reactiveCoarse2((xn-1)/cellx,yn/celly), sedx
  REAL(4) :: sed(xn), sed1(xn), sed2(xn), sed3(xn)

  ! EXTRA BULLSHIT
  REAL(4) :: permeable(xn,4), permeability(xn,yn)
  REAL(4) :: mask(xn,yn), maskLong((xn-2)*(yn-2)), maskLongT((xn-2)*(yn-2)), h0Long((xn-2)*(yn-2)), h0LongT((xn-2)*(yn-2)), h0T(yn,xn)
  REAL(4) :: maskPLongFull(xn*yn)
  REAL(4) :: maskLongU((xn-2)*(yn-0)), maskLongTV((xn-0)*(yn-2))
  REAL(4) :: maskP(xn,yn), maskPLong((xn-2)*((yn/2)-2)), maskPLongT((xn-2)*((yn/2)-2))
  REAL(4) :: outerBand((xn-2)*((yn/2)-2),2*((yn/2)-2) + 1), bigBand((xn-2)*((yn/2)-2),4*((yn/2)-2) + 3)
  REAL(4) :: stretch(xn,yn), stretchLong((xn-2)*(yn-2)), stretchT(xn,yn), stretchLongT((xn-2)*(yn-2))
  REAL(4) :: u_inter(xn,yn), v_inter(xn,yn), u_step_coarse((xn-1)/cellx,yn/celly), v_step_coarse((xn-1)/cellx,yn/celly)

  ! 05/06 INPUT STUFF
  REAL(4) :: primary((xn-1)/cellx,yn/(2*celly),g_pri), primaryMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_pri)
  REAL(4) :: secondary((xn-1)/cellx,yn/(2*celly),g_sec), secondaryMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sec)
  REAL(4) :: solute((xn-1)/cellx,yn/(2*celly),g_sol), soluteMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sol)
  REAL(4) :: solute_inter((xn-1)/cellx,yn/(2*celly)), solute_inter_long(((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: medium((xn-1)/cellx,yn/(2*celly),g_med), mediumMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_med)
  REAL(4) :: saturation((xn-1)/cellx,yn/(2*celly),g_sec/2), saturationMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sec/2)
  REAL(4) :: checkMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly)), checkMatF((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly))
  REAL(4) :: soluteOcean(g_sol), soluteOcean_a(g_sol), soluteOcean_b(g_sol), sea(g_sol)
  REAL(4) :: inter, slope, inter2, slope2, boundary, buffer, edge, edge2
  REAL(4) :: bit_thing((xn-1)/cellx,yn/(2*celly)), bit_thing_t(yn/(2*celly),(xn-1)/cellx), bit_thing_t1((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: solute_fine(xn,yn,g_sol), solute_fine_a(xn,yn,g_sol), solute_fine_b(xn,yn,g_sol)
  ! trial 04/29
  REAL(4) :: bit_thing_a((xn-1)/cellx,yn/(2*celly)), bit_thing_b((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: phiCalc((xn-1)/cellx,yn/(2*celly)), phiCalc_a((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: phiCalcMat((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly)), phiCalcMat_a((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly))

  ! CHAMBER INPUT STUFF
  REAL(4) :: primary_a((xn-1)/cellx,yn/(2*celly),g_pri), primaryMat_a((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_pri)
  REAL(4) :: secondary_a((xn-1)/cellx,yn/(2*celly),g_sec), secondaryMat_a((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sec)
  REAL(4) :: solute_a((xn-1)/cellx,yn/(2*celly),g_sol), soluteMat_a((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sol)
  REAL(4) :: medium_a((xn-1)/cellx,yn/(2*celly),g_med), mediumMat_a((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_med)

  REAL(4) :: primary_b((xn-1)/cellx,yn/(2*celly),g_pri), primaryMat_b((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_pri)
  REAL(4) :: secondary_b((xn-1)/cellx,yn/(2*celly),g_sec), secondaryMat_b((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sec)
  REAL(4) :: solute_b((xn-1)/cellx,yn/(2*celly),g_sol), soluteMat_b((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sol)
  REAL(4) :: medium_b((xn-1)/cellx,yn/(2*celly),g_med), mediumMat_b((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_med)

  REAL(4) :: primaryMat_d((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_pri)
  REAL(4) :: secondaryMat_d((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sec)
  REAL(4) :: soluteMat_d((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_sol)
  REAL(4) :: mediumMat_d((xn-1)*tn/(cellx*(mstep*ar)),yn/(2*celly),g_med)

  REAL(4) :: t_vol_s, t_vol_a, t_vol_b

  REAL(4) :: sec_density(g_sec/2), sec_molar(g_sec/2)
  REAL(4) :: pri_density(g_pri), pri_molar(g_pri)

  ! ! coarse/mean arrays added for super duper coarse grid 03/25/17
  REAL(4) :: h_coarse((xn-1)/cellx,yn/(2*celly))
  ! real(4) :: pri_coarse((xn-1)/cellx,yn/celly,g_pri)
  ! real(4) :: sec_coarse((xn-1)/cellx,yn/celly,g_sec)
  ! real(4) :: sol_coarse((xn-1)/cellx,yn/celly,g_sol)
  ! real(4) :: med_coarse((xn-1)/cellx,yn/celly,g_med)
  !
  ! real(4) :: pri_coarse_a((xn-1)/cellx,yn/celly,g_pri)
  ! real(4) :: sec_coarse_a((xn-1)/cellx,yn/celly,g_sec)
  ! real(4) :: sol_coarse_a((xn-1)/cellx,yn/celly,g_sol)
  ! real(4) :: med_coarse_a((xn-1)/cellx,yn/celly,g_med)
  !
  ! real(4) :: pri_coarse_b((xn-1)/cellx,yn/celly,g_pri)
  ! real(4) :: sec_coarse_b((xn-1)/cellx,yn/celly,g_sec)
  ! real(4) :: sol_coarse_b((xn-1)/cellx,yn/celly,g_sol)
  ! real(4) :: med_coarse_b((xn-1)/cellx,yn/celly,g_med)

  REAL(4) :: coarse_mask((xn-1)/cellx,yn/(2*celly)), coarse_mask_long(((xn-1)/cellx)*(yn/(2*celly)))

  REAL(4) :: volume_ratio, mix_ratio, volume_ratio_i, volume_ratio_u
  REAL(4) :: vol_i, vol_i_a, vol_i_b

  ! coordinates for optimization
  REAL(4) :: fives(2,yn), tens(2,yn), twentyfives(2,xn), fifties(2,xn)

  ! shell file parameters
  CHARACTER(len=300) :: path, path2, path_final, crashstring, restartstring, iso_path
  CHARACTER(len=300) :: param_o_string, param_w_string, param_w_rhs_string, param_h_string, param_o_rhs_string, param_tsw_string
  CHARACTER(len=300) :: param_dic_string, param_scope_string, param_trace_string, param_ch_string, param_f_dx_string, param_f_k_string
  CHARACTER(len=300) :: param_paq_string, param_ch_rhs_string, param_f_freq_string, param_f_por_string
  CHARACTER(len=300) :: param_t_diff_string
  INTEGER :: in, crashstep, restart, param_trace
  REAL(4):: param_o, param_w, param_w_rhs, param_h, param_o_rhs, param_tsw, param_dic, param_scope, param_ch
  REAL(4) :: param_paq, param_ch_rhs, param_f_dx, param_f_k, param_f_freq, param_f_por
  REAL(4) :: param_t_diff


  ! TRANSPOSED
  REAL(4) :: hTrans(yn,xn), psiTrans(yn,xn), permeabilityTrans(yn,xn), phiTrans(yn,xn)
  REAL(4) :: primaryTrans(yn/celly,(xn-1)/cellx,g_pri), secondaryTrans(yn/celly,(xn-1)/cellx,g_sec), soluteTrans(yn/celly,(xn-1)/cellx,g_sol), mediumTrans(yn/celly,(xn-1)/cellx,g_med)
  REAL(4) :: primaryTrans_a(yn/celly,(xn-1)/cellx,g_pri), secondaryTrans_a(yn/celly,(xn-1)/cellx,g_sec), soluteTrans_a(yn/celly,(xn-1)/cellx,g_sol), mediumTrans_a(yn/celly,(xn-1)/cellx,g_med)
  REAL(4) :: primaryTrans_b(yn/celly,(xn-1)/cellx,g_pri), secondaryTrans_b(yn/celly,(xn-1)/cellx,g_sec), soluteTrans_b(yn/celly,(xn-1)/cellx,g_sol), mediumTrans_b(yn/celly,(xn-1)/cellx,g_med)
  REAL(4) :: genTrans(yn,xn)


  ! D2 STUFF
  INTEGER :: xn2, yn2, long2, m2
  REAL(4) :: perm2 = 1e-12


  REAL(4) :: frac6(yn,2), frac6_last(yn,2), temp6(yn,2), temp6_last(yn,2), temp6_mid(yn,2)
  INTEGER :: f_index1 = xn-60, iter = 0, spinup = 20!50000


  REAL(4) :: temp6_a(yn), temp6_b(yn), temp6_c(yn), temp6_rhs(yn)

  ! SOLUTE ADVECTION DISTRIBUTION
  INTEGER :: sol_index(11)





CONTAINS

  SUBROUTINE init_mini ()
    USE globals

    sec_density = (/2.65, 2.3, 3.05, 2.17, 5.01, 2.5, 3.8, & ! 7
         & 2.7, 2.71, 2.56, 2.3, 2.28, 2.28, 3.05, 2.28, & ! 8
         & 2.25, 5.3, 2.5, 2.55, 2.27, 2.2, 2.5, 2.26, & ! 8
         & 2.55, 2.25, 2.75, 2.7, 2.87, 2.9, 2.275, 2.8, & !8
         & 2.8, 2.3, 2.55, 4.61, 2.3, 2.3/) ! 6

    pri_density = (/1.0, 2.7, 3.0, 3.0, 3.0/)


    sec_molar = (/258.156, 480.19, 429.02, 2742.13, 119.98, 549.07, 88.851, &
         & 549.07, 100.0869, 287.327, 480.19, 495.90, 495.90, 429.02, 495.90, &
         & 380.22, 159.6882, 549.07, 504.19, 220.15, 649.86, 549.07, 649.86, &
         & 504.19, 380.22, 379.259, 549.07, 395.38, 64.448, 392.34, 64.448, &
         & 64.448, 480.19, 504.19, 85.12, 480.19, 480.19/)

    ! pri_molar = (/1.0, 110.0, 153.0, 158.81, 277.0/)
    pri_molar = (/1.0, 277.0, 153.0, 158.81, 110.0/)


    dx = ( x_max - x_min ) / REAL ( xn - 1, kind = 4 )
    x = linspace ( xn, x_min,x_max )
    dy = ( y_max - y_min ) / REAL ( yn - 1, kind = 4 )
    y = linspace ( yn, y_min, y_max )
    dt = ( t_max - t_min ) / REAL ( tn - 1, kind = 4 )
    t = linspace ( tn, t_min, t_max)


    sea = (/8.2, 0.00243, 0.0, 0.0021, 0.01028, 0.0528, 0.460, 0.00995, 0.0, 0.028, 0.0, 0.540, 0.0, 0.00245, 0.0/)
    !         1    2       3     4         5      6       7       8      9     10    11   12      13      14      15

    sol_index = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/)
    !            alk dic ca mg na k  fe s   si  cl  al

    RETURN

  END SUBROUTINE init_mini

  ! ----------------------------------------------------------------------------------%%
  !
  ! SUBROUTINE TO INITIALIZE, JUST CALL IT
  !
  ! ----------------------------------------------------------------------------------%%

  SUBROUTINE init ()
    USE globals
    INTEGER :: m,n
    !integer :: i,ii


    in = iargc()
    CALL getarg(1,restartstring)
    CALL getarg(2,path)
    CALL getarg(3,path_final)
    CALL getarg(4,crashstring)
    CALL getarg(5,param_o_string)
    CALL getarg(6,param_w_string)
    CALL getarg(7,param_w_rhs_string)
    CALL getarg(8,param_h_string)
    CALL getarg(9,param_o_rhs_string)
    CALL getarg(10,param_tsw_string)
    CALL getarg(11,param_dic_string)
    CALL getarg(12,param_scope_string)
    CALL getarg(13,param_trace_string)
    CALL getarg(14,param_ch_string)
    CALL getarg(15,param_paq_string)
    CALL getarg(16,param_ch_rhs_string)
    CALL getarg(17,param_f_dx_string)
    CALL getarg(18,param_f_k_string)
    CALL getarg(19,param_f_freq_string)
    CALL getarg(20,param_f_por_string)
    CALL getarg(21,iso_path)
    CALL getarg(22,param_t_diff_string)

    READ (crashstring, *) crashstep
    READ (restartstring, *) restart
    READ (param_o_string, *) param_o
    param_o = param_o
    READ (param_w_string, *) param_w
    READ (param_w_rhs_string, *) param_w_rhs
    READ (param_h_string, *) param_h
    param_h = param_h
    READ (param_o_rhs_string, *) param_o_rhs
    param_o_rhs = param_o_rhs
    READ (param_tsw_string, *) param_tsw
    param_tsw = param_tsw/10.0
    READ (param_dic_string, *) param_dic
    param_dic = param_dic/100000.0
    READ (param_scope_string, *) param_scope
    param_scope = 1.0*param_scope*(1e-10)
    scope = param_scope

    scope = 1.0*scope
    param_scope = 1.0*param_scope

    READ (param_trace_string, *) param_trace
    READ (param_ch_string, *) param_ch
    param_ch = -1.0*param_ch

    READ (param_paq_string, *) param_paq
    READ (param_ch_rhs_string, *) param_ch_rhs

    READ (param_f_dx_string, *) param_f_dx

    param_f_dx = 10.0**(param_f_dx)


    READ (param_f_k_string, *) param_f_k
    READ (param_f_freq_string, *) param_f_freq
    param_f_freq = param_f_dx/param_f_freq

    READ (param_f_por_string, *) param_f_por
    READ (param_t_diff_string, *) param_t_diff



    !permf = rho_fluid*grav*4.0*param_f_dx*param_f_dx/(12.0*viscosity)
    permf = param_f_dx*param_f_dx/3.0

    ! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
    dx = ( x_max - x_min ) / REAL ( xn - 1, kind = 4 )
    x = linspace ( xn, x_min,x_max )
    dy = ( y_max - y_min ) / REAL ( yn - 1, kind = 4 )
    y = linspace ( yn, y_min, y_max )
    dt = ( t_max - t_min ) / REAL ( tn - 1, kind = 4 )
    t = linspace ( tn, t_min, t_max)

    ! BOUNDARY CONDITIONS
    ic0(:,:) = 273.16 ! IC
    bcx0(:,1) = 273.16 ! bottom
    bcx0(:,2) = 273.16 ! top
    bcy0(1,:) = 273.16 ! left
    bcy0(2,:) = 273.16 ! right
    bcyPsi(1,1:yn) = 0.0 ! left
    bcyPsi(2,1:yn) = 0.0 ! right
    bcxPsi(1:xn,1) = 0.0 ! bottom
    bcxPsi(1:xn,2) = 0.0 ! top

    ! PERMEABILITY
    lambdaMat = 1.8
    permeability = 1e-15

    ! HEAT TRANSFER PARAMETERS
    kMat = 2.0/(1000.0*4186.0)
    ki=2.0/(1000.0*4186.0)


    !-Primary initial conditions

    ! t_vol_s = 0.30
    ! t_vol_a = 0.28
    ! t_vol_b = 0.02

    ! t_vol_s = 0.01156
    ! t_vol_a = 0.01156 - 0.001156
    ! t_vol_b = 0.001156

    ! t_vol_s = 0.0112
    ! t_vol_a = 0.0112 - 0.00112
    ! t_vol_b = 0.00112

    ! t_vol_s = 0.0112
    ! t_vol_a = 0.0112 - 0.00112
    ! t_vol_b = 0.00112

    ! good alk values
    ! t_vol_s = 0.012
    ! t_vol_a = 0.008! - 0.0112/2.0
    ! t_vol_b = 0.004!/2.0

    t_vol_s = 0.012*1.5
    t_vol_a = 0.009*1.5! - 0.0112/2.0
    t_vol_b = 0.003*1.5!/2.0

    ! primary minerals [mol]
    primary(:,:,:) = 0.0
    primary(:,:,1) = 0.0 !1.29600 ! feldspar
    primary(:,:,2) = 0.1 !.69600 ! plag
    primary(:,:,3) = 0.1 !.12600 ! pyr
    primary(:,:,4) = 0.95068 !.04000 ! ol
    primary(:,:,5) = 1.0 !9.67700 ! basaltic glass

    primary_a(:,:,:) = 0.0
    primary_a(:,:,1) = 0.0   ! feldspar
    primary_a(:,:,2) = 0.0   ! plag
    primary_a(:,:,3) = 0.0   ! pyr
    primary_a(:,:,4) = 0.0   ! ol
    primary_a(:,:,5) = 1.00  ! basaltic glass

    primary_b(:,:,:) = 0.0
    primary_b(:,:,1) = 0.0  ! feldspar
    primary_b(:,:,2) = 0.1 ! plag
    primary_b(:,:,3) = 0.1 ! pyr
    primary_b(:,:,4) = 0.95068 ! ol
    primary_b(:,:,5) = 0.0  ! basaltic glass

    ! ! speed test
    ! primary(:,:,:) = 0.0
    ! primary_a(:,:,:) = 0.0
    ! primary_b(:,:,:) = 0.0

    ! secondary minerals [mol]
    secondary(:,:,:) = 0.0
    secondary_a(:,:,:) = 0.0
    secondary_b(:,:,:) = 0.0

    sec_density = (/2.65, 2.3, 3.05, 2.17, 5.01, 2.5, 3.8, & ! 7
         & 2.7, 2.71, 2.56, 2.3, 2.28, 2.28, 3.05, 2.28, & ! 8
         & 2.25, 5.3, 2.5, 2.55, 2.27, 2.2, 2.5, 2.26, & ! 8
         & 2.55, 2.25, 2.75, 2.7, 2.87, 2.9, 2.275, 2.8, & !8
         & 2.8, 2.3, 2.55, 4.61, 2.3, 2.3/) ! 6

    pri_density = (/1.0, 2.7, 3.0, 3.0, 3.0/)


    sec_molar = (/258.156, 480.19, 429.02, 2742.13, 119.98, 549.07, 88.851, &
         & 549.07, 100.0869, 287.327, 480.19, 495.90, 495.90, 429.02, 495.90, &
         & 380.22, 159.6882, 549.07, 504.19, 220.15, 649.86, 549.07, 649.86, &
         & 504.19, 380.22, 379.259, 549.07, 395.38, 64.448, 392.34, 64.448, &
         & 64.448, 480.19, 504.19, 85.12, 480.19, 480.19/)

    ! pri_molar = (/1.0, 110.0, 153.0, 158.81, 277.0/)
    pri_molar = (/1.0, 277.0, 153.0, 158.81, 110.0/)
    ! saturation
    saturation(:,:,:) = 0.0

    ! ! SURFACE SEAWATER SITE 858-ish, JUAN DE FUCA AREA
    ! from elderfield 1999, and other places
    solute(:,:,1) = 8.2      ! ph
    solute(:,:,2) = .00243   ! Alk 1.6e-3
    solute(:,:,3) = t_vol_s     ! water mass
    solute(:,:,4) = .002100  ! param_dic , TOTAL C
    solute(:,:,5) = .01028   ! Ca
    ! solute(1:xn/(4*cellx),:,5) = .01428
    solute(:,:,6) = .0528    ! Mg
    solute(:,:,7) = .460     ! Na
    solute(:,:,8) = .00995   ! K
    solute(:,:,9) = 0.0      ! Fe
    solute(:,:,10) = .028    ! S(6)
    solute(:,:,11) = 0.0     ! Si
    solute(:,:,12) = .540    ! Cl
    solute(:,:,13) =  0.0 ! Al
    solute(:,:,14) = .00245  ! inert
    solute(:,:,15) = 0.0     ! CO3-2

    solute_a(:,:,1) = 8.2      ! ph
    solute_a(:,:,2) = .00243   ! Alk 1.6e-3
    solute_a(:,:,3) = t_vol_a   ! water mass
    solute_a(:,:,4) = .002100  ! param_dic , TOTAL C
    solute_a(:,:,5) = .01028   ! Ca
    ! solute_a(1:xn/(4*cellx),:,5) = .01428
    solute_a(:,:,6) = .0528    ! Mg
    solute_a(:,:,7) = .460     ! Na
    solute_a(:,:,8) = .00995   ! K
    solute_a(:,:,9) = 0.0      ! Fe
    solute_a(:,:,10) = .028    ! S(6)
    solute_a(:,:,11) = 0.0     ! Si
    solute_a(:,:,12) = .540    ! Cl
    solute_a(:,:,13) =  0.0 ! Al
    solute_a(:,:,14) = .00245  ! inert
    solute_a(:,:,15) = 0.0     ! CO3-2

    solute_b(:,:,1) = 8.2      ! ph
    solute_b(:,:,2) = .00243   ! Alk 1.6e-3
    solute_b(:,:,3) = t_vol_b  ! water mass
    solute_b(:,:,4) = .002100  ! param_dic , TOTAL C
    solute_b(:,:,5) = .01028   ! Ca
    ! solute_b(1:xn/(4*cellx),:,5) = .01428
    solute_b(:,:,6) = .0528    ! Mg
    solute_b(:,:,7) = .460     ! Na
    solute_b(:,:,8) = .00995   ! K
    solute_b(:,:,9) = 0.0      ! Fe
    solute_b(:,:,10) = .028    ! S(6)
    solute_b(:,:,11) = 0.0    ! Si
    solute_b(:,:,12) = .540    ! Cl
    solute_b(:,:,13) =  0.0 ! Al
    solute_b(:,:,14) = .00245  ! inert
    solute_b(:,:,15) = 0.0     ! CO3-2

    ! seawater solute concentrations [mol/kgw]
    soluteOcean = (/ solute(1,1,1), solute(1,1,2), solute(1,1,3), solute(1,1,4), solute(1,1,5), &
         & solute(1,1,6), solute(1,1,7), solute(1,1,8), solute(1,1,9), solute(1,1,10), &
         & solute(1,1,11), solute(1,1,12), solute(1,1,13), solute(1,1,14), solute(1,1,15) /)

    soluteOcean_a = (/ solute_a(1,1,1), solute_a(1,1,2), solute_a(1,1,3), solute_a(1,1,4), solute_a(1,1,5), &
         & solute_a(1,1,6), solute_a(1,1,7), solute_a(1,1,8), solute_a(1,1,9), solute_a(1,1,10), &
         & solute_a(1,1,11), solute_a(1,1,12), solute_a(1,1,13), solute_a(1,1,14), solute_a(1,1,15) /)

    soluteOcean_b = (/ solute_b(1,1,1), solute_b(1,1,2), solute_b(1,1,3), solute_b(1,1,4), solute_b(1,1,5), &
         & solute_b(1,1,6), solute_b(1,1,7), solute_b(1,1,8), solute_b(1,1,9), solute_b(1,1,10), &
         & solute_b(1,1,11), solute_b(1,1,12), solute_b(1,1,13), solute_b(1,1,14), solute_b(1,1,15) /)

    DO g=1,g_sol
       solute_fine(:,:,g) = soluteOcean(g)
       solute_fine_a(:,:,g) = soluteOcean_a(g)
       solute_fine_b(:,:,g) = soluteOcean_b(g)
    END DO

    sol_index = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/)
    !            alk dic ca mg na k  fe s   si  cl  al

    volume_ratio = t_vol_a/t_vol_b
    mix_ratio = dt*mstep/param_t_diff

    !soluteOcean(5) = 0.011
    !solute_fine(1:xn/4,:,5) = 0.011
    !solute(1:xn/(4*cellx),:,5) = 0.011

    ! write(*,*) "fine" , solute_fine(:,:,5)
    ! write(*,*) " "
    ! write(*,*) "coarse" , solute(:,:,5)

    medium(:,:,1) = .1          ! phiCoarse
    medium(:,:,2) = 0.0         ! precip (saturation state)
    medium(:,:,3) = t_vol_s       ! water_volume
    vol_i = medium(1,1,3)
    medium(:,:,4) = 0.01        ! reactive fraction now!
    medium(:,:,5) = 1.0         ! rxn toggle
    medium(:,:,6) = 0.0         ! x-coord
    medium(:,:,7) = 0.0         ! y-coord


    medium_a(:,:,1) = .1          ! phiCoarse
    medium_a(:,:,2) = 0.0         ! precip (saturation state)
    medium_a(:,:,3) = t_vol_a       ! water_volume
    vol_i_a = medium_a(1,1,3)
    medium_a(:,:,4) = 0.01        ! reactive fraction now!
    medium_a(:,:,5) = 1.0         ! rxn toggle
    medium_a(:,:,6) = 0.0         ! x-coord
    medium_a(:,:,7) = 0.0         ! y-coord

    medium_b(:,:,1) = .1          ! phiCoarse
    medium_b(:,:,2) = 0.0         ! precip (saturation state)
    medium_b(:,:,3) = t_vol_b     ! water_volume
    vol_i_b = medium_b(1,1,3)
    medium_b(:,:,4) = 0.01        ! reactive fraction now!
    medium_b(:,:,5) = 1.0         ! rxn toggle
    medium_b(:,:,6) = 0.0         ! x-coord
    medium_b(:,:,7) = 0.0         ! y-coord

    sea = soluteOcean


    !-PERMEABILITY SET UP

    slope = param_w
    slope2 = x_max-param_w_rhs
    buffer = 1000.0
    edge = param_w
    edge2 = x_max-param_w_rhs

    sed = -150.0
    sed1 = -300.0



    sed = -50.0



    sed1 = -200.0

    !sed1 = sum(sed1)/xn

    sed1((param_w/dx)+1:) = sed1(:xn-(param_w/dx))
    sed1(1:param_w/dx) = sed1((param_w/dx)+1)

    sed((param_w/dx)+1:) = sed(:xn-(param_w/dx))
    sed(1:param_w/dx) = sed((param_w/dx)+1)

    sed1(f_index1-2:f_index1+2) = SUM(sed1(f_index1-3:f_index1+2))/6.0
    sed(f_index1-2:f_index1+2) = SUM(sed(f_index1-3:f_index1+2))/6.0

    sed2 = sed-sed1!-sed1! - sed1


    sed2((param_w/dx)+1:) = sed2(:xn-(param_w/dx))
    sed2(1:param_w/dx) = sed2((param_w/dx)+1)

    !sed2 = (-sum(sed2-sed1)/xn)*(sed2-sed1)

    !sed = sum(sed)/xn
    !sed = -100.0

    IF (param_o_rhs .GT. param_o) THEN
       sed = sed-(param_o_rhs)
    END IF

    IF (param_o .GE. param_o_rhs) THEN
       sed = sed-(param_o)
    END IF


    IF (param_o_rhs .GT. param_o) THEN
       sed1 = sed1-(param_o_rhs)
    END IF

    IF (param_o .GE. param_o_rhs) THEN
       sed1 = sed1-(param_o)
    END IF

    DO gg=1,yn
       DO g=1,xn
          IF ((y(gg) .GE. sed1(g)) .AND. (x(g) .GT. edge) .AND. (x(g) .LT. edge2)) THEN
             lambdaMat(g,gg) = 1.2
          END IF
       END DO
    END DO


    ! with sediment cap
    sed3 = sed1 - (param_h)


    ! the mask
    mask = 1.0
    DO gg=1,yn
       DO g =2,xn
          IF ((x(g) .GE. edge) .AND. (x(g) .LE. edge2+5000.0) .AND. (y(gg) .GE. sed(g))) THEN
             mask(g,gg) = 0.0
          END IF
       END DO
    END DO

    DO g =1,xn
       ! left outcrop top
       IF ((x(g) .LT. edge)) THEN
          IF (param_o_rhs .GT. param_o) THEN
             mask(g,yn-(ABS(param_o-param_o_rhs)/dy)-1) = 25.0
             mask(g,yn-(ABS(param_o-param_o_rhs)/dy):yn) = 0.0
          END IF
          IF (param_o_rhs .LT. param_o) THEN
             mask(g,yn-1) = 25.0
             mask(g,yn:yn) = 0.0
          END IF
          IF (param_o_rhs .EQ. param_o) THEN
             mask(g,yn-1) = 25.0
             mask(g,yn-0:yn) = 0.0
          END IF
       END IF
    END DO

    ! inner vertical edges
    DO gg=1,yn
       DO g =2,xn-1
          IF ((mask(g,gg) - mask(g-1,gg) .EQ. -1.0) .AND. (x(g) .LE. edge+dx)) THEN
             mask(g-1,gg) = 5.0
          END IF
          IF ((mask(g,gg) - mask(g-1,gg) .EQ. 1.0) .AND. (x(g) .GE. edge2-dx)) THEN
             mask(g,gg) = 10.0
          END IF
       END DO
    END DO

    DO gg=2,yn-1
       DO g =1,xn
          IF ((mask(g,gg) .EQ. 1.0) .AND. (mask(g,gg+1) .EQ. 0.0)) THEN
             mask(g,gg) = 50.0
          END IF
       END DO
    END DO

    DO gg=2,yn-1
       DO g =2,xn-1
          ! left upper corner
          IF ((mask(g,gg) .EQ. 25.0) .AND. (mask(g+1,gg-1) .EQ. 5.0)) THEN
             mask(g+1,gg) = 12.5
          END IF
          ! right upper corner
          IF ((mask(g,gg) .EQ. 25.0) .AND. (mask(g-1,gg-1) .EQ. 10.0)) THEN
             mask(g-1,gg) = 17.5
          END IF
       END DO
    END DO

    DO gg=2,yn-1
       DO g =2,xn-1
          ! left bottom corner
          IF ((mask(g,gg) .EQ. 5.0) .AND. (mask(g+1,gg-1) .EQ. 50.0)) THEN
             mask(g,gg-1) = 2.5
          END IF

          ! right bottom corner
          IF ((mask(g,gg) .EQ. 10.0) .AND. (mask(g-1,gg-1) .EQ. 50.0)) THEN
             mask(g,gg-1) = 7.5
          END IF
       END DO
    END DO



    maskP = mask

    DO gg=1,yn-1
       DO g =1,xn
          IF (y(gg) .LE. sed3(g)) THEN
             maskP(g,gg) = 0.0
          END IF
       END DO
    END DO

    DO gg=2,yn-3
       DO g =1,xn
          IF (gg .EQ. yn/2 + 3) THEN
             maskP(:,gg) = 100.0
             code = gg
          END IF
       END DO
    END DO

    DO gg=2,yn-3
       DO g =1,xn
          IF ((y(gg) .LE. sed3(g)) .AND. (gg .GT. code)) THEN
             maskP(:,gg) = 1.0
          END IF
       END DO
    END DO



    ! coarse mask
    coarse_mask = 0.0
    DO gg=1,yn/(2*celly)
       DO g =1,(xn-1)/cellx
          IF (MAXVAL(maskP((g-1)*cellx+1:g*cellx,   (gg-1)*celly+1+yn/(2*celly):(gg*celly)+yn/(2*celly)  )) .GT. 0.0) THEN
             coarse_mask(g,gg) = 1.0
          END IF
          IF (coarse_mask(g,gg) .EQ. 0.0) THEN
             medium(g,gg,5) = 0.0    ! cell toggle
             primary(g,gg,:) = 0.0   ! basaltic glass
             ! solute(g,gg,3) = 0.0    ! solute water
             ! medium(g,gg,3) = 0.0    ! medium water

             medium_a(g,gg,5) = 0.0    ! cell toggle
             primary_a(g,gg,:) = 0.0   ! basaltic glass
             ! solute_a(g,gg,3) = 0.0    ! solute water
             ! medium_a(g,gg,3) = 0.0    ! medium water

             medium_b(g,gg,5) = 0.0    ! cell toggle
             primary_b(g,gg,:) = 0.0   ! basaltic glass
             ! solute_b(g,gg,3) = 0.0    ! solute water
             ! medium_b(g,gg,3) = 0.0    ! medium water
          END IF

          ! 			if (x((g-1)*cellx) .ge. x(f_index1-1)) then
          ! 				medium(g,gg,5) = 0.0    ! cell toggle
          ! 				primary(g,gg,:) = 0.0   ! basaltic glass
          ! 				solute(g,gg,3) = 0.0    ! solute water
          ! 				medium(g,gg,3) = 0.0    ! medium water
          !
          ! 				medium_a(g,gg,5) = 0.0    ! cell toggle
          ! 				primary_a(g,gg,:) = 0.0   ! basaltic glass
          ! 				solute_a(g,gg,3) = 0.0    ! solute water
          ! 				medium_a(g,gg,3) = 0.0    ! medium water
          !
          ! 				medium_b(g,gg,5) = 0.0    ! cell toggle
          ! 				primary_b(g,gg,:) = 0.0   ! basaltic glass
          ! 				solute_b(g,gg,3) = 0.0    ! solute water
          ! 				medium_b(g,gg,3) = 0.0    ! medium water
          ! 			end if

       END DO
    END DO

    !-QUICK FIX
    ! medium(:,14:,5) = 0.0
    ! medium_a(:,14:,5) = 0.0
    ! medium_b(:,14:,5) = 0.0
    ! primary(:,14:,:) = 0.0
    ! primary_a(:,14:,:) = 0.0
    ! primary_b(:,14:,:) = 0.0

    !-QUICK FIX 2 CELLS HIGH ONLY
    ! medium(:,8:,5) = 0.0
    ! medium_a(:,8:,5) = 0.0
    ! medium_b(:,8:,5) = 0.0
    ! primary(:,8:,:) = 0.0
    ! primary_a(:,8:,:) = 0.0
    ! primary_b(:,8:,:) = 0.0

    ! 3 cells high
    medium(:,9:,5) = 0.0
    medium_a(:,9:,5) = 0.0
    medium_b(:,9:,5) = 0.0
    primary(:,9:,:) = 0.0
    primary_a(:,9:,:) = 0.0
    primary_b(:,9:,:) = 0.0


    ! quick fix full column, early august
    ! medium(:2,8:,5) = 0.0
    ! medium_a(:2,8:,5) = 0.0
    ! medium_b(:2,8:,5) = 0.0
    ! primary(:2,8:,:) = 0.0
    ! primary_a(:2,8:,:) = 0.0
    ! primary_b(:2,8:,:) = 0.0

    ! !-QUICK FIX to fix 3950, -150
    ! medium(:2,14:,5) = 0.0
    ! medium_a(:2,14:,5) = 0.0
    ! medium_b(:2,14:,5) = 0.0
    ! primary(:2,14:,:) = 0.0
    ! primary_a(:2,14:,:) = 0.0
    ! primary_b(:2,14:,:) = 0.0



    ! 	do gg=1,yn/celly
    ! 		do g =1,(xn-1)/cellx-1
    ! 			if ((coarse_mask(g,gg) .eq. 1.0) .and. (coarse_mask(g+1,gg) .eq. 0.0)) then
    ! 		end do
    ! 	end do






    fives = 1.0
    tens = 1.0
    twentyfives = 1.0
    fifties = 1.0

    DO gg=1,yn
       DO g=1,xn
          IF ((maskP(g,gg) .EQ. 5.0) .OR. (maskP(g,gg) .EQ. 12.5)) THEN
             fives(1,gg) = g
             fives(2,gg) = gg
          END IF
          IF ((maskP(g,gg) .EQ. 10.0) .OR. (maskP(g,gg) .EQ. 17.5)) THEN
             tens(1,gg) = g
             tens(2,gg) = gg
          END IF
       END DO
    END DO


    DO g=1,xn
       DO gg=1,yn
          IF ((maskP(g,gg) .EQ. 25.0) .OR. (maskP(g,gg) .EQ. 12.5) .OR. (maskP(g,gg) .EQ. 17.5)) THEN
             twentyfives(1,g) = g
             twentyfives(2,g) = gg
          END IF
          IF (maskP(g,gg) .EQ. 50.0) THEN
             fifties(1,g) = g
             fifties(2,g) = gg
          END IF
       END DO
    END DO

    ! TOO SIMPLE
    permeability = param_paq
    DO gg=1,yn-2
       DO g=1,xn
          IF ((ANY(maskP(g,:) .EQ. 50.0)) .AND. (y(gg) .GE. sed1(g)) .AND. (y(gg) .LE. sed(g)) .AND. (x(g) .LE. x_max-param_w_rhs+5000.0)) THEN
             permeability(g:g,gg:gg) = 1e-17
          END IF
          IF ((y(gg) .LE. sed3(g))) THEN
             permeability(g,gg) = 1e-16
          END IF
       END DO
    END DO


    ! 369 fracture goes HERE
    DO gg=2,yn-3
       DO g =1,xn
          IF ((maskP(g,gg) .EQ. 1.0) .AND. (g .EQ. f_index1) .AND. (y(gg) .GE. sed1(g))) THEN
             maskP(g,gg) = 6.0
             maskP(g-1,gg) = 3.0
             mask(g,gg) = 6.0
             mask(g-1,gg) = 3.0
          END IF
          !if ((maskP(g,gg) .eq. 100.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
          IF ((maskP(g,gg-1) .EQ. 1.0) .AND. (maskP(g,gg) .EQ. 6.0)) THEN
             maskP(g,code) = 6.01
             maskP(g,code+1:gg-1) = 6.05
             maskP(g,gg-1) = 6.1
             maskP(g-1,code) = 3.01
             maskP(g-1,code+1:gg-1) = 3.05
             maskP(g-1,gg-1) = 3.1
             mask(g,2:gg-1) = 6.05
             mask(g,gg-1) = 6.1
             mask(g-1,2:gg-1) = 3.05
             mask(g-1,gg-1) = 3.1
          END IF
          IF ((maskP(g,gg) .EQ. 50.0) .AND. (g .EQ. f_index1) .AND. (y(gg) .GE. sed1(g))) THEN
             maskP(g,gg) = 6.5
             maskP(g-1,gg) = 3.5
             mask(g,gg) = 6.5
             mask(g-1,gg) = 3.5
          END IF
       END DO
    END DO


    !
    ! coarse mask 2
    DO gg=1,yn/(2*celly)
       DO g =1,(xn-1)/cellx
          IF (MAXVAL(permeability((g-1)*cellx+1:g*cellx,(gg-1)*celly+1+yn/(2*celly):(gg*celly)+yn/(2*celly))) .LT. param_paq) THEN
             coarse_mask(g,gg) = 0.0
          END IF
          IF (coarse_mask(g,gg) .EQ. 0.0) THEN
             medium(g,gg,5) = 0.0    ! cell toggle
             primary(g,gg,:) = 0.0   ! basaltic glass
             ! solute(g,gg,3) = 0.0    ! solute water
             ! medium(g,gg,3) = 0.0    ! medium water

             medium_a(g,gg,5) = 0.0    ! cell toggle
             primary_a(g,gg,:) = 0.0   ! basaltic glass
             ! solute_a(g,gg,3) = 0.0    ! solute water
             ! medium_a(g,gg,3) = 0.0    ! medium water

             medium_b(g,gg,5) = 0.0    ! cell toggle
             primary_b(g,gg,:) = 0.0   ! basaltic glass
             ! solute_b(g,gg,3) = 0.0    ! solute water
             ! medium_b(g,gg,3) = 0.0    ! medium water
          END IF
       END DO
    END DO

    !- COARSE MASK QUICK FIX
    coarse_mask(:,9:) = 0.0
    coarse_mask_long = RESHAPE(coarse_mask,(/((xn-1)/cellx)*(yn/(2*celly))/))

    ! high lambda in deep basalt
    DO gg=1,yn
       DO g=1,xn
          IF ((y(gg) .LT. sed1(g)) .AND. (permeability(g,gg) .EQ. 1e-18)) THEN
             lambdaMat(g,gg) = 1.8
          END IF
       END DO
    END DO

    active_cells = 0
    DO gg=1,yn
       DO g=1,xn
          IF (maskP(g,gg) .NE. 0.0) THEN
             active_cells = active_cells + 1
          END IF
       END DO
    END DO

    long = (xn-2)*(yn-2)
    longP = (xn-2)*((yn/2)-2)

    maskLong = RESHAPE(mask(2:xn-1,2:yn-1), (/long/))
    maskLongT = RESHAPE(TRANSPOSE(mask(2:xn-1,2:yn-1)), (/long/))

    maskLongU = RESHAPE(mask(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
    maskLongTV = RESHAPE(TRANSPOSE(mask(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))

    maskPLong = RESHAPE(maskP(2:xn-1,(yn/2)+2:yn-1), (/longP/))
    maskPLongT = RESHAPE(TRANSPOSE(maskP(2:xn-1,(yn/2)+2:yn-1)), (/longP/))



    RETURN

  END SUBROUTINE init






  FUNCTION h_bc(h_in)

    USE globals
    REAL(4) :: h_in(xn,yn), h_bc(xn,yn), rip_lith_y(xn)
    INTEGER :: p, pp

    rip_lith_y = 0.48
    ! 	rip_lith_y = (/0.602087135445, &
    ! & 0.60116684559, 0.600246555734, 0.599326265879, 0.598405976024, 0.597485686169, 0.596565396314, &
    ! & 0.595645106459, 0.594724816604, 0.593804526749, 0.592884236894, 0.591963947039, 0.591043657184, &
    ! & 0.590223660487, 0.589514072973, 0.588804485459, 0.588094897945, 0.587385310431, 0.586675722917, &
    ! & 0.585966135403, 0.585256547889, 0.584546960375, 0.583837372862, 0.583127785348, 0.582418197834, &
    ! & 0.58170861032, 0.580999022806, 0.580289435292, 0.579579847778, 0.578870260264, 0.57816067275, &
    ! & 0.577451085236, 0.576741497722, 0.576031910208, 0.575322322694, 0.57461273518, 0.573903147666, &
    ! & 0.573193560152, 0.572483972639, 0.571774385125, 0.571064797611, 0.570355210097, 0.569645622583, &
    ! & 0.568936035069, 0.568226447555, 0.567516860041, 0.566810372198, 0.566178703989, 0.565547035781, &
    ! & 0.564915367572, 0.564283699364, 0.563652031155, 0.563020362947, 0.562388694739, 0.56175702653, &
    ! & 0.561125358322, 0.560493690113, 0.559862021905, 0.559230353696, 0.558598685488, 0.557967017279, &
    ! & 0.557335349071, 0.556703680862, 0.556072012654, 0.555440344446, 0.554808676237, 0.554177008029, &
    ! & 0.55354533982, 0.552913671612, 0.552282003403, 0.551650335195, 0.551018666986, 0.550386998778, &
    ! & 0.54975533057, 0.549123662361, 0.548491994153, 0.547860325944, 0.547228657736, 0.546596989527, &
    ! & 0.545965321319, 0.54533365311, 0.544701984902, 0.544070316694, 0.543438648485, 0.542806980277, &
    ! & 0.542175312068, 0.54154364386, 0.540911975651, 0.540280307443, 0.539648639234, 0.539016971026, &
    ! & 0.538385302817, 0.537753634609, 0.537121966401, 0.536490298192, 0.535858629984, 0.535226961775, &
    ! & 0.534595293567, 0.533963625358, 0.53333195715, 0.532700288941, 0.532068620733, 0.531436952525, &
    ! & 0.530805284316, 0.530173616108, 0.529541947899, 0.528910279691, 0.528278611482, 0.527646943274, &
    ! & 0.527015275065, 0.526383606857, 0.525751938648, 0.52512027044, 0.524488602232, 0.523856934023, &
    ! & 0.523225265815, 0.522593597606, 0.521961929398, 0.521330261189, 0.520698592981, 0.520066924772, &
    ! & 0.519435256564, 0.518803588356, 0.518171920147, 0.517540251939, 0.51690858373, 0.516276915522, &
    ! & 0.515681666368, 0.515130941731, 0.514580217093, 0.514029492455, 0.513478767818, 0.51292804318, &
    ! & 0.512377318542, 0.511826593905, 0.511275869267, 0.510725144629, 0.510174419992, 0.509623695354, &
    ! & 0.509072970716, 0.508522246079, 0.507971521441, 0.507420796803, 0.506870072166, 0.506319347528, &
    ! & 0.50576862289, 0.505217898253, 0.504667173615, 0.504116448977, 0.503565724339, 0.503014999702, &
    ! & 0.502464275064, 0.501913550426, 0.501362825789, 0.500812101151, 0.500261376513, 0.499710651876, &
    ! & 0.499159927238, 0.498679752216, 0.498249474112, 0.497819196008, 0.497388917903, 0.496958639799, &
    ! & 0.496528361695, 0.496098083591, 0.495667805487, 0.495237527382, 0.494807249278, 0.494376971174, &
    ! & 0.49394669307, 0.493516414966, 0.493086136861, 0.492655858757, 0.492225580653, 0.491795302549, &
    ! & 0.491365024445, 0.490934746341, 0.490504468236, 0.490074190132, 0.489643912028, 0.489213633924, &
    ! & 0.48878335582, 0.488353077715, 0.487922799611, 0.487492521507, 0.487062243403, 0.486631965299, &
    ! & 0.486201687194, 0.48577140909, 0.485341130986, 0.484873859367, 0.48440077448, 0.483927689594, &
    ! & 0.483454604708, 0.482981519822, 0.482508434936, 0.48203535005, 0.481562265164, 0.481089180278, &
    ! & 0.480616095391, 0.480143010505, 0.479669925619, 0.479196840733, 0.478723755847, 0.478250670961, &
    ! & 0.477777586075, 0.477304501189, 0.476831416302, 0.476358331416, 0.47588524653, 0.475412161644, &
    ! & 0.474939076758, 0.474465991872, 0.473992906986, 0.4735198221, 0.473046737213, 0.472573652327, &
    ! & 0.472100567441, 0.471627482555, 0.471154397669, 0.470681312783, 0.470208227897, 0.469735143011, &
    ! & 0.469262058124, 0.468788973238, 0.468315888352, 0.467944418472, 0.46759004727, 0.467235676069, &
    ! & 0.466881304867, 0.466526933666, 0.466172562464, 0.465818191263, 0.465463820061, 0.46510944886, &
    ! & 0.464755077658, 0.464400706457, 0.464046335255, 0.463691964054, 0.463337592852, 0.462983221651, &
    ! & 0.462628850449, 0.462274479248, 0.461920108046, 0.461565736845, 0.461211365643, 0.460856994442, &
    ! & 0.46050262324, 0.460148252039, 0.459793880837, 0.459439509636, 0.459085138434, 0.458730767233, &
    ! & 0.458376396031, 0.45802202483, 0.457667653628, 0.457313282427, 0.456958911225, 0.456604540024, &
    ! & 0.456250168822, 0.455895797621, 0.455541426419, 0.455187055218, 0.454832684016, 0.454478312815, &
    ! & 0.454123941613, 0.453769570412, 0.45341519921, 0.453060828009, 0.452706456807, 0.452352085606, &
    ! & 0.451997714404, 0.451643343203, 0.451288972001, 0.4509346008, 0.450580229598, 0.450225858397, &
    ! & 0.449871487195, 0.449517115994, 0.449162744792, 0.448808373591, 0.448454002389, 0.448099631188, &
    ! & 0.447745259986, 0.447390888785, 0.447036517584, 0.446682146382, 0.446327775181, 0.445973403979, &
    ! & 0.445619032778, 0.445264661576, 0.444910290375, 0.444555919173, 0.444201547972, 0.44384717677, &
    ! & 0.443492805569, 0.443138434367, 0.442784063166, 0.442429691964, 0.442075320763, 0.441720949561, &
    ! & 0.44136657836, 0.441012207158, 0.440657835957, 0.440303464755, 0.439949093554, 0.439594722352, &
    ! & 0.439240351151, 0.438885979949, 0.438531608748, 0.438177237546, 0.437822866345, 0.437468495143, &
    ! & 0.437114123942, 0.43675975274, 0.436405381539, 0.436051010337, 0.435696639136, 0.435342267934, &
    ! & 0.434987896733, 0.434633525531, 0.43427915433, 0.433924783128, 0.433570411927, 0.433216040725, &
    ! & 0.432861669524, 0.432507298322, 0.432152927121, 0.431798555919, 0.431444184718, 0.431089813516, &
    ! & 0.430735442315, 0.430422684669, 0.430144535282, 0.429866385895, 0.429588236509, 0.429310087122, &
    ! & 0.429031937735, 0.428753788348, 0.428475638961, 0.428197489574, 0.427919340187, 0.427641190801, &
    ! & 0.427363041414, 0.427084892027, 0.42680674264, 0.426528593253, 0.426250443866, 0.42597229448, &
    ! & 0.425694145093, 0.425415995706, 0.425137846319, 0.424859696932, 0.424581547545, 0.424303398158, &
    ! & 0.424025248772, 0.423747099385, 0.423468949998, 0.423190800611, 0.422912651224, 0.422634501837, &
    ! & 0.422356352451, 0.422078203064, 0.421800053677, 0.42152190429, 0.421243754903, 0.420965605516, &
    ! & 0.420687456129, 0.420409306743, 0.420131157356, 0.419853007969, 0.419574858582, 0.419296709195, &
    ! & 0.419018559808, 0.418740410422, 0.418462261035, 0.418184111648, 0.417905962261, 0.417627812874, &
    ! & 0.417349663487, 0.4170715141, 0.416793364714, 0.416515215327, 0.41623706594, 0.415958916553, &
    ! & 0.415680767166, 0.415402617779, 0.415124468393, 0.414878071636, 0.414638695561, 0.414399319487, &
    ! & 0.414159943412, 0.413920567337, 0.413681191262, 0.413441815188, 0.413202439113, 0.412963063038, &
    ! & 0.412723686964, 0.412484310889, 0.412244934814, 0.41200555874, 0.411766182665, 0.41152680659, &
    ! & 0.411287430516, 0.411048054441, 0.410808678366, 0.410569302292, 0.410329926217, 0.410090550142, &
    ! & 0.409851174068, 0.409611797993, 0.409372421918, 0.409133045844, 0.408893669769, 0.408654293694, &
    ! & 0.40841491762, 0.408175541545, 0.40793616547, 0.407696789396, 0.407457413321, 0.407218037246, &
    ! & 0.406978661172, 0.406739285097, 0.406499909022, 0.406260532948, 0.406021156873, 0.405781780798, &
    ! & 0.405542404724, 0.405303028649, 0.405063652574, 0.4048242765, 0.404584900425, 0.40434552435, &
    ! & 0.404106148276, 0.403866772201, 0.403627396126, 0.403388020051, 0.403148643977, 0.402914577103, &
    ! & 0.402697703397, 0.402480829691, 0.402263955985, 0.402047082279, 0.401830208573, 0.401613334867, &
    ! & 0.401396461161, 0.401179587455, 0.400962713749, 0.400745840043, 0.400528966337, 0.400312092631, &
    ! & 0.400095218925, 0.399878345219, 0.399661471513, 0.399444597807, 0.399227724101, 0.399010850395, &
    ! & 0.398793976689, 0.398577102983, 0.398360229277, 0.398143355571, 0.397926481865, 0.397709608159, &
    ! & 0.397492734453, 0.397275860747, 0.397058987041, 0.396842113335, 0.396625239629, 0.396408365923, &
    ! & 0.396191492217, 0.395974618511, 0.395757744805, 0.395540871099, 0.395323997393, 0.395107123687, &
    ! & 0.394890249981, 0.394673376275, 0.394456502569, 0.394239628863, 0.394022755157, 0.393805881451, &
    ! & 0.393589007745, 0.393372134039, 0.393155260333, 0.392938386627, 0.392721512921, 0.3925131885, &
    ! & 0.392334317409, 0.392155446318, 0.391976575228, 0.391797704137, 0.391618833046, 0.391439961955, &
    ! & 0.391261090865, 0.391082219774, 0.390903348683, 0.390724477592, 0.390545606501, 0.390366735411, &
    ! & 0.39018786432, 0.390008993229, 0.389830122138, 0.389651251048, 0.389472379957, 0.389293508866, &
    ! & 0.389114637775, 0.388935766685, 0.388756895594, 0.388578024503, 0.388399153412, 0.388220282321, &
    ! & 0.388041411231, 0.38786254014, 0.387683669049, 0.387504797958, 0.387325926868, 0.387147055777, &
    ! & 0.386968184686, 0.386789313595, 0.386610442505, 0.386431571414, 0.386252700323, 0.386073829232, &
    ! & 0.385894958141, 0.385716087051, 0.38553721596, 0.385358344869, 0.385179473778, 0.385000602688, &
    ! & 0.384821731597, 0.384642860506, 0.384463989415, 0.384285118325, 0.384106247234, 0.383927376143, &
    ! & 0.383763195563, 0.383625514403, 0.383487833244, 0.383350152084, 0.383212470925, 0.383074789766, &
    ! & 0.382937108606, 0.382799427447, 0.382661746287, 0.382524065128, 0.382386383969, 0.382248702809, &
    ! & 0.38211102165, 0.38197334049, 0.381835659331, 0.381697978171, 0.381560297012, 0.381422615853, &
    ! & 0.381284934693, 0.381147253534, 0.381009572374, 0.380871891215, 0.380734210055, 0.380596528896, &
    ! & 0.380458847737, 0.380321166577, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
    ! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
    ! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
    ! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787/)

    h_bc = h_in


    ! top of outcrops
    DO pp=1,yn
       DO p=1,xn
          IF (mask(p,pp) .EQ. 0.0) THEN
             h_bc(p,pp) = param_tsw
          END IF
       END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    VERTICAL OUTER   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO pp=1,yn
       IF (mask(xn,pp) .NE. 0.0) THEN
          h_bc(xn,pp) = (4.0/3.0)*h_in(xn-1,pp) - (1.0/3.0)*h_in(xn-2,pp) ! right
       END IF
    END DO

    DO pp=1,yn
       IF (mask(1,pp) .NE. 0.0) THEN
          h_bc(1,pp) = (4.0/3.0)*h_in(2,pp) - (1.0/3.0)*h_in(3,pp) ! left
       END IF
    END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   HORIZONTAL OUTER  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! bottom
    DO p = 1,xn
       h_bc(p,1) = h_in(p,2) + ( rip_lith_y(p)) * dy/1.8
    END DO

    ! two lines recent...
    h_bc(xn,1) = h_bc(xn,2)
    h_bc(1,1) = h_bc(1,2)

    ! top of outcrops
    DO p=1,xn
       DO pp=2,yn-1
          IF (mask(p,pp) .EQ. 25.0) THEN
             h_bc(p,pp+1) = param_tsw
          END IF

          IF (mask(p,pp) .EQ. 12.5) THEN
             h_bc(p,pp+1) = param_tsw
          END IF
          IF (mask(p,pp) .EQ. 17.5) THEN
             h_bc(p,pp+1) = param_tsw
          END IF
       END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   HORIZONTAL INNER  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! top of sediment
    DO p=2,xn-1
       DO pp=2,yn-1
          IF ((mask(p,pp) .EQ. 50.0) .OR. (mask(p,pp) .EQ. 3.5) .OR. (mask(p,pp) .EQ. 6.5)) THEN
             h_bc(p,pp+1) = param_tsw
          END IF
       END DO
    END DO

    DO p=1,xn
       DO pp=2,yn-1
          IF ((mask(p,pp) .EQ. 12.5) .OR. (mask(p,pp) .EQ. 17.5)) THEN
             h_bc(p,pp+1) = param_tsw
          END IF
       END DO
    END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    VERTICAL INNER   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO p=2,xn-1
       DO pp=2,yn
          ! left outcrop
          IF ((mask(p,pp) .EQ. 12.5)) THEN
             h_bc(p+1,pp) = param_tsw
          END IF
          IF ((mask(p,pp) .EQ. 5.0)) THEN
             h_bc(p+1,pp) = param_tsw
          END IF
          ! right outcrop
          IF ((mask(p,pp) .EQ. 17.5)) THEN
             h_bc(p-1,pp) = param_tsw
          END IF
          IF ((mask(p,pp) .EQ. 10.0)) THEN
             h_bc(p-1,pp) = param_tsw
          END IF
       END DO
    END DO


    RETURN
  END FUNCTION h_bc










  FUNCTION psi_bc(psi_in)

    USE globals
    REAL(4) :: psi_in(xn,yn), psi_bc(xn,yn)
    INTEGER :: p,pp

    psi_bc = psi_in


    DO pp=1,yn
       DO p=1,xn
          IF ((maskP(p,pp) .EQ. 0.0)) THEN
             psi_bc(p,pp) = 0.0
          END IF
       END DO
    END DO

    DO p=1,xn
       DO pp=2,yn-1
          ! top of outcrops
          IF ((maskP(p,pp) .EQ. 25.0) .OR. (maskP(p,pp) .EQ. 12.5) .OR. (maskP(p,pp) .EQ. 17.5)) THEN
             psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
          END IF
       END DO
    END DO


    DO pp=1,yn
       ! left
       IF (maskP(1,pp) .NE. 0.0) THEN
          psi_bc(1,pp) = 0.0
       END IF
       ! right
       IF (maskP(xn,pp) .NE. 0.0) THEN
          psi_bc(xn,pp) =0.0
       END IF
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   HORIZONTAL OUTER  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO p=1,xn
       DO pp=2,yn-1
          IF ((maskP(p,pp) .EQ. 100.0) .OR. (maskP(p,pp) .EQ. 3.01) .OR. (maskP(p,pp) .EQ. 6.01)) THEN
             psi_bc(p,pp-1) = 0.0
          END IF
       END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   HORIZONTAL INNER  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO p=2,xn-1
       DO pp=2,yn-1
          ! right outcrop (left boundary)
          IF ((maskP(p,pp) .EQ. 10.0)) THEN
             psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
          END IF

          IF ((maskP(p,pp) .EQ. 17.5)) THEN
             psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
          END IF

          ! left outcrop (right boundary)
          IF ((maskP(p,pp) .EQ. 5.0)) THEN
             psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
          END IF

          IF ((maskP(p,pp) .EQ. 12.5)) THEN
             psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
          END IF
       END DO
    END DO

    DO p=2,xn
       DO pp=2,yn-1
          ! top of sediment
          IF ((maskP(p,pp) .EQ. 50.0) .OR. (maskP(p,pp) .EQ. 3.5) .OR. (maskP(p,pp) .EQ. 6.5)) THEN
             psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
          END IF
       END DO
    END DO

    psi_bc(f_index1,:) = 0.0


    RETURN
  END FUNCTION psi_bc









  FUNCTION psi_mod(psi_in)

    USE globals
    REAL(4) :: psi_in(xn,yn), psi_mod(xn,yn)
    INTEGER :: p,pp

    psi_mod = psi_in

    DO pp=1,yn
       DO p=1,xn
          IF (mask(p,pp) .EQ. 0.0) THEN
             psi_mod(p,pp) = dx*(5.0e-12)*1000.0 - dx*(5.0e-12)*0.0
          END IF
       END DO
    END DO




    RETURN
  END FUNCTION psi_mod

END MODULE initialize
