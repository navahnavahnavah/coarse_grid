module initialize

use globals
implicit none

save
 
integer :: g, gg, f, ff, long, longP, maskBit
real(4) :: x(xn), y(yn), t(tn)
real(4) :: rho0(xn,yn)
real(4) :: bcx0(xn,2), bcy0(2,yn), bcxPsi(xn,2), bcyPsi(2,yn), ic0(xn,yn)
real(4) :: kMat(xn,yn), lambdaMat(xn,yn), lambdaMatW(xn,yn), porosity(xn,yn), reactive(xn,yn), reactiveCoarse(xn/cellx,yn/celly)
real(4) :: reactiveCoarse1(xn/cellx,yn/celly), reactiveCoarse2(xn/cellx,yn/celly), sedx
real(4) :: sed(xn), sed1(xn), sed2(xn), sed3(xn)

! EXTRA BULLSHIT
real(4) :: permeable(xn,4), permeability(xn,yn)
real(4) :: mask(xn,yn), maskLong((xn-2)*(yn-2)), maskLongT((xn-2)*(yn-2)), h0Long((xn-2)*(yn-2)), h0LongT((xn-2)*(yn-2)), h0T(yn,xn)
real(4) :: maskPLongFull(xn*yn)
real(4) :: maskLongU((xn-2)*(yn-0)), maskLongTV((xn-0)*(yn-2))
real(4) :: maskP(xn,yn), maskPLong((xn-2)*((yn/2)-2)), maskPLongT((xn-2)*((yn/2)-2))
real(4) :: outerBand((xn-2)*((yn/2)-2),2*((yn/2)-2) + 1), bigBand((xn-2)*((yn/2)-2),4*((yn/2)-2) + 3)
real(4) :: stretch(xn,yn), stretchLong((xn-2)*(yn-2)), stretchT(xn,yn), stretchLongT((xn-2)*(yn-2))
real(4) :: u_inter(xn,yn), v_inter(xn,yn), u_step_coarse(xn/cellx,yn/celly), v_step_coarse(xn/cellx,yn/celly)

! 05/06 INPUT STUFF
real(4) :: primary(xn/cellx,yn/celly,g_pri), primaryMat(xn*tn/(cellx*(mstep*ar)),yn/celly,g_pri)
real(4) :: secondary(xn/cellx,yn/celly,g_sec), secondaryMat(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sec)
real(4) :: solute(xn/cellx,yn/celly,g_sol), soluteMat(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sol)
real(4) :: solute_inter(xn/cellx,yn/celly)
real(4) :: medium(xn/cellx,yn/celly,g_med), mediumMat(xn*tn/(cellx*(mstep*ar)),yn/celly,g_med)
real(4) :: saturation(xn/cellx,yn/celly,g_sec/2), saturationMat(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sec/2)
real(4) :: checkMat(xn*tn/(cellx*(mstep*ar)),yn/celly), checkMatF(xn*tn/(cellx*(mstep*ar)),yn/celly)
real(4) :: soluteOcean(g_sol), soluteOcean_a(g_sol), soluteOcean_b(g_sol), sea(g_sol)
real(4) :: inter, slope, inter2, slope2, boundary, buffer, edge, edge2
real(4) :: bit_thing(xn/cellx,yn/(2*celly)), bit_thing_t(yn/(2*celly),xn/cellx), bit_thing_t1(xn/cellx,yn/(2*celly))
real(4) :: solute_fine(xn,yn,g_sol), solute_fine_a(xn,yn,g_sol), solute_fine_b(xn,yn,g_sol)

! CHAMBER INPUT STUFF

real(4) :: primary_a(xn/cellx,yn/celly,g_pri), primaryMat_a(xn*tn/(cellx*(mstep*ar)),yn/celly,g_pri)
real(4) :: secondary_a(xn/cellx,yn/celly,g_sec), secondaryMat_a(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sec)
real(4) :: solute_a(xn/cellx,yn/celly,g_sol), soluteMat_a(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sol)
real(4) :: medium_a(xn/cellx,yn/celly,g_med), mediumMat_a(xn*tn/(cellx*(mstep*ar)),yn/celly,g_med)

real(4) :: primary_b(xn/cellx,yn/celly,g_pri), primaryMat_b(xn*tn/(cellx*(mstep*ar)),yn/celly,g_pri)
real(4) :: secondary_b(xn/cellx,yn/celly,g_sec), secondaryMat_b(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sec)
real(4) :: solute_b(xn/cellx,yn/celly,g_sol), soluteMat_b(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sol)
real(4) :: medium_b(xn/cellx,yn/celly,g_med), mediumMat_b(xn*tn/(cellx*(mstep*ar)),yn/celly,g_med)

real(4) :: primaryMat_d(xn*tn/(cellx*(mstep*ar)),yn/celly,g_pri)
real(4) :: secondaryMat_d(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sec)
real(4) :: soluteMat_d(xn*tn/(cellx*(mstep*ar)),yn/celly,g_sol)
real(4) :: mediumMat_d(xn*tn/(cellx*(mstep*ar)),yn/celly,g_med)

! coarse/mean arrays added for super duper coarse grid 03/25/17
real(4) :: h_coarse(xn/cellx,yn/celly)
real(4) :: pri_coarse(xn/cellx,yn/celly,g_pri)
real(4) :: sec_coarse(xn/cellx,yn/celly,g_sec)
real(4) :: sol_coarse(xn/cellx,yn/celly,g_sol)
real(4) :: med_coarse(xn/cellx,yn/celly,g_med)

real(4) :: pri_coarse_a(xn/cellx,yn/celly,g_pri)
real(4) :: sec_coarse_a(xn/cellx,yn/celly,g_sec)
real(4) :: sol_coarse_a(xn/cellx,yn/celly,g_sol)
real(4) :: med_coarse_a(xn/cellx,yn/celly,g_med)

real(4) :: pri_coarse_b(xn/cellx,yn/celly,g_pri)
real(4) :: sec_coarse_b(xn/cellx,yn/celly,g_sec)
real(4) :: sol_coarse_b(xn/cellx,yn/celly,g_sol)
real(4) :: med_coarse_b(xn/cellx,yn/celly,g_med)

real(4) :: coarse_mask(xn/cellx,yn/celly)

real(4) :: volume_ratio, mix_ratio, volume_ratio_i
real(4) :: vol_i, vol_i_a, vol_i_b

! coordinates for optimization
real(4) :: fives(2,yn), tens(2,yn), twentyfives(2,xn), fifties(2,xn)

! shell file parameters
character(len=300) :: path, path2, path_final, crashstring, restartstring, iso_path
character(len=300) :: param_o_string, param_w_string, param_w_rhs_string, param_h_string, param_o_rhs_string, param_tsw_string
character(len=300) :: param_dic_string, param_scope_string, param_trace_string, param_ch_string, param_f_dx_string, param_f_k_string
character(len=300) :: param_paq_string, param_ch_rhs_string, param_f_freq_string, param_f_por_string
integer :: in, crashstep, restart, param_trace
real(4):: param_o, param_w, param_w_rhs, param_h, param_o_rhs, param_tsw, param_dic, param_scope, param_ch
real(4) :: param_paq, param_ch_rhs, param_f_dx, param_f_k, param_f_freq, param_f_por


! TRANSPOSED
real(4) :: hTrans(yn,xn), psiTrans(yn,xn), permeabilityTrans(yn,xn), phiTrans(yn,xn)
real(4) :: primaryTrans(yn/celly,xn/cellx,g_pri), secondaryTrans(yn/celly,xn/cellx,g_sec), soluteTrans(yn/celly,xn/cellx,g_sol), mediumTrans(yn/celly,xn/cellx,g_med)
real(4) :: primaryTrans_a(yn/celly,xn/cellx,g_pri), secondaryTrans_a(yn/celly,xn/cellx,g_sec), soluteTrans_a(yn/celly,xn/cellx,g_sol), mediumTrans_a(yn/celly,xn/cellx,g_med)
real(4) :: primaryTrans_b(yn/celly,xn/cellx,g_pri), secondaryTrans_b(yn/celly,xn/cellx,g_sec), soluteTrans_b(yn/celly,xn/cellx,g_sol), mediumTrans_b(yn/celly,xn/cellx,g_med)
real(4) :: genTrans(yn,xn)


! D2 STUFF
integer :: xn2, yn2, long2, m2
real(4) :: perm2 = 1e-12


real(4) :: frac6(yn,2), frac6_last(yn,2), temp6(yn,2), temp6_last(yn,2), temp6_mid(yn,2)
integer :: f_index1 = xn-60, iter = 0, spinup = 20!50000


real(4) :: temp6_a(yn), temp6_b(yn), temp6_c(yn), temp6_rhs(yn)

! SOLUTE ADVECTION DISTRIBUTION
integer :: sol_index(11)





contains
	
subroutine init_mini ()
use globals

dx = ( x_max - x_min ) / real ( xn - 1, kind = 4 ) 
x = linspace ( xn, x_min,x_max )
dy = ( y_max - y_min ) / real ( yn - 1, kind = 4 ) 
y = linspace ( yn, y_min, y_max )
dt = ( t_max - t_min ) / real ( tn - 1, kind = 4 ) 
t = linspace ( tn, t_min, t_max)


sea = (/8.2, 0.00243, 0.0, 0.0021, 0.01028, 0.0528, 0.460, 0.00995, 0.0, 0.028, 0.0, 0.540, 1.0e-8, 0.00245, 0.0/)


sol_index = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/)


return

end subroutine init_mini

! ----------------------------------------------------------------------------------%%
!
! SUBROUTINE TO INITIALIZE, JUST CALL IT
!
! ----------------------------------------------------------------------------------%%
  
subroutine init ()
use globals
integer :: m,n
!integer :: i,ii


in = iargc()
call getarg(1,restartstring)
call getarg(2,path)
call getarg(3,path_final)
call getarg(4,crashstring)
call getarg(5,param_o_string)
call getarg(6,param_w_string)
call getarg(7,param_w_rhs_string)
call getarg(8,param_h_string)
call getarg(9,param_o_rhs_string)
call getarg(10,param_tsw_string)
call getarg(11,param_dic_string)
call getarg(12,param_scope_string)
call getarg(13,param_trace_string)
call getarg(14,param_ch_string)
call getarg(15,param_paq_string)
call getarg(16,param_ch_rhs_string)
call getarg(17,param_f_dx_string)
call getarg(18,param_f_k_string)
call getarg(19,param_f_freq_string)
call getarg(20,param_f_por_string)
call getarg(21,iso_path)

read (crashstring, *) crashstep
read (restartstring, *) restart
read (param_o_string, *) param_o
param_o = param_o
read (param_w_string, *) param_w
read (param_w_rhs_string, *) param_w_rhs
read (param_h_string, *) param_h
param_h = param_h
read (param_o_rhs_string, *) param_o_rhs
param_o_rhs = param_o_rhs
read (param_tsw_string, *) param_tsw
param_tsw = param_tsw/10.0
read (param_dic_string, *) param_dic
param_dic = param_dic/100000.0
read (param_scope_string, *) param_scope
param_scope = 1.0*param_scope*(1e-10)
scope = param_scope

scope = 1.0*scope
param_scope = 1.0*param_scope

read (param_trace_string, *) param_trace
read (param_ch_string, *) param_ch
param_ch = -1.0*param_ch

read (param_paq_string, *) param_paq
read (param_ch_rhs_string, *) param_ch_rhs

read (param_f_dx_string, *) param_f_dx

param_f_dx = 10.0**(param_f_dx)


read (param_f_k_string, *) param_f_k
read (param_f_freq_string, *) param_f_freq
param_f_freq = param_f_dx/param_f_freq

read (param_f_por_string, *) param_f_por



!permf = rho_fluid*grav*4.0*param_f_dx*param_f_dx/(12.0*viscosity)
permf = param_f_dx*param_f_dx/3.0

! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
dx = ( x_max - x_min ) / real ( xn - 1, kind = 4 ) 
x = linspace ( xn, x_min,x_max )
dy = ( y_max - y_min ) / real ( yn - 1, kind = 4 ) 
y = linspace ( yn, y_min, y_max )
dt = ( t_max - t_min ) / real ( tn - 1, kind = 4 ) 
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


!--------------GEOCHEMICAL INITIAL CONDITIONS

! primary minerals [mol]
primary(:,:,:) = 0.0
primary(:,:,1) = 0.0 !1.29600 ! feldspar
primary(:,:,2) = 0.1 !.69600 ! plag
primary(:,:,3) = 0.1 !.12600 ! pyr
primary(:,:,4) = 0.1 !.04000 ! ol
primary(:,:,5) = 1.0 !9.67700 ! basaltic glass


primary_a(:,:,:) = 0.0
primary_a(:,:,1) = 0.0   ! feldspar
primary_a(:,:,2) = 0.1   ! plag
primary_a(:,:,3) = 0.1   ! pyr
primary_a(:,:,4) = 0.1   ! ol
primary_a(:,:,5) = 1.0  ! basaltic glass

primary_b(:,:,:) = 0.0
primary_b(:,:,1) = 0.0  ! feldspar
primary_b(:,:,2) = 0.0 ! plag
primary_b(:,:,3) = 0.0 ! pyr
primary_b(:,:,4) = 0.0 ! ol
primary_b(:,:,5) = 0.0  ! basaltic glass

! ! speed test
! primary(:,:,:) = 0.0
! primary_a(:,:,:) = 0.0
! primary_b(:,:,:) = 0.0

! secondary minerals [mol]
secondary(:,:,:) = 0.0
secondary_a(:,:,:) = 0.0
secondary_b(:,:,:) = 0.0

! saturation
saturation(:,:,:) = 0.0

! ! SURFACE SEAWATER SITE 858-ish, JUAN DE FUCA AREA
! from elderfield 1999, and other places
solute(:,:,1) = 8.2      ! ph
solute(:,:,2) = .00243   ! Alk 1.6e-3
solute(:,:,3) = .266     ! water mass
solute(:,:,4) = .002100  ! param_dic , TOTAL C
solute(:,:,5) = .01028   ! Ca
solute(1:xn/(2*cellx),:,5) = .01428
solute(:,:,6) = .0528    ! Mg
solute(:,:,7) = .460     ! Na
solute(:,:,8) = .00995   ! K
solute(:,:,9) = 0.0      ! Fe
solute(:,:,10) = .028    ! S(6)
solute(:,:,11) = 0.0     ! Si
solute(:,:,12) = .540    ! Cl
solute(:,:,13) =  1.0e-8 ! Al
solute(:,:,14) = .00245  ! inert
solute(:,:,15) = 0.0     ! CO3-2

solute_a(:,:,1) = 8.2      ! ph
solute_a(:,:,2) = .00243   ! Alk 1.6e-3
solute_a(:,:,3) = .266     ! water mass
solute_a(:,:,4) = .002100  ! param_dic , TOTAL C
solute_a(:,:,5) = .01028   ! Ca
solute_a(1:xn/(2*cellx),:,5) = .01428
solute_a(:,:,6) = .0528    ! Mg
solute_a(:,:,7) = .460     ! Na
solute_a(:,:,8) = .00995   ! K
solute_a(:,:,9) = 0.0      ! Fe
solute_a(:,:,10) = .028    ! S(6)
solute_a(:,:,11) = 0.0     ! Si
solute_a(:,:,12) = .540    ! Cl
solute_a(:,:,13) =  1.0e-8 ! Al
solute_a(:,:,14) = .00245  ! inert
solute_a(:,:,15) = 0.0     ! CO3-2

solute_b(:,:,1) = 8.2      ! ph
solute_b(:,:,2) = .00243   ! Alk 1.6e-3
solute_b(:,:,3) = 0.0266   ! water mass
solute_b(:,:,4) = .002100  ! param_dic , TOTAL C
solute_b(:,:,5) = .01028   ! Ca
solute_b(1:xn/(2*cellx),:,5) = .01428
solute_b(:,:,6) = .0528    ! Mg
solute_b(:,:,7) = .460     ! Na
solute_b(:,:,8) = .00995   ! K
solute_b(:,:,9) = 0.0      ! Fe
solute_b(:,:,10) = .028    ! S(6)
solute_b(:,:,11) = 0.0     ! Si
solute_b(:,:,12) = .540    ! Cl
solute_b(:,:,13) =  1.0e-8 ! Al
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
			  
do g=1,g_sol
	solute_fine(:,:,g) = soluteOcean(g)
	solute_fine_a(:,:,g) = soluteOcean_a(g)
	solute_fine_b(:,:,g) = soluteOcean_b(g)
end do

sol_index = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/)


volume_ratio = 10.0
mix_ratio = 0.01

!soluteOcean(5) = 0.011
!solute_fine(1:xn/4,:,5) = 0.011
!solute(1:xn/(4*cellx),:,5) = 0.011

! write(*,*) "fine" , solute_fine(:,:,5)
! write(*,*) " "
! write(*,*) "coarse" , solute(:,:,5)

medium(:,:,1) = .1          ! phiCoarse
medium(:,:,2) = 0.0         ! precip
medium(:,:,3) = 0.266       ! water_volume
vol_i = medium(1,1,3)
medium(:,:,4) = 0.01        ! reactive fraction now!
medium(:,:,5) = 0.0         ! rxn toggle
! medium(:,:,5) = 0.0         ! rxn toggle
medium(xn,yn/cell-3:,5) = 1.0         ! rxn toggle
medium(:,:,6) = 0.0         ! x-coord
medium(:,:,7) = 0.0         ! y-coord


medium_a(:,:,1) = .1          ! phiCoarse
medium_a(:,:,2) = 0.0         ! precip toggle
medium_a(:,:,3) = 0.266       ! water_volume
vol_i_a = medium_a(1,1,3)
medium_a(:,:,4) = 0.01        ! reactive fraction now!
medium_a(:,:,5) = 1.0         ! rxn toggle
medium_a(:,:,6) = 0.0         ! x-coord
medium_a(:,:,7) = 0.0         ! y-coord

medium_b(:,:,1) = .1          ! phiCoarse
medium_b(:,:,2) = 0.0         ! precip toggle
medium_b(:,:,3) = 0.0266      ! water_volume
vol_i_b = medium_b(1,1,3)
medium_b(:,:,4) = 0.01        ! reactive fraction now!
medium_b(:,:,5) = 1.0         ! rxn toggle
medium_b(:,:,6) = 0.0         ! x-coord
medium_b(:,:,7) = 0.0         ! y-coord

sea = soluteOcean


!-----------------PERMEABILITY SET UP

slope = param_w
slope2 = x_max-param_w_rhs
buffer = 1000.0
edge = param_w
edge2 = x_max-param_w_rhs

sed = -150.0
sed1 = -300.0

!
!
! sed = (/-77.9798996719, &
! & -77.9798996719, -78.1602449181, -78.5534729091, -78.9467009001, -79.339928891, -79.733156882, &
! & -80.126384873, -80.5196128639, -80.9128408549, -81.3060688459, -81.6992968368, -82.0925248278, &
! & -82.4857528188, -82.8789808098, -83.2722088007, -83.6654367917, -84.0586647827, -84.4518927736, &
! & -84.8451207646, -85.2383487556, -85.6315767466, -86.0248047375, -86.4180327285, -86.8112607195, &
! & -87.2044887104, -87.5977167014, -87.9909446924, -88.3841726834, -88.7774006743, -89.1706286653, &
! & -89.5638566563, -89.9570846472, -90.3503126382, -90.7435406292, -91.1367686201, -91.5299966111, &
! & -91.9232246021, -92.3164525931, -92.709680584, -93.1011873234, -93.3394730885, -93.5777588536, &
! & -93.8160446187, -94.0543303838, -94.2926161489, -94.530901914, -94.7691876791, -95.0074734442, &
! & -95.2457592093, -95.4840449744, -95.7223307394, -95.9606165045, -96.1989022696, -96.4371880347, &
! & -96.6754737998, -96.9137595649, -97.15204533, -97.3903310951, -97.6286168602, -97.8669026253, &
! & -98.1051883904, -98.3434741555, -98.5817599206, -98.8200456857, -99.0583314507, -99.2966172158, &
! & -99.5349029809, -99.773188746, -100.011474511, -100.249760276, -100.488046041, -100.726331806, &
! & -100.964617571, -101.202903337, -101.441189102, -101.679474867, -101.917760632, -102.156046397, &
! & -102.394332162, -102.789284177, -103.189669428, -103.590054678, -103.990439929, -104.390825179, &
! & -104.791210429, -105.19159568, -105.59198093, -105.992366181, -106.392751431, -106.793136682, &
! & -107.193521932, -107.573818396, -107.916943848, -108.260069299, -108.603194751, -108.946320202, &
! & -109.289445654, -109.632571105, -109.975696557, -110.318822008, -110.66194746, -111.005072911, &
! & -111.348198363, -111.691323814, -112.034449265, -112.377574717, -112.720700168, -113.06382562, &
! & -113.406951071, -113.750076523, -114.093201974, -114.436327426, -114.779452877, -115.122578329, &
! & -115.46570378, -115.808829232, -116.151954683, -116.495080135, -116.838205586, -117.181331038, &
! & -117.512288351, -117.831835093, -118.151381835, -118.470928578, -118.79047532, -119.110022062, &
! & -119.429568804, -119.749115547, -120.068662289, -120.388209031, -120.707755773, -121.027302515, &
! & -121.346849258, -121.666396, -121.985942742, -122.305489484, -122.625036226, -122.944582969, &
! & -123.264129711, -123.583676453, -123.903223195, -124.222769937, -124.54231668, -124.861863422, &
! & -125.181410164, -125.500956906, -125.820503649, -126.140050391, -126.459597133, -126.779143875, &
! & -127.098690617, -127.086337378, -127.028533909, -126.970730441, -126.912926973, -126.855123505, &
! & -126.797320037, -126.739516568, -126.6817131, -126.623909632, -126.566106164, -126.508302695, &
! & -126.450499227, -126.392695759, -126.334892291, -126.277088823, -126.219285354, -126.161481886, &
! & -126.103678418, -126.04587495, -125.988071482, -125.930268013, -125.872464545, -125.814661077, &
! & -125.756857609, -125.699054141, -125.092619163, -124.197472067, -123.302324971, -122.407177875, &
! & -121.512030779, -120.616883684, -119.721736588, -118.826589492, -117.931442396, -117.0362953, &
! & -116.141148205, -115.246001109, -114.350854013, -113.455706917, -115.310736614, -117.655613958, &
! & -120.000491303, -122.345368648, -124.690245993, -127.035123338, -129.380000683, -131.724878028, &
! & -134.069755373, -135.279218424, -132.592376319, -129.905534214, -127.218692108, -124.531850003, &
! & -121.845007898, -119.158165793, -118.250469091, -118.389063099, -118.527657108, -118.666251116, &
! & -118.804845124, -118.943439132, -119.08203314, -119.220627148, -119.359221156, -119.497815164, &
! & -119.636409172, -119.77500318, -119.913597188, -120.052191196, -120.190785204, -120.329379212, &
! & -120.46797322, -120.606567228, -120.745161236, -120.883755244, -121.022349252, -121.16094326, &
! & -121.299537268, -121.438131276, -121.576725284, -121.715319292, -121.8539133, -121.992507308, &
! & -122.131101316, -122.278435148, -122.461212204, -122.643989259, -122.826766314, -123.00954337, &
! & -123.192320425, -123.37509748, -123.557874536, -123.740651591, -123.923428646, -124.106205702, &
! & -124.288982757, -124.471759812, -124.654536868, -124.837313923, -125.020090978, -125.202868034, &
! & -125.385645089, -125.568422144, -125.7511992, -125.933976255, -126.11675331, -126.299530366, &
! & -126.482307421, -126.560540475, -126.132477322, -125.704414169, -125.276351015, -124.848287862, &
! & -124.420224709, -123.992161556, -123.564098403, -123.136035249, -122.707972096, -122.279908943, &
! & -121.85184579, -121.423782637, -120.995719483, -120.56765633, -120.139593177, -119.91870732, &
! & -119.860903852, -119.803100383, -119.745296915, -119.687493447, -119.629689979, -119.571886511, &
! & -119.514083042, -119.456279574, -119.398476106, -119.340672638, -119.282869169, -119.225065701, &
! & -119.167262233, -119.109458765, -119.051655297, -118.993851828, -118.93604836, -118.878244892, &
! & -118.820441424, -118.762637956, -118.704834487, -118.647031019, -118.589227551, -118.531424083, &
! & -118.473620615, -118.415817146, -118.358013678, -118.30021021, -118.031844928, -116.368392662, &
! & -114.704940395, -113.041488129, -111.378035863, -109.714583597, -108.051131331, -106.387679065, &
! & -106.219285354, -106.161481886, -106.103678418, -106.04587495, -105.988071482, -105.930268013, &
! & -105.872464545, -105.814661077, -105.756857609, -105.699054141, -105.641250672, -105.583447204, &
! & -105.525643736, -105.467840268, -105.4100368, -105.352233331, -105.294429863, -105.236626395, &
! & -105.178822927, -105.121019458, -105.06321599, -105.005412522, -104.947609054, -104.889805586, &
! & -104.832002117, -104.774198649, -104.716395181, -104.658591713, -104.600788245, -104.542984776, &
! & -104.485181308, -104.42737784, -104.369574372, -104.311770904, -104.253967435, -104.196163967, &
! & -104.138360499, -104.080557031, -104.022753563, -103.964950094, -103.907146626, -103.849343158, &
! & -103.79153969, -103.733736221, -103.675932753, -103.618129285, -103.560325817, -103.502522349, &
! & -103.44471888, -103.386915412, -103.329111944, -103.271308476, -103.213505008, -103.155701539, &
! & -103.097898071, -103.040094603, -102.982291135, -102.924487667, -102.866684198, -102.80888073, &
! & -102.751077262, -102.693273794, -102.635470326, -102.577666857, -102.600548391, -103.720435223, &
! & -104.840322056, -105.960208888, -107.080095721, -108.199982553, -109.319869385, -110.439756218, &
! & -111.55964305, -112.679529883, -113.799416715, -114.919303548, -116.03919038, -117.159077213, &
! & -118.278964045, -119.398850878, -120.51873771, -121.638624543, -122.758511375, -123.878398208, &
! & -124.99828504, -126.118171873, -127.238058705, -128.357945538, -129.47783237, -130.597719203, &
! & -131.717606035, -132.837492868, -133.9573797, -135.077266533, -134.984348317, -134.145805022, &
! & -133.307261727, -132.468718432, -131.630175138, -130.791631843, -129.953088548, -129.114545253, &
! & -128.276001958, -127.437458664, -126.598915369, -125.760372074, -124.921828779, -124.083285485, &
! & -123.24474219, -122.406198895, -121.5676556, -120.729112305, -119.890569011, -119.052025716, &
! & -118.213482421, -117.374939126, -117.622998358, -119.446143834, -121.269289309, -123.092434784, &
! & -124.915580259, -126.738725734, -128.561871209, -130.385016684, -132.208162159, -134.031307634, &
! & -135.854453109, -137.677598584, -139.500744059, -141.323889535, -143.14703501, -144.970180485, &
! & -146.79332596, -148.616471435, -150.43961691, -152.262762385, -154.08590786, -155.909053335, &
! & -157.73219881, -159.555344285, -161.37848976, -163.201635236, -165.024780711, -166.847926186, &
! & -167.331959519, -167.370396604, -167.408833688, -167.447270773, -167.485707857, -167.524144942, &
! & -167.562582026, -167.601019111, -167.639456196, -167.67789328, -167.716330365, -167.754767449, &
! & -167.793204534, -167.831641618, -167.870078703, -167.908515788, -167.946952872, -167.985389957, &
! & -168.023827041, -168.062264126, -168.10070121, -168.139138295, -168.177575379, -168.216012464, &
! & -168.254449549, -168.292886633, -168.331323718, -168.369760802, -168.408197887, -168.446634971, &
! & -168.485072056, -168.52350914, -168.561946225, -168.60038331, -168.638820394, -168.677257479, &
! & -168.715694563, -168.754131648, -168.792568732, -168.831005817, -168.869442902, -168.907879986, &
! & -168.946317071, -168.984754155, -169.02319124, -169.061628324, -169.100065409, -169.138502493, &
! & -169.176939578, -169.215376663, -169.253813747, -169.292250832, -169.330687916, -169.369125001, &
! & -169.407562085, -169.44599917, -169.484436255, -169.522873339, -169.561310424, -169.599747508, &
! & -169.571886511, -169.514083042, -169.456279574, -169.398476106, -169.340672638, -169.282869169, &
! & -169.225065701, -169.167262233, -169.109458765, -169.051655297, -168.993851828, -168.93604836, &
! & -168.878244892, -168.820441424, -168.762637956, -168.704834487, -168.647031019, -168.589227551, &
! & -168.531424083, -168.473620615, -168.415817146, -168.358013678, -168.30021021, -168.242406742, &
! & -168.184603274, -168.126799805, -168.068996337, -168.011192869, -167.953389401, -167.895585932, &
! & -167.837782464, -167.779978996, -167.722175528, -167.66437206, -167.606568591, -167.548765123, &
! & -167.490961655, -167.433158187, -167.375354719, -167.31755125, -167.259747782, -167.201944314, &
! & -167.144140846, -167.086337378, -167.090825388, -167.286262651, -167.481699913, -167.677137175, &
! & -167.872574438, -168.0680117, -168.263448962, -168.458886225, -168.654323487, -168.84976075, &
! & -169.045198012, -169.240635274, -169.436072537, -169.631509799, -169.826947061, -170.022384324, &
! & -170.217821586, -170.413258848, -170.608696111, -170.804133373/)

sed = -50.0

!
! ! USED FOR K AUGUST 18
! sed1 = (/-193.617538317, &
! & -193.617538317, -192.879601041, -192.134041355, -191.388481669, -190.642921983, -189.897362297, &
! & -189.151802611, -188.406242925, -187.660683239, -190.340575274, -195.449174564, -200.557773855, &
! & -205.666373146, -210.774972437, -215.883571728, -220.992171018, -226.100770309, -228.524748033, &
! & -230.069502604, -231.614257176, -233.159011748, -234.703766319, -236.248520891, -237.793275463, &
! & -240.657128024, -244.647403247, -248.637678471, -252.627953695, -256.618228918, -260.608504142, &
! & -263.174969362, -263.117165894, -263.059362426, -263.001558958, -262.943755489, -262.885952021, &
! & -262.828148553, -262.770345085, -262.712541616, -262.654738148, -262.59693468, -267.379886543, &
! & -272.55537109, -277.730855637, -282.906340183, -285.459126794, -282.50892413, -279.558721467, &
! & -276.608518803, -274.119092676, -274.061289208, -274.00348574, -273.945682272, -273.887878803, &
! & -273.830075335, -273.772271867, -283.018879848, -293.26192006, -303.504960272, -313.748000484, &
! & -323.991040696, -330.340596605, -331.885351177, -333.430105748, -332.455119295, -331.327227264, &
! & -330.199335233, -328.351851698, -323.467793727, -322.958644261, -326.358505821, -329.758367381, &
! & -333.158228941, -336.5580905, -339.95795206, -343.35781362, -346.75767518, -350.15753674, &
! & -353.557398299, -356.454063823, -353.539631054, -350.625198285, -347.710765517, -344.796332748, &
! & -341.88189998, -338.967467211, -336.053034442, -333.138601674, -330.224168905, -327.309736137, &
! & -324.395303368, -323.875033588, -322.997584387, -318.113526416, -313.229468445, -308.345410474, &
! & -303.461352503, -298.577294532, -293.693236561, -288.80917859, -283.925120619, -279.041062648, &
! & -275.390486461, -277.029413091, -278.668339721, -280.307266352, -281.946192982, -283.585119613, &
! & -285.224046243, -286.982674685, -289.927139164, -292.871603642, -295.816068121, -298.7605326, &
! & -301.704997079, -304.649461558, -308.971847309, -315.305820438, -321.639793566, -327.973766695, &
! & -330.226992484, -330.169189015, -330.111385547, -323.871028269, -315.213830077, -306.556631885, &
! & -297.899433692, -289.2422355, -280.585037307, -271.927839115, -265.237504932, -264.439017388, &
! & -263.640529843, -262.842042299, -262.043554754, -261.24506721, -265.616041277, -271.03990959, &
! & -276.463777902, -281.887646215, -287.311514528, -292.735382841, -301.177339147, -303.522216492, &
! & -305.867093837, -308.211971181, -308.653453627, -308.595650159, -308.53784669, -308.480043222, &
! & -308.422239754, -308.364436286, -308.306632818, -308.248829349, -308.156820852, -307.085283562, &
! & -306.013746271, -304.94220898, -303.870671689, -302.799134399, -301.727597108, -300.656059817, &
! & -299.584522527, -298.512985236, -297.441447945, -296.369910654, -295.298373364, -294.226836073, &
! & -293.155298782, -292.083761492, -291.026678248, -290.029358457, -289.032038665, -288.034718874, &
! & -287.037399083, -286.040079292, -285.0427595, -284.045439709, -283.048119918, -282.050800127, &
! & -281.053480336, -280.056160544, -279.058840753, -278.061520962, -277.064201171, -276.066881379, &
! & -275.069561588, -274.261907991, -275.078609672, -275.895311353, -276.712013033, -277.528714714, &
! & -277.550046449, -276.820487546, -276.090928644, -275.361369741, -274.631810838, -273.902251936, &
! & -273.172693033, -272.443134131, -271.713575228, -270.984016325, -270.254457423, -269.52489852, &
! & -268.795339617, -268.065780715, -267.336221812, -266.606662909, -265.877104007, -265.147545104, &
! & -264.812734295, -264.754930827, -264.697127358, -264.63932389, -264.878646769, -265.533472363, &
! & -266.188297957, -266.84312355, -267.497949144, -268.152774738, -268.807600332, -269.462425926, &
! & -270.11725152, -270.772077113, -271.426902707, -272.181733237, -273.497739838, -274.813746438, &
! & -275.396575339, -273.515905962, -271.635236584, -269.754567207, -267.873897829, -265.993228452, &
! & -264.112559074, -262.231889697, -260.351220319, -258.470550942, -256.589881564, -254.709212187, &
! & -252.828542809, -250.947873431, -249.067204054, -247.186534676, -248.748286718, -251.093164063, &
! & -253.438041407, -255.782918752, -258.127796097, -258.895302676, -259.274866638, -259.654430599, &
! & -260.03399456, -260.413558521, -260.793122482, -261.172686443, -261.552250404, -261.931814365, &
! & -261.075323668, -259.476156809, -257.876989949, -256.27782309, -254.67865623, -253.079489371, &
! & -251.480322511, -249.881155652, -248.281988792, -246.682821933, -245.49872587, -245.89911112, &
! & -246.299496371, -246.699881621, -247.100266872, -247.500652122, -247.901037373, -248.301422623, &
! & -248.701807873, -248.212876259, -246.622475207, -245.032074155, -243.441673103, -241.851272051, &
! & -240.260871, -238.670469948, -237.080068896, -235.489667844, -233.899266792, -232.30886574, &
! & -230.718464688, -229.128063636, -227.537662584, -225.947261532, -224.35686048, -222.766459428, &
! & -221.176058376, -219.931292901, -221.412004381, -222.892715862, -224.373427343, -225.854138824, &
! & -227.334850304, -228.815561785, -230.296273266, -231.776984746, -233.257696227, -234.738407708, &
! & -235.504122155, -235.927410222, -236.350698289, -236.773986356, -237.197274423, -237.62056249, &
! & -238.043850557, -238.467138624, -239.335072401, -241.096152301, -242.857232201, -244.618312101, &
! & -246.379392001, -248.140471901, -249.901551801, -251.662631701, -253.423711602, -255.184791502, &
! & -256.945871402, -258.706951302, -260.468031202, -262.229111102, -263.990191002, -265.751270902, &
! & -266.274376965, -266.634929729, -266.995482493, -267.356035257, -267.716588021, -268.077140784, &
! & -268.437693548, -268.798246312, -269.158799076, -269.000884849, -267.567000552, -266.133116254, &
! & -264.699231957, -263.265347659, -261.831463362, -260.397579064, -258.963694766, -257.529810469, &
! & -256.095926171, -254.662041874, -253.228157576, -251.794273279, -250.360388981, -248.926504684, &
! & -247.492620386, -246.058736089, -244.624851791, -243.992818985, -243.57847209, -243.164125195, &
! & -242.7497783, -242.335431405, -241.92108451, -241.506737615, -241.09239072, -240.678043825, &
! & -240.263696929, -239.849350034, -237.713360139, -234.038486858, -230.363613578, -226.688740298, &
! & -223.013867018, -219.338993737, -215.664120457, -214.76997019, -214.293600155, -213.81723012, &
! & -213.340860085, -212.86449005, -212.388120014, -211.911749979, -211.435379944, -210.959009909, &
! & -210.83305529, -211.611754109, -212.390452928, -213.169151748, -213.947850567, -214.726549387, &
! & -215.505248206, -216.283947025, -217.062645845, -217.841344664, -218.459700326, -219.003217841, &
! & -219.546735357, -220.090252872, -220.633770387, -221.177287903, -221.720805418, -222.264322933, &
! & -222.807840449, -223.351357964, -223.894875479, -224.438392995, -224.98191051, -224.99703777, &
! & -223.942989106, -222.888940442, -221.834891778, -220.780843115, -219.726794451, -218.672745787, &
! & -217.618697124, -216.56464846, -215.510599796, -214.456551132, -213.402502469, -213.454771956, &
! & -216.954024252, -220.453276549, -223.952528845, -227.451781142, -230.951033438, -234.450285735, &
! & -237.949538031, -241.448790328, -244.948042624, -248.447294921, -251.946547218, -253.492677604, &
! & -254.953163926, -256.413650248, -257.874136571, -259.334622893, -260.795109215, -262.255595537, &
! & -263.739468376, -266.08434572, -268.429223065, -270.77410041, -273.118977755, -275.4638551, &
! & -277.808732445, -281.961968645, -289.568412833, -297.174857022, -304.781301211, -312.387745399, &
! & -319.994189588, -327.432343793, -332.172977799, -336.913611804, -341.65424581, -346.394879816, &
! & -346.932194796, -346.874391327, -346.816587859, -346.758784391, -346.700980923, -346.060020539, &
! & -344.717945722, -343.375870905, -342.033796088, -340.691721271, -339.349646455, -338.062452511, &
! & -336.990915221, -335.91937793, -334.847840639, -333.776303348, -332.704766058, -331.633228767, &
! & -330.561691476, -329.920906789, -330.369507763, -330.818108736, -331.266709709, -331.715310682, &
! & -332.163911656, -332.612512629, -333.061113602, -333.934228912, -335.573155543, -337.212082173, &
! & -338.851008804, -340.489935434, -342.128862064, -343.767788695, -346.458679173, -352.155683781, &
! & -357.852688389, -363.549692998, -368.772897044, -368.073205448, -367.373513852, -366.673822256, &
! & -365.97413066, -365.274439064, -364.574747468, -362.204707337, -359.469101032, -356.733494728, &
! & -353.997888424, -351.26228212, -348.526675815, -345.791069511, -342.356682616, -337.472624645, &
! & -332.588566674, -327.704508703, -323.423962666, -322.162212847, -320.900463028, -319.638713208, &
! & -318.376963389, -317.11521357, -315.853463751, -316.2181465, -318.220208356, -320.222270212, &
! & -322.224332067, -324.226393923, -326.228455779, -328.362032987, -330.565768132, -332.769503276, &
! & -334.97323842, -337.176973565, -339.380708709, -341.584443853, -343.088699797, -344.025950793, &
! & -344.963201789, -345.900452785, -346.837703781, -347.774954777, -348.712205773, -349.649456769, &
! & -350.586707765, -351.523958761, -352.461209757, -353.398460753, -353.738482114, -353.348729323, &
! & -352.958976532, -352.56922374, -352.179470949, -351.789718158, -351.399965367, -351.010212576, &
! & -350.620459785, -350.230706994, -349.840954203, -349.451201412, -349.180107448, -349.12230398, &
! & -349.064500512, -349.006697044, -348.948893575, -348.891090107, -348.833286639, -348.775483171, &
! & -348.884247832, -359.188221025, -369.492194218, -379.796167411, -390.100140604, -400.404113796, &
! & -404.236682923, -407.936532472, -411.63638202, -415.336231569, -419.036081117, -422.735930666, &
! & -426.435780215, -430.135629763, -433.835479312, -436.086602814/)

sed1 = -200.0

!sed1 = sum(sed1)/xn

sed1((param_w/dx)+1:) = sed1(:xn-(param_w/dx))
sed1(1:param_w/dx) = sed1((param_w/dx)+1)

sed((param_w/dx)+1:) = sed(:xn-(param_w/dx))
sed(1:param_w/dx) = sed((param_w/dx)+1)

sed1(f_index1-2:f_index1+2) = sum(sed1(f_index1-3:f_index1+2))/6.0
sed(f_index1-2:f_index1+2) = sum(sed(f_index1-3:f_index1+2))/6.0

sed2 = sed-sed1!-sed1! - sed1


sed2((param_w/dx)+1:) = sed2(:xn-(param_w/dx))
sed2(1:param_w/dx) = sed2((param_w/dx)+1)

!sed2 = (-sum(sed2-sed1)/xn)*(sed2-sed1)

!sed = sum(sed)/xn
!sed = -100.0
 
if (param_o_rhs .gt. param_o) then
	sed = sed-(param_o_rhs)
end if

if (param_o .ge. param_o_rhs) then
	sed = sed-(param_o)
end if


if (param_o_rhs .gt. param_o) then
	sed1 = sed1-(param_o_rhs)
end if

if (param_o .ge. param_o_rhs) then
	sed1 = sed1-(param_o)
end if

do gg=1,yn
	do g=1,xn
		if ((y(gg) .ge. sed1(g)) .and. (x(g) .gt. edge) .and. (x(g) .lt. edge2)) then
			 lambdaMat(g,gg) = 1.2
		end if
	end do
end do


! with sediment cap
sed3 = sed1 - (param_h)

	
	! the mask
	mask = 1.0
	do gg=1,yn
		do g =2,xn
			if ((x(g) .ge. edge) .and. (x(g) .le. edge2+5000.0) .and. (y(gg) .ge. sed(g))) then
				mask(g,gg) = 0.0
			end if
		end do
	end do

	do g =1,xn
			! left outcrop top
			if ((x(g) .lt. edge)) then
				if (param_o_rhs .gt. param_o) then
					mask(g,yn-(abs(param_o-param_o_rhs)/dy)-1) = 25.0
					mask(g,yn-(abs(param_o-param_o_rhs)/dy):yn) = 0.0
				end if
				if (param_o_rhs .lt. param_o) then
					mask(g,yn-1) = 25.0
					mask(g,yn:yn) = 0.0
				end if
				if (param_o_rhs .eq. param_o) then
					mask(g,yn-1) = 25.0
					mask(g,yn-0:yn) = 0.0
				end if
			end if
	end do

	! inner vertical edges
	do gg=1,yn
		do g =2,xn-1
			if ((mask(g,gg) - mask(g-1,gg) .eq. -1.0) .and. (x(g) .le. edge+dx)) then
				mask(g-1,gg) = 5.0
			end if
			if ((mask(g,gg) - mask(g-1,gg) .eq. 1.0) .and. (x(g) .ge. edge2-dx)) then
				mask(g,gg) = 10.0
			end if
		end do
	end do

	do gg=2,yn-1
		do g =1,xn
			if ((mask(g,gg) .eq. 1.0) .and. (mask(g,gg+1) .eq. 0.0)) then
				mask(g,gg) = 50.0
			end if
		end do
	end do
	
	do gg=2,yn-1
		do g =2,xn-1
			! left upper corner
			if ((mask(g,gg) .eq. 25.0) .and. (mask(g+1,gg-1) .eq. 5.0)) then
				mask(g+1,gg) = 12.5
			end if
			! right upper corner
			if ((mask(g,gg) .eq. 25.0) .and. (mask(g-1,gg-1) .eq. 10.0)) then
				mask(g-1,gg) = 17.5
			end if
		end do
	end do
		
	do gg=2,yn-1
		do g =2,xn-1 
			! left bottom corner
			if ((mask(g,gg) .eq. 5.0) .and. (mask(g+1,gg-1) .eq. 50.0)) then
				mask(g,gg-1) = 2.5
			end if
		
			! right bottom corner
			if ((mask(g,gg) .eq. 10.0) .and. (mask(g-1,gg-1) .eq. 50.0)) then
				mask(g,gg-1) = 7.5
			end if
		end do
	end do
	


	maskP = mask

	do gg=1,yn-1
		do g =1,xn
			if (y(gg) .le. sed3(g)) then
				maskP(g,gg) = 0.0
			end if
		end do
	end do
	
	do gg=2,yn-3
		do g =1,xn
			if (gg .eq. yn/2 + 3) then
				maskP(:,gg) = 100.0
				code = gg
			end if
		end do
	end do
	
	do gg=2,yn-3
		do g =1,xn
			if ((y(gg) .le. sed3(g)) .and. (gg .gt. code)) then
				maskP(:,gg) = 1.0
			end if
		end do
	end do
	
	

	! coarse mask
	coarse_mask = 0.0
	do gg=1,yn/celly
		do g =1,xn/cellx
			if (maxval(maskP((g-1)*cellx+1:g*cellx,(gg-1)*celly+1:gg*celly)) .gt. 0.0) then
				coarse_mask(g,gg) = 1.0
			end if
			if (coarse_mask(g,gg) .eq. 0.0) then
				medium(g,gg,5) = 0.0    ! cell toggle
				primary(g,gg,:) = 0.0   ! basaltic glass
				solute(g,gg,3) = 0.0    ! solute water
				medium(g,gg,3) = 0.0    ! medium water
			
				medium_a(g,gg,5) = 0.0    ! cell toggle
				primary_a(g,gg,:) = 0.0   ! basaltic glass
				solute_a(g,gg,3) = 0.0    ! solute water
				medium_a(g,gg,3) = 0.0    ! medium water
			
				medium_b(g,gg,5) = 0.0    ! cell toggle
				primary_b(g,gg,:) = 0.0   ! basaltic glass
				solute_b(g,gg,3) = 0.0    ! solute water
				medium_b(g,gg,3) = 0.0    ! medium water
			end if
			
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
			
		end do
	end do
	
! 	do gg=1,yn/celly
! 		do g =1,xn/cellx-1
! 			if ((coarse_mask(g,gg) .eq. 1.0) .and. (coarse_mask(g+1,gg) .eq. 0.0)) then
! 		end do
! 	end do
	
	
	
	
	

	fives = 1.0
	tens = 1.0
	twentyfives = 1.0
	fifties = 1.0
	
	do gg=1,yn
		do g=1,xn
			if ((maskP(g,gg) .eq. 5.0) .or. (maskP(g,gg) .eq. 12.5)) then
				fives(1,gg) = g
				fives(2,gg) = gg
			end if
			if ((maskP(g,gg) .eq. 10.0) .or. (maskP(g,gg) .eq. 17.5)) then
				tens(1,gg) = g
				tens(2,gg) = gg
			end if
		end do
	end do
		
	
	do g=1,xn
		do gg=1,yn
			if ((maskP(g,gg) .eq. 25.0) .or. (maskP(g,gg) .eq. 12.5) .or. (maskP(g,gg) .eq. 17.5)) then
				twentyfives(1,g) = g
				twentyfives(2,g) = gg
			end if
			if (maskP(g,gg) .eq. 50.0) then
				fifties(1,g) = g
				fifties(2,g) = gg
			end if
		end do
	end do

	! TOO SIMPLE
	permeability = param_paq
	do gg=1,yn-2
		do g=1,xn
			if ((any(maskP(g,:) .eq. 50.0)) .and. (y(gg) .ge. sed1(g)) .and. (y(gg) .le. sed(g)) .and. (x(g) .le. x_max-param_w_rhs+5000.0)) then
				permeability(g:g,gg:gg) = 1e-17
			end if
			if ((y(gg) .le. sed3(g))) then
				permeability(g,gg) = 1e-16
			end if
		end do
	end do
	
	
	! 369 fracture goes HERE
	do gg=2,yn-3
		do g =1,xn
			if ((maskP(g,gg) .eq. 1.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
				maskP(g,gg) = 6.0
				maskP(g-1,gg) = 3.0
				mask(g,gg) = 6.0
				mask(g-1,gg) = 3.0
			end if
			!if ((maskP(g,gg) .eq. 100.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
			if ((maskP(g,gg-1) .eq. 1.0) .and. (maskP(g,gg) .eq. 6.0)) then
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
			end if
			if ((maskP(g,gg) .eq. 50.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
				maskP(g,gg) = 6.5
				maskP(g-1,gg) = 3.5
				mask(g,gg) = 6.5
				mask(g-1,gg) = 3.5
			end if
		end do
	end do
	
	
	
	! coarse mask 2
	do gg=1,yn/celly
		do g =1,xn/cellx
			if (minval(permeability((g-1)*cellx+1:g*cellx,(gg-1)*celly+1:gg*celly)) .lt. param_paq) then
				coarse_mask(g,gg) = 0.0
			end if
			if (coarse_mask(g,gg) .eq. 0.0) then
				medium(g,gg,5) = 0.0    ! cell toggle
				primary(g,gg,:) = 0.0   ! basaltic glass
				solute(g,gg,3) = 0.0    ! solute water
				medium(g,gg,3) = 0.0    ! medium water
			
				medium_a(g,gg,5) = 0.0    ! cell toggle
				primary_a(g,gg,:) = 0.0   ! basaltic glass
				solute_a(g,gg,3) = 0.0    ! solute water
				medium_a(g,gg,3) = 0.0    ! medium water
			
				medium_b(g,gg,5) = 0.0    ! cell toggle
				primary_b(g,gg,:) = 0.0   ! basaltic glass
				solute_b(g,gg,3) = 0.0    ! solute water
				medium_b(g,gg,3) = 0.0    ! medium water
			end if
		end do
	end do

	! high lambda in deep basalt
	do gg=1,yn
		do g=1,xn
			if ((y(gg) .lt. sed1(g)) .and. (permeability(g,gg) .eq. 1e-18)) then
				lambdaMat(g,gg) = 1.8
			end if
		end do
	end do
	
	active_cells = 0
	do gg=1,yn
		do g=1,xn
			if (maskP(g,gg) .ne. 0.0) then
				active_cells = active_cells + 1
			end if
		end do
	end do

	long = (xn-2)*(yn-2)
	longP = (xn-2)*((yn/2)-2)

	maskLong = reshape(mask(2:xn-1,2:yn-1), (/long/))
	maskLongT = reshape(transpose(mask(2:xn-1,2:yn-1)), (/long/))
	
	maskLongU = reshape(mask(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
	maskLongTV = reshape(transpose(mask(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))
	
	maskPLong = reshape(maskP(2:xn-1,(yn/2)+2:yn-1), (/longP/))
	maskPLongT = reshape(transpose(maskP(2:xn-1,(yn/2)+2:yn-1)), (/longP/))
	


return

end subroutine init






function h_bc(h_in)
	
	use globals
	real(4) :: h_in(xn,yn), h_bc(xn,yn), rip_lith_y(xn)
	integer :: p, pp
	
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
	do pp=1,yn
		do p=1,xn
			if (mask(p,pp) .eq. 0.0) then
				 h_bc(p,pp) = param_tsw
			end if
		end do
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!    VERTICAL OUTER   !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	do pp=1,yn	
		if (mask(xn,pp) .ne. 0.0) then
			h_bc(xn,pp) = (4.0/3.0)*h_in(xn-1,pp) - (1.0/3.0)*h_in(xn-2,pp) ! right
		end if
	end do
	
	do pp=1,yn	
		if (mask(1,pp) .ne. 0.0) then
			h_bc(1,pp) = (4.0/3.0)*h_in(2,pp) - (1.0/3.0)*h_in(3,pp) ! left 
		end if
	end do


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL OUTER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 
	! bottom
	do p = 1,xn
		h_bc(p,1) = h_in(p,2) + ( rip_lith_y(p)) * dy/1.8
	end do
	
	! two lines recent...
 	h_bc(xn,1) = h_bc(xn,2)
 	h_bc(1,1) = h_bc(1,2)

	! top of outcrops
	do p=1,xn
	 	do pp=2,yn-1
			if (mask(p,pp) .eq. 25.0) then
				 h_bc(p,pp+1) = param_tsw
			end if

			if (mask(p,pp) .eq. 12.5) then
				 h_bc(p,pp+1) = param_tsw
			end if
			if (mask(p,pp) .eq. 17.5) then
				 h_bc(p,pp+1) = param_tsw
			end if
		end do
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL INNER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! top of sediment
	do p=2,xn-1
	 	do pp=2,yn-1
			if ((mask(p,pp) .eq. 50.0) .or. (mask(p,pp) .eq. 3.5) .or. (mask(p,pp) .eq. 6.5)) then
				 h_bc(p,pp+1) = param_tsw
			end if
		end do
	end do

	do p=1,xn
	 	do pp=2,yn-1
			if ((mask(p,pp) .eq. 12.5) .or. (mask(p,pp) .eq. 17.5)) then
				 h_bc(p,pp+1) = param_tsw
			end if
		end do
	end do

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!    VERTICAL INNER   !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do p=2,xn-1
	 	do pp=2,yn
			! left outcrop
			if ((mask(p,pp) .eq. 12.5)) then
				h_bc(p+1,pp) = param_tsw
			end if
			if ((mask(p,pp) .eq. 5.0)) then
				h_bc(p+1,pp) = param_tsw
			end if
			! right outcrop
			if ((mask(p,pp) .eq. 17.5)) then
				h_bc(p-1,pp) = param_tsw
			end if
			if ((mask(p,pp) .eq. 10.0)) then
				h_bc(p-1,pp) = param_tsw
			end if
		end do
	end do

		
return
end function h_bc










function psi_bc(psi_in)
	
	use globals
	real(4) :: psi_in(xn,yn), psi_bc(xn,yn)
	integer :: p,pp
	
	psi_bc = psi_in
	
	
	do pp=1,yn
	    do p=1,xn
			if ((maskP(p,pp) .eq. 0.0)) then
				psi_bc(p,pp) = 0.0
			end if
		end do
	end do

    do p=1,xn
    	do pp=2,yn-1
    		! top of outcrops
    		if ((maskP(p,pp) .eq. 25.0) .or. (maskP(p,pp) .eq. 12.5) .or. (maskP(p,pp) .eq. 17.5)) then
   				psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
    		end if
    	end do
    end do
	
	
	do pp=1,yn
		! left
		if (maskP(1,pp) .ne. 0.0) then
			psi_bc(1,pp) = 0.0
		end if
		! right
		if (maskP(xn,pp) .ne. 0.0) then
			psi_bc(xn,pp) =0.0
		end if
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL OUTER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do p=1,xn
		do pp=2,yn-1
			if ((maskP(p,pp) .eq. 100.0) .or. (maskP(p,pp) .eq. 3.01) .or. (maskP(p,pp) .eq. 6.01)) then
				psi_bc(p,pp-1) = 0.0
			end if
		end do
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL INNER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do p=2,xn-1
		do pp=2,yn-1
			! right outcrop (left boundary)
			if ((maskP(p,pp) .eq. 10.0)) then
				psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
			end if

			if ((maskP(p,pp) .eq. 17.5)) then
				psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
			end if

			! left outcrop (right boundary)
			if ((maskP(p,pp) .eq. 5.0)) then
				psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
			end if

			if ((maskP(p,pp) .eq. 12.5)) then
				psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
			end if
		end do
	end do
	
	do p=2,xn
		do pp=2,yn-1
			! top of sediment
			if ((maskP(p,pp) .eq. 50.0) .or. (maskP(p,pp) .eq. 3.5) .or. (maskP(p,pp) .eq. 6.5)) then
				psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
			end if
		end do
	end do		

	psi_bc(f_index1,:) = 0.0
	

	return
end function psi_bc









function psi_mod(psi_in)
	
	use globals
	real(4) :: psi_in(xn,yn), psi_mod(xn,yn)
	integer :: p,pp
	
	psi_mod = psi_in
	
	do pp=1,yn
    do p=1,xn
		if (mask(p,pp) .eq. 0.0) then
			psi_mod(p,pp) = dx*(5.0e-12)*1000.0 - dx*(5.0e-12)*0.0
		end if
	end do
	end do
	



	return
end function psi_mod

end module initialize