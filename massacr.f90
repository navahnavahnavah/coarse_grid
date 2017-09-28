! ----------------------------------------------------------------------------------%%
!
! MASSACR
!
! SUMMARY: main method runs fluid dynamic simulation coupled to geochemical
!          simulation and writes selected model output to file
!
! TO RUN: make -f berMakeFile
!		  qsub kinsub.csh
!
! ----------------------------------------------------------------------------------%%

PROGRAM main

  USE globals
  USE initialize
  !use alteration
  !use netcdf

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE "IPhreeqc.f90.inc"
  SAVE

  ! functions within massacr.f90
  INTERFACE

     ! solves thermal energy equation
     FUNCTION h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, n, ii, m=3
       REAL(4) :: sx, sy, qx, qy, rho_in(xn,yn), phi_in(xn,yn)
       REAL(4) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
       REAL(4) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
       REAL(4) ::  velocities0(xn,2*yn)
       REAL(4) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
       REAL(4) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
       REAL(4) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
       REAL(4) :: kMatLong((xn-2)*(yn-2))
       REAL(4) :: mn(xn,yn)
       REAL(4) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
       REAL(4) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
       REAL(4) :: frac6_in(yn,2), temp6_in(yn,2), dt_in
     END FUNCTION h_next

     ! solves streamfunction vorticity equation
     FUNCTION psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage,frac6_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, ii, n, m, stage
       REAL(4) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
       REAL(4) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
       REAL(4) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
       REAL(4) :: psi_next(xn,yn)
       REAL(4) :: mn(xn,yn)
       REAL(4) :: aBand0(longP,4*((yn/2)-2) + 3), band_in(longP,4*((yn/2)-2) + 3)
       REAL(4) :: rhoLong(longP)
       REAL(4) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
     END FUNCTION psi_next

     FUNCTION make_band(perm_in,phi_in,permx,permy,rho_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, ii, n, m
       REAL(4) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
       REAL(4) :: permx(xn,yn), permy(xn,yn), permLong(longP)
       REAL(4) :: permxLong(longP), permyLong(longP)
       REAL(4) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
       REAL(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
       REAL(4) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
       REAL(4) :: perm_long(longP)
     END FUNCTION make_band

     FUNCTION particles_next (trace_in, uTransport, vTransport, inval, num, num_sat)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, ii, n, m, mm, nn, num, num_sat
       REAL(4) :: trace_in(5,num), particles_next(5,num)
       REAL(4) :: uTransport(xn,yn), vTransport(xn,yn)
       REAL(4) :: u_wt, v_wt, inval
       REAL(4) :: rando
     END FUNCTION particles_next

     FUNCTION solute_next(sol, uTransport, vTransport, seaw)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, ii, n, m
       REAL(4) :: sol(xn,yn), sol0(xn,yn)
       REAL(4) :: uTransport(xn,yn), vTransport(xn,yn)
       REAL(4) :: uLong(((xn)-2)*((yn)-0)), vLong(((xn)-0)*((yn)-2))
       REAL(4) :: aBand(((xn)-2)*((yn)-0),5), bBand(((xn)-0)*((yn)-2),5)
       REAL(4) :: qx, qy, solute_next(xn,yn), vec(((xn)-2)*((yn)-0))
       REAL(4) :: sol_nextRow(((xn)-2)*((yn)-0)), sol_nextRowB(((xn)-0)*((yn)-2))
       REAL(4) :: seaw
       REAL(4) :: bm1(xn,yn), b0(xn,yn), bp1(xn,yn), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
       REAL(4) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6
     END FUNCTION solute_next

     FUNCTION solute_next_coarse (sol, uTransport, vTransport, phiTransport, seaw)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i, j, ii, n, m
       REAL(4) :: sol((xn-1)/cellx,yn/(2*celly)), sol0((xn-1)/cellx,yn/(2*celly))
       REAL(4) :: uTransport((xn-1)/cellx,yn/(2*celly)), vTransport((xn-1)/cellx,yn/(2*celly))
       real(4) :: phiTransport((xn-1)/cellx,yn/(2*celly))
       REAL(4) :: uLong((((xn-1)/cellx)-2)*((yn/(2*celly))-0)), vLong((((xn-1)/cellx)-0)*((yn/(2*celly))-2))
       REAL(4) :: aBand((((xn-1)/cellx)-2)*((yn/(2*celly))-0),5), bBand((((xn-1)/cellx)-0)*((yn/(2*celly))-2),5)
       REAL(4) :: qx, qy, solute_next_coarse((xn-1)/cellx,yn/(2*celly)), vec((((xn-1)/cellx)-2)*((yn/(2*celly))-0))
       REAL(4) :: sol_nextRow((((xn-1)/cellx)-2)*((yn/(2*celly))-0)), sol_nextRowB((((xn-1)/cellx)-0)*((yn/(2*celly))-2))
       REAL(4) :: seaw
       REAL(4) :: bm1((xn-1)/cellx,yn/(2*celly)), b0((xn-1)/cellx,yn/(2*celly)), bp1((xn-1)/cellx,yn/(2*celly)), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
       REAL(4) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6
     END FUNCTION solute_next_coarse

     ! calculates fluid density
     FUNCTION rho_next (h_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i,j
       REAL(4) :: h_in(xn,yn), rho_next(xn,yn)
     END FUNCTION rho_next

     ! calculates fluid density
     FUNCTION rho_one(h_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i,j
       REAL(4) :: h_in, rho_one
     END FUNCTION rho_one

     ! calculates viscosity
     FUNCTION visc_next(h_in)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: i,j
       REAL(4) :: h_in(xn,yn), visc_next(xn,yn)
     END FUNCTION visc_next

     ! calculates velocities from streamfunction values
     FUNCTION velocities(psi)
       USE globals
       USE initialize
       IMPLICIT NONE
       REAL(4) :: velocities(xn,2*yn), psi(xn,yn)
       REAL(4) :: u0(xn,yn), v0(xn,yn)
     END FUNCTION velocities

     ! calculates velocities from COARSE streamfunction values
     FUNCTION velocities_coarse(psi_coarse)
       USE globals
       USE initialize
       IMPLICIT NONE
       REAL(4) :: velocities_coarse((xn-1)/cellx,yn/celly), psi_coarse((xn-1)/cellx,yn/(2*celly))
       REAL(4) :: u0((xn-1)/cellx,yn/(2*celly)), v0((xn-1)/cellx,yn/(2*celly))
     END FUNCTION velocities_coarse

     ! calculates partial derivative of any 1D or 2D array
     FUNCTION partial(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial(rows,cols)
     END FUNCTION partial

     FUNCTION partial_coarse(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_coarse(rows,cols)
     END FUNCTION partial_coarse

     ! calculates partial derivative of any 1D or 2D array
     FUNCTION partial_edge(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge(rows,cols)
     END FUNCTION partial_edge

     ! calculates partial derivative of any 1D or 2D array
     FUNCTION partial_edge_p(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge_p(rows,cols)
     END FUNCTION partial_edge_p

     ! writes 2D array to file
     FUNCTION write_matrix ( m, n, table, filename )
       USE globals
       IMPLICIT NONE
       INTEGER :: m, n, j, output_status, unit0, reclen
       CHARACTER ( len = * ) filename
       CHARACTER ( len = 30 ) string
       REAL(4)  :: table(m,n) , write_matrix
     END FUNCTION write_matrix

     ! writes 1D array to file
     FUNCTION write_vec ( n, vector, filename )
       USE globals
       IMPLICIT NONE
       INTEGER :: n, j, output_status, unit0
       CHARACTER ( len = * ) filename
       REAL(4)  :: vector(n), write_vec
     END FUNCTION write_vec

  END INTERFACE


  ! dependent variable arrays
  REAL(4) :: h(xn,yn), psi(xn,yn), pside(xn,yn) ! xn rows deep & yn columns wide
  REAL(4) :: hmat((xn*tn/(mstep*ar)),yn), psimat((xn*tn/(mstep*ar)),yn), rhomat((xn*tn/(mstep*ar)),yn)
  REAL(4) :: velocities0(xn,2*yn)
  REAL(4) :: phimat((xn*tn/(mstep*ar)),yn), umat((xn*tn/(mstep*ar)),yn), vmat((xn*tn/(mstep*ar)),yn), rhsmat((xn*tn/(mstep*ar)),yn)
  REAL(4) :: permxMat((xn*tn/(mstep*ar)),yn), permyMat((xn*tn/(mstep*ar)),yn)
  REAL(4) :: psiCoarseMat((xn*tn/(mstep*ar)),yn), uCoarseMat((xn*tn/(mstep*ar)),yn), vCoarseMat((xn*tn/(mstep*ar)),yn)
  REAL(4) :: permmat((xn*tn/(mstep*ar)),yn)
  REAL(4) :: u(xn,yn), v(xn,yn),  permx(xn,yn), permy(xn,yn), uTest(xn,yn), vTest(xn,yn), psiTest(xn,yn), hTest(xn,yn)


  ! autumn performance review
  INTEGER :: counti, countf, count_rate, count_max
  REAL :: timeBit

  ! material properties
  REAL(4) :: rho(xn,yn), visc(xn,yn)
  REAL(4) :: rhs0(xn,yn)
  INTEGER :: unit
  REAL(4) :: phi_coarse((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: phi0(xn,yn), phi(xn,yn)

  ! netCDF & output stuff
  INTEGER :: xInt, yInt, tInt, hInt, uInt, vInt
  INTEGER :: ncid
  INTEGER :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
  INTEGER :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid
  INTEGER :: i, j, ii, m, n, jj, jjj
  REAL(4) :: yep

  ! benchmark stuff
  REAL(4) :: nusseltLocalv(xn,1), nuBar

  ! geochemical alteration stuff
  REAL(4) :: alt0(1,altnum), alt_mat(3*((xn-1)/cellx)*(yn/(2*celly)),altnum)
  REAL(4) :: primaryShift((xn-1)/cellx,yn/celly,g_pri), secondaryShift((xn-1)/cellx,yn/celly,g_sec)

  ! solute transport stuff
  REAL(4) :: uTransport((xn-1)/cellx,yn/celly), vTransport((xn-1)/cellx,yn/celly)
  REAL(4) :: u_coarse((xn-1)/cellx,yn/(2*celly)), v_coarse((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: psi_coarse((xn-1)/cellx,yn/(2*celly)), velocities_coarse0((xn-1)/cellx,yn/celly)
  REAL(4) :: permeability0(xn,yn)

  ! message passing stuff
  INTEGER, PARAMETER :: max_rows = 10000000
  INTEGER, PARAMETER :: send_data_tag = 2001, return_data_tag = 2002
  INTEGER :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
  INTEGER :: num_procs, an_id, num_rows_to_receive, an_id0
  INTEGER :: avg_rows_per_process, num_rows, num_rows_to_send
  INTEGER :: end_row, sender, start_row, num_rows_received
  REAL(4) :: vector(max_rows), vector2(max_rows), partial_sum, sum
  REAL(4) :: local_mean, global_mean
  REAL(4) :: hLocal(((xn-1)/cellx)*(yn/(2*celly))), dt_local
  INTEGER :: order

  ! formatted message passing arrays
  REAL(4) :: hLong(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: priLong(((xn-1)/cellx)*(yn/(2*celly)),g_pri), priLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: secLong(((xn-1)/cellx)*(yn/(2*celly)),g_sec), secLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_sec)
  REAL(4) :: solLong(((xn-1)/cellx)*(yn/(2*celly)),g_sol), solLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_sol)
  REAL(4) :: medLong(((xn-1)/cellx)*(yn/(2*celly)),g_med), medLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_med)
  REAL(4) :: ph_fix_Local(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: priLongBit(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: secLongBit(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: solLongBit(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: medLongBit(3*((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: dpriLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: dsecLocal(3*((xn-1)/cellx)*(yn/(2*celly)),g_sec)


  REAL(4) :: sol_coarse_long(((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: sol_coarse_long_local(((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: sol_coarse_local((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: u_coarse_long(((xn-1)/cellx)*(yn/(2*celly))), v_coarse_long(((xn-1)/cellx)*(yn/(2*celly))), phi_coarse_long(((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: u_coarse_local((xn-1)/cellx,yn/(2*celly)), v_coarse_local((xn-1)/cellx,yn/(2*celly)), phi_coarse_local((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: u_coarse_long_local(((xn-1)/cellx)*(yn/(2*celly))), v_coarse_long_local(((xn-1)/cellx)*(yn/(2*celly))), phi_coarse_long_local(((xn-1)/cellx)*(yn/(2*celly)))
  INTEGER :: an_id_local

  REAL(4) :: priLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: secLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_sec)
  REAL(4) :: solLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_sol)
  REAL(4) :: volLongBitFull(((xn-1)/cellx)*(yn/(2*celly)))
  REAL(4) :: medLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_med)
  REAL(4) :: phiLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)))

  REAL(4) :: dpriLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: dsecLongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)),g_sec)

  REAL(4) :: ph_fix((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: ph_fix_a((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: ph_fix_b((xn-1)/cellx,yn/(2*celly))

  REAL(4) :: ph_fix_LongBitFull(3*((xn-1)/cellx)*(yn/(2*celly)))

  ! chamber stuff
  REAL(4) :: priLong_a(((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: secLong_a(((xn-1)/cellx)*(yn/(2*celly)),g_sec)
  REAL(4) :: solLong_a(((xn-1)/cellx)*(yn/(2*celly)),g_sol)
  REAL(4) :: medLong_a(((xn-1)/cellx)*(yn/(2*celly)),g_med)

  REAL(4) :: priLong_b(((xn-1)/cellx)*(yn/(2*celly)),g_pri)
  REAL(4) :: secLong_b(((xn-1)/cellx)*(yn/(2*celly)),g_sec)
  REAL(4) :: solLong_b(((xn-1)/cellx)*(yn/(2*celly)),g_sol)
  REAL(4) :: medLong_b(((xn-1)/cellx)*(yn/(2*celly)),g_med)


  REAL(4) :: cml3((3*(xn-1)/cellx)*(yn/(2*celly)))



  INTEGER(KIND=4) :: id, all=136
  CHARACTER(LEN=5000) :: line
  CHARACTER(len=20000) :: inputz0
  REAL(4) :: outmat(4,136)
  REAL(4) :: temp3, timestep3, primary3(5), secondary3(80), solute3(15), medium3(7), ph_fix3 ! important information
  REAL(4) :: dsecondary3(80), dprimary3(5)
  REAL(4) :: water

  ! STRINGS
  CHARACTER(len=25) :: s_verm_ca, s_analcime, s_phillipsite, s_clinozoisite, s_verm_na
  CHARACTER(len=25) :: s_diopside, s_epidote, s_minnesotaite, s_ferrite_ca, s_foshagite
  CHARACTER(len=25) :: s_gismondine, s_gyrolite, s_hedenbergite, s_chalcedony, s_verm_mg
  CHARACTER(len=25) :: s_ferrihydrite, s_lawsonite, s_merwinite, s_monticellite, s_natrolite
  CHARACTER(len=25) :: s_talc, s_smectite_low, s_prehnite, s_chlorite, s_rankinite
  CHARACTER(len=25) :: s_scolecite, s_tobermorite_9a, s_tremolite, s_chamosite7a, s_clinochlore14a
  CHARACTER(len=25) :: s_clinochlore7a, s_andradite
  CHARACTER(len=25) :: s_saponite_ca, s_troilite, s_pyrrhotite, s_lepidocrocite, s_daphnite_7a
  CHARACTER(len=25) :: s_daphnite_14a, s_verm_k, s_greenalite, s_aragonite
  CHARACTER(len=25) :: s_siderite, s_kaolinite, s_goethite, s_dolomite, s_celadonite ! secondary
  CHARACTER(len=25) :: s_sio2, s_albite, s_calcite, s_mont_na, s_smectite, s_saponite ! secondary
  CHARACTER(len=25) :: s_stilbite, s_saponite_k, s_anhydrite, s_clinoptilolite, s_pyrite ! secondary
  CHARACTER(len=25) :: s_quartz, s_kspar, s_saponite_na, s_nont_na, s_nont_mg, s_nont_k ! secondary
  CHARACTER(len=25) :: s_fe_celadonite, s_nont_ca, s_muscovite, s_mesolite, s_hematite, s_diaspore !
  CHARACTER(len=25) :: s_feldspar, s_pigeonite, s_augite, s_glass, s_magnetite ! primary
  CHARACTER(len=25) :: s_laumontite, s_mont_k, s_mont_mg, s_mont_ca
  ! mid may 2017 minerals
  CHARACTER(len=25) :: s_fe_saponite_ca, s_fe_saponite_mg
  CHARACTER(len=25) :: s_temp, s_timestep ! important information
  CHARACTER(len=25) :: s_ph, s_ca, s_mg, s_na, s_k, s_fe, s_s, s_si, s_cl, s_al, s_alk, s_co2 ! solutes
  CHARACTER(len=25) :: s_hco3, s_co3, s_pe
  CHARACTER(len=25) :: s_water, s_w, s_reactive ! medium
  CHARACTER(len=25) :: s_basalt1, s_basalt2, s_basalt3, s_precip, s_clay, s_chlor, s_ph_fix

  ! staggered equilibrium
  CHARACTER(len=25) :: sd_dbasalt1, sd_dbasalt2, sd_dbasalt3, sd_dglass

  CHARACTER(len=25) :: sd_verm_ca, sd_analcime, sd_phillipsite, sd_clinozoisite, sd_verm_na
  CHARACTER(len=25) :: sd_diopside, sd_epidote, sd_minnesotaite, sd_ferrite_ca, sd_foshagite
  CHARACTER(len=25) :: sd_gismondine, sd_gyrolite, sd_hedenbergite, sd_chalcedony, sd_verm_mg
  CHARACTER(len=25) :: sd_ferrihydrite, sd_lawsonite, sd_merwinite, sd_monticellite, sd_natrolite
  CHARACTER(len=25) :: sd_talc, sd_smectite_low, sd_prehnite, sd_chlorite, sd_rankinite
  CHARACTER(len=25) :: sd_scolecite, sd_tobermorite_9a, sd_tremolite, sd_chamosite7a, sd_clinochlore14a
  CHARACTER(len=25) :: sd_clinochlore7a, sd_andradite
  CHARACTER(len=25) :: sd_saponite_ca, sd_troilite, sd_pyrrhotite, sd_lepidocrocite, sd_daphnite_7a
  CHARACTER(len=25) :: sd_daphnite_14a, sd_verm_k, sd_greenalite, sd_aragonite
  CHARACTER(len=25) :: sd_siderite, sd_kaolinite, sd_goethite, sd_dolomite, sd_celadonite ! secondary
  CHARACTER(len=25) :: sd_sio2, sd_albite, sd_calcite, sd_mont_na, sd_smectite, sd_saponite ! secondary
  CHARACTER(len=25) :: sd_stilbite, sd_saponite_k, sd_anhydrite, sd_clinoptilolite, sd_pyrite ! secondary
  CHARACTER(len=25) :: sd_quartz, sd_kspar, sd_saponite_na, sd_nont_na, sd_nont_mg, sd_nont_k ! secondary
  CHARACTER(len=25) :: sd_fe_celadonite, sd_nont_ca, sd_muscovite, sd_mesolite, sd_hematite, sd_diaspore !
  CHARACTER(len=25) :: sd_laumontite, sd_mont_k, sd_mont_mg, sd_mont_ca
  CHARACTER(len=25) :: sd_fe_saponite_ca, sd_fe_saponite_mg


  ! NEWER GEOCHEM
  ! UPDATED BASALT_BOX CHEM STUFF 03/26/17
  CHARACTER(len=300) :: param_tra_string, param_xb_string, param_exp_string, param_exp1_string, param_ol_string, param_pyr_string, param_plag_string
  REAL(4) :: param_temp, param_tra, param_xb, param_exp, param_exp1
  CHARACTER(len=25) :: ol_k1, ol_e1, ol_n1, ol_k2, ol_e2, ol_k3, ol_e3, ol_n3
  CHARACTER(len=25) :: pyr_k1, pyr_e1, pyr_n1, pyr_k2, pyr_e2, pyr_k3, pyr_e3, pyr_n3
  CHARACTER(len=25) :: plag_k1, plag_e1, plag_n1, plag_k2, plag_e2, plag_k3, plag_e3, plag_n3
  CHARACTER(len=25) :: exp_ol1, exp_ol2, exp_ol3, exp_ol
  CHARACTER(len=25) :: exp_pyr1, exp_pyr2, exp_pyr3, exp_pyr
  CHARACTER(len=25) :: exp_plag1, exp_plag2, exp_plag3, exp_plag




  CHARACTER(len=60000) :: L5

  param_tra = 10.0**11
  param_xb = 10.0**(-2.0)


  !path = "/home/navah/input/"
  path2 = ""


  !L5 = " "
  L5 = "#  $Id: llnl.dat 4023 2010-02-09 21:02:42Z dlpark $" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"LLNL_AQUEOUS_MODEL_PARAMETERS" //NEW_LINE('')// &
       &"-temperatures" //NEW_LINE('')// &
       &"         0.0100   25.0000   60.0000  100.0000" //NEW_LINE('')// &
       &"       150.0000  200.0000  250.0000  300.0000" //NEW_LINE('')// &
     !  &"#debye huckel a (adh)" //NEW_LINE('')// &
       &"-dh_a" //NEW_LINE('')// &
       &"         0.4939    0.5114    0.5465    0.5995" //NEW_LINE('')// &
       &"         0.6855    0.7994    0.9593    1.2180" //NEW_LINE('')// &
     !  &"#debye huckel b (bdh)" //NEW_LINE('')// &
       &"-dh_b" //NEW_LINE('')// &
       &"         0.3253    0.3288    0.3346    0.3421" //NEW_LINE('')// &
       &"         0.3525    0.3639    0.3766    0.3925" //NEW_LINE('')// &
       &"-bdot" //NEW_LINE('')// &
       &"         0.0374    0.0410    0.0438    0.0460" //NEW_LINE('')// &
       &"         0.0470    0.0470    0.0340    0.0000" //NEW_LINE('')// &
     !  &"#cco2   (coefficients for the Drummond (1981) polynomial)" //NEW_LINE('')// &
       &"-co2_coefs" //NEW_LINE('')// &
       &"        -1.0312              0.0012806" //NEW_LINE('')// &
       &"          255.9                 0.4445" //NEW_LINE('')// &
       &"      -0.001606" //NEW_LINE('')// &
       &"NAMED_EXPRESSIONS" //NEW_LINE('')// &
     !  &"#" //NEW_LINE('')// &
     !  &"# formation of O2 from H2O " //NEW_LINE('')// &
     !  &"# 2H2O =  O2 + 4H+ + 4e-  " //NEW_LINE('')// &
     !  &"#" //NEW_LINE('')// &
       &"	Log_K_O2" //NEW_LINE('')// &
       &"	 	log_k      -85.9951" //NEW_LINE('')// &
       &"		-delta_H	559.543	kJ/mol	# 	O2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-2.9 kcal/mol" //NEW_LINE('')// &
       &"	        -analytic   38.0229    7.99407E-03   -2.7655e+004  -1.4506e+001  199838.45" //NEW_LINE('')// &
     !  &"#	Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"SOLUTION_MASTER_SPECIES" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
     !  &"#element species        alk     gfw_formula     element_gfw" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Al       Al+3           0.0     Al              26.9815" //NEW_LINE('')// &
                                !&"Alkalinity HCO3-        1.0     Ca0.5(CO3)0.5   50.05" //NEW_LINE('')// &
       &"Alkalinity	CO3-2	1.0			Ca0.5(CO3)0.5	50.05" //NEW_LINE('')// &
       ! &"C(-4)	CH4		0.0	CH4" //NEW_LINE('')// &
       ! &"C(-3)	C2H6		0.0	C2H6" //NEW_LINE('')// &
       ! &"C(-2)	C2H4		0.0	C2H4" //NEW_LINE('')// &
       &"C        HCO3-          1.0     HCO3            12.0110" //NEW_LINE('')// &
       !&"C        CO3-2          1.0     HCO3            12.0110" //NEW_LINE('')// &
       ! &"C(+2)	 CO		0	C" //NEW_LINE('')// &
       &"C(+4)    HCO3-          1.0     HCO3" //NEW_LINE('')// &
       !&"C(+4)    CO3-2          1.0     HCO3" //NEW_LINE('')// &
       &"Ca       Ca+2           0.0     Ca              40.078" //NEW_LINE('')// &
       &"Cl       Cl-            0.0     Cl              35.4527" //NEW_LINE('')// &
       &"Cl(-1)	 Cl-		0	Cl" //NEW_LINE('')// &
       &"Cl(1)	 ClO-		0	Cl" //NEW_LINE('')// &
       &"Cl(3)	 ClO2-		0	Cl" //NEW_LINE('')// &
       &"Cl(5)	 ClO3-		0	Cl" //NEW_LINE('')// &
       &"Cl(7)	 ClO4-		0	Cl" //NEW_LINE('')// &
       &"E        e-             0.0     0.0             0.0" //NEW_LINE('')// &
       &"Fe       Fe+2           0.0     Fe              55.847" //NEW_LINE('')// &
       &"Fe(+2)   Fe+2           0.0     Fe" //NEW_LINE('')// &
       &"Fe(+3)   Fe+3           -2.0    Fe" //NEW_LINE('')// &
       &"H        H+             -1.     H               1.0079" //NEW_LINE('')// &
       &"H(0)     H2             0.0     H" //NEW_LINE('')// &
       &"H(+1)    H+             -1.     0.0" //NEW_LINE('')// &
       &"K        K+             0.0     K               39.0983" //NEW_LINE('')// &
       &"Mg       Mg+2           0.0     Mg              24.305" //NEW_LINE('')// &
       &"Na       Na+            0.0     Na              22.9898" //NEW_LINE('')// &
       &"O        H2O            0.0     O               15.994" //NEW_LINE('')// &
       &"O(-2)    H2O            0.0     0.0" //NEW_LINE('')// &
       &"O(0)     O2             0.0     O" //NEW_LINE('')// &
       &"S	 SO4-2          0.0     SO4             32.066" //NEW_LINE('')// &
       &"S(-2)	 HS-            1.0     S" //NEW_LINE('')// &
       &"S(+2)	 S2O3-2		0.0	S" //NEW_LINE('')// &
       &"S(+3)	 S2O4-2		0.0	S" //NEW_LINE('')// &
       &"S(+4)	 SO3-2		0.0	S" //NEW_LINE('')// &
       &"S(+5)	 S2O5-2		0.0	S" //NEW_LINE('')// &
       &"S(+6)	 SO4-2          0.0     SO4" //NEW_LINE('')// &
       &"S(+7)	 S2O8-2		0.0	S" //NEW_LINE('')// &
       &"S(+8)	 HSO5-		0.0	S" //NEW_LINE('')// &
       &"Si       SiO2         0.0     SiO2            28.0855" //NEW_LINE('')// &
       !&"Si		H4SiO4	0.0	SiO2		28.0843" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       !-DB SOLUTION SPECIES
       &"SOLUTION_SPECIES" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &

       &"H2O + 0.01e- = H2O-0.01" //NEW_LINE('')// &
       &"	log_k -9" //NEW_LINE('')// &

       &"Al+3 =  Al+3 " //NEW_LINE('')// &
       &"	-llnl_gamma	9.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Al+3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-128.681 kcal/mol" //NEW_LINE('')// &
       &"-Vm  -3.3404  -17.1108  14.9917  -2.0716  2.8711 9 # supcrt" //NEW_LINE('')// &


       &"Ca+2 =  Ca+2 " //NEW_LINE('')// &
       &"	-llnl_gamma	6.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Ca+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-129.8 kcal/mol" //NEW_LINE('')// &
       !&"	-millero -19.69 0.1058 -0.001256 1.617 -0.075 0.0008262" //NEW_LINE('')// &
       &"-Vm  -0.3456  -7.252  6.149  -2.479  1.239  5  1.60  -57.1  -6.12e-3  1 # supcrt modified" //NEW_LINE('')// &


       &"Cl- =  Cl- " //NEW_LINE('')// &
       &"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Cl-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-39.933 kcal/mol" //NEW_LINE('')// &
       !&"	-millero 16.37 0.0896 -0.001264 -1.494 0.034 -0.000621" //NEW_LINE('')// &
       &"-Vm  4.465  4.801  4.325  -2.847  1.748  0  -0.331  20.16  0  1 # supcrt modified" //NEW_LINE('')// &


       &"e- =  e- " //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol		e-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kJ/mol" //NEW_LINE('')// &

       &"Fe+2 =  Fe+2 " //NEW_LINE('')// &
       &"	-llnl_gamma	6.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Fe+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-22.05 kcal/mol" //NEW_LINE('')// &
       &"-Vm  -0.3255  -9.687  1.536  -2.379  0.3033  5.5  -4.21e-2  37.96  0  1 # supcrt modified" //NEW_LINE('')// &


       &"H+ =  H+ " //NEW_LINE('')// &
       &"	-llnl_gamma	9.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	H+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kJ/mol" //NEW_LINE('')// &

       &"HCO3- =  HCO3- " //NEW_LINE('')// &
       &"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	HCO3-" //NEW_LINE('')// &
       !&"-Vm  5.224  0  0  -5.85  5.126  0  0.404  110  -5.74e-3  1 # supcrt modified" //NEW_LINE('')// &
       ! not entirely sure about carbonate vs. bicarbonate decision here...

       ! &"CO3-2 =  CO3-2 " //NEW_LINE('')// &
       ! &"	-llnl_gamma	5.4	" //NEW_LINE('')// &
       ! &"	log_k 0" //NEW_LINE('')// &
       ! &"	-delta_H	0	kJ/mol	# 	CO3-2" //NEW_LINE('')// &
       ! &"-Vm  5.224  0  0  -5.85  5.126  0  0.404  110  -5.74e-3  1 # supcrt modified" //NEW_LINE('')// &

       ! &"#	Enthalpy of formation:	-164.898 kcal/mol" //NEW_LINE('')// &
       &"K+ =  K+ " //NEW_LINE('')// &
       &"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	K+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-60.27 kcal/mol" //NEW_LINE('')// &
       !&"	-millero 7.26 0.0892 -0.000736 2.722 -0.101 0.00151" //NEW_LINE('')// &
       &"-Vm  3.322  -1.473  6.534  -2.712  9.06e-2  3.5  0  29.7  0  1 # supcrt modified" //NEW_LINE('')// &

       &"Mg+2 =  Mg+2 " //NEW_LINE('')// &
       &"	-llnl_gamma	8.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Mg+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-111.367 kcal/mol" //NEW_LINE('')// &
       !&"	-millero -22.32 0.0868 -0.0016 2.017 -0.125 0.001457" //NEW_LINE('')// &
       &"-Vm  -1.410  -8.6  11.13  -2.39  1.332  5.5  1.29  -32.9  -5.86e-3  1 # supcrt modified" //NEW_LINE('')// &

       &"Na+ =  Na+ " //NEW_LINE('')// &
       &"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	Na+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-57.433 kcal/mol" //NEW_LINE('')// &
       !&"	-millero -3.46 0.1092 -0.000768 2.698 -0.106 0.001651" //NEW_LINE('')// &
       &"-Vm  1.403  -2.285  4.419  -2.726  -5.125e-5  4.0  0.162  47.67  -3.09e-3  0.725 # sup" //NEW_LINE('')// &

       &"H2O =  H2O " //NEW_LINE('')// &
       &"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
       &"        log_k   0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	H2O" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-68.317 kcal/mol" //NEW_LINE('')// &
       &"SO4-2 =  SO4-2 " //NEW_LINE('')// &
       &"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	SO4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-217.4 kcal/mol" //NEW_LINE('')// &
       !&"	-millero 9.26 0.284 -0.003808 0.4348 -0.0099143 -8.4762e-05" //NEW_LINE('')// &
       &"-Vm  8.0  2.51  -42.5 5.41  4.23 0 0 0 0 1 # supcrt modified" //NEW_LINE('')// &


       &"SiO2 =  SiO2 " //NEW_LINE('')// &
       &"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
       &"	log_k 0" //NEW_LINE('')// &
       &"	-delta_H	0	kJ/mol	# 	SiO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-209.775 kcal/mol" //NEW_LINE('')// &
       !&"-Vm  10.5  1.7  20  -2.7  0.1291 # supcrt + 2*H2O in a1" //NEW_LINE('')// &
       ! not sure of ion vs. sio2 here

       &"2H2O =  O2 + 4H+ + 4e-  " //NEW_LINE('')// &
       &"	-CO2_llnl_gamma" //NEW_LINE('')// &
       &" 	log_k      -85.9951" //NEW_LINE('')// &
       &"	-delta_H	559.543	kJ/mol	# 	O2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-2.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic   38.0229    7.99407E-03   -2.7655e+004  -1.4506e+001  199838.45" //NEW_LINE('')// &
     !  &"#	Range:  0-300" //NEW_LINE('')// &
       &"-Vm  5.7889  6.3536  3.2528  -3.0417  -0.3943 # supcrt" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.0000 SO4-- + 1.0000 H+  =  HS- +2.0000 O2  " //NEW_LINE('')// &
       &"        -llnl_gamma           3.5    " //NEW_LINE('')// &
       &"        log_k           -138.3169" //NEW_LINE('')// &
       &"	-delta_H	869.226	kJ/mol	# 	HS-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-3.85 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 2.6251e+001 3.9525e-002 -4.5443e+004 -1.1107e+001 3.1843e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       !&"-Vm 8.2 9.2590  2.1108  -3.1618 1.1748  0 -0.3 15 0 1 # supcrt modified" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" .5000 O2 + 2.0000 HS-  = S2--  + H2O" //NEW_LINE('')// &
     !  &"#2 HS- = S2-- +2 H+ + 2e-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           33.2673" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.21730E+02   -0.12307E-02    0.10098E+05   -0.88813E+01    0.15757E+03" //NEW_LINE('')// &
       &"	-mass_balance	S(-2)2" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
     !  &"#	-add_logk	Log_K_O2	0.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"2.0000 H+  + 2.0000 SO3--  = S2O3--  + O2  + H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -40.2906" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O3-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic  0.77679E+02    0.65761E-01   -0.15438E+05   -0.34651E+02   -0.24092E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       ! ! new
       ! &"CO3-2 + H+ = HCO3-" //NEW_LINE('')// &
       ! &" 	-log_k	10.329" //NEW_LINE('')// &
       ! &"	-delta_H	-3.561	kcal" //NEW_LINE('')// &
       ! &"  -analytic	107.8871	0.03252849	-5151.79	-38.92561	563713.9" //NEW_LINE('')// &
       ! &"-Vm  8.615  0  -12.21 0  1.667  0  0  264  0  1 # supcrt modified" //NEW_LINE('')// &
       !
       ! ! new
       ! &"CO3-2 + 2 H+ = CO2 + H2O" //NEW_LINE('')// &
       ! &" 	-log_k	16.681" //NEW_LINE('')// &
       ! &"	-delta_h -5.738	kcal" //NEW_LINE('')// &
       ! &"  -analytic	464.1965	0.09344813	-26986.16	-165.75951	2248628.9" //NEW_LINE('')// &
       ! &"-Vm  21.78  -49.4  -91.7  31.96 # supcrt modified" //NEW_LINE('')// &

       ! &" H+  + HCO3-  + H2O  = CH4 + 2.0000 O2" //NEW_LINE('')// &
       ! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       ! &"        log_k            -144.1412" //NEW_LINE('')// &
       ! &"	-delta_H	863.599	kJ/mol	# 	CH4" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-21.01 kcal/mol" //NEW_LINE('')// &
       ! &"	-analytic    -0.41698E+02    0.36584E-01   -0.40675E+05    0.93479E+01   -0.63468E+03" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       !&"-Vm 7.7" //NEW_LINE('')// &
       ! unsure but unimportant


       ! &" CO3-2 + 10 H+ + 8 e- = CH4 + 3 H2O" //NEW_LINE('')// &
       ! &"        -log_k	41.071" //NEW_LINE('')// &
       ! &"	-delta_h -61.039 kcal" //NEW_LINE('')// &
       ! &"	-Vm 7.7" //NEW_LINE('')// &

       ! &"" //NEW_LINE('')// &
       ! &" 2.0000 H+  + 2.0000 HCO3-  + H2O  = C2H6  + 3.5000 O2" //NEW_LINE('')// &
       ! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       ! &"        log_k            -228.6072" //NEW_LINE('')// &
       ! &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	C2H6" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic    0.10777E+02    0.72105E-01   -0.67489E+05   -0.13915E+02   -0.10531E+04" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &

       ! &"" //NEW_LINE('')// &
       ! &" 2.000 H+  + 2.0000 HCO3-  = C2H4 + 3.0000 O2" //NEW_LINE('')// &
       ! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       ! &"        log_k            -254.5034" //NEW_LINE('')// &
       ! &"	-delta_H	1446.6	kJ/mol	# 	C2H4" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	24.65 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic    -0.30329E+02    0.71187E-01   -0.73140E+05    0.00000E+00    0.00000E+00" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       ! &" 1.0000 HCO3- + 1.0000 H+  =  CO +1.0000 H2O +0.5000 O2 " //NEW_LINE('')// &
       ! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       ! &"        log_k           -41.7002" //NEW_LINE('')// &
       ! &"	-delta_H	277.069	kJ/mol	# 	CO" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-28.91 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.0028e+002 4.6877e-002 -1.8062e+004 -4.0263e+001 3.8031e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &

       &" 1.0000 Cl- + 0.5000 O2  =  ClO-   " //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -15.1014" //NEW_LINE('')// &
       &"	-delta_H	66.0361	kJ/mol	# 	ClO-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-25.6 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 6.1314e+001 3.4812e-003 -6.0952e+003 -2.3043e+001 -9.5128e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.0000 O2 + 1.0000 Cl-  =  ClO2-   " //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -23.108" //NEW_LINE('')// &
       &"	-delta_H	112.688	kJ/mol	# 	ClO2-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-15.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 3.3638e+000 -6.1675e-003 -4.9726e+003 -2.0467e+000 -2.5769e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.5000 O2 + 1.0000 Cl-  =  ClO3-   " //NEW_LINE('')// &
       &"        -llnl_gamma           3.5    " //NEW_LINE('')// &
       &"        log_k           -17.2608" //NEW_LINE('')// &
       &"	-delta_H	81.3077	kJ/mol	# 	ClO3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-24.85 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 2.8852e+001 -4.8281e-003 -4.6779e+003 -1.0772e+001 -2.0783e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 2.0000 O2 + 1.0000 Cl-  =  ClO4-   " //NEW_LINE('')// &
       &"        -llnl_gamma           3.5    " //NEW_LINE('')// &
       &"        log_k           -15.7091" //NEW_LINE('')// &
       &"	-delta_H	62.0194	kJ/mol	# 	ClO4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-30.91 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 7.0280e+001 -6.8927e-005 -5.5690e+003 -2.6446e+001 -1.6596e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &" 1.0000 H+ + 1.0000 Fe++ + 0.2500 O2  =  Fe+++ +0.5000 H2O " //NEW_LINE('')// &
       &"        -llnl_gamma           9.0    " //NEW_LINE('')// &
       &"        log_k           +8.4899" //NEW_LINE('')// &
       &"	-delta_H	-97.209	kJ/mol	# 	Fe+3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-11.85 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.7808e+001 -1.1753e-002 4.7609e+003 5.5866e+000 7.4295e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.0000 H2O  =  H2 +0.5000 O2   " //NEW_LINE('')// &
       &"	-CO2_llnl_gamma" //NEW_LINE('')// &
       &"        log_k           -46.1066" //NEW_LINE('')// &
       &"	-delta_H	275.588	kJ/mol	# 	H2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 6.6835e+001 1.7172e-002 -1.8849e+004 -2.4092e+001 4.2501e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.0000 SO4-- + 1.0000 H+ + 0.5000 O2  =  HSO5-  " //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -17.2865" //NEW_LINE('')// &
       &"	-delta_H	140.038	kJ/mol	# 	HSO5-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-185.38 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 5.9944e+001 3.0904e-002 -7.7494e+003 -2.4420e+001 -1.2094e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &" 2.0000 H+  + 2.0000 SO3--  = S2O4--  + .500 O2  + H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           5.0    " //NEW_LINE('')// &
     !  &"#        log_k           -25.2075" //NEW_LINE('')// &
       &"        log_k           -25.2076" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
     !  &"#        -analytic  -0.15158E+05   -0.31356E+01    0.47072E+06    0.58544E+04    0.73497E+04" //NEW_LINE('')// &
       &"	-analytic	-2.3172e2	2.0393e-3	-7.1011e0	8.3239e1	9.4155e-1" //NEW_LINE('')// &
     !  &"#	changed 3/23/04, corrected to supcrt temperature dependence, GMA" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
     !  &"# 2.0000 SO3--  + .500 O2  + 2.0000 H+  = S2O6--  + H2O" //NEW_LINE('')// &
     !  &"#  H2O = .5 O2 + 2H+ + 2e- " //NEW_LINE('')// &
       &"2SO3-- = S2O6-- + 2e-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           41.8289" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O6-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.14458E+03    0.61449E-01    0.71877E+04   -0.58657E+02    0.11211E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"	-add_logk  Log_K_O2	0.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 2.0000 SO3--  + 1.500 O2  + 2.0000 H+  = S2O8--  + H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           70.7489" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O8-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.18394E+03    0.60414E-01    0.13864E+05   -0.71804E+02    0.21628E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"O2 + H+ + 3.0000 HS-  = S3--  + 2.0000 H2O" //NEW_LINE('')// &
     !  &"# 2H2O = O2 + 4H+ + 4e-" //NEW_LINE('')// &
     !  &"#3HS- = S3-- + 3H+ + 4e-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           79.3915" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S3-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -0.51626E+02    0.70208E-02    0.31797E+05    0.11927E+02   -0.64249E+06" //NEW_LINE('')// &
       &"	-mass_balance	S(-2)3" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
     !  &"#	-add_logk  Log_K_O2	1.0" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
     !  &"# 3.0000 SO3--  + 4.0000 H+  = S3O6-- + .500 O2 + 2.0000 H2O" //NEW_LINE('')// &
     !  &"# .5 O2 + 2H+ + 2e- = H2O" //NEW_LINE('')// &
       &"3SO3-- + 6 H+ + 2e- = S3O6-- + 3H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -6.2316" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S3O6-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.23664E+03    0.12702E+00   -0.10110E+05   -0.99715E+02   -0.15783E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"	-add_logk	Log_K_O2	-0.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.5000 O2 + 2.0000 H+ + 4.0000 HS-  = S4--  + 3.0000 H2O" //NEW_LINE('')// &
     !  &"#4 HS- = S4-- + 4H+ + 6e-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           125.2958" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.20875E+03    0.58133E-01    0.33278E+05   -0.85833E+02    0.51921E+03" //NEW_LINE('')// &
       &"	-mass_balance	S(-2)4" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
     !  &"#	-add_logk	Log_K_O2	1.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
     !  &"# 4.0000 SO3-- + 6.0000 H+  = S4O6-- + 1.500 O2 + 3.0000 H2O" //NEW_LINE('')// &
       &"4 SO3-- + 12 H+ + 6e- = S4O6-- + 6H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -38.3859" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S4O6-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.32239E+03    0.19555E+00   -0.23617E+05   -0.13729E+03   -0.36862E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"	-add_logk	Log_K_O2	-1.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"2.0000 O2 + 3.0000 H+  + 5.0000 HS-  = S5--  + 4.0000 H2O" //NEW_LINE('')// &
     !  &"#5 HS- = S5-- + 5H+ + 8e-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           170.9802" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S5-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.30329E+03    0.88033E-01    0.44739E+05   -0.12471E+03    0.69803E+03" //NEW_LINE('')// &
       &"	-mass_balance	S(-2)5" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
     !  &"#	-add_logk	Log_K_O2	2" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
     !  &"# 5.0000 SO3-- + 8.0000 H+  = S5O6-- + 2.5000 O2 + 4.0000 H2O" //NEW_LINE('')// &
     !  &"# 2.5O2 + 10 H+ + 10e- = 5H2O" //NEW_LINE('')// &
       &"5SO3-- + 18H+ + 10e- = S5O6-- + 9H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -99.4206" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S5O6-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.42074E+03    0.25833E+00   -0.43878E+05   -0.18178E+03   -0.68480E+03" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"	-add_logk	Log_K_O2	-2.5" //NEW_LINE('')// &

       ! &"" //NEW_LINE('')// &
       ! &"# 1.0000 H+  + HCO3-  + HS-  + NH3 = SCN-  + 3.0000 H2O" //NEW_LINE('')// &
       ! &"#        -llnl_gamma           3.5    " //NEW_LINE('')// &
       ! &"#        log_k            3.0070" //NEW_LINE('')// &
       ! &"#	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	SCN-" //NEW_LINE('')// &
       ! &"##	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       ! &"#        -analytic  0.16539E+03    0.49623E-01   -0.44624E+04   -0.65544E+02   -0.69680E+02" //NEW_LINE('')// &
       ! &"##       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &" 1.0000 SO4--  =  SO3-- +0.5000 O2   " //NEW_LINE('')// &
       &"        -llnl_gamma           4.5    " //NEW_LINE('')// &
       &"        log_k           -46.6244" //NEW_LINE('')// &
       &"	-delta_H	267.985	kJ/mol	# 	SO3-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-151.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.3771e+001 6.5102e-004 -1.3330e+004 4.7164e+000 -2.0800e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 Al+++  =  Al(OH)2+ +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -10.5945" //NEW_LINE('')// &
       &"	-delta_H	98.2822	kJ/mol	# 	Al(OH)2+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-241.825 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 4.4036e+001 2.0168e-002 -5.5455e+003 -1.6987e+001 -8.6545e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &


       &"" //NEW_LINE('')// &
       &"2.0000 SO4-- + 1.0000 Al+++  =  Al(SO4)2-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +4.9000" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al(SO4)2-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

       &" " //NEW_LINE('')// &
       &"28.0000 H2O + 13.0000 Al+++  =  Al13O4(OH)24+7 +32.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           6.0    " //NEW_LINE('')// &
       &"        log_k           -98.73" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al13O4(OH)24+7" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 2.0000 Al+++  =  Al2(OH)2++++ +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           5.5    " //NEW_LINE('')// &
       &"        log_k           -7.6902" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al2(OH)2+4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

       &" " //NEW_LINE('')// &
       &"4.0000 H2O + 3.0000 Al+++  =  Al3(OH)4+5 +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           6.0    " //NEW_LINE('')// &
       &"        log_k           -13.8803" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al3(OH)4+5" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 Al+++  =  AlO2- +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -22.8833" //NEW_LINE('')// &
       &"	-delta_H	180.899	kJ/mol	# 	AlO2-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-222.079 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.0803e+001 -3.4379e-003 -9.7391e+003 0.0000e+000 0.0000e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.0000 H2O + 1.0000 Al+++  =  AlOH++ +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.5    " //NEW_LINE('')// &
       &"        log_k           -4.9571" //NEW_LINE('')// &
       &"	-delta_H	49.798	kJ/mol	# 	AlOH+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-185.096 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.6224e-001 8.8816e-003 -1.8686e+003 -4.3195e-001 -2.9158e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 Al+++  =  AlSO4+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +3.0100" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	AlSO4+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

       &" " //NEW_LINE('')// &
       &"1.0000 HCO3- + 1.0000 H+  =  CO2 +1.0000 H2O" //NEW_LINE('')// &
       &"        -CO2_llnl_gamma" //NEW_LINE('')// &
       &"        log_k           +6.3447" //NEW_LINE('')// &
       &"	-delta_H	-9.7027	kJ/mol	# 	CO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-98.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.0534e+001 2.1746e-002 2.5216e+003 7.9125e-001 3.9351e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       ! &"CO3-2 + 2 H+ = CO2 + H2O" //NEW_LINE('')// &
       ! &"        -log_k	16.681" //NEW_LINE('')// &
       ! &"	-delta_h -5.738	kcal" //NEW_LINE('')// &
       ! &"        -analytic	464.1965	0.09344813	-26986.16	-165.75951	2248628.9" //NEW_LINE('')// &
       ! &"        -Vm  21.78  -49.4  -91.7  31.96 # supcrt modified" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.0000 HCO3-  =  CO3-- +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.5    " //NEW_LINE('')// &
       &"        log_k           -10.3288" //NEW_LINE('')// &
       &"	-delta_H	14.6984	kJ/mol	# 	CO3-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-161.385 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -6.9958e+001 -3.3526e-002 -7.0846e+001 2.8224e+001 -1.0849e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"CO3-2 + H+ = HCO3-" //NEW_LINE('')// &
       &"        -llnl_gamma           5.4    " //NEW_LINE('')// &
       &"        log_k           10.3288" //NEW_LINE('')// &
       &"	-delta_h -3.561	kcal" //NEW_LINE('')// &
       &"        -analytic	107.8871	0.03252849	-5151.79	-38.92561	563713.9" //NEW_LINE('')// &
     !  &"#       -Vm  8.615  0  -12.21 0  1.667  0  0  264  0  1 # supcrt modified" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.0000 HCO3- + 1.0000 Ca++  =  CaCO3 +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -7.0017" //NEW_LINE('')// &
       &"	-delta_H	30.5767	kJ/mol	# 	CaCO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-287.39 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 2.3045e+002 5.5350e-002 -8.5056e+003 -9.1096e+001 -1.3279e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       !&"-Vm  -.2430  -8.3748  9.0417  -2.4328  -.0300 # supcrt" //NEW_LINE('')// &
       ! again

       ! &"" //NEW_LINE('')// &
       ! &"Ca+2 + CO3-2 = CaCO3" //NEW_LINE('')// &
       ! &"        -log_k	3.224" //NEW_LINE('')// &
       ! &"	-delta_h 3.545	kcal" //NEW_LINE('')// &
       ! &"        -analytic	-1228.732	-0.299440	35512.75	485.818" //NEW_LINE('')// &
       ! &"#       -Vm  -.2430  -8.3748  9.0417  -2.4328  -.0300 # supcrt" //NEW_LINE('')// &
       !

       &"" //NEW_LINE('')// &
       &"1.0000 Cl- + 1.0000 Ca++  =  CaCl+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -0.6956" //NEW_LINE('')// &
       &"	-delta_H	2.02087	kJ/mol	# 	CaCl+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-169.25 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 8.1498e+001 3.8387e-002 -1.3763e+003 -3.5968e+001 -2.1501e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 Cl- + 1.0000 Ca++  =  CaCl2" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -0.6436" //NEW_LINE('')// &
       &"	-delta_H	-5.8325	kJ/mol	# 	CaCl2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-211.06 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.8178e+002 7.6910e-002 -3.1088e+003 -7.8760e+001 -4.8563e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"1.0000 HCO3- + 1.0000 Ca++  =  CaHCO3+" //NEW_LINE('')// &
       ! &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       ! &"        log_k           +1.0467" //NEW_LINE('')// &
       ! &"	-delta_H	1.45603	kJ/mol	# 	CaHCO3+" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-294.35 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 5.5985e+001 3.4639e-002 -3.6972e+002 -2.5864e+001 -5.7859e+000" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"-Vm  3.1911  .0104  5.7459  -2.7794  .3084 5.4 # supcrt" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"1.0000 H2O + 1.0000 Ca++  =  CaOH+ +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -12.85" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	CaOH+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &

     !  &"#Al(OH)4-            82" //NEW_LINE('')// &
       &"        Al+3 + 4H2O = Al(OH)4- + 4H+ " //NEW_LINE('')// &
       &"        log_k           -22.7" //NEW_LINE('')// &
       &"        delta_h 42.3 kcal" //NEW_LINE('')// &
       &"        -analytical     51.578          0.0     -11168.9        -14.865         0.0" //NEW_LINE('')// &



       &"1.0000 SO4-- + 1.0000 Ca++  =  CaSO4" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +2.1111" //NEW_LINE('')// &
       &"	-delta_H	5.4392	kJ/mol	# 	CaSO4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-345.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 2.8618e+002 8.4084e-002 -7.6880e+003 -1.1449e+002 -1.2005e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"-Vm  2.7910  -.9666  6.1300  -2.7390  -.0010 # supcrt" //NEW_LINE('')// &

       &"2.0000 H2O + 1.0000 Fe++  =  Fe(OH)2 +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -20.6" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 Fe+++  =  Fe(OH)2+ +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -5.67" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)2+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"3.0000 H2O + 1.0000 Fe+++  =  Fe(OH)3 +3.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -12" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"3.0000 H2O + 1.0000 Fe++  =  Fe(OH)3- +3.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -31" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"4.0000 H2O + 1.0000 Fe+++  =  Fe(OH)4- +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -21.6" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"4.0000 H2O + 1.0000 Fe++  =  Fe(OH)4-- +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -46" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 SO4-- + 1.0000 Fe+++  =  Fe(SO4)2-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +3.2137" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(SO4)2-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 2.0000 Fe+++  =  Fe2(OH)2++++ +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           5.5    " //NEW_LINE('')// &
       &"        log_k           -2.95" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe2(OH)2+4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"4.0000 H2O + 3.0000 Fe+++  =  Fe3(OH)4+5 +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           6.0    " //NEW_LINE('')// &
       &"        log_k           -6.3" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe3(OH)4+5" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 HCO3- + 1.0000 Fe++  =  FeCO3 +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -5.5988" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 HCO3- + 1.0000 Fe+++  =  FeCO3+ +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -0.6088" //NEW_LINE('')// &
       &"	-delta_H	-50.208	kJ/mol	# 	FeCO3+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-188.748 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.7100e+002 8.0413e-002 -4.3217e+002 -7.8449e+001 -6.7948e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Fe++ + 1.0000 Cl-  =  FeCl+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -0.1605" //NEW_LINE('')// &
       &"	-delta_H	3.02503	kJ/mol	# 	FeCl+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-61.26 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 8.2435e+001 3.7755e-002 -1.4765e+003 -3.5918e+001 -2.3064e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Fe+++ + 1.0000 Cl-  =  FeCl++" //NEW_LINE('')// &
       &"        -llnl_gamma           4.5    " //NEW_LINE('')// &
       &"        log_k           -0.8108" //NEW_LINE('')// &
       &"	-delta_H	36.6421	kJ/mol	# 	FeCl+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-180.018 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 1.6186e+002 5.9436e-002 -5.1913e+003 -6.5852e+001 -8.1053e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 Cl- + 1.0000 Fe++  =  FeCl2" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -2.4541" //NEW_LINE('')// &
       &"	-delta_H	6.46846	kJ/mol	# 	FeCl2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-100.37 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.9171e+002 7.8070e-002 -4.1048e+003 -8.2292e+001 -6.4108e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 Cl- + 1.0000 Fe+++  =  FeCl2+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +2.1300" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl2+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"4.0000 Cl- + 1.0000 Fe+++  =  FeCl4-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -0.79" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"4.0000 Cl- + 1.0000 Fe++  =  FeCl4--" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -1.9" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.4108e+002 -6.0086e-003 9.7979e+003 8.4084e+001 1.5296e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 HCO3- + 1.0000 Fe++  =  FeHCO3+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +2.7200" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeHCO3+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 H2O + 1.0000 Fe++  =  FeOH+ +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -9.5" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeOH+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 H2O + 1.0000 Fe+++  =  FeOH++ +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.5    " //NEW_LINE('')// &
       &"        log_k           -2.19" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeOH+2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 Fe++  =  FeSO4" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +2.2000" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeSO4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 Fe+++  =  FeSO4+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +1.9276" //NEW_LINE('')// &
       &"	-delta_H	27.181	kJ/mol	# 	FeSO4+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-932.001 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 2.5178e+002 1.0080e-001 -6.0977e+003 -1.0483e+002 -9.5223e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
     !  &"#1.0000 HS- + 1.0000 H+  =  H2S" //NEW_LINE('')// &
     !  &"#        -llnl_gamma           3.0    " //NEW_LINE('')// &
     !  &"#        log_k           +6.99" //NEW_LINE('')// &
     !  &"#        -analytic 1.2833e+002 5.1641e-002 -1.1681e+003 -5.3665e+001 -1.8266e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"# these (above) H2S values are from " //NEW_LINE('')// &
       ! &"# Suleimenov & Seward, Geochim. Cosmochim. Acta, v. 61, p. 5187-5198." //NEW_LINE('')// &
       ! &"# values below are the original Thermo.com.v8.r6.230 data from somewhere" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 HS- + 1.0000 H+  =  H2S" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +6.9877" //NEW_LINE('')// &
       &"	-delta_H	-21.5518	kJ/mol	# 	H2S" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-9.001 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 3.9283e+001 2.8727e-002  1.3477e+003 -1.8331e+001  2.1018e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 H+ + 1.0000 SO3--  =  H2SO3" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +9.2132" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 H+ + 1.0000 SO4--  =  H2SO4" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -1.0209" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SO4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 SiO2  =  H2SiO4-- +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -22.96" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SiO4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"8.0000 H2O + 4.0000 SiO2  =  H4(H2SiO4)4---- +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -35.94" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H4(H2SiO4)4-4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"8.0000 H2O + 4.0000 SiO2  =  H6(H2SiO4)4-- +2.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -13.64" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H6(H2SiO4)4-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 Al+++  =  HAlO2 +3.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -16.4329" //NEW_LINE('')// &
       &"	-delta_H	144.704	kJ/mol	# 	HAlO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-230.73 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 4.2012e+001 1.9980e-002 -7.7847e+003 -1.5470e+001 -1.2149e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 H+ + 1.0000 Cl-  =  HCl" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -0.67" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HCl" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 4.1893e+002 1.1103e-001 -1.1784e+004 -1.6697e+002 -1.8400e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 H+ + 1.0000 ClO-  =  HClO" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +7.5692" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HClO" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 H+ + 1.0000 ClO2-  =  HClO2" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +3.1698" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HClO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"  " //NEW_LINE('')// &
       &"1.0000 H+ + 1.0000 S2O3--  =  HS2O3-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k            1.0139" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HS2O3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SO3-- + 1.0000 H+  =  HSO3-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +7.2054" //NEW_LINE('')// &
       &"	-delta_H	9.33032	kJ/mol	# 	HSO3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-149.67 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 5.5899e+001 3.3623e-002 -5.0120e+002 -2.3040e+001 -7.8373e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 H+  =  HSO4-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +1.9791" //NEW_LINE('')// &
       &"	-delta_H	20.5016	kJ/mol	# 	HSO4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-212.5 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 4.9619e+001 3.0368e-002 -1.1558e+003 -2.1335e+001 -1.8051e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 SiO2 + 1.0000 H2O  =  HSiO3- +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -9.9525" //NEW_LINE('')// &
       &"	-delta_H	25.991	kJ/mol	# 	HSiO3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-271.88 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 6.4211e+001 -2.4872e-002 -1.2707e+004 -1.4681e+001 1.0853e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 K+ + 1.0000 Cl-  =  KCl" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -1.4946" //NEW_LINE('')// &
       &"	-delta_H	14.1963	kJ/mol	# 	KCl" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-96.81 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.3650e+002 3.8405e-002 -4.4014e+003 -5.4421e+001 -6.8721e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 K+ + 1.0000 H+  =  KHSO4" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +0.8136" //NEW_LINE('')// &
       &"	-delta_H	29.8319	kJ/mol	# 	KHSO4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-270.54 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.2620e+002 5.7349e-002 -3.3670e+003 -5.3003e+001 -5.2576e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 K+ + 1.0000 H2O  =  KOH +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -14.46" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	KOH" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 K+  =  KSO4-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +0.8796" //NEW_LINE('')// &
       &"	-delta_H	2.88696	kJ/mol	# 	KSO4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-276.98 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 9.9073e+001 3.7817e-002 -2.1628e+003 -4.1297e+001 -3.3779e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"4.0000 Mg++ + 4.0000 H2O  =  Mg4(OH)4++++ +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           5.5    " //NEW_LINE('')// &
       &"        log_k           -39.75" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Mg4(OH)4+4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 Mg++ + 1.0000 HCO3-  =  MgCO3 +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -7.3499" //NEW_LINE('')// &
       &"	-delta_H	23.8279	kJ/mol	# 	MgCO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-270.57 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 2.3465e+002 5.5538e-002 -8.3947e+003 -9.3104e+001 -1.3106e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Mg++ + 1.0000 Cl-  =  MgCl+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -0.1349" //NEW_LINE('')// &
       &"	-delta_H	-0.58576	kJ/mol	# 	MgCl+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-151.44 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 4.3363e+001 3.2858e-002 1.1878e+002 -2.1688e+001 1.8403e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Mg++ + 1.0000 HCO3-  =  MgHCO3+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +1.0357" //NEW_LINE('')// &
       &"	-delta_H	2.15476	kJ/mol	# 	MgHCO3+" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-275.75 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 3.8459e+001 3.0076e-002 9.8068e+001 -1.8869e+001 1.5187e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 Mg++  =  MgSO4" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +2.4117" //NEW_LINE('')// &
       &"	-delta_H	19.6051	kJ/mol	# 	MgSO4" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1355.96 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 1.7994e+002 6.4715e-002 -4.7314e+003 -7.3123e+001 -8.0408e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-200" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 H2O + 1.0000 Na+ + 1.0000 Al+++  =  NaAlO2 +4.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -23.6266" //NEW_LINE('')// &
       &"	-delta_H	190.326	kJ/mol	# 	NaAlO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-277.259 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.2288e+002 3.4921e-002 -1.2808e+004 -4.6046e+001 -1.9990e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Na+ + 1.0000 HCO3-  =  NaCO3- +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           -9.8144" //NEW_LINE('')// &
       &"	-delta_H	-5.6521	kJ/mol	# 	NaCO3-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-935.885 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 1.6939e+002 5.3122e-004 -7.6768e+003 -6.2078e+001 -1.1984e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Na+ + 1.0000 Cl-  =  NaCl" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -0.777" //NEW_LINE('')// &
       &"	-delta_H	5.21326	kJ/mol	# 	NaCl" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-96.12 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.1398e+002 3.6386e-002 -3.0847e+003 -4.6571e+001 -4.8167e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Na+ + 1.0000 HCO3-  =  NaHCO3" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +0.1541" //NEW_LINE('')// &
       &"	-delta_H	-13.7741	kJ/mol	# 	NaHCO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-944.007 kJ/mol" //NEW_LINE('')// &
       &"        -analytic -9.0668e+001 -2.9866e-002 2.7947e+003 3.6515e+001 4.7489e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-200" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 SiO2 + 1.0000 Na+ + 1.0000 H2O  =  NaHSiO3 +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -8.304" //NEW_LINE('')// &
       &"	-delta_H	11.6524	kJ/mol	# 	NaHSiO3" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-332.74 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 3.6045e+001 -9.0411e-003 -6.6605e+003 -1.0447e+001 5.8415e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 Na+ + 1.0000 H2O  =  NaOH +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           -14.7948" //NEW_LINE('')// &
       &"	-delta_H	53.6514	kJ/mol	# 	NaOH" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-112.927 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 8.7326e+001 2.3555e-002 -5.4770e+003 -3.6678e+001 -8.5489e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"1.0000 SO4-- + 1.0000 Na+  =  NaSO4-" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           +0.8200" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	NaSO4-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 H2O  =  OH- +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           3.5    " //NEW_LINE('')// &
       &"        log_k           -13.9951" //NEW_LINE('')// &
       &"	-delta_H	55.8146	kJ/mol	# 	OH-" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-54.977 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -6.7506e+001 -3.0619e-002 -1.9901e+003 2.8004e+001 -3.1033e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       &"1.0000 HS-  =  S-- +1.0000 H+" //NEW_LINE('')// &
       &"        -llnl_gamma           5.0    " //NEW_LINE('')// &
       &"        log_k           -12.9351" //NEW_LINE('')// &
       &"	-delta_H	49.0364	kJ/mol	# 	S-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	32.928 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 9.7756e+001 3.2913e-002 -5.0784e+003 -4.1812e+001 -7.9273e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 H+  + 2.0000 SO3--  = S2O5--  + H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
       &"        log_k           9.5934" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O5-2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 0.12262E+03    0.62883E-01   -0.18005E+04   -0.50798E+02   -0.28132E+02" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"2.0000 H+ + 1.0000 SO3--  =  SO2 +1.0000 H2O" //NEW_LINE('')// &
       &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
       &"        log_k           +9.0656" //NEW_LINE('')// &
       &"	-delta_H	26.7316	kJ/mol	# 	SO2" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-77.194 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 9.4048e+001 6.2127e-002 -1.1072e+003 -4.0310e+001 -1.7305e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &" " //NEW_LINE('')// &
       !-DB PHASES
       &"PHASES" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
     !  &"#  1122 minerals" //NEW_LINE('')// &
       ! &"Afwillite" //NEW_LINE('')// &
       ! &"        Ca3Si2O4(OH)6 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Ca++ + 6.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           60.0452" //NEW_LINE('')// &
       ! &"	-delta_H	-316.059	kJ/mol	# 	Afwillite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1143.31 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.8353e+001 1.9014e-003 1.8478e+004 -6.6311e+000 -4.0227e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"Akermanite" //NEW_LINE('')// &
       ! &"        Ca2MgSi2O7 +6.0000 H+  =  + 1.0000 Mg++ + 2.0000 Ca++ + 2.0000 SiO2 + 3.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           45.3190" //NEW_LINE('')// &
       ! &"	-delta_H	-288.575	kJ/mol	# 	Akermanite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-926.497 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -4.8295e+001 -8.5613e-003 2.0880e+004 1.3798e+001 -7.1975e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       &"Al" //NEW_LINE('')// &
       &"       Al +3.0000 H+ +0.7500 O2  =  + 1.0000 Al+++ + 1.5000 H2O" //NEW_LINE('')// &
       &"        log_k           149.9292" //NEW_LINE('')// &
       &"	-delta_H	-958.059	kJ/mol	# 	Al" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	0 kJ/mol" //NEW_LINE('')// &
       &"        -analytic -1.8752e+002 -4.6187e-002 5.7127e+004 6.6270e+001 -3.8952e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Al2(SO4)3" //NEW_LINE('')// &
       ! &"       Al2(SO4)3  =  + 2.0000 Al+++ + 3.0000 SO4--" //NEW_LINE('')// &
       ! &"        log_k           19.0535" //NEW_LINE('')// &
       ! &"	-delta_H	-364.566	kJ/mol	# 	Al2(SO4)3" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-3441.04 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -6.1001e+002 -2.4268e-001 2.9194e+004 2.4383e+002 4.5573e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Al2(SO4)3:6H2O" //NEW_LINE('')// &
       ! &"       Al2(SO4)3:6H2O  =  + 2.0000 Al+++ + 3.0000 SO4-- + 6.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           1.6849" //NEW_LINE('')// &
       ! &"	-delta_H	-208.575	kJ/mol	# 	Al2(SO4)3:6H2O" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-5312.06 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -7.1642e+002 -2.4552e-001 2.6064e+004 2.8441e+002 4.0691e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"################################" //NEW_LINE('')// &
       ! &"#      ADDITIONS BY NAVAH      #" //NEW_LINE('')// &
       ! &"################################" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Augite" //NEW_LINE('')// &
       ! &"        Ca.7Fe.6Mg.7Si2O6 +4.0000 H+  =  2.0 H2O + 2.0 SiO2 + .7Ca+2 + 0.6Fe+2 + 0.7Mg+2" //NEW_LINE('')// &
       ! &"        log_k           21.00" //NEW_LINE('')// &
       ! &"	-delta_H	-51.8523	kJ/mol	# 	Augite" //NEW_LINE('')// &
       ! &"	-analytic 7.84710902e+00   7.21674649e-03   1.25039649e+04  -8.82692820e+00  -8.09786954e+05" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Pigeonite" //NEW_LINE('')// &
       ! &"        Ca1.14Fe.64Mg.22Si2O6 +4.0000 H+  =  2.0 H2O + 2.0 SiO2 + 1.14Ca+2 + 0.64Fe+2 + 0.22Mg+2" //NEW_LINE('')// &
       ! &"        log_k           21.40" //NEW_LINE('')// &
       ! &"	-delta_H	-51.8523	kJ/mol	# 	Pigeonite" //NEW_LINE('')// &
       ! &"	-analytic 3.92773074e+01   1.11617261e-02   1.07613145e+04  -1.98006851e+01  -7.39527557e+05" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Plagioclase" //NEW_LINE('')// &
       ! &"        Ca.5Na.5Al1.5Si2.5O8 +6.0000 H+  =  3.0000 H2O + 2.5 SiO2 + .5 Ca+2 + 1.5 Al+3 + 0.5Na+" //NEW_LINE('')// &
       ! &"        log_k           14.20" //NEW_LINE('')// &
       ! &"	-delta_H	-51.8523	kJ/mol	# 	Plagioclase" //NEW_LINE('')// &
       ! &"	-analytic -3.80343385e+01  -7.37083665e-03   1.59944487e+04   4.95599390e+00  -1.01574822e+06" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"################################" //NEW_LINE('')// &
       ! &"#        END NAVAHBLOCK        #" //NEW_LINE('')// &
       ! &"################################" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Albite" //NEW_LINE('')// &
       ! &"        NaAlSi3O8 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 Na+ + 2.0000 H2O + 3.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           2.7645" //NEW_LINE('')// &
       ! &"	-delta_H	-51.8523	kJ/mol	# 	Albite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-939.68 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.1694e+001 1.4429e-002 1.3784e+004 -7.2866e+000 -1.6136e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Analcime" //NEW_LINE('')// &
       &"        Na.96Al.96Si2.04O6:H2O +3.8400 H+  =  + 0.9600 Al+++ + 0.9600 Na+ + 2.0400 SiO2 + 2.9200 H2O" //NEW_LINE('')// &
       &"        log_k           6.1396" //NEW_LINE('')// &
       &"	-delta_H	-75.844	kJ/mol	# 	Analcime" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-3296.86 kJ/mol" //NEW_LINE('')// &
       &"        -analytic -6.8694e+000 6.6052e-003 9.8260e+003 -4.8540e+000 -8.8780e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Andradite" //NEW_LINE('')// &
       ! &"        Ca3Fe2(SiO4)3 +12.0000 H+  =  + 2.0000 Fe+++ + 3.0000 Ca++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           33.3352" //NEW_LINE('')// &
       ! &"	-delta_H	-301.173	kJ/mol	# 	Andradite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1380.35 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.3884e+001 -2.3886e-002 1.5314e+004 -8.1606e+000 -4.2193e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Saponite-Ca" //NEW_LINE('')// &
       &"        Ca.165Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.1650 Ca++ + 0.3300 Al+++ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           26.2900" //NEW_LINE('')// &
       &"	-delta_H	-207.971	kJ/mol	# Saponite-Ca" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1436.51 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -4.6904e+001 6.2555e-003 2.2572e+004 5.3198e+000 -1.5725e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Anhydrite" //NEW_LINE('')// &
       ! &"        CaSO4  =  + 1.0000 Ca++ + 1.0000 SO4--" //NEW_LINE('')// &
       ! &"        log_k           -4.3064" //NEW_LINE('')// &
       ! &"	-delta_H	-18.577	kJ/mol	# 	Anhydrite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-342.76 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.0986e+002 -7.8823e-002 5.0969e+003 8.5642e+001 7.9594e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"-Vm 46.1 # 136.14 / 2.95" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"Phillipsite" //NEW_LINE('')// &
       ! &"        Na0.5K0.5AlSi3O8:H2O + 4H+ = 0.5Na+ +0.5K+ + 3H2O + Al+3 + 3SiO2" //NEW_LINE('')// &
       &"        Na0.5K0.5AlSi3O8:H2O + 7H2O = 0.5Na+ +0.5K+ + Al(OH)4- + 6H2O + 3SiO2" //NEW_LINE('')// &
       &"        log_k           -19.874" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &

       ! &"Aragonite" //NEW_LINE('')// &
       ! &"        CaCO3 +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 HCO3-" //NEW_LINE('')// &
       ! &"        log_k           1.9931" //NEW_LINE('')// &
       ! &"	-delta_H	-25.8027	kJ/mol	# 	Aragonite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-288.531 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.4934e+002 -4.8043e-002 4.9089e+003 6.0284e+001 7.6644e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! !&"-Vm 34.04" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"Calcite" //NEW_LINE('')// &
       &"        CaCO3 +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 HCO3-" //NEW_LINE('')// &
       &"        log_k           1.8487" //NEW_LINE('')// &
       &"	-delta_H	-25.7149	kJ/mol	# 	Calcite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-288.552 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.4978e+002 -4.8370e-002 4.8974e+003 6.0458e+001 7.6464e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       !&"-Vm 36.9 cm3/mol # MW (100.09 g/mol) / rho (2.71 g/cm3)" //NEW_LINE('')// &

       ! &"" //NEW_LINE('')// &
       !  &"Calcite" //NEW_LINE('')// &
       !  &"	CaCO3 = CO3-2 + Ca+2" //NEW_LINE('')// &
       !  &"	-log_k	-8.48" //NEW_LINE('')// &
       !  &"	-delta_h -2.297 kcal" //NEW_LINE('')// &
       !  &"	-analytic	-171.9065	-0.077993	2839.319	71.595" //NEW_LINE('')// &
       !  &"	-Vm 36.9 cm3/mol # MW (100.09 g/mol) / rho (2.71 g/cm3)" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       !  &"Aragonite" //NEW_LINE('')// &
       !  &"	CaCO3 = CO3-2 + Ca+2" //NEW_LINE('')// &
       !  &"	-log_k	-8.336" //NEW_LINE('')// &
       !  &"	-delta_h -2.589 kcal" //NEW_LINE('')// &
       !  &"	-analytic	-171.9773	-0.077993	2903.293	71.595" //NEW_LINE('')// &
       !  &"	-Vm 34.04" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"Celadonite" //NEW_LINE('')// &
       &"        KMgAlSi4O10(OH)2 +6.0000 H+  =  + 1.0000 Al+++ + 1.0000 K+ + 1.0000 Mg++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           7.4575" //NEW_LINE('')// &
       &"	-delta_H	-74.3957	kJ/mol	# 	Celadonite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1394.9 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -3.3097e+001 1.7989e-002 1.8919e+004 -2.1219e+000 -2.0588e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Clinoptilolite-Ca" //NEW_LINE('')// &
       &"        Ca1.7335Al3.45Fe.017Si14.533O36:10.922H2O +13.8680 H+  =  + 0.0170 Fe+++ + 1.7335 Ca++ + 3.4500 Al+++ + 14.5330 SiO2 + 17.8560 H2O" //NEW_LINE('')// &
       &"        log_k           -7.0095" //NEW_LINE('')// &
       &"	-delta_H	-74.6745	kJ/mol	# 	Clinoptilolite-Ca" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-4919.84 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -4.4820e+001 5.3696e-002 5.4878e+004 -3.1459e+001 -7.5491e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Clinozoisite" //NEW_LINE('')// &
       ! &"        Ca2Al3Si3O12(OH) +13.0000 H+  =  + 2.0000 Ca++ + 3.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           43.2569" //NEW_LINE('')// &
       ! &"	-delta_H	-457.755	kJ/mol	# 	Clinozoisite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1643.78 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.8690e+001 -3.7056e-002 2.2770e+004 3.7880e+000 -2.5834e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Cronstedtite-7A" //NEW_LINE('')// &
       ! &"        Fe2Fe2SiO5(OH)4 +10.0000 H+  =  + 1.0000 SiO2 + 2.0000 Fe++ + 2.0000 Fe+++ + 7.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           16.2603" //NEW_LINE('')// &
       ! &"	-delta_H	-244.266	kJ/mol	# 	Cronstedtite-7A" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-697.413 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.3783e+002 -7.1026e-002 1.7752e+004 8.7147e+001 2.7707e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Daphnite-14A" //NEW_LINE('')// &
       &"        Fe5AlAlSi3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Fe++ + 12.0000 H2O" //NEW_LINE('')// &
       &"        log_k           52.2821" //NEW_LINE('')// &
       &"	-delta_H	-517.561	kJ/mol	# 	Daphnite-14A" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1693.04 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.5261e+002 -6.1392e-002 2.8283e+004 5.1788e+001 4.4137e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Daphnite-7A" //NEW_LINE('')// &
       &"        Fe5AlAlSi3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Fe++ + 12.0000 H2O" //NEW_LINE('')// &
       &"        log_k           55.6554" //NEW_LINE('')// &
       &"	-delta_H	-532.326	kJ/mol	# 	Daphnite-7A" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1689.51 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.6430e+002 -6.3160e-002 2.9499e+004 5.6442e+001 4.6035e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Dawsonite" //NEW_LINE('')// &
       ! &"        NaAlCO3(OH)2 +3.0000 H+  =  + 1.0000 Al+++ + 1.0000 HCO3- + 1.0000 Na+ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           4.3464" //NEW_LINE('')// &
       ! &"	-delta_H	-76.3549	kJ/mol	# 	Dawsonite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1963.96 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.1393e+002 -2.3487e-002 7.1758e+003 4.0900e+001 1.2189e+002" //NEW_LINE('')// &
       ! ! &"#       -Range:  0-200" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Diaspore" //NEW_LINE('')// &
       ! &"        AlHO2 +3.0000 H+  =  + 1.0000 Al+++ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           7.1603" //NEW_LINE('')// &
       ! &"	-delta_H	-110.42	kJ/mol	# 	Diaspore" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-238.924 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.2618e+002 -3.1671e-002 8.8737e+003 4.5669e+001 1.3850e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Dicalcium_silicate" //NEW_LINE('')// &
       ! &"        Ca2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           37.1725" //NEW_LINE('')// &
       ! &"	-delta_H	-217.642	kJ/mol	# 	Dicalcium_silicate" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-2317.9 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -5.9723e+001 -1.3682e-002 1.5461e+004 2.1547e+001 -3.7732e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Diopside" //NEW_LINE('')// &
       ! &"        CaMgSi2O6 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 2.0000 H2O + 2.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           20.9643" //NEW_LINE('')// &
       ! &"	-delta_H	-133.775	kJ/mol	# 	Diopside" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-765.378 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 7.1240e+001 1.5514e-002 8.1437e+003 -3.0672e+001 -5.6880e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Dolomite" //NEW_LINE('')// &
       ! &"        CaMg(CO3)2 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 2.0000 HCO3-" //NEW_LINE('')// &
       ! &"        log_k           2.5135" //NEW_LINE('')// &
       ! &"	-delta_H	-59.9651	kJ/mol	# 	Dolomite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-556.631 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -3.1782e+002 -9.8179e-002 1.0845e+004 1.2657e+002 1.6932e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"-Vm 64.5" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"Epidote" //NEW_LINE('')// &
       &"        Ca2FeAl2Si3O12OH +13.0000 H+  =  + 1.0000 Fe+++ + 2.0000 Al+++ + 2.0000 Ca++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
       &"        log_k           32.9296" //NEW_LINE('')// &
       &"	-delta_H	-386.451	kJ/mol	# 	Epidote" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1543.99 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.6187e+001 -3.6436e-002 1.9351e+004 3.3671e+000 -3.0319e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Ettringite" //NEW_LINE('')// &
       ! &"        Ca6Al2(SO4)3(OH)12:26H2O +12.0000 H+  =  + 2.0000 Al+++ + 3.0000 SO4-- + 6.0000 Ca++ + 38.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           62.5362" //NEW_LINE('')// &
       ! &"	-delta_H	-382.451	kJ/mol	# 	Ettringite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-4193 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.0576e+003 -1.1585e-001 5.9580e+004 3.8585e+002 1.0121e+003" //NEW_LINE('')// &
       ! ! &"#       -Range:  0-200" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Fayalite" //NEW_LINE('')// &
       ! &"        Fe2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Fe++ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           19.1113" //NEW_LINE('')// &
       ! &"	-delta_H	-152.256	kJ/mol	# 	Fayalite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-354.119 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.3853e+001 -3.5501e-003 7.1496e+003 -6.8710e+000 -6.3310e+004" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Ferrite-Ca" //NEW_LINE('')// &
       ! &"        CaFe2O4 +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Fe+++ + 4.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           21.5217" //NEW_LINE('')// &
       ! &"	-delta_H	-264.738	kJ/mol	# 	Ferrite-Ca" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-363.494 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.8472e+002 -7.5870e-002 2.0688e+004 1.0485e+002 3.2289e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Foshagite" //NEW_LINE('')// &
       ! &"        Ca4Si3O9(OH)2:0.5H2O +8.0000 H+  =  + 3.0000 SiO2 + 4.0000 Ca++ + 5.5000 H2O" //NEW_LINE('')// &
       ! &"        log_k           65.9210" //NEW_LINE('')// &
       ! &"	-delta_H	-359.839	kJ/mol	# 	Foshagite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1438.27 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 2.9983e+001 5.5272e-003 2.3427e+004 -1.3879e+001 -8.9461e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Gismondine" //NEW_LINE('')// &
       &"        Ca2Al4Si4O16:9H2O +16.0000 H+  =  + 2.0000 Ca++ + 4.0000 Al+++ + 4.0000 SiO2 + 17.0000 H2O" //NEW_LINE('')// &
       &"        log_k           41.7170" //NEW_LINE('')// &
       &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Gismondine" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	0 kcal/mol" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Goethite" //NEW_LINE('')// &
       &"        FeOOH +3.0000 H+  =  + 1.0000 Fe+++ + 2.0000 H2O" //NEW_LINE('')// &
       &"        log_k           0.5345" //NEW_LINE('')// &
       &"	-delta_H	-61.9291	kJ/mol	# 	Goethite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-559.328 kJ/mol" //NEW_LINE('')// &
       &"        -analytic -6.0331e+001 -1.0847e-002 4.7759e+003 1.9429e+001 8.1122e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-200" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Greenalite" //NEW_LINE('')// &
       ! &"        Fe3Si2O5(OH)4 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Fe++ + 5.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           22.6701" //NEW_LINE('')// &
       ! &"	-delta_H	-165.297	kJ/mol	# 	Greenalite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-787.778 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.4187e+001 -3.8377e-003 1.1710e+004 1.6442e+000 -4.8290e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Gyrolite" //NEW_LINE('')// &
       ! &"        Ca2Si3O7(OH)2:1.5H2O +4.0000 H+  =  + 2.0000 Ca++ + 3.0000 SiO2 + 4.5000 H2O" //NEW_LINE('')// &
       ! &"        log_k           22.9099" //NEW_LINE('')// &
       ! &"	-delta_H	-82.862	kJ/mol	# 	Gyrolite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1176.55 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.4416e+001 1.4646e-002 1.6181e+004 2.3723e+000 -1.5369e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Hedenbergite" //NEW_LINE('')// &
       ! &"        CaFe(SiO3)2 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Fe++ + 2.0000 H2O + 2.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           19.6060" //NEW_LINE('')// &
       ! &"	-delta_H	-124.507	kJ/mol	# 	Hedenbergite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-678.276 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.9473e+001 1.5288e-003 1.2910e+004 2.1729e+000 -9.0058e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Hematite" //NEW_LINE('')// &
       &"        Fe2O3 +6.0000 H+  =  + 2.0000 Fe+++ + 3.0000 H2O" //NEW_LINE('')// &
       &"        log_k           0.1086" //NEW_LINE('')// &
       &"	-delta_H	-129.415	kJ/mol	# 	Hematite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-197.72 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.2015e+002 -6.0290e-002 1.1812e+004 8.0253e+001 1.8438e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Hillebrandite" //NEW_LINE('')// &
       ! &"        Ca2SiO3(OH)2:0.17H2O +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 3.1700 H2O" //NEW_LINE('')// &
       ! &"        log_k           36.8190" //NEW_LINE('')// &
       ! &"	-delta_H	-203.074	kJ/mol	# 	Hillebrandite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-637.404 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.9360e+001 -7.5176e-003 1.1947e+004 8.0558e+000 -1.4504e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"K-Feldspar" //NEW_LINE('')// &
       &"        KAlSi3O8 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 K+ + 2.0000 H2O + 3.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           -0.2753" //NEW_LINE('')// &
       &"	-delta_H	-23.9408	kJ/mol	# 	K-Feldspar" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-949.188 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.0684e+000 1.3111e-002 1.1671e+004 -9.9129e+000 -1.5855e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Kaolinite" //NEW_LINE('')// &
       &"        Al2Si2O5(OH)4 +6.0000 H+  =  + 2.0000 Al+++ + 2.0000 SiO2 + 5.0000 H2O" //NEW_LINE('')// &
       &"        log_k           6.8101" //NEW_LINE('')// &
       &"	-delta_H	-151.779	kJ/mol	# 	Kaolinite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-982.221 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.6835e+001 -7.8939e-003 7.7636e+003 -1.2190e+001 -3.2354e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Larnite" //NEW_LINE('')// &
       ! &"        Ca2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           38.4665" //NEW_LINE('')// &
       ! &"	-delta_H	-227.061	kJ/mol	# 	Larnite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-551.74 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 2.6900e+001 -2.1833e-003 1.0900e+004 -9.5257e+000 -7.2537e+004" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Laumontite" //NEW_LINE('')// &
       ! &"        CaAl2Si4O12:4H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 4.0000 SiO2 + 8.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           13.6667" //NEW_LINE('')// &
       ! &"	-delta_H	-184.657	kJ/mol	# 	Laumontite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1728.66 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.1904e+000 8.1763e-003 1.9005e+004 -1.4561e+001 -1.5851e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Lawsonite" //NEW_LINE('')// &
       ! &"        CaAl2Si2O7(OH)2:H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 2.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           22.2132" //NEW_LINE('')// &
       ! &"	-delta_H	-244.806	kJ/mol	# 	Lawsonite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1158.1 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.3995e+001 -1.7668e-002 1.0119e+004 -8.3100e+000 1.5789e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &

       ! &"Magnesite" //NEW_LINE('')// &
       ! &"        MgCO3 +1.0000 H+  =  + 1.0000 HCO3- + 1.0000 Mg++" //NEW_LINE('')// &
       ! &"        log_k           2.2936" //NEW_LINE('')// &
       ! &"	-delta_H	-44.4968	kJ/mol	# 	Magnesite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-265.63 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.6665e+002 -4.9469e-002 6.4344e+003 6.5506e+001 1.0045e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Magnetite" //NEW_LINE('')// &
       ! &"        Fe3O4 +8.0000 H+  =  + 1.0000 Fe++ + 2.0000 Fe+++ + 4.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           10.4724" //NEW_LINE('')// &
       ! &"	-delta_H	-216.597	kJ/mol	# 	Magnetite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-267.25 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -3.0510e+002 -7.9919e-002 1.8709e+004 1.1178e+002 2.9203e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Merwinite" //NEW_LINE('')// &
       ! &"        MgCa3(SiO4)2 +8.0000 H+  =  + 1.0000 Mg++ + 2.0000 SiO2 + 3.0000 Ca++ + 4.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           68.5140" //NEW_LINE('')// &
       ! &"	-delta_H	-430.069	kJ/mol	# 	Merwinite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1090.8 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.2524e+002 -4.2525e-002 3.5619e+004 7.9984e+001 -9.8259e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Mesolite" //NEW_LINE('')// &
       &"        Na.676Ca.657Al1.99Si3.01O10:2.647H2O +7.9600 H+  =  + 0.6570 Ca++ + 0.6760 Na+ + 1.9900 Al+++ + 3.0100 SiO2 + 6.6270 H2O" //NEW_LINE('')// &
       &"        log_k           13.6191" //NEW_LINE('')// &
       &"	-delta_H	-179.744	kJ/mol	# 	Mesolite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-5947.05 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 7.1993e+000 5.9356e-003 1.4717e+004 -1.3627e+001 -9.8863e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Minnesotaite" //NEW_LINE('')// &
       ! &"        Fe3Si4O10(OH)2 +6.0000 H+  =  + 3.0000 Fe++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           13.9805" //NEW_LINE('')// &
       ! &"	-delta_H	-105.211	kJ/mol	# 	Minnesotaite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1153.37 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.8812e+001 1.7261e-002 1.9804e+004 -6.4410e+000 -2.0433e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Monticellite" //NEW_LINE('')// &
       ! &"        CaMgSiO4 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 1.0000 SiO2 + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           29.5852" //NEW_LINE('')// &
       ! &"	-delta_H	-195.711	kJ/mol	# 	Monticellite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-540.8 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.5730e+001 -3.5567e-003 9.0789e+003 -6.3007e+000 1.4166e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Montmor-Na" //NEW_LINE('')// &
       &"        Na.33Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.3300 Mg++ + 0.3300 Na+ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           2.4844" //NEW_LINE('')// &
       &"	-delta_H	-93.2165	kJ/mol	# 	Montmor-Na" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1360.69 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.9601e+000 1.1342e-002 1.6051e+004 -1.4718e+001 -1.8160e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &



       ! &"Montmor-k" //NEW_LINE('')// &
       ! &"        K.33Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.3300 K+ + 0.3300 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           2.1423" //NEW_LINE('')// &
       ! &"	-delta_H	-88.184	kJ/mol	# Calculated enthalpy of reaction	Montmor-K" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1362.83 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 8.4757e+000 1.1219e-002 1.5654e+004 -1.6833e+001 -1.8386e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &

       &"Montmor-Mg" //NEW_LINE('')// &
       &"        Mg.495Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.4950 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           2.3879" //NEW_LINE('')// &
       &"	-delta_H	-102.608	kJ/mol	# Calculated enthalpy of reaction	Montmor-Mg" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1357.87 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -6.8505e+000 9.0710e-003 1.6817e+004 -1.1887e+001 -1.8323e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &

       &"Montmor-Ca" //NEW_LINE('')// &
       &"        Ca.165Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.1650 Ca++ + 0.3300 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           2.4952" //NEW_LINE('')// &
       &"	-delta_H	-100.154	kJ/mol	# Calculated enthalpy of reaction	Montmor-Ca" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1361.5 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 6.0725e+000 1.0644e-002 1.6024e+004 -1.6334e+001 -1.7982e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &



       ! &"Muscovite" //NEW_LINE('')// &
       ! &"        KAl3Si3O10(OH)2 +10.0000 H+  =  + 1.0000 K+ + 3.0000 Al+++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           13.5858" //NEW_LINE('')// &
       ! &"	-delta_H	-243.224	kJ/mol	# 	Muscovite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1427.41 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 3.3085e+001 -1.2425e-002 1.2477e+004 -2.0865e+001 -5.4692e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Natrolite" //NEW_LINE('')// &
       &"        Na2Al2Si3O10:2H2O +8.0000 H+  =  + 2.0000 Al+++ + 2.0000 Na+ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
       &"        log_k           18.5204" //NEW_LINE('')// &
       &"	-delta_H	-186.971	kJ/mol	# 	Natrolite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-5718.56 kJ/mol" //NEW_LINE('')// &
       &"        -analytic -2.7712e+001 -2.7963e-003 1.6075e+004 1.5332e+000 -9.5765e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Nontronite-Ca" //NEW_LINE('')// &
       &"        Ca.165Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.1650 Ca++ + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           -11.5822" //NEW_LINE('')// &
       &"	-delta_H	-38.138	kJ/mol	# 	Nontronite-Ca" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1166.7 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.6291e+001 4.3557e-003 1.0221e+004 -1.8690e+001 -1.5427e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Nontronite-H" //NEW_LINE('')// &
       ! &"        H.33Fe2Al.33Si3.67H2O12 +6.9900 H+  =  + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       ! &"        log_k           -12.5401" //NEW_LINE('')// &
       ! &"	-delta_H	-30.452	kJ/mol	# 	Nontronite-H" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1147.12 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 9.7794e+001 1.4055e-002 4.7440e+003 -4.7272e+001 -1.2103e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Nontronite-K" //NEW_LINE('')// &
       &"        K.33Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 K+ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           -11.8648" //NEW_LINE('')// &
       &"	-delta_H	-26.5822	kJ/mol	# 	Nontronite-K" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1167.93 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.3630e+001 4.7708e-003 1.0073e+004 -1.7407e+001 -1.5803e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Nontronite-Mg" //NEW_LINE('')// &
       &"        Mg.165Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.1650 Mg++ + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           -11.6200" //NEW_LINE('')// &
       &"	-delta_H	-41.1779	kJ/mol	# 	Nontronite-Mg" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1162.93 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 5.5961e+001 1.0139e-002 8.0777e+003 -3.3164e+001 -1.4031e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Nontronite-Na" //NEW_LINE('')// &
       &"        Na.33Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 Na+ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           -11.5263" //NEW_LINE('')// &
       &"	-delta_H	-31.5687	kJ/mol	# 	Nontronite-Na" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1165.8 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 6.7915e+001 1.2851e-002 7.1218e+003 -3.7112e+001 -1.3758e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Fe-Celadonite" //NEW_LINE('')// &
       &"        KFeAlSi4O10(OH)2 = +1.000K+     +1.000Fe+2     +1.000Al+3     -6.000H+     +8.0H2O + 4SiO2    -4.000H2O     " //NEW_LINE('')// &
       &"        log_k           2.73" //NEW_LINE('')// &
       &"	-delta_H	-83.838	kJ/mol	# 	Fe-Celadonite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-5498.159 kJ/mol" //NEW_LINE('')// &
       ! &"Okenite" //NEW_LINE('')// &
       ! &"        CaSi2O4(OH)2:H2O +2.0000 H+  =  + 1.0000 Ca++ + 2.0000 SiO2 + 3.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           10.3816" //NEW_LINE('')// &
       ! &"	-delta_H	-19.4974	kJ/mol	# 	Okenite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-749.641 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -7.7353e+001 1.5091e-002 1.3023e+004 2.1337e+001 -1.1831e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Smectite-low-Fe-Mg" //NEW_LINE('')// &
       &"        Ca.02Na.15K.2Fe.29Fe.16Mg.9Al1.25Si3.75H2O12 +7.0000 H+  =  + 0.0200 Ca++ + 0.1500 Na+ + 0.1600 Fe+++ + 0.2000" // &
       &" K+ + 0.2900 Fe++ + 0.9000 Mg++ + 1.2500 Al+++ + 3.7500 SiO2 + 4.5000 H2O" //NEW_LINE('')// &
       &"        log_k           11.0405" //NEW_LINE('')// &
       &"	-delta_H	-144.774	kJ/mol	# Smectite-low-Fe-Mg" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1352.12 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.7003e+001 6.9848e-003 1.8359e+004 -6.8896e+000 -1.6637e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Prehnite" //NEW_LINE('')// &
       &"        Ca2Al2Si3O10(OH)2 +10.0000 H+  =  + 2.0000 Al+++ + 2.0000 Ca++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
       &"        log_k           32.9305" //NEW_LINE('')// &
       &"	-delta_H	-311.875	kJ/mol	# 	Prehnite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1481.65 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -3.5763e+001 -2.1396e-002 2.0167e+004 6.3554e+000 -7.4967e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Pseudowollastonite" //NEW_LINE('')// &
       ! &"        CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           13.9997" //NEW_LINE('')// &
       ! &"	-delta_H	-79.4625	kJ/mol	# 	Pseudowollastonite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-388.9 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 2.6691e+001 6.3323e-003 5.5723e+003 -1.1822e+001 -3.6038e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Pyrite" //NEW_LINE('')// &
       &"        FeS2 +1.0000 H2O  =  + 0.2500 H+ + 0.2500 SO4-- + 1.0000 Fe++ + 1.7500 HS-" //NEW_LINE('')// &
       &"        log_k           -24.6534" //NEW_LINE('')// &
       &"	-delta_H	109.535	kJ/mol	# 	Pyrite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-41 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.4195e+002 -8.7948e-002 -6.2911e+002 9.9248e+001 -9.7454e+000" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Pyrrhotite" //NEW_LINE('')// &
       &"        FeS +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HS-" //NEW_LINE('')// &
       &"        log_k           -3.7193" //NEW_LINE('')// &
       &"	-delta_H	-7.9496	kJ/mol	# 	Pyrrhotite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-24 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -1.5785e+002 -5.2258e-002 3.9711e+003 6.3195e+001 6.2012e+001" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Quartz" //NEW_LINE('')// &
       ! &"        SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           -3.9993" //NEW_LINE('')// &
       ! &"	-delta_H	32.949	kJ/mol	# 	Quartz" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-217.65 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 7.7698e-002 1.0612e-002 3.4651e+003 -4.3551e+000 -7.2138e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"-Vm 22.67" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       ! &"Rankinite" //NEW_LINE('')// &
       ! &"        Ca3Si2O7 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Ca++ + 3.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           51.9078" //NEW_LINE('')// &
       ! &"	-delta_H	-302.089	kJ/mol	# 	Rankinite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-941.7 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -9.6393e+001 -1.6592e-002 2.4832e+004 3.2541e+001 -9.4630e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       &"Saponite-Mg" //NEW_LINE('')// &
       !&"        Mg3Ca.165Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + .1650 Ca ++ + 3 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       !&"" //NEW_LINE('')// &
       &"	Mg3.165Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 3.1650 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           26.2523" //NEW_LINE('')// &
       &"	-delta_H	-210.822	kJ/mol	# 	Saponite-Mg" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1432.79 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 9.8888e+000 1.4320e-002 1.9418e+004 -1.5259e+001 -1.3716e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Saponite-Na" //NEW_LINE('')// &
       &"        Na.33Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 Na+ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
       &"        log_k           26.3459" //NEW_LINE('')// &
       &"	-delta_H	-201.401	kJ/mol	# 	Saponite-Na" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1435.61 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -6.7611e+001 4.7327e-003 2.3586e+004 1.2868e+001 -1.6493e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       &"Scolecite" //NEW_LINE('')// &
       &"        CaAl2Si3O10:3H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
       &"        log_k           15.8767" //NEW_LINE('')// &
       &"	-delta_H	-204.93	kJ/mol	# 	Scolecite" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-6048.92 kJ/mol" //NEW_LINE('')// &
       &"        -analytic 5.0656e+001 -3.1485e-003 1.0574e+004 -2.5663e+001 -5.2769e+005" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"SiO2(am)" //NEW_LINE('')// &
       ! &"       SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           -2.7136" //NEW_LINE('')// &
       ! &"	-delta_H	20.0539	kJ/mol	# 	SiO2(am)" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-214.568 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.2109e+000 7.0767e-003 2.3634e+003 -3.4449e+000 -4.8591e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Siderite" //NEW_LINE('')// &
       ! &"        FeCO3 +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HCO3-" //NEW_LINE('')// &
       ! &"        log_k           -0.1920" //NEW_LINE('')// &
       ! &"	-delta_H	-32.5306	kJ/mol	# 	Siderite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-179.173 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.5990e+002 -4.9361e-002 5.4947e+003 6.3032e+001 8.5787e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"-Vm 29.2" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       &"Smectite-high-Fe-Mg" //NEW_LINE('')// &
       ! &"#        Ca.025Na.1K.2Fe++.5Fe+++.2Mg1.15Al1.25Si3.5H2O12 +8.0000 H+  =  + 0.0250 Ca++ + 0.1000 Na+ + 0.2000 Fe+++ + 0.2000 K+ " //&
       ! &"+ 0.5000 Fe++ + 1.1500 Mg++ + 1.2500 Al+++ + 3.5000 SiO2 + 5.0000 H2O" //NEW_LINE('')// &
       &"        Ca.025Na.1K.2Fe.5Fe.2Mg1.15Al1.25Si3.5H2O12 +8.0000 H+  =  + 0.0250 Ca++ + 0.1000 Na+ + 0.2000 Fe+++ + 0.2000 K+ + 0.5000 Fe++ " //&
       &"+ 1.1500 Mg++ + 1.2500 Al+++ + 3.5000 SiO2 + 5.0000 H2O " //NEW_LINE('')// &
       &"        log_k           17.4200" //NEW_LINE('')// &
       &"	-delta_H	-199.841	kJ/mol	# 	Smectite-high-Fe-Mg" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1351.39 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -9.6102e+000 1.2551e-003 1.8157e+004 -7.9862e+000 -1.3005e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &
       &"" //NEW_LINE('')// &
       ! &"Stilbite" //NEW_LINE('')// &
       ! &"        Ca1.019Na.136K.006Al2.18Si6.82O18:7.33H2O +8.7200 H+  =  + 0.0060 K+ + 0.1360 Na+ " //&
       ! &"+ 1.0190 Ca++ + 2.1800 Al+++ + 6.8200 SiO2 + 11.6900 H2O" //NEW_LINE('')// &
       ! &"        log_k           1.0545" //NEW_LINE('')// &
       ! &"	-delta_H	-83.0019	kJ/mol	# 	Stilbite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-11005.7 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.4483e+001 3.0987e-002 2.8013e+004 -1.5802e+001 -3.4491e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Tobermorite-9A" //NEW_LINE('')// &
       ! &"        Ca5Si6H6O20 +10.0000 H+  =  + 5.0000 Ca++ + 6.0000 SiO2 + 8.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           69.0798" //NEW_LINE('')// &
       ! &"	-delta_H	-329.557	kJ/mol	# 	Tobermorite-9A" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-2375.42 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -6.3384e+001 1.1722e-002 3.8954e+004 1.2268e+001 -2.8681e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Tremolite" //NEW_LINE('')// &
       ! &"        Ca2Mg5Si8O22(OH)2 +14.0000 H+  =  + 2.0000 Ca++ + 5.0000 Mg++ + 8.0000 H2O + 8.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           61.2367" //NEW_LINE('')// &
       ! &"	-delta_H	-406.404	kJ/mol	# 	Tremolite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-2944.04 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 8.5291e+001 4.6337e-002 3.9465e+004 -5.4414e+001 -3.1913e+006" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &

       ! &"Fe" //NEW_LINE('')// &
       ! &"        Fe +2.0000 H+ +0.5000 O2  =  + 1.0000 Fe++ + 1.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           59.0325" //NEW_LINE('')// &
       ! &"	-delta_H	-372.029	kJ/mol		Fe" //NEW_LINE('')// &
       ! &"# Enthalpy of formation:	0 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -6.2882e+001 -2.0379e-002 2.0690e+004 2.3673e+001 3.2287e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &


       ! &"Ferrihydrite" //NEW_LINE('')// &
       ! &"        Fe(OH)3 + 3H+ = Fe+3 + 3H2O" //NEW_LINE('')// &
       ! &"        log_k	3.191" //NEW_LINE('')// &
       ! &"	delta_h	-73.374	kJ" //NEW_LINE('')// &

       ! &"Fe3(OH)8" //NEW_LINE('')// &
       ! &"        Fe3(OH)8 + 8H+ = 2Fe+3 + Fe+2 + 8H2O" //NEW_LINE('')// &
       ! &"        log_k   20.222" //NEW_LINE('')// &
       ! &"	delta_h -0      kcal" //NEW_LINE('')// &

       &"Talc" //NEW_LINE('')// &
       &"        Mg3Si4O10(OH)2 +6.0000 H+  =  + 3.0000 Mg++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
       &"        log_k           21.1383" //NEW_LINE('')// &
       &"-delta_H	-148.737	kJ/mol	# Calculated enthalpy of reaction	Talc" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-1410.92 kcal/mol" //NEW_LINE('')// &
       &"        -analytic 1.1164e+001 2.4724e-002 1.9810e+004 -1.7568e+001 -1.8241e+006" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &


       &"Chlorite(14A)" //NEW_LINE('')// &
       &"        Mg5Al2Si3O10(OH)8 + 16H+ = 5Mg+2 + 2Al+3 + 3.0 SiO2 + 12H2O" //NEW_LINE('')// &
       &"        log_k           68.38" //NEW_LINE('')// &
       &"delta_h -151.494 kcal" //NEW_LINE('')// &

       &"" //NEW_LINE('')// &
       ! &"Fe(OH)2" //NEW_LINE('')// &
       ! &"        Fe(OH)2 +2.0000 H+  =  + 1.0000 Fe++ + 2.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           13.9045" //NEW_LINE('')// &
       ! &"	-delta_H	-95.4089	kJ/mol		Fe(OH)2" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-568.525 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -8.6666e+001 -1.8440e-002 7.5723e+003 3.2597e+001 1.1818e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &

       ! CHLORITE MINERALS
       !
       ! &"Chamosite-7A" //NEW_LINE('')// &
       ! &"        Fe2Al2SiO5(OH)4 +10.0000 H+  =  + 1.0000 SiO2 + 2.0000 Al+++ + 2.0000 Fe++ + 7.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           32.8416" //NEW_LINE('')// &
       ! &"-delta_H	-364.213	kJ/mol	# Calculated enthalpy of reaction	Chamosite-7A" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-902.407 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -2.5581e+002 -7.0890e-002 2.4619e+004 9.1789e+001 3.8424e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       !
       !
       &"Clinochlore-14A" //NEW_LINE('')// &
       &"        Mg5Al2Si3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Mg++ + 12.0000 H2O" //NEW_LINE('')// &
       &"        log_k           67.2391" //NEW_LINE('')// &
       &"-delta_H	-612.379	kJ/mol	# Calculated enthalpy of reaction	Clinochlore-14A" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-2116.96 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.0441e+002 -6.2268e-002 3.5388e+004 6.9239e+001 5.5225e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &


       &"Clinochlore-7A" //NEW_LINE('')// &
       &"        Mg5Al2Si3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Mg++ + 12.0000 H2O" //NEW_LINE('')// &
       &"        log_k           70.6124" //NEW_LINE('')// &
       &"-delta_H	-628.14	kJ/mol	# Calculated enthalpy of reaction	Clinochlore-7A" //NEW_LINE('')// &
       ! &"#	Enthalpy of formation:	-2113.2 kcal/mol" //NEW_LINE('')// &
       &"        -analytic -2.1644e+002 -6.4187e-002 3.6548e+004 7.4123e+001 5.7037e+002" //NEW_LINE('')// &
     !  &"#       -Range:  0-300" //NEW_LINE('')// &


       ! &"Ripidolite-14A" //NEW_LINE('')// &
       ! &"        Mg3Fe2Al2Si3O10(OH)8 +16.00 H+  =  + 2.00 Al+++ + 2.00 Fe++ + 3.00 Mg++ + 3.00 SiO2 + 12.00 H2O" //NEW_LINE('')// &
       ! &"        log_k           60.9638" //NEW_LINE('')// &
       ! &"-delta_H	-572.472	kJ/mol	# Calculated enthalpy of reaction	Ripidolite-14A" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1947.87 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.8376e+002 -6.1934e-002 3.2458e+004 6.2290e+001 5.0653e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       !
       ! &"Ripidolite-7A" //NEW_LINE('')// &
       ! &"        Mg3Fe2Al2Si3O10(OH)8 +16.00 H+  =  + 2.00 Al+++ + 2.00 Fe++ + 3.00 Mg++ + 3.00 SiO2 + 12.00 H2O" //NEW_LINE('')// &
       ! &"        log_k           64.3371" //NEW_LINE('')// &
       ! &"-delta_H	-586.325	kJ/mol	# Calculated enthalpy of reaction	Ripidolite-7A" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1944.56 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.9557e+002 -6.3779e-002 3.3634e+004 6.7057e+001 5.2489e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &

       ! &"Fe(OH)3" //NEW_LINE('')// &
       ! &"        Fe(OH)3 +3.0000 H+  =  + 1.0000 Fe+++ + 3.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           5.6556" //NEW_LINE('')// &
       ! &"	-delta_H	-84.0824	kJ/mol		Fe(OH)3" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-823.013 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.3316e+002 -3.1284e-002 7.9753e+003 4.9052e+001 1.2449e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &

       ! &"Troilite" //NEW_LINE('')// &
       ! &"        FeS +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HS-" //NEW_LINE('')// &
       ! &"        log_k           -3.8184" //NEW_LINE('')// &
       ! &"	-delta_H	-7.3296	kJ/mol	# 	Troilite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-101.036 kJ/mol" //NEW_LINE('')// &
       ! &"        -analytic -1.6146e+002 -5.3170e-002 4.0461e+003 6.4620e+001 6.3183e+001" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Wollastonite" //NEW_LINE('')// &
       ! &"        CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           13.7605" //NEW_LINE('')// &
       ! &"	-delta_H	-76.5756	kJ/mol	# 	Wollastonite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-389.59 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 3.0931e+001 6.7466e-003 5.1749e+003 -1.3209e+001 -3.4579e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Xonotlite" //NEW_LINE('')// &
       ! &"        Ca6Si6O17(OH)2 +12.0000 H+  =  + 6.0000 Ca++ + 6.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           91.8267" //NEW_LINE('')// &
       ! &"	-delta_H	-495.457	kJ/mol	# 	Xonotlite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-2397.25 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 1.6080e+003 3.7309e-001 -2.2548e+004 -6.2716e+002 -3.8346e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-200" //NEW_LINE('')// &
       ! &"" //NEW_LINE('')// &
       ! &"Zoisite" //NEW_LINE('')// &
       ! &"        Ca2Al3(SiO4)3OH +13.0000 H+  =  + 2.0000 Ca++ + 3.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
       ! &"        log_k           43.3017" //NEW_LINE('')// &
       ! &"	-delta_H	-458.131	kJ/mol	# 	Zoisite" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-1643.69 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic 2.5321e+000 -3.5886e-002 1.9902e+004 -6.2443e+000 3.1055e+002" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &

       ! MINS FROM ANOTHER DATABASE
       &"Vermiculite-Na" //NEW_LINE('')// &
       ! &"Na0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850Na+  +0.850Al+3  -9.400H+  +3.150H4(SiO4)  -0.600H2O" //NEW_LINE('')// &
       &"Na0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850Na+  +0.850Al+3  -9.400H+  +6.30H2O + 3.15SiO2  -0.600H2O" //NEW_LINE('')// &
       &"  log_k  40.17  #" //NEW_LINE('')// &
       &"  delta_h  -354.987   kJ/mol  #" //NEW_LINE('')// &
       !&"  # Enthalpy of formation:    -6139.206  kJ/mol  07VIE" //NEW_LINE('')// &

       ! FROM SIT.dat
       &"Fe-Saponite-Ca" //NEW_LINE('')// &
       ! &"Ca0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Ca+2     +3.000Fe+2     +0.340Al+3     -7.360H+     +3.660H4(SiO4)     -2.640H2O" //NEW_LINE('')// &
       ! &"Ca0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Ca+2     +3.000Fe+2     +0.340Al+3     -7.360H+   +7.32H2O  +3.66SiO2     -2.640H2O" //NEW_LINE('')// &
       &"Ca0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Ca+2     +3.000Fe+2     +0.340Al+3     -7.360H+   +4.68H2O  +3.66SiO2" //NEW_LINE('')// &
       &"  log_k     22.43        #" //NEW_LINE('')// &
       &"  delta_h  -222.166      kJ/mol        #" //NEW_LINE('')// &
       !&"  # Enthalpy of formation:             -4916.58      kJ/mol        07VIE" //NEW_LINE('')// &



       ! FROM SIT.dat
       &"Fe-Saponite-Mg" //NEW_LINE('')// &
       ! &"Mg0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Mg+2     +3.000Fe+2     +0.340Al+3     -7.360H+     +3.660H4(SiO4)     -2.640H2O" //NEW_LINE('')// &
       ! &"Mg0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Mg+2     +3.000Fe+2     +0.340Al+3     -7.360H+   +7.32H2O  +3.66SiO2     -2.640H2O" //NEW_LINE('')// &
       &"Mg0.17Fe3Al0.34Si3.66O10(OH)2 = +0.170Mg+2     +3.000Fe+2     +0.340Al+3     -7.360H+   +4.68H2O  +3.66SiO2" //NEW_LINE('')// &
       &"  log_k     21.73        #" //NEW_LINE('')// &
       &"  delta_h  -222.096      kJ/mol        #" //NEW_LINE('')// &
       !&"  # Enthalpy of formation:             -4903.73      kJ/mol        07VIE" //NEW_LINE('')// &
       !
       ! &"Lepidocrocite" //NEW_LINE('')// &
       ! &"FeOOH = +1.000Fe+3  -3.000H+  +2.000H2O  " //NEW_LINE('')// &
       ! &"  log_k  0.75   #98DIA in 98CHI" //NEW_LINE('')// &
       ! &"  delta_h  -64.26  kJ/mol  #" //NEW_LINE('')// &
       ! &"  # Enthalpy of formation:    -556.4  kJ/mol  " //NEW_LINE('')// &
       ! !
       ! &"Vermiculite-K" //NEW_LINE('')// &
       ! &"K0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850K+  +0.850Al+3  -9.400H+  +6.30H2O + 3.15SiO2  -0.600H2O " //NEW_LINE('')// &
       ! &"  log_k  36.86  #" //NEW_LINE('')// &
       ! &"  delta_h  -331.639   kJ/mol  #" //NEW_LINE('')// &
       ! &"  # Enthalpy of formation:    -6172.584  kJ/mol  07VIE" //NEW_LINE('')// &
       !
       ! &"Saponite-K" //NEW_LINE('')// &
       ! &"K0.33Mg3Al0.33Si3.67O10(OH)2 = +3.000Mg+2  +0.330K+  +0.330Al+3  -7.320H+  +7.34H2O + 3.67SiO2  -2.680H2O" //NEW_LINE('')// &
       ! &"  log_k  28.1   #" //NEW_LINE('')// &
       ! &"  delta_h  -252.497   kJ/mol  #" //NEW_LINE('')// &
       ! &"  # Enthalpy of formation:    -6005.94   kJ/mol  07VIE" //NEW_LINE('')// &
       !
       &"Vermiculite-Ca" //NEW_LINE('')// &
       &"Ca0.43Mg3Si3.14Al0.86O10(OH)2 = +0.430Ca+2  +3.000Mg+2  +0.860Al+3  -9.440H+  +6.28H2O + 3.14SiO2  -0.560H2O" //NEW_LINE('')// &
       &"  log_k  40.68  #" //NEW_LINE('')// &
       &"  delta_h  -378.219   kJ/mol  #" //NEW_LINE('')// &
       !&"  # Enthalpy of formation:    -6147.38   kJ/mol  07VIE" //NEW_LINE('')// &


       &"Vermiculite-Mg" //NEW_LINE('')// &
       &"Mg0.43Mg3Si3.14Al0.86O10(OH)2 = +3.430Mg+2  +0.860Al+3  -9.440H+  +6.28H2O + 3.14SiO2  -0.560H2O" //NEW_LINE('')// &
       &"  log_k  38.8   #" //NEW_LINE('')// &
       &"  delta_h  -377.469   kJ/mol  #" //NEW_LINE('')// &
       !&"  # Enthalpy of formation:    -6115.45   kJ/mol  07VIE" //NEW_LINE('')// &

       ! &"Chalcedony" //NEW_LINE('')// &
       ! &"        SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
       ! &"        log_k           -3.7281" //NEW_LINE('')// &
       ! &"	-delta_H	31.4093	kJ/mol	# Calculated enthalpy of reaction	Chalcedony" //NEW_LINE('')// &
       ! ! &"#	Enthalpy of formation:	-217.282 kcal/mol" //NEW_LINE('')// &
       ! &"        -analytic -9.0068e+000 9.3241e-003 4.0535e+003 -1.0830e+000 -7.5077e+005" //NEW_LINE('')// &
       ! &"#       -Range:  0-300" //NEW_LINE('')// &


       &""! //NEW_LINE('')// &



  WRITE(*,*) "testing..."

  ! ! initialize domain geometry
  ! call init()

  !#initialize all processors


  ! process #0 is the root process
  root_process = 0

  ! initialize a process
  CALL MPI_INIT ( ierr )

  ! find out the process ID and how many processes were started so far
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  ! print out current processor id
  WRITE(*,*) "my_id:", my_id
  WRITE(*,*) " "

  ! what to do if you are the master processor
  IF (my_id .EQ. root_process) THEN

     !#master processor does stuff

     ! initialize domain geometry
     CALL init()

     cml3 = coarse_mask_long3

     if (param_d_only .EQ. 1) then
         cml3(1:((xn-1)/cellx)*(yn/(2*celly))) = 0.0
     end if

     !cml3(2*((xn-1)/cellx)*(yn/(2*celly)):) = 0.0

     write(*,*) "before bcast" , SUM(cml3)

     do an_id = 1 , num_procs - 1
         CALL MPI_SEND( cml3, (3*(xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
     end do

    !  CALL MPI_Bcast( cml3, (3*(xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, 0, MPI_COMM_WORLD, status, ierr)
    !  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

     write(*,*) "after bcast"
     active_coarse = SUM(cml3)
     write(*,*) "active_coarse" , active_coarse


     do an_id = 1,num_procs - 1
        par_rounds = active_coarse / (num_procs - 1)
        end_loop = par_rounds
        if (an_id .LE. mod(active_coarse,num_procs-1)) then
            end_loop = end_loop + 1
        end if
        ! send end_loop to processor an_id
        CALL MPI_SEND( end_loop, 1, MPI_INTEGER, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
    end do




     ! diss_toggle
     diss_toggle = 0

     permeability0 = permeability
     phi_coarse = .1
     phi = .1
     phi = 0.1
     phi0 = phi


     ! fill coordinate arrays
     DO i = 1,((xn-1)/cellx)
	DO ii = 1,(yn/(2*celly))
           medium(i,ii,6) = x(i*cellx)
           medium(i,ii,7) = y(ii*celly+yn/(2*cell))
           medium_a(i,ii,6) = x(i*cellx)
           medium_a(i,ii,7) = y(ii*celly+yn/(2*cell))
           medium_b(i,ii,6) = x(i*cellx)
           medium_b(i,ii,7) = y(ii*celly+yn/(2*cell))
	END DO
     END DO

     ! boundary & initial condtions for flow equations
     psi=0.0
     frac6 = 0.0
     ! frac6(:,1) = 1.0e-3
     ! frac6(:,2) = -1.0e-3
     temp6 = 280.0
     psi(1,1:yn) = bcyPsi(1,1:yn)
     psi(xn,1:yn) = bcyPsi(2,1:yn)
     psi(1:xn,1) = bcxPsi(1:xn,1)
     psi(1:xn,yn) = bcxPsi(1:xn,2)
     psi = 0.0!dx*scope*1000.0 - dx*scope*0.0
     rho = rho_fluid
     visc = viscosity

     h = param_tsw

     !no slant IC, 11/15/15
     DO i=1,xn
	DO ii=1,yn

    ! (.68-(0.00133*i*1.0))

           h(i,1) = param_tsw!480.0 + ( first-(0.0016*i*factor)) * dy/1.8

           h(i,ii) = param_tsw!480.0*(h(i,1)/h(1,1)) + (param_tsw-(480.0*(h(i,1)/h(1,1)) ))*((-y_min)+y(ii))/((-y_min))
           !h(i,ii) = h(i,ii) - (400.0 + (param_tsw-400.0)*((-y_min)-max(param_o_rhs,param_o))/((max(param_o_rhs,param_o))))
	END DO
     END DO

     ! do i=1,xn
     ! 	do ii=1,yn
     ! 		if (y(ii) .ge. -max(param_o_rhs,param_o)) then
     ! 			h(i,ii) = 275.0
     ! 		end if
     ! 	end do
     ! end do


     ! put initial values into array for file
     !hmat(1:xn,1:yn) = h
     !psimat(1:xn,1:yn) = psi
     !umat(1:xn,1:yn) = u
     !vmat(1:xn,1:yn) = v

     uTransport = 0.0
     vTransport = 0.0
     u_coarse = 0.0
     v_coarse = 0.0



     IF (restart .EQ. 1) THEN


        ! FLUID DYNAMIC TIME SLICES
        !
        ! 		OPEN(UNIT=10, FILE='/data/navah/ic_saturday/h_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')
        ! 		OPEN(UNIT=11, FILE='/data/navah/ic_saturday/psi_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')


        ! THIS IS SAVEBLOCK!!!
        OPEN(UNIT=10, FILE=TRIM(path) // 'h.txt')
        OPEN(UNIT=11, FILE=TRIM(path) // 'psi.txt')

        DO i = 1,yn
           READ(10,*) (hTrans(i,ii),ii=1,xn)
           READ(11,*) (psiTrans(i,ii),ii=1,xn)
        END DO
        h = TRANSPOSE(hTrans)
        psi = TRANSPOSE(psiTrans)

        ! end restart loop
     END IF




     permx = partial((phi/(permeability)),xn,yn,dx,dy,1)
     permy = partial((phi/(permeability)),xn,yn,dx,dy,2)

     outerBand = make_band(permeability,phi,permx,permy,rho)
     outerBand = band(outerBand,2*((yn/2)-2) + 1,longP)

     ! permx = 0.0
     ! permy = 0.0






     !-dynamics loop

     ! this is the main loop that does all the solving for tn timesteps
     DO j = crashstep, tn
	!do j = 2, 50

        IF (j .EQ. crashstep) THEN
           WRITE(*,*) "STARTING STEP:" , j
           WRITE(*,*) " "
           WRITE(*,*) " "
           WRITE(*,*) " "
           OPEN(UNIT=8, status = 'replace', FILE=TRIM(path_final) // 'dynamicStep.txt')
           WRITE(*,*) "opened"
           WRITE(8,*) 0
           CLOSE ( 8 )
        END IF

        ! write(*,*) "j = ", j


        IF (MOD(j,mstep*10) .EQ. 0) THEN
           WRITE(*,*) "STARTING STEP:" , j
           WRITE(*,*) " "
           WRITE(*,*) " "
           WRITE(*,*) " "
           OPEN(UNIT=8, status = 'replace', FILE=TRIM(path_final) // 'dynamicStep.txt')
           WRITE(*,*) "opened"
           WRITE(8,*) "j/mstep:" , j/mstep
           WRITE(8,*) "j:" , j
           CLOSE ( 8 )
        END IF


        ! write(*,*) "error string"

        ! ! IF (MOD(j,mstep/10) .EQ. 0) THEN
        ! IF (MOD(j,mstep) .EQ. 0) THEN
        !    WRITE(*,*) "WRITING TO DYNAMIC SUB STEP"
        !    WRITE(*,*) " "
        !    WRITE(*,*) " "
        !    WRITE(*,*) " "
        !    ! 		OPEN(UNIT=88, status = 'replace', FILE=trim(path_final) // 'dynamicSubStep.txt')
        !    ! 		write(88,*) mod(j,mstep)
        !    ! 		close ( 88 )
        ! END IF

        ! write(*,*) "error string"
        !write(*,*) "error string"


        IF (restart .NE. 1) THEN

           dt_bit = dt

           h = h_next(h, psi,rho,phi,u,v,frac6,temp6,dt_bit)
           h = h_bc(h)

           ! short outcrop outflow condition
           DO jj=2,yn-1
              DO i=1,xn

                 IF ((v(i,jj) .GE. 0.0) .AND. (mask(i,jj) .EQ. 25.0)  ) THEN
                    h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
                 END IF

                 IF ((v(i,jj) .GE. 0.0) .AND. (mask(i,jj) .EQ. 12.5) ) THEN
                    h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
                 END IF
                 IF ((v(i,jj) .GE. 0.0) .AND. (mask(i,jj) .EQ. 17.5) ) THEN
                    h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
                 END IF

              END DO
           END DO

           ! inner boundaries outflow condtion
           DO jj=2,yn-1
              DO i=1,xn

                 IF ((mask(i,jj) .EQ. 5.0) .AND. (mask(i,jj) .EQ. 5.0) .AND. (u(i,jj) .GT. 0.0)) THEN
                    h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
                 END IF
                 IF ((mask(i,jj) .EQ. 12.5) .AND. (u(i,jj) .GT. 0.0)) THEN
                    h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
                 END IF

                 ! since these are also a 50.0 bc, they maybe aren't wilcock's but isothermal??
                 IF ((mask(i,jj) .EQ. 10.0) .AND. (mask(i,jj) .EQ. 10.0).AND. (u(i,jj) .LT. 0.0)) THEN
                    h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
                 END IF
                 IF ((mask(i,jj) .EQ. 17.5) .AND. (u(i,jj) .LT. 0.0)) THEN
                    h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
                 END IF

              END DO
           END DO

           ! insulating boundaries during spinup phase
           IF (iter .EQ. 0) THEN
              DO ii=2,yn-1
                 IF ((mask(f_index1,ii) .EQ. 6.0) .OR. (mask(f_index1,ii) .EQ. 6.5) .OR. (mask(f_index1,ii) .EQ. 6.1) .OR. (mask(f_index1,ii) .EQ. 6.05)) THEN
                    temp6(ii,2) = (4.0/3.0)*h(f_index1,ii) - (1.0/3.0)*h(f_index1+1,ii)
                 END IF
                 IF ((mask(f_index1-1,ii) .EQ. 3.0) .OR. (mask(f_index1-1,ii) .EQ. 3.5) .OR. (mask(f_index1-1,ii) .EQ. 3.1) .OR. (mask(f_index1-1,ii) .EQ. 3.05)) THEN
                    temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
                 END IF

                 IF ((mask(f_index1,ii) .EQ. 6.5)) THEN
                    temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
                 END IF

                 IF ((mask(f_index1-1,ii) .EQ. 3.5)) THEN
                    temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
                 END IF
              END DO
           END IF

           ! fracture temperature set by vertical heat advection after spinup phase
           IF ((iter .EQ. 1)) THEN

              ! 		temp6_mid = temp6
              !
              DO ii=2,yn-1
                 IF ((mask(f_index1-1,ii) .EQ. 3.05)) THEN
                    temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
                 END IF
              END DO

              DO ii=2,yn-1
                 IF ((mask(f_index1-1,ii) .EQ. 3.1) .OR. (mask(f_index1-1,ii) .EQ. 3.0) .OR. (mask(f_index1-1,ii) .EQ. 3.5)) THEN
                    temp6(ii,1) = temp6(ii-1,1) - (dy*lambdaMat(f_index1-1,ii)/(dx*4179.0*frac6(ii,1))) * (h(f_index1-1,ii) - h(f_index1-2,ii))
                 END IF
              END DO

              ! 				temp6 = temp6_mid
              !
              ! 				do ii=2,yn-1
              ! 					if ((mask(f_index1-1,ii) .eq. 3.1)) then
              ! 						temp6_a(ii) = 0.0
              ! 						temp6_b(ii) = 1.0 + frac6(ii,1)
              ! 					end if
              ! 				end do
              !

              ! 				if (mod(j,50) .eq. 0) then
              !
              ! 				do jj=1,10000000
              !
              ! 					temp6 = temp6_mid
              ! 					do ii=2,yn-1
              !
              ! 						if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
              ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
              ! 						end if
              ! ! 						if (mask(f_index1-1,ii) .eq. 3.0) then
              ! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
              ! ! 						end if
              ! ! 						if (mask(f_index1-1,ii) .eq. 3.5) then
              ! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
              ! ! 						end if
              ! 					end do
              !
              ! 				end do
              ! ! 				write(*,*) "mod 50"
              ! ! 				write(*,*) temp6(:,1)
              ! ! 				write(*,*) " "
              !
              ! ! ! TESTING SOMETHING
              ! ! do ii=2,yn-1
              ! ! 	if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
              ! ! 		temp6_mid(ii,1) = 328.0
              ! ! 	end if
              ! ! end do
              !
              ! 				end if



              ! top of fracture is no heat flow out? wilcox condition?
              DO ii=yn/2,yn-1
                 IF ((mask(f_index1,ii) .EQ. 6.5)) THEN
                    temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
                 END IF

                 IF ((mask(f_index1-1,ii) .EQ. 3.5)) THEN
                    temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
                 END IF
              END DO

           END IF


           DO ii=1,yn-1
              temp6(ii,2) = temp6(ii,1)
           END DO



           rho = rho_next(h)
           visc = visc_next(h)




           ! solve streamfunction-vorticity equation

           rhs0 = -1.0*partial((rho/rho_fluid)-1.0,xn,yn,dx,dy,1)

           frac6_last = frac6

           ! find properties at top and base of fracture
           DO jj=yn/2,yn-1
              IF (maskP(f_index1-1,jj) .EQ. 3.1) THEN
                 h_base = temp6(jj,1) ! h(f_index1-1,jj)!
                 y_base = y(jj)
                 jj_base = jj
              END IF
              IF (maskP(f_index1-1,jj) .EQ. 3.5) THEN
                 h_top = param_tsw ! temp6(jj,1)!
                 y_top = y(jj)
                 jj_top = jj
              END IF
           END DO

           h_adjacent = SUM(temp6(jj_base:jj_top,1))/(jj_top-jj_base+1)

           IF ((j .GT. spinup-1)) THEN
              iter = 1
           END IF

           ! 		if ((j .eq. spinup-1)) then
           ! 			h_adjacent = 273.0
           ! 		end if


           IF (iter .EQ. 1) THEN
              !write(*,*) "loop of fracj 1"
              DO jj=yn/2,yn-1
                 IF ((maskP(f_index1,jj) .EQ. 6.0) .OR. (maskP(f_index1,jj) .EQ. 6.5) .OR. (maskP(f_index1,jj) .EQ. 6.1)) THEN
                    !frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_fluid*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*y_top - rho_one(h_base)*y_base - rho_fluid*(y_top-y_base))
                    frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_one(h_base)*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*(y_top-y_top) - rho_one(h_base)*(y_base-y_top) - rho_one(param_tsw)*(y_top-y_base))
                    !frac6(jj,1) = param_f_por*psi(f_index1-1,jj) + (dx*(12.0*(1.0e-16)*1.0)/(param_f_dx*param_f_dx*param_f_dx))
                 END IF
              END DO
              !write(*,*) " "

              psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
              psi = psi_bc(psi)

              DO jj=yn/2,yn-1
                 DO i=1,xn
                    IF ((maskP(i,jj) .EQ. 50.0) .AND. (i .LT. f_index1)) THEN
                       psi(i,jj+1) = MAXVAL(frac6(:,1))
                    END IF
                    ! 					if ((maskP(i,jj) .eq. 3.5)) then
                    ! 						psi(i,jj+1) = maxval(frac6(:,1))
                    ! 					end if
                 END DO
              END DO


           END IF


           ! run with frac6 = 0 during spinup phase
           IF (iter .EQ. 0) THEN

              psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
              psi = psi_bc(psi)

           END IF



        END IF ! end if restart .ne. 1

        !write(*,*) "error string 2"

        IF (j .EQ. 3) THEN

           !-INITIALIZE FOR GEOCHEM

           leng = ((yn/(2*celly))*((xn-1)/cellx))

           ! get velocities from streamfunction

           DO ii = yn/(2*celly)+1,yn/(celly)
              DO i = 1,(xn-1)/cellx
                 h_coarse(i,ii-yn/(2*celly)) = SUM(h((i-1)*cellx+1:i*cellx,(ii-1)*celly+1:ii*celly))/(cellx*celly)
                 psi_coarse(i,ii-yn/(2*celly)) = psi(i*cellx,ii*celly)
              END DO
           END DO

           DO i = 1,yn/(2*celly)
              h_coarse(:,i) = SUM(h_coarse(:,i))/((xn-1)/cellx)
           END DO



           !h_coarse((xn-1)/cellx,:) = h_coarse(((xn-1)/cellx)-1,:)

        !    DO ii = yn/(2*celly)+1,yn/(celly)
        !       DO i = 1,((xn-1)/cellx)/2
        !           h_coarse(2*i-1:2*i,ii-yn/(2*celly)) = SUM(h((i-1)*cellx*2+1:i*cellx*2,(ii-1)*celly+1:ii*celly))/(2*cellx*celly)
        !       END DO
        !    END DO

           !h_coarse = 360.0


           hLong = (/ RESHAPE(TRANSPOSE(h_coarse), (/ leng /)), RESHAPE(TRANSPOSE(h_coarse), (/ leng /)), RESHAPE(TRANSPOSE(h_coarse), (/ leng /)) /) ! for cell = 1
        !hLong = (/ RESHAPE(h_coarse, (/ leng /)), RESHAPE(h_coarse, (/ leng /)), RESHAPE(h_coarse, (/ leng /)) /)

           ! 			psi_coarse(2,:) = psi_coarse(1,:)
           psi_coarse((xn-1)/cellx,:) = 0.0
        !    psi_coarse((xn-1)/cellx-1,:) = 0.0
        !    psi_coarse((xn-1)/cellx-2,:) = 0.0

           velocities_coarse0 = velocities_coarse(psi_coarse)
        !    u_coarse = phi_coarse*velocities_coarse0(1:(xn-1)/cellx,1:yn/(2*celly))/(rho_fluid)
        !    v_coarse = phi_coarse*velocities_coarse0(1:(xn-1)/cellx,yn/(2*celly)+1:2*yn/(2*celly))/(rho_fluid)
        u_coarse = velocities_coarse0(1:(xn-1)/cellx,1:yn/(2*celly))/(rho_fluid)
        v_coarse = velocities_coarse0(1:(xn-1)/cellx,yn/(2*celly)+1:2*yn/(2*celly))/(rho_fluid)

           ! EQUALIZE U_COARSE

           DO i = 1,yn/(2*celly)
              u_coarse(:,i) = SUM(u_coarse(:,i))/((xn-1)/cellx)
           END DO
           !u_coarse((xn-1)/cellx,:) = 0.0
        !    u_coarse((xn-1)/cellx-1,:) = 0.0
        !    u_coarse((xn-1)/cellx-2,:) = 0.0

        ! SAME VELOCITY ON ALL 3 ROWS
        DO i = 3,yn/(2*celly)-1
            if ((u_coarse(5,i) .gt. 0.0) .and. (u_coarse(5,i+1) .eq. 0.0)) then
            u_coarse(:,i-1) = u_coarse(5,i)
            u_coarse(:,i-2) = u_coarse(5,i)
            u_coarse(:,i-3) = u_coarse(5,i)
            u_coarse(:,i-4) = u_coarse(5,i)
            u_coarse(:,i-5) = u_coarse(5,i)
            u_coarse(:,i-6) = u_coarse(5,i)
            end if
        END DO
        !u_coarse = u_coarse/1.5

           velocities0 = velocities(psi)
           u = phi*velocities0(1:xn,1:yn)/(rho_fluid)
           v = phi*velocities0(1:xn,yn+1:2*yn)/(rho_fluid)

           u_coarse_long = RESHAPE(u_coarse, (/ leng /))
           v_coarse_long = RESHAPE(v_coarse, (/ leng /))
           phi_coarse_long = RESHAPE(phi_coarse, (/ leng /))

           ! 			u(f_index1-1:,:) = 0.0
           ! 			v(f_index1-1:,:) = 0.0



           WRITE(*,*) "BEGIN INITIAL STRETCHING EVERYTHING OUT FOR GEOCHEM"


           ! stretch everything out
           !hLong = reshape(h(1:xn-1:cell,1:yn-1:cell), (/(xn/cell)*(yn/cell)/)) ! for cell > 1
           DO i = 1,g_pri
              priLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(primary(:,:,i)), (/ leng /))
              priLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(primary_a(:,:,i)), (/ leng /))
              priLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(primary_b(:,:,i)), (/ leng /))
            !   priLongBitFull(:leng,i) = RESHAPE(primary(:,:,i), (/ leng /))
            !   priLongBitFull(leng+1:2*leng,i) = RESHAPE(primary_a(:,:,i), (/ leng /))
            !   priLongBitFull(2*leng+1:,i) = RESHAPE(primary_b(:,:,i), (/ leng /))
           END DO

           DO i = 1,g_sec/2
              secLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(secondary(:,:,i)), (/ leng /))
              secLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(secondary_a(:,:,i)), (/ leng /))
              secLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(secondary_b(:,:,i)), (/ leng /))
            ! secLongBitFull(:leng,i) = RESHAPE(secondary(:,:,i), (/ leng /))
            ! secLongBitFull(leng+1:2*leng,i) = RESHAPE(secondary_a(:,:,i), (/ leng /))
            ! secLongBitFull(2*leng+1:,i) = RESHAPE(secondary_b(:,:,i), (/ leng /))
           END DO

           DO i = 1,g_sol
              solLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(solute(:,:,i)), (/ leng /))
              solLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(solute_a(:,:,i)), (/ leng /))
              solLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(solute_b(:,:,i)), (/ leng /))
            ! solLongBitFull(:leng,i) = RESHAPE(solute(:,:,i), (/ leng /))
            ! solLongBitFull(leng+1:2*leng,i) = RESHAPE(solute_a(:,:,i), (/ leng /))
            ! solLongBitFull(2*leng+1:,i) = RESHAPE(solute_b(:,:,i), (/ leng /))
           END DO

           DO i = 1,g_med
              medLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(medium(:,:,i)), (/ leng /))
              medLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(medium_a(:,:,i)), (/ leng /))
              medLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(medium_b(:,:,i)), (/ leng /))
            ! medLongBitFull(:leng,i) = RESHAPE(medium(:,:,i), (/ leng /))
            ! medLongBitFull(leng+1:2*leng,i) = RESHAPE(medium_a(:,:,i), (/ leng /))
            ! medLongBitFull(2*leng+1:,i) = RESHAPE(medium_b(:,:,i), (/ leng /))
           END DO






        END IF ! if j == 5

        ! 			do ii=1,yn
        ! 				do i=2,xn-2
        ! 					if ((maskP(i,ii) .eq. 6.0)) then
        ! 						v(i,ii) = -1.0*(frac6(ii,2) - frac6(ii,1))/(rho(i,ii)*dx)
        ! 					end if
        ! 				end do
        ! 			end do


        ! 	write(*,*) "about to start mstep loop"

	! things only done every mth timestep go here
	IF (MOD(j,mstep) .EQ. 0) THEN

           WRITE(*,*) " "
           WRITE(*,*) "STEP" , j , "STUFF"



           ! 		! OUTER BAND THING
           ! 		outerBand = make_band(permeability,phi)
           ! 		outerBand = band(outerBand,2*(yn-2) + 1,(xn-2)*(yn-2))
           ! 		permx = partial_edge((phi/(permeability)),xn,yn,dx,dy,1)
           ! 		permy = partial_edge((phi/(permeability)),xn,yn,dx,dy,2)

           ! make coarse grid average velocities
           ! 		uTransport = (uTransport/(mstep*wscale))
           ! 		vTransport = (vTransport/(mstep*wscale))
           !uTransport(1,:) = 0.0
           !uTransport(xn/cell,:) = 0.0

           ! 		write(*,*) "SOLUTE COURANT NUMBERS"
           ! 		!write(*,*) (dt*mstep/(cstep*dx*cell))*maxval(abs(uTransport))
           ! 		!write(*,*) (dt*mstep/(cstep*dy*cell))*maxval(abs(vTransport))
           ! 		write(*,*) dt*mstep/(cstep*dy*dy*cell*cell)

           !NOT RUNNING TRANSPORT RIGHT NOW, JUST CHEMISTRY
           ! if (mod(j,tn/res_factor) .ne. 0) then
           ! 	medium(:,:,2) = 1000.0
           ! 	medium_a(:,:,2) = 1000.0
           ! 	medium_b(:,:,2) = 1000.0
           ! end if
           !
           ! if (mod(j,tn/res_factor) .eq. 0) then
           medium(:,:,2) = 0.0
           medium_a(:,:,2) = 0.0
           medium_b(:,:,2) = 0.0
           ! end if





           !-phi_coarse_calc_long (s + a)
           DO i = 1,leng
              phi_calc_denom = 0.0
              phi_calc_denom = solLongBitFull(i,3)*1000.0
              phi_calc_denom = phi_calc_denom + (priLongBitFull(i,2)*pri_molar(2)/pri_density(2)) + (priLongBitFull(i,3)*pri_molar(3)/pri_density(3))
              phi_calc_denom = phi_calc_denom + (priLongBitFull(i,4)*pri_molar(4)/pri_density(4)) + (priLongBitFull(i,5)*pri_molar(5)/pri_density(5))
              DO ii=1,g_sec/2
                 phi_calc_denom = phi_calc_denom + (secLongBitFull(i,ii)*sec_molar(ii)/sec_density(ii))
              END DO
              phi_coarse_long(i) = solLongBitFull(i,3)*1000.0 / phi_calc_denom


              ! stop at 5% porosity
              if (phi_coarse_long(i) .lt. 0.08) then
                  medLongBitFull(i,2) = precip_th!1000.0
              else
                  medLongBitFull(i,2) = precip_th
              end if

           END DO
           phiLongBitFull(:leng) = phi_coarse_long

           !WRITE(*,*) "done with phi_coarse_long (s)"

           DO i = 1,leng
              phi_calc_denom = 0.0
              phi_calc_denom = (solLongBitFull(leng+i,3)*1000.0) + (solLongBitFull(2*leng+i,3)*1000.0)
              phi_calc_denom = phi_calc_denom + (priLongBitFull(leng+i,2)*pri_molar(2)/pri_density(2)) + (priLongBitFull(leng+i,3)*pri_molar(3)/pri_density(3))
              phi_calc_denom = phi_calc_denom + (priLongBitFull(leng+i,4)*pri_molar(4)/pri_density(4)) + (priLongBitFull(leng+i,5)*pri_molar(5)/pri_density(5))
              ! adding porosity of b
              phi_calc_denom = phi_calc_denom + (priLongBitFull(2*leng+i,2)*pri_molar(2)/pri_density(2)) + (priLongBitFull(2*leng+i,3)*pri_molar(3)/pri_density(3))
              phi_calc_denom = phi_calc_denom + (priLongBitFull(2*leng+i,4)*pri_molar(4)/pri_density(4)) + (priLongBitFull(2*leng+i,5)*pri_molar(5)/pri_density(5))
              DO ii=1,g_sec/2
                 phi_calc_denom = phi_calc_denom + (secLongBitFull(leng+i,ii)*sec_molar(ii)/sec_density(ii))
                 phi_calc_denom = phi_calc_denom + (secLongBitFull(2*leng+i,ii)*sec_molar(ii)/sec_density(ii))
              END DO
              phi_coarse_long(i) = ((solLongBitFull(leng+i,3)*1000.0)+(solLongBitFull(2*leng+i,3)*1000.0)) / phi_calc_denom

              ! stop at 5% porosity
              if (phi_coarse_long(i) .lt. 0.08) then
                  medLongBitFull(leng+i,2) = precip_th!1000.0
              else
                  medLongBitFull(leng+i,2) = precip_th
              end if
           END DO
           phiLongBitFull(leng+1:2*leng) = phi_coarse_long


           !WRITE(*,*) "done with phi_coarse_long (a)"

           call system_clock(counti, count_rate, count_max)
           !write(*,*) "about to send advection..."

           DO an_id = 1, 11

              !IF (an_id .LE. 11) THEN

                                ! aug 7th wrong element glitch
                 			! 	do i=1,leng
                 			! 		if (coarse_mask_long(i) .eq. 0.0) then
                 			! 			solLongBitFull(i,sol_index(an_id)) = sea(sol_index(an_id))
                 			! 		end if
                 			! 	end do
                 !
                 sol_coarse_long = solLongBitFull(:leng,sol_index(an_id))
                 phi_coarse_long = phiLongBitFull(:leng)

            !    ! reshape them all
            !    sol_coarse_local = TRANSPOSE(RESHAPE(sol_coarse_long,(/yn/(2*celly),(xn-1)/cellx/)))
            !    u_coarse_local = RESHAPE(u_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
            !    v_coarse_local = RESHAPE(v_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
            !    phi_coarse_local = TRANSPOSE(RESHAPE(phi_coarse_long,(/yn/(2*celly),(xn-1)/cellx/)))
            ! ! sol_coarse_local = RESHAPE(sol_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
            ! ! u_coarse_local = RESHAPE(u_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
            ! ! v_coarse_local = RESHAPE(v_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
            ! ! phi_coarse_local = RESHAPE(phi_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))

            !    DO ii = 1,cstep
            !      sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id)))
            !   END DO
              !
            !   sol_coarse_long_local = RESHAPE(TRANSPOSE(sol_coarse_local),(/((xn-1)/cellx)*(yn/(2*celly))/))
            !   !sol_coarse_long_local = RESHAPE(sol_coarse_local,(/((xn-1)/cellx)*(yn/(2*celly))/))
              !
            !   IF (an_id .LE. 11) THEN
            !      solLongBitFull(:leng,sol_index(an_id)) = sol_coarse_long_local
            !   END IF




                 !#ADVECTION: send from master to slaves (s)

                 ! send an_id name
                 CALL MPI_SEND( an_id, 1, MPI_INTEGER, &
                      an_id, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long sol coarse
                 CALL MPI_SEND( sol_coarse_long, leng, MPI_REAL4, &
                      an_id, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long u coarse
                 CALL MPI_SEND( u_coarse_long, leng, MPI_REAL4, &
                      an_id, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long v coarse
                 CALL MPI_SEND( v_coarse_long, leng, MPI_REAL4, &
                      an_id, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long phi coarse
                 CALL MPI_SEND( phi_coarse_long, leng, MPI_REAL4, &
                      an_id, send_data_tag, MPI_COMM_WORLD, ierr)


              !END IF

               !write(*,*) "sent to first set" , an_id

            END DO



            DO an_id = 1, 11

               CALL MPI_RECV( sol_coarse_long, leng, MPI_REAL4, &
                    an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

               IF (an_id .LE. 11) THEN
                  solLongBitFull(:leng,sol_index(an_id)) = sol_coarse_long
               END IF

            !    IF (an_id .GT. 11) THEN
            !       solLongBitFull(leng+1:2*leng,sol_index(an_id-11)) = sol_coarse_long
            !    END IF

            END DO



            DO an_id = 12, 22

              !IF (an_id .GT. 11) THEN
                                ! aug 7th wrong element glitch
                 			! 	do i=1,leng
                 			! 		if (coarse_mask_long(i) .eq. 0.0) then
                 			! 			solLongBitFull(leng+1+i,sol_index(an_id-11)) = sea(sol_index(an_id-11))
                 			! 		end if
                 			! 	end do

                 sol_coarse_long = solLongBitFull(leng+1:2*leng,sol_index(an_id-11))
                 phi_coarse_long = phiLongBitFull(leng+1:2*leng)


                !  ! reshape them all
                !  sol_coarse_local = TRANSPOSE(RESHAPE(sol_coarse_long,(/yn/(2*celly),(xn-1)/cellx/)))
                !  u_coarse_local = RESHAPE(u_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
                !  v_coarse_local = RESHAPE(v_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
                !  phi_coarse_local = TRANSPOSE(RESHAPE(phi_coarse_long,(/yn/(2*celly),(xn-1)/cellx/)))
                ! ! sol_coarse_local = RESHAPE(sol_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
                ! ! u_coarse_local = RESHAPE(u_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
                ! ! v_coarse_local = RESHAPE(v_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))
                ! ! phi_coarse_local = RESHAPE(phi_coarse_long,(/(xn-1)/cellx,yn/(2*celly)/))


            !      DO ii = 1,cstep
            !      sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id-11)))
            !  END DO
            !
            !   sol_coarse_long_local = RESHAPE(TRANSPOSE(sol_coarse_local),(/((xn-1)/cellx)*(yn/(2*celly))/))
            ! !sol_coarse_long_local = RESHAPE(sol_coarse_local,(/((xn-1)/cellx)*(yn/(2*celly))/))
            !
            !
            !   IF (an_id .GT. 11) THEN
            !      solLongBitFull(leng+1:2*leng,sol_index(an_id-11)) = sol_coarse_long_local
            !   END IF



                !  !#ADVECTION: send from master to slaves (a)
                 !
                 ! send an_id name
                 CALL MPI_SEND( an_id-11, 1, MPI_INTEGER, &
                      an_id-11, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long sol coarse
                 CALL MPI_SEND( sol_coarse_long, leng, MPI_REAL4, &
                      an_id-11, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long u coarse
                 CALL MPI_SEND( u_coarse_long, leng, MPI_REAL4, &
                      an_id-11, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long v coarse
                 CALL MPI_SEND( v_coarse_long, leng, MPI_REAL4, &
                      an_id-11, send_data_tag, MPI_COMM_WORLD, ierr)

                 ! send long phi coarse
                 CALL MPI_SEND( phi_coarse_long, leng, MPI_REAL4, &
                      an_id-11, send_data_tag, MPI_COMM_WORLD, ierr)


              !END IF

               !write(*,*) "sent to second set" , an_id

              !write(*,*) "DONE SENDING SOLUTE TO PROCESSOR", an_id

              END DO


              DO an_id = 12, 22

                 CALL MPI_RECV( sol_coarse_long, leng, MPI_REAL4, &
                      an_id-11, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

                !  IF (an_id .LE. 11) THEN
                !     solLongBitFull(:leng,sol_index(an_id)) = sol_coarse_long
                !  END IF

                 IF (an_id .GT. 11) THEN
                    solLongBitFull(leng+1:2*leng,sol_index(an_id-11)) = sol_coarse_long
                 END IF

              END DO

           !#ADVECTION: master receives from slaves

           ! call system_clock(counti, count_rate, count_max)
           !write(*,*) "BEGIN RECEIVING ADVECTED SOLUTES"

        !    DO an_id = 1, 22
           !
        !       CALL MPI_RECV( sol_coarse_long, leng, MPI_REAL4, &
        !            an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           !
        !       IF (an_id .LE. 11) THEN
        !          solLongBitFull(:leng,sol_index(an_id)) = sol_coarse_long
        !       END IF
           !
        !       IF (an_id .GT. 11) THEN
        !          solLongBitFull(leng+1:2*leng,sol_index(an_id-11)) = sol_coarse_long
        !       END IF
           !
        !    END DO

           call system_clock(countf, count_rate, count_max)
           write(*,*) "advection time: " , countf - counti

           !write(*,*) "END RECEIVING ADVECTED SOLUTES"

           ! old chamber mixing
           ! solute_inter_long = solLongBitFull(leng+1:2*leng,n)
           ! solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volume_ratio) + solLongBitFull(2*leng+1:,n)*mix_ratio/volume_ratio ! a mix
           ! solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio
           !
           ! do n=4,13 ! solutes
           ! 	solute_inter_long = solLongBitFull(leng+1:2*leng,n)
           ! 	solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volume_ratio) + solLongBitFull(2*leng+1:,n)*mix_ratio/volume_ratio ! a mix
           ! 	solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio
           ! end do

           !-mixing between chambers
        !    DO i=1,leng
        !       volLongBitFull(i) = MAX(solLongBitFull(leng+i,3)/MAX(solLongBitFull(2*leng+i,3),1e-10),1e-10)
        !    END DO
        call system_clock(counti, count_rate, count_max)
        DO i=1,leng
           volLongBitFull(i) = solLongBitFull(leng+i,3)/solLongBitFull(2*leng+i,3)
        END DO
        !    WRITE(*,*) "made volLongBitFull"
        !    n=2 ! alk
        !    solute_inter_long = solLongBitFull(leng+1:2*leng,n)
        !    solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volLongBitFull) + solLongBitFull(2*leng+1:,n)*mix_ratio/volLongBitFull ! a mix
        !    solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio
           !
        !    DO n=4,13 ! solutes
        !       solute_inter_long = solLongBitFull(leng+1:2*leng,n)
        !       solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volLongBitFull) + solLongBitFull(2*leng+1:,n)*mix_ratio/volLongBitFull ! a mix
        !       solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio
        !    END DO

        ! !-new mixing between chambers
        ! DO i=1,leng
        !    volLongBitFull(i) = t_vol_a/t_vol_b
        ! END DO
        WRITE(*,*) "made volLongBitFull"
        n=2 ! alk
        solute_inter_long = solLongBitFull(leng+1:2*leng,n)
        solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volLongBitFull) + solLongBitFull(2*leng+1:,n)*mix_ratio/volLongBitFull
        solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio

        DO n=4,13 ! solutes
           solute_inter_long = solLongBitFull(leng+1:2*leng,n)
           solLongBitFull(leng+1:2*leng,n) = solLongBitFull(leng+1:2*leng,n)*(1.0-mix_ratio/volLongBitFull) + solLongBitFull(2*leng+1:,n)*mix_ratio/volLongBitFull
           solLongBitFull(2*leng+1:,n) = solLongBitFull(2*leng+1:,n)*(1.0-mix_ratio) + solute_inter_long*mix_ratio
        END DO
        call system_clock(countf, count_rate, count_max)
        write(*,*) "mixing time: " , countf - counti


           !-move aged cells
           num_rows = 3*((xn-1)/cellx)*(yn/(2*celly))


           age_toggle = 0
           ! if (t(j) .gt. 2.512e13) then
           IF (t(j) .LT. 8.478e13) THEN
              !WRITE(*,*) "past 0.8 myr"
              !for x cell blocks = 1000m
              ! if (floor((t(j-mstep)-2.512e13)/9.42e11) .lt. floor((t(j)-2.512e13)/9.42e11)) then

            !   ! for x cell blocks = 2000m
            !   IF (FLOOR((t(j-mstep)-2.512e13)/18.84e11) .LT. FLOOR((t(j)-2.512e13)/18.84e11)) THEN
            !      WRITE(*,*) "moving cells now..."

            ! for x cell blocks = 4500m (cellx = 90)
            ! stop moving after 2.7 ma, still for 0.8 ma
            IF (FLOOR((t(j-mstep+1)-2.512e13)/4.239e12) .LT. FLOOR((t(j+1)-2.512e13)/4.239e12)) THEN
            ! if (floor((t(j-1)-2.512e13)/4.239e12) .gt. floor((t(j)-2.512e13)/4.239e12)) then
               call system_clock(counti, count_rate, count_max)
               WRITE(*,*) "moving cells now..."

                 DO i = 1,g_pri
                    bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    priLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                    bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    priLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                    bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    priLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                    ! bit_thing_t1 = RESHAPE(priLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! priLongBitFull(:leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                    !
                    ! bit_thing_t1 = RESHAPE(priLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! priLongBitFull(leng+1:2*leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                    !
                    ! bit_thing_t1 = RESHAPE(priLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! priLongBitFull(2*leng+1:,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                 END DO

                 DO i = 1,g_sec
                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    secLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    secLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    secLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                    ! bit_thing_t1 = RESHAPE(secLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! secLongBitFull(:leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                    !
                    ! bit_thing_t1 = RESHAPE(secLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! secLongBitFull(leng+1:2*leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                    !
                    ! bit_thing_t1 = RESHAPE(secLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                    ! secLongBitFull(2*leng+1:,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                 END DO

                !  bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(1:leng,5),(/yn/(2*celly), (xn-1)/cellx/)))
                !  bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !  medLongBitFull(:leng,5) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                 !
                !  bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(leng+1:2*leng,5),(/yn/(2*celly), (xn-1)/cellx/)))
                !  bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !  medLongBitFull(leng+1:2*leng,5) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                 !
                !  bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(2*leng+1:,5),(/yn/(2*celly), (xn-1)/cellx/)))
                !  bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !  medLongBitFull(2*leng+1:,5) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                !  DO i = 1,g_sol
                !     bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                !     bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     solLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                 !
                !     bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                !     bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     solLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                 !
                !     bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                !     bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     solLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                !     ! bit_thing_t1 = RESHAPE(solLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                !     ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     ! solLongBitFull(:leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                !     !
                !     ! bit_thing_t1 = RESHAPE(solLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                !     ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     ! solLongBitFull(leng+1:2*leng,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                !     !
                !     ! bit_thing_t1 = RESHAPE(solLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                !     ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                !     ! solLongBitFull(2*leng+1:,i) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                !  END DO


                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(1:leng,2),(/yn/(2*celly), (xn-1)/cellx/)))
                 bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                 medLongBitFull(:leng,2) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(leng+1:2*leng,2),(/yn/(2*celly), (xn-1)/cellx/)))
                 bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                 medLongBitFull(leng+1:2*leng,2) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))

                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(2*leng+1:,2),(/yn/(2*celly), (xn-1)/cellx/)))
                 bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                 medLongBitFull(2*leng+1:,2) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
                ! bit_thing_t1 = RESHAPE(medLongBitFull(1:leng,2),(/(xn-1)/cellx,yn/(2*celly)/))
                ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                ! medLongBitFull(:leng,2) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                !
                ! bit_thing_t1 = RESHAPE(medLongBitFull(leng+1:2*leng,2),(/(xn-1)/cellx,yn/(2*celly)/))
                ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                ! medLongBitFull(leng+1:2*leng,2) = RESHAPE(bit_thing_t1(:,:), (/ leng /))
                !
                ! bit_thing_t1 = RESHAPE(medLongBitFull(2*leng+1:,2),(/(xn-1)/cellx,yn/(2*celly)/))
                ! bit_thing_t1(2:,:) = bit_thing_t1(:(xn-1)/cellx-1,:)
                ! medLongBitFull(2*leng+1:,2) = RESHAPE(bit_thing_t1(:,:), (/ leng /))

                call system_clock(countf, count_rate, count_max)
                write(*,*) "moving time: " , countf - counti
                age_toggle = 1



              END IF
           END IF


        !    !-adjust RHS solutes (avoid drift?)
        !    DO i = 1,g_sol
           !
        !       bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
        !       !DO ii = (xn-1)/cellx-2,(xn-1)/cellx
        !          ii = (xn-1)/cellx
        !          bit_thing_t1(ii,:) = bit_thing_t1((xn-1)/cellx-1,:)
        !       !END DO
        !       solLongBitFull(:leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
           !
           !
        !       bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
        !       !DO ii = (xn-1)/cellx-2,(xn-1)/cellx
        !          ii = (xn-1)/cellx
        !          bit_thing_t1(ii,:) = bit_thing_t1((xn-1)/cellx-1,:)
        !       !END DO
        !       solLongBitFull(leng+1:2*leng,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
           !
           !
        !       bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
        !       !DO ii = (xn-1)/cellx-2,(xn-1)/cellx
        !          ii = (xn-1)/cellx
        !          bit_thing_t1(ii,:) = bit_thing_t1((xn-1)/cellx-1,:)
        !       !END DO
        !       solLongBitFull(2*leng+1:,i) = RESHAPE(TRANSPOSE(bit_thing_t1(:,:)), (/ leng /))
           !
        !    END DO


        ! calculate ph_fix (entire domain)

        ph_fix_LongBitFull = 0.0


        ph_fix_LongBitFull(:leng) = RESHAPE(TRANSPOSE(ph_fix(:,:)), (/ leng /))
        ph_fix_LongBitFull(leng+1:2*leng) = RESHAPE(TRANSPOSE(ph_fix_a(:,:)), (/ leng /))
        ph_fix_LongBitFull(2*leng+1:) = RESHAPE(TRANSPOSE(ph_fix_b(:,:)), (/ leng /))
        ! ph_fix_LongBitFull(:leng) = RESHAPE(ph_fix(:,:), (/ leng /))
        ! ph_fix_LongBitFull(leng+1:2*leng) = RESHAPE(ph_fix_a(:,:), (/ leng /))
        ! ph_fix_LongBitFull(2*leng+1:) = RESHAPE(ph_fix_b(:,:), (/ leng /))





        !#GEOCHEM: send from master to slaves
        DO an_id = 1 , num_procs - 1
            CALL MPI_SEND( j, 1, MPI_INTEGER, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            !CALL MPI_SEND( age_toggle, 1, MPI_INTEGER, an_id, send_data_tag, MPI_COMM_WORLD, ierr)

            CALL MPI_SEND( hLong, 3*leng, MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)

            DO ii = 1 , g_pri
                CALL MPI_SEND( priLongBitFull(:,ii), 3*leng, MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            END DO

            DO ii = 1 , g_sec/2
                CALL MPI_SEND( secLongBitFull(:,ii), 3*leng, MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            END DO

            DO ii = 1 , g_sol
                CALL MPI_SEND( solLongBitFull(:,ii), 3*leng, MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            END DO

            DO ii = 1 , g_med
                CALL MPI_SEND( medLongBitFull(:,ii), 3*leng, MPI_REAL4, an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            END DO

        END DO




        !#GEOCHEM: receive from slaves

        DO an_id = 1,num_procs - 1
            call system_clock(counti, count_rate, count_max)

            CALL MPI_RECV( end_loop, 1, MPI_INTEGER, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            CALL MPI_RECV( slave_vector, end_loop, MPI_INTEGER, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)


                DO ii = 1 , g_pri
                    CALL MPI_RECV( priLocal(1,ii), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                    priLongBitFull(slave_vector(1:end_loop),ii) = priLocal(1:end_loop,ii)
                END DO

                DO ii = 1 , g_sec/2
                    CALL MPI_RECV( secLocal(1,ii), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                    secLongBitFull(slave_vector(1:end_loop),ii) = secLocal(1:end_loop,ii)
                END DO

                ! CALL MPI_RECV( solLocal(1:end_loop,3), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                !
                ! CALL MPI_RECV( solLocal(1:end_loop,2), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                ! solLongBitFull(slave_vector(1:end_loop),2) = solLocal(1:end_loop,2)*solLocal(1:end_loop,3)/(solLongBitFull(slave_vector(1:end_loop),3))
                !
                ! CALL MPI_RECV( solLocal(1:end_loop,1), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                ! solLongBitFull(slave_vector(1:end_loop),1) = -1.0*LOG10(10.0**(-1.0*solLocal(1:end_loop,1))*solLocal(1:end_loop,3)/(solLongBitFull(slave_vector(1:end_loop),3)))
                !
                ! DO ii = 4 , g_sol
                !     CALL MPI_RECV( solLocal(1:end_loop,ii), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                !     solLongBitFull(slave_vector(1:end_loop),ii) = solLocal(1:end_loop,ii)*solLocal(1:end_loop,3)/(solLongBitFull(slave_vector(1:end_loop),3))
                ! END DO

                DO ii = 1 , g_sol
                    CALL MPI_RECV( solLocal(1:end_loop,ii), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                    solLongBitFull(slave_vector(1:end_loop),ii) = solLocal(1:end_loop,ii)
                END DO



                DO ii = 1 , g_med
                    CALL MPI_RECV( medLocal(1,ii), end_loop, MPI_REAL4, an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                    medLongBitFull(slave_vector(1:end_loop),ii) = medLocal(1:end_loop,ii)
                END DO

                call system_clock(countf, count_rate, count_max)
                write(*,*) "receive from" , an_id , countf - counti


        END DO






           write(*,*) "done receiving geochem from slaves"



           !-water volume correction here

           ! 		medium(:,:,3) = vol_i
           ! 		medium_a(:,:,3) = vol_i_a
           ! 		medium_b(:,:,3) = vol_i_b
           !
           ! 		solute(:,:,3) = vol_i
           ! 		solute_a(:,:,3) = vol_i_a
           ! 		solute_b(:,:,3) = vol_i_b

           ! 		i=2
           ! 		solute(:,:,i) = solute(:,:,i)*solute(:,:,3)/vol_i
           ! 		solute_a(:,:,i) = solute_a(:,:,i)*solute_a(:,:,3)/vol_i_a
           ! 		solute_b(:,:,i) = solute_b(:,:,i)*solute_b(:,:,3)/vol_i_b
           !
           ! 		do i=4,13
           ! 			solute(:,:,i) = solute(:,:,i)*solute(:,:,3)/vol_i
           ! 			solute_a(:,:,i) = solute_a(:,:,i)*solute_a(:,:,3)/vol_i_a
           ! 			solute_b(:,:,i) = solute_b(:,:,i)*solute_b(:,:,3)/vol_i_b
           ! 		end do


           !-add to full output arrays
           IF (MOD(j,mstep*ar) .EQ. 0) THEN

               call system_clock(counti, count_rate, count_max)

              !leng = (yn/(2*celly))*((xn-1)/cellx)

              DO i = 1,g_pri
                 bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 primary(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 primary_a(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(priLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 primary_b(:,:,i) = bit_thing_t1
                ! bit_thing_t1 = RESHAPE(priLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! primary(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(priLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! primary_a(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(priLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! primary_b(:,:,i) = bit_thing_t1
              END DO

              IF (MAXVAL(medium(:,:,2)) .EQ. 0.0) THEN
                 DO i = 1,g_sec

                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    secondary(:,:,i) = bit_thing_t1

                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    secondary_a(:,:,i) = bit_thing_t1

                    bit_thing_t1 = TRANSPOSE(RESHAPE(secLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                    secondary_b(:,:,i) = bit_thing_t1

                    ! bit_thing_t1 = RESHAPE(secLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! secondary(:,:,i) = bit_thing_t1
                    !
                    ! bit_thing_t1 = RESHAPE(secLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! secondary_a(:,:,i) = bit_thing_t1
                    !
                    ! bit_thing_t1 = RESHAPE(secLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                    ! secondary_b(:,:,i) = bit_thing_t1
                 END DO
              END IF

              DO i = 1,g_sol

                 bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 solute(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 solute_a(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(solLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 solute_b(:,:,i) = bit_thing_t1
                ! bit_thing_t1 = RESHAPE(solLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! solute(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(solLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! solute_a(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(solLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! solute_b(:,:,i) = bit_thing_t1

              END DO

              DO i = 1,g_med

                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(1:leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 medium(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(leng+1:2*leng,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 medium_a(:,:,i) = bit_thing_t1

                 bit_thing_t1 = TRANSPOSE(RESHAPE(medLongBitFull(2*leng+1:,i),(/yn/(2*celly), (xn-1)/cellx/)))
                 medium_b(:,:,i) = bit_thing_t1
                ! bit_thing_t1 = RESHAPE(medLongBitFull(1:leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! medium(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(medLongBitFull(leng+1:2*leng,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! medium_a(:,:,i) = bit_thing_t1
                !
                ! bit_thing_t1 = RESHAPE(medLongBitFull(2*leng+1:,i),(/(xn-1)/cellx,yn/(2*celly)/))
                ! medium_b(:,:,i) = bit_thing_t1

              END DO

            !   ! phi_calc
              bit_thing_t1 = TRANSPOSE(RESHAPE(phiLongBitFull(1:leng),(/yn/(2*celly), (xn-1)/cellx/)))
              phiCalc(:,:) = bit_thing_t1
              bit_thing_t1 = TRANSPOSE(RESHAPE(phiLongBitFull(leng+1:2*leng),(/yn/(2*celly), (xn-1)/cellx/)))
              phiCalc_a(:,:) = bit_thing_t1
            ! phi_calc
            ! bit_thing_t1 = RESHAPE(phiLongBitFull(1:leng),(/(xn-1)/cellx,yn/(2*celly)/))
            ! phiCalc(:,:) = bit_thing_t1
            ! bit_thing_t1 = RESHAPE(phiLongBitFull(leng+1:2*leng),(/(xn-1)/cellx,yn/(2*celly)/))
            ! phiCalc_a(:,:) = bit_thing_t1


              ! pH fix
              bit_thing_t1 = TRANSPOSE(RESHAPE(ph_fix_LongBitFull(1:leng),(/yn/(2*celly), (xn-1)/cellx/)))
              ph_fix(:,:) = bit_thing_t1
              bit_thing_t1 = TRANSPOSE(RESHAPE(ph_fix_LongBitFull(leng+1:2*leng),(/yn/(2*celly), (xn-1)/cellx/)))
              ph_fix_a(:,:) = bit_thing_t1
              bit_thing_t1 = TRANSPOSE(RESHAPE(ph_fix_LongBitFull(2*leng+1:),(/yn/(2*celly), (xn-1)/cellx/)))
              ph_fix_b(:,:) = bit_thing_t1
            ! bit_thing_t1 = RESHAPE(ph_fix_LongBitFull(1:leng),(/(xn-1)/cellx,yn/(2*celly)/))
            ! ph_fix(:,:) = bit_thing_t1
            ! bit_thing_t1 = RESHAPE(ph_fix_LongBitFull(leng+1:2*leng),(/(xn-1)/cellx,yn/(2*celly)/))
            ! ph_fix_a(:,:) = bit_thing_t1
            ! bit_thing_t1 = RESHAPE(ph_fix_LongBitFull(2*leng+1:),(/(xn-1)/cellx,yn/(2*celly)/))
            ! ph_fix_b(:,:) = bit_thing_t1

            call system_clock(countf, count_rate, count_max)
            write(*,*) "1d to 2d matrix (ar) time: " , countf - counti


              WRITE(*,*) "BEGIN UPDATING _MAT ARRAYS"
              rhsmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rhs0
              rhomat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rho
              hmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = h
              psimat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = psi
              umat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = u
              vmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = v
              permmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permeability
              permxMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permx
              permyMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permy


              psiCoarseMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = psi_coarse
              ! 			 uCoarseMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = uTransport
              ! 			 vCoarseMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = vTransport
              ! reset coarse grid velocities for next timestep
              uTransport = 0.0
              vTransport = 0.0

              primaryMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = primary
              secondaryMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:)= secondary
              soluteMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = solute
              mediumMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = medium
              saturationMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = saturation

              primaryMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = primary_a
              secondaryMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:)= secondary_a
              soluteMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = solute_a
              mediumMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = medium_a

              primaryMat_b(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = primary_b
              secondaryMat_b(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:)= secondary_b
              soluteMat_b(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = solute_b
              mediumMat_b(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = medium_b

              primaryMat_d(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = primary_a + primary_b
              secondaryMat_d(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = secondary_a + secondary_b
              soluteMat_d(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly),:) = solute_a*(volume_ratio/(1.0+volume_ratio)) + solute_b*(1.0/(1.0+volume_ratio))

              phiCalcMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = phiCalc
              phiCalcMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = phiCalc_a

              ph_fixMat(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = ph_fix
              ph_fixMat_a(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = ph_fix_a
              ph_fixMat_b(1+((xn-1)/cellx)*(j/(mstep*ar)-1):((xn-1)/cellx)*(j/(mstep*ar)),1:yn/(2*celly)) = ph_fix_b


              WRITE(*,*) "...DONE UPDATING _MAT ARRAYS"


              ! 		 ! get new porosity
              phi_coarse = 0.1 !medium(:,:,1) ! 1.0


              ! only write to file 5 times total
              IF (MOD(j,tn/write_factor) .EQ. 0) THEN

                 WRITE(*,*) "BEGIN WRITING TO FILE"

                 !-write everything to file

                 yep = write_matrix ( xn, yn, REAL(psi,kind=4), TRIM(path) // 'psi.txt' )
                 yep = write_matrix ( xn, yn, REAL(h,kind=4), TRIM(path) // 'h.txt' )

                 yep = write_matrix ( xn, yn, REAL(u,kind=4), TRIM(path) // 'u.txt' )
                 yep = write_matrix ( xn, yn, REAL(v,kind=4), TRIM(path) // 'v.txt' )

                 yep = write_matrix ( (xn-1)/cellx, yn/(2*celly), REAL(h_coarse,kind=4), TRIM(path) // 'h_coarse.txt' )
                 yep = write_matrix ( (xn-1)/cellx, yn/(2*celly), REAL(u_coarse,kind=4), TRIM(path) // 'u_coarse.txt' )
                 yep = write_matrix ( (xn-1)/cellx, yn/(2*celly), REAL(v_coarse,kind=4), TRIM(path) // 'v_coarse.txt' )
                 yep = write_matrix ( (xn-1)/cellx, yn/(2*celly), REAL(psi_coarse,kind=4), TRIM(path) // 'psi_coarse.txt' )

                 yep = write_matrix ( xn, yn/2,REAL(permeability(:,(yn/2)+1:),kind=4), TRIM(path) // 'permeability.txt' )
                 yep = write_vec ( xn, REAL(x,kind=4), TRIM(path) // 'x.txt' )
                 yep = write_vec ( yn/2, REAL(y(yn/2:),kind=4), TRIM(path) // 'y.txt' )
                 yep = write_matrix ( xn, yn/2, REAL(mask(:,(yn/2)+1:), kind = 4), TRIM(path) // 'mask.txt' )
                 yep = write_matrix ( xn, yn/2, REAL(maskP(:,(yn/2)+1:), kind = 4), TRIM(path) // 'maskP.txt' )
                 yep = write_matrix ( (xn-1)/(cellx), yn/(2*celly), REAL(coarse_mask,kind=4), TRIM(path) // 'mask_coarse.txt' )


                 IF (restart .NE. 1) THEN

                    yep = write_matrix ( yn, 2, REAL(frac6, kind = 4), TRIM(path) // 'frac6.txt' )
                    yep = write_matrix ( yn, 2, REAL(temp6, kind = 4), TRIM(path) // 'temp6.txt' )

                    yep = write_matrix (xn*tn/(mstep*ar), yn/2, REAL(hmat(:,(yn/2)+1:), kind = 4), TRIM(path) // 'hMat.txt' )
                    yep = write_matrix (xn*tn/(mstep*ar), yn/2, REAL(psimat(:,(yn/2)+1:),kind=4), TRIM(path) // 'psiMat.txt' )
                    yep = write_matrix (xn*tn/(mstep*ar), yn/2, REAL(rhomat(:,(yn/2)+1:),kind=4), TRIM(path) // 'rhoMat.txt' )
                    yep = write_matrix (xn*tn/(mstep*ar), yn/2, REAL(umat(:,(yn/2)+1:), kind = 4), TRIM(path) // 'uMat.txt' )
                    yep = write_matrix (xn*tn/(mstep*ar), yn/2, REAL(vmat(:,(yn/2)+1:),kind=4), TRIM(path) // 'vMat.txt' )
                    yep = write_matrix ( xn, yn/2,REAL(lambdaMat(:,(yn/2)+1:),kind=4), TRIM(path) // 'lambdaMat.txt' )
                    yep = write_matrix ( xn, yn/2,REAL(phi(:,(yn/2)+1:),kind=4), TRIM(path) // 'phi.txt' )

                 END IF ! end write only if restart ne 1

                 write(*,*) "done writing basics"


                 ! solute concentrations
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,1),kind=4), TRIM(path) // 'ch_s/z_sol_ph.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,3),kind=4), TRIM(path) // 'ch_s/z_sol_w.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,2),kind=4), TRIM(path) // 'ch_s/z_sol_alk.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,4),kind=4), TRIM(path) // 'ch_s/z_sol_c.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,5),kind=4), TRIM(path) // 'ch_s/z_sol_ca.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,6),kind=4), TRIM(path) // 'ch_s/z_sol_mg.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,7),kind=4), TRIM(path) // 'ch_s/z_sol_na.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,8),kind=4), TRIM(path) // 'ch_s/z_sol_k.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,9),kind=4), TRIM(path) // 'ch_s/z_sol_fe.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,10),kind=4), TRIM(path) // 'ch_s/z_sol_s.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,11),kind=4), TRIM(path) // 'ch_s/z_sol_si.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,12),kind=4), TRIM(path) // 'ch_s/z_sol_cl.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,13),kind=4), TRIM(path) // 'ch_s/z_sol_al.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat(:,:,14),kind=4), TRIM(path) // 'ch_s/z_sol_inert.txt' )

                 ! primary minerals
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat(:,:,2),kind=4), TRIM(path) // 'ch_s/z_pri_plag.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat(:,:,3),kind=4), TRIM(path) // 'ch_s/z_pri_pyr.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat(:,:,4),kind=4), TRIM(path) // 'ch_s/z_pri_ol.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat(:,:,5),kind=4), TRIM(path) // 'ch_s/z_pri_glass.txt' )

                 ! medium properties
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat(:,:,1),kind=4), TRIM(path) // 'ch_s/z_med_phi.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat(:,:,2),kind=4), TRIM(path) // 'ch_s/z_med_precip.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat(:,:,3),kind=4), TRIM(path) // 'ch_s/z_med_v_water.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat(:,:,4),kind=4), TRIM(path) // 'ch_s/z_med_reactive.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat(:,:,5),kind=4), TRIM(path) // 'ch_s/z_med_cell_toggle.txt' )

                 ! phi
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(phiCalcMat(:,:),kind=4), TRIM(path) // 'ch_s/z_phiCalc.txt' )

                 ! ph_fix
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(ph_fixMat(:,:),kind=4), TRIM(path) // 'ch_s/z_ph_fix.txt' )

                 write(*,*) "done writing ch_s sol, med"



                 ! solute concentrations
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,1),kind=4), TRIM(path) // 'ch_a/z_sol_ph.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,3),kind=4), TRIM(path) // 'ch_a/z_sol_w.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,2),kind=4), TRIM(path) // 'ch_a/z_sol_alk.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,4),kind=4), TRIM(path) // 'ch_a/z_sol_c.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,5),kind=4), TRIM(path) // 'ch_a/z_sol_ca.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,6),kind=4), TRIM(path) // 'ch_a/z_sol_mg.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,7),kind=4), TRIM(path) // 'ch_a/z_sol_na.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,8),kind=4), TRIM(path) // 'ch_a/z_sol_k.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,9),kind=4), TRIM(path) // 'ch_a/z_sol_fe.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,10),kind=4), TRIM(path) // 'ch_a/z_sol_s.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,11),kind=4), TRIM(path) // 'ch_a/z_sol_si.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,12),kind=4), TRIM(path) // 'ch_a/z_sol_cl.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,13),kind=4), TRIM(path) // 'ch_a/z_sol_al.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_a(:,:,14),kind=4), TRIM(path) // 'ch_a/z_sol_inert.txt' )

                 ! primary minerals
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_a(:,:,2),kind=4), TRIM(path) // 'ch_a/z_pri_plag.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_a(:,:,3),kind=4), TRIM(path) // 'ch_a/z_pri_pyr.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_a(:,:,4),kind=4), TRIM(path) // 'ch_a/z_pri_ol.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_a(:,:,5),kind=4), TRIM(path) // 'ch_a/z_pri_glass.txt' )

                 ! medium properties
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_a(:,:,1),kind=4), TRIM(path) // 'ch_a/z_med_phi.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_a(:,:,2),kind=4), TRIM(path) // 'ch_a/z_med_precip.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_a(:,:,3),kind=4), TRIM(path) // 'ch_a/z_med_v_water.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_a(:,:,4),kind=4), TRIM(path) // 'ch_a/z_med_reactive.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_a(:,:,5),kind=4), TRIM(path) // 'ch_a/z_med_cell_toggle.txt' )

                 ! phi
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(phiCalcMat_a(:,:),kind=4), TRIM(path) // 'ch_a/z_phiCalc.txt' )

                 ! ph_fix
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(ph_fixMat_a(:,:),kind=4), TRIM(path) // 'ch_a/z_ph_fix.txt' )


                 write(*,*) "done writing ch_a sol, med"





                 ! solute concentrations
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,1),kind=4), TRIM(path) // 'ch_b/z_sol_ph.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,3),kind=4), TRIM(path) // 'ch_b/z_sol_w.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,2),kind=4), TRIM(path) // 'ch_b/z_sol_alk.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,4),kind=4), TRIM(path) // 'ch_b/z_sol_c.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,5),kind=4), TRIM(path) // 'ch_b/z_sol_ca.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,6),kind=4), TRIM(path) // 'ch_b/z_sol_mg.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,7),kind=4), TRIM(path) // 'ch_b/z_sol_na.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,8),kind=4), TRIM(path) // 'ch_b/z_sol_k.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,9),kind=4), TRIM(path) // 'ch_b/z_sol_fe.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,10),kind=4), TRIM(path) // 'ch_b/z_sol_s.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,11),kind=4), TRIM(path) // 'ch_b/z_sol_si.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,12),kind=4), TRIM(path) // 'ch_b/z_sol_cl.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,13),kind=4), TRIM(path) // 'ch_b/z_sol_al.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_b(:,:,14),kind=4), TRIM(path) // 'ch_b/z_sol_inert.txt' )

                 ! primary minerals
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_b(:,:,2),kind=4), TRIM(path) // 'ch_b/z_pri_plag.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_b(:,:,3),kind=4), TRIM(path) // 'ch_b/z_pri_pyr.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_b(:,:,4),kind=4), TRIM(path) // 'ch_b/z_pri_ol.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_b(:,:,5),kind=4), TRIM(path) // 'ch_b/z_pri_glass.txt' )

                 ! medium properties
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_b(:,:,1),kind=4), TRIM(path) // 'ch_b/z_med_phi.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_b(:,:,2),kind=4), TRIM(path) // 'ch_b/z_med_precip.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_b(:,:,3),kind=4), TRIM(path) // 'ch_b/z_med_v_water.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_b(:,:,4),kind=4), TRIM(path) // 'ch_b/z_med_reactive.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(mediumMat_b(:,:,5),kind=4), TRIM(path) // 'ch_b/z_med_cell_toggle.txt' )

                 ! ph_fix
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(ph_fixMat_b(:,:),kind=4), TRIM(path) // 'ch_b/z_ph_fix.txt' )


                 write(*,*) "done writing ch_b sol, med"

                 ! solute concentrations

                !  DO ii = (yn/(2*celly))+1,yn/(2*celly)
                DO ii = 1,yn/(2*celly)
                    DO i = 1,(xn-1)*tn/(cellx*mstep*ar)
                       soluteMat_d(i,ii,1) = -1.0*LOG10((volume_ratio/(1.0+volume_ratio))*10.0**(-1.0*soluteMat_a(i,ii,1)) + (1.0/(1.0+volume_ratio))*10.0**(-1.0*soluteMat_b(i,ii,1)))
                    END DO
		         END DO


                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,1),kind=4), TRIM(path) // 'ch_d/z_sol_ph.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,3),kind=4), TRIM(path) // 'ch_d/z_sol_w.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,2),kind=4), TRIM(path) // 'ch_d/z_sol_alk.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,4),kind=4), TRIM(path) // 'ch_d/z_sol_c.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,5),kind=4), TRIM(path) // 'ch_d/z_sol_ca.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,6),kind=4), TRIM(path) // 'ch_d/z_sol_mg.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,7),kind=4), TRIM(path) // 'ch_d/z_sol_na.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,8),kind=4), TRIM(path) // 'ch_d/z_sol_k.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,9),kind=4), TRIM(path) // 'ch_d/z_sol_fe.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,10),kind=4), TRIM(path) // 'ch_d/z_sol_s.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,11),kind=4), TRIM(path) // 'ch_d/z_sol_si.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,12),kind=4), TRIM(path) // 'ch_d/z_sol_cl.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,13),kind=4), TRIM(path) // 'ch_d/z_sol_al.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(soluteMat_d(:,:,14),kind=4), TRIM(path) // 'ch_d/z_sol_inert.txt' )

                 ! primary minerals
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_d(:,:,2),kind=4), TRIM(path) // 'ch_d/z_pri_plag.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_d(:,:,3),kind=4), TRIM(path) // 'ch_d/z_pri_pyr.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_d(:,:,4),kind=4), TRIM(path) // 'ch_d/z_pri_ol.txt' )
                 yep = write_matrix ( (xn-1)*tn/(cellx*mstep*ar), yn/(2*celly), REAL(primaryMat_d(:,:,5),kind=4), TRIM(path) // 'ch_d/z_pri_glass.txt' )

                 write(*,*) "done writing ch_d sol, med"







                 WRITE(*,*) "precipitated"
                 DO i = 1,g_sec/2
                    IF (MAXVAL(secondaryMat(:,:,i)) .GT. 0.0) THEN
                       WRITE(*,*) i

                       IF (i < 10) THEN
                          WRITE(s_i,'(i1)') i
                       ELSE
                          WRITE(s_i,'(i2)') i
                       END IF
                       yep = write_matrix((xn-1)*tn/(cellx*mstep*ar),yn/(2*celly),REAL(secondaryMat(:,:,i),kind=4),TRIM(path)//'ch_s/z_sec'//TRIM(s_i)//'.txt')
                    END IF

                    !write(*,*) "done writing ch_s sec"

                    IF (MAXVAL(secondaryMat_a(:,:,i)) .GT. 0.0) THEN
                       WRITE(*,*) i

                       IF (i < 10) THEN
                          WRITE(s_i,'(i1)') i
                       ELSE
                          WRITE(s_i,'(i2)') i
                       END IF
                       yep = write_matrix((xn-1)*tn/(cellx*mstep*ar),yn/(2*celly),REAL(secondaryMat_a(:,:,i),kind=4),TRIM(path)//'ch_a/z_sec'//TRIM(s_i)//'.txt')
                    END IF

                    !write(*,*) "done writing ch_a sec"

                    IF (MAXVAL(secondaryMat_b(:,:,i)) .GT. 0.0) THEN
                       WRITE(*,*) i

                       IF (i < 10) THEN
                          WRITE(s_i,'(i1)') i
                       ELSE
                          WRITE(s_i,'(i2)') i
                       END IF
                       yep = write_matrix((xn-1)*tn/(cellx*mstep*ar),yn/(2*celly),REAL(secondaryMat_b(:,:,i),kind=4),TRIM(path)//'ch_b/z_sec'//TRIM(s_i)//'.txt')
                    END IF

                    !write(*,*) "done writing ch_b sec"

                    IF (MAXVAL(secondaryMat_d(:,:,i)) .GT. 0.0) THEN
                       WRITE(*,*) i

                       IF (i < 10) THEN
                          WRITE(s_i,'(i1)') i
                       ELSE
                          WRITE(s_i,'(i2)') i
                       END IF
                       yep = write_matrix((xn-1)*tn/(cellx*mstep*ar),yn/(2*celly),REAL(secondaryMat_d(:,:,i),kind=4),TRIM(path)//'ch_d/z_sec'//TRIM(s_i)//'.txt')
                    END IF

                    !write(*,*) "done writing ch_d sec"

                 END DO


                 WRITE(*,*) "...DONE WRITING TO FILE"
             END IF ! end only write to file 50 times total


           END IF ! end if (mod(j,mstep*ar) .eq. 0)


        END IF
        ! end mstep timestep loop, finally


     END DO
     ! end all timestep loop




     WRITE(*,*) " "
     WRITE(*,*) "ALL DONE!"
     WRITE(*,*) tn
     WRITE(*,*) "steps"
     WRITE(*,*) tn/mstep
     WRITE(*,*) "msteps"




     ! what to do if you are a slave processor
  ELSE






     CALL init_mini()
    !  CALL MPI_Bcast( cml3, (3*(xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, 0, MPI_COMM_WORLD, status, ierr)
    !  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

    CALL MPI_RECV ( cml3, (3*(xn-1)/cellx)*(yn/(2*celly)) , MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

     write(*,*) my_id, sum(cml3)

     ! receive timestep size
     CALL MPI_RECV ( end_loop, 1 , MPI_INTEGER, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
     write(*,*) "my_id: " , my_id , "end loop is: " , end_loop




     leng = (yn/(2*celly))*((xn-1)/cellx)
     ! message receiving has to happen every mth step

     !-slave_vector!!!!!!
    !  write(*,*) "my_id: " , my_id , "sum: " , SUM(coarse_mask_long3)
    slave_count = 1
    index_count = 1
    do jj = 1 , 3*leng
        if (cml3(jj) .EQ. 1.0) then
            ! if (slave_count .EQ. (slave_count/(num_procs-1)) + my_id) then
            if (my_id .EQ. mod(slave_count-((num_procs-1)*(index_count-1)),num_procs)) then
                slave_vector(index_count) = jj
                index_count = index_count + 1
            end if
            slave_count = slave_count + 1
        end if
    end do
    write(*,*) "my_id: " , my_id , "slave_vector(end_loop)" , slave_vector(1:end_loop)

     !-primary compositions + amounts

     param_ol_string ='-f MgO 1.0 FeO 1.0 SiO2 1.0'
     !param_ol_string ='-f MgO 2.0 SiO2 1.0'
     !param_ol_string ='-f FeO 2.0 SiO2 1.0'


    !param_pyr_string='-f CaO 1.0 MgO 1.0 SiO2 2.0'
    param_pyr_string='-f MgO 2.0 SiO2 2.0'

    ! param_plag_string='-f NaAlSi3O8 0.5 CaAl2Si2O8 0.5'
    param_plag_string='-f NaAlSi3O8 1.0'

     ! 		&"-f CaO 1.0 FeO 1.0 SiO2 2.0 " //NEW_LINE('')// & ! hedenbergite
     ! 		&"-f CaO 1.0 MgO 1.0 SiO2 2.0 " //NEW_LINE('')// & ! diopside
     ! 		&"-f FeO 1.0 MgO 1.0 SiO2 2.0 " //NEW_LINE('')// & ! fer mag
     ! 		&"-f MgO 2.0 SiO2 2.0 " //NEW_LINE('')// & ! enstatite
     ! 		&"-f FeO 2.0 SiO2 2.0 " //NEW_LINE('')// & ! ferrosilite
     ! 		&"-f CaO 2.0 SiO2 2.0 " //NEW_LINE('')// & ! wollastonite


     ! 			param_exp_string = '0.0025'
     ! 			param_exp1_string = '0.0025'
     !
     ! !			exp_ol1 = "2.0"
     ! 			exp_ol1 = "0.0025"
     ! ! 			exp_ol2 = "0.1"
     ! ! 			exp_ol3 = "0.1"
     !
     ! 			exp_ol = exp_ol1
     !
     ! 			! exp_pyr1 = "2.0"
     ! 			exp_pyr1 = "0.0025"
     ! ! 			exp_pyr2 = "0.1"
     ! ! 			exp_pyr3 = "0.1"
     !
     ! 			exp_pyr = exp_pyr1
     !
     ! 			exp_plag1 = "0.25"
     ! ! 			exp_plag2 = "0.1"
     ! ! 			exp_plag3 = "0.1"
     !
     ! 			exp_plag = exp_plag1

     ol_k1 = "10.0^(-4.8)"
     ol_e1 = "94.4"
     ol_n1 = "1.0"
     ol_k2 = "10.0^(-12.8)"
     ol_e2 = "94.4"
     ol_k3 = ""
     ol_e3 = ""
     ol_n3 = ""

     ! 			ol_k1 = "10.0^(-6.85)"
     ! 			ol_e1 = "67.2"
     ! 			ol_n1 = "0.470"
     ! 			ol_k2 = "10.0^(-10.64)"
     ! 			ol_e2 = "79.0"
     ! 			ol_k3 = ""
     ! 			ol_e3 = ""
     ! 			ol_n3 = ""

     pyr_k1 = "10.0^(-6.82)"
     pyr_e1 = "78.0"
     pyr_n1 = "0.7"
     pyr_k2 = "10.0^(-11.97)"
     pyr_e2 = "78.0"
     pyr_k3 = ""
     pyr_e3 = ""
     pyr_n3 = ""

     plag_k1 = "10.0^(-7.87)"
     plag_e1 = "42.1"
     plag_n1 = "0.626"
     plag_k2 = "10.0^(-10.91)"
     plag_e2 = "45.2"
     plag_k3 = ""
     plag_e3 = ""
     plag_n3 = ""


     kinetics = " precipitate_only"
     !kinetics = " "


     DO jj = 1, tn/mstep

        IF (my_id .LE. 11) THEN

           !#ADVECTION: slave receives message

           ! receive an_id
           CALL MPI_RECV ( an_id_local, 1 , MPI_INTEGER, &
                root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

           ! receive solute long for advection
           CALL MPI_RECV ( sol_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

           ! receive u long for advection
           CALL MPI_RECV ( u_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

           ! receive v long for advection
           CALL MPI_RECV ( v_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

           ! receive phi long for advection
           CALL MPI_RECV ( phi_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

           !write(*,*) "advection proc my_id" , my_id , an_id_local

           ! reshape them all
           sol_coarse_local = TRANSPOSE(RESHAPE(sol_coarse_long_local,(/yn/(2*celly),(xn-1)/cellx/)))
           u_coarse_local = RESHAPE(u_coarse_long_local,(/(xn-1)/cellx,yn/(2*celly)/))
           v_coarse_local = RESHAPE(v_coarse_long_local,(/(xn-1)/cellx,yn/(2*celly)/))
           phi_coarse_local = TRANSPOSE(RESHAPE(phi_coarse_long_local,(/yn/(2*celly),(xn-1)/cellx/))) ! this is where the phi transpose bug was
           !phi_coarse_local = 0.5
           !write(*,*) maxval(phi_coarse_local)

           ! do advection
           IF (an_id_local .LE. 11) THEN
              DO ii = 1,cstep
                 sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id_local)))
              END DO
              !-proc 34 cell fix
            !   sol_coarse_local(2,:) = sol_coarse_local(3,:)
           END IF

           IF (an_id_local .GT. 11) THEN
              DO ii = 1,cstep
                 sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id_local-11)))
              END DO
              !-proc 34 cell fix
            !   sol_coarse_local(2,:) = sol_coarse_local(3,:)
           END IF

           sol_coarse_long_local = RESHAPE(TRANSPOSE(sol_coarse_local),(/((xn-1)/cellx)*(yn/(2*celly))/))

           ! send advected solutes back :)
           CALL MPI_SEND( sol_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, root_process, &
                return_data_tag, MPI_COMM_WORLD, ierr)




            ! receive an_id
            CALL MPI_RECV ( an_id_local, 1 , MPI_INTEGER, &
                 root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

            ! receive solute long for advection
            CALL MPI_RECV ( sol_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                 root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

            ! receive u long for advection
            CALL MPI_RECV ( u_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                 root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

            ! receive v long for advection
            CALL MPI_RECV ( v_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                 root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

            ! receive phi long for advection
            CALL MPI_RECV ( phi_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, &
                 root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

            !write(*,*) "advection proc my_id" , my_id , an_id_local

            ! reshape them all
            sol_coarse_local = TRANSPOSE(RESHAPE(sol_coarse_long_local,(/yn/(2*celly),(xn-1)/cellx/)))
            u_coarse_local = RESHAPE(u_coarse_long_local,(/(xn-1)/cellx,yn/(2*celly)/))
            v_coarse_local = RESHAPE(v_coarse_long_local,(/(xn-1)/cellx,yn/(2*celly)/))
            phi_coarse_local = TRANSPOSE(RESHAPE(phi_coarse_long_local,(/yn/(2*celly),(xn-1)/cellx/))) ! this is where the phi transpose bug was
            !phi_coarse_local = 0.5
            !write(*,*) maxval(phi_coarse_local)

            ! do advection
            IF (an_id_local .LE. 11) THEN
               DO ii = 1,cstep
                  sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id_local)))
               END DO
            END IF

            IF (an_id_local .GT. 11) THEN
               DO ii = 1,cstep
                  sol_coarse_local = solute_next_coarse(sol_coarse_local,u_coarse_local,v_coarse_local,phi_coarse_local,sea(sol_index(an_id_local-11)))
               END DO
            END IF

            sol_coarse_long_local = RESHAPE(TRANSPOSE(sol_coarse_local),(/((xn-1)/cellx)*(yn/(2*celly))/))

            ! send advected solutes back :)
            CALL MPI_SEND( sol_coarse_long_local, ((xn-1)/cellx)*(yn/(2*celly)), MPI_REAL4, root_process, &
                 return_data_tag, MPI_COMM_WORLD, ierr)

        END IF ! end if my_id .le. 11


        !#GEOCHEM: slave receives from root

        CALL MPI_RECV ( j_root, 1, MPI_INTEGER, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

        !CALL MPI_RECV ( age_root, 1, MPI_INTEGER, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

        se_toggle = 0
        !-set SE_TOGGLE
        IF ((MOD(j_root,mstep*se_factor) .EQ. 0) .OR. (j_root .LE. mstep*se_spinup) .OR. (MOD(j_root-mstep,mstep*se_factor) .EQ. 0) .OR. (MOD(j_root-(2*mstep),mstep*se_factor) .EQ. 0) .OR. (MOD(j_root-(3*mstep),mstep*se_factor) .EQ. 0)) THEN
            se_toggle = 1
        END IF

        if (my_id .eq. 10) then
            write(*,*) "j_root:" , j_root , "se_toggle:" , se_toggle
        end if

        CALL MPI_RECV ( hLong, 3*leng, MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

        DO ii = 1 , g_pri
            CALL MPI_RECV ( priLongBitFull(:,ii), 3*leng, MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        END DO

        DO ii = 1 , g_sec/2
            CALL MPI_RECV ( secLongBitFull(:,ii), 3*leng, MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        END DO

        DO ii = 1 , g_sol
            CALL MPI_RECV ( solLongBitFull(:,ii), 3*leng, MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        END DO

        DO ii = 1 , g_med
            CALL MPI_RECV ( medLongBitFull(:,ii), 3*leng, MPI_REAL4, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        END DO

        dt_local = dt



        call system_clock(counti, count_rate, count_max)

        DO jjj = 1,end_loop



        timestep3 = dt_local*mstep
        WRITE(s_timestep,'(F25.10)') timestep3

        ! slave processor loops through each coarse cell
        m=slave_vector(jjj)

           IF (medLongBitFull(m,5) .EQ. 1.0) THEN

               !-primary3, secondar 3, etc.
              primary3 = priLongBitFull(m,:)
              secondary3 = secLongBitFull(m,:)
              solute3 = solLongBitFull(m,:)
              medium3 = medLongBitFull(m,:)
             !  dprimary3 = dpriLongBitFull(m,:)
              dprimary3 = dpriLocal(m,:)
            !   dsecondary3 = dsecLongBitFull(m,:)
              dsecondary3 = dsecLocal(m,:)

              temp3 = hLong(m)-273.0
              IF (temp3 .GE. 300.0) THEN
                 temp3 = 299.0
              END IF


              ! SOLUTES TO STRINGS
              WRITE(s_ph,'(F25.10)') solute3(1)
              WRITE(s_alk,'(F25.10)') solute3(2)
              WRITE(s_water,'(F25.10)') solute3(3)
              WRITE(s_co2,'(F25.10)') solute3(4)
              WRITE(s_ca,'(F25.10)') solute3(5)
              WRITE(s_mg,'(F25.10)') solute3(6)
              WRITE(s_na,'(F25.10)') solute3(7)
              WRITE(s_k,'(F25.10)') solute3(8)
              WRITE(s_fe,'(F25.10)') solute3(9)
              WRITE(s_s,'(F25.10)') solute3(10)
              WRITE(s_si,'(F25.10)') solute3(11)
              WRITE(s_cl,'(F25.10)') solute3(12)
              WRITE(s_al,'(F25.10)') solute3(13)
              WRITE(s_hco3,'(F25.10)') solute3(14)
              WRITE(s_co3,'(F25.10)') solute3(15)

              ! MEDIUM TO STRINGS
              WRITE(s_w,'(F25.10)') medium3(3) !solute3(3)

              WRITE(s_basalt3,'(F25.10)') primary3(2)
              WRITE(s_basalt2,'(F25.10)') primary3(3)
              WRITE(s_basalt1,'(F25.10)') primary3(4)
              WRITE(s_glass,'(F25.10)') primary3(5)

              !-rate constants!
              exp_ol = "0.0375"
              exp_pyr = "0.005"
              exp_plag = "0.05"
              param_exp_string = "0.000125"

              WRITE(s_kaolinite,'(F25.10)') secondary3(1)
              WRITE(s_saponite,'(F25.10)') secondary3(2)
              WRITE(s_celadonite,'(F25.10)') secondary3(3)
              WRITE(s_clinoptilolite,'(F25.10)') secondary3(4)
              WRITE(s_pyrite,'(F25.10)') secondary3(5)
              WRITE(s_mont_na,'(F25.10)') secondary3(6)
              WRITE(s_goethite,'(F25.10)') secondary3(7)
              WRITE(s_smectite,'(F25.10)') secondary3(8)
              WRITE(s_calcite,'(F25.10)') secondary3(9)
              WRITE(s_kspar,'(F25.10)') secondary3(10)
              WRITE(s_saponite_na,'(F25.10)') secondary3(11) !!!!
              WRITE(s_nont_na,'(F25.10)') secondary3(12)
              WRITE(s_nont_mg,'(F25.10)') secondary3(13)
              WRITE(s_fe_celadonite,'(F25.10)') secondary3(14)
              WRITE(s_nont_ca,'(F25.10)') secondary3(15)
              WRITE(s_mesolite,'(F25.10)') secondary3(16)
              WRITE(s_hematite,'(F25.10)') secondary3(17)
              WRITE(s_mont_ca,'(F25.10)') secondary3(18)
              WRITE(s_verm_ca,'(F25.10)') secondary3(19)
              WRITE(s_analcime,'(F25.10)') secondary3(20)
              WRITE(s_phillipsite,'(F25.10)') secondary3(21)
              WRITE(s_mont_mg,'(F25.10)') secondary3(22)
              WRITE(s_gismondine,'(F25.10)') secondary3(23)
              WRITE(s_verm_mg,'(F25.10)') secondary3(24)
              WRITE(s_natrolite,'(F25.10)') secondary3(25)
              WRITE(s_talc,'(F25.10)') secondary3(26) !!!!!!!!!
              WRITE(s_smectite_low,'(F25.10)') secondary3(27)
              WRITE(s_prehnite,'(F25.10)') secondary3(28)
              WRITE(s_chlorite,'(F25.10)') secondary3(29) !!!!!!!!
              WRITE(s_scolecite,'(F25.10)') secondary3(30)
              WRITE(s_clinochlore14a,'(F25.10)') secondary3(31)
              WRITE(s_clinochlore7a,'(F25.10)') secondary3(32)
              WRITE(s_saponite_ca,'(F25.10)') secondary3(33)
              WRITE(s_verm_na,'(F25.10)') secondary3(34)
              WRITE(s_pyrrhotite,'(F25.10)') secondary3(35)
              WRITE(s_fe_saponite_ca,'(F25.10)') secondary3(36) !!!
              WRITE(s_fe_saponite_mg,'(F25.10)') secondary3(37) !!!
              WRITE(s_daphnite_7a,'(F25.10)') secondary3(38) !!!
              WRITE(s_daphnite_14a,'(F25.10)') secondary3(39) !!!
              WRITE(s_epidote,'(F25.10)') secondary3(40) !!!


              !-deltas to strings
              WRITE(sd_dbasalt3,'(F25.10)') -1.0*dprimary3(2)
              WRITE(sd_dbasalt2,'(F25.10)') -1.0*dprimary3(3)
              WRITE(sd_dbasalt1,'(F25.10)') -1.0*dprimary3(4)
              WRITE(sd_dglass,'(F25.10)') -1.0*dprimary3(5)

            !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
            !       write(*,*) "plag sd_dbasalt3" , sd_dbasalt3
            !       write(*,*) "pyr sd_dbasalt2" , sd_dbasalt2
            !       write(*,*) "ol sd_dbasalt1" , sd_dbasalt1
            !       write(*,*) "sd_dglass" , sd_dglass
            !   end if

            !   WRITE(sd_kaolinite,'(F25.10)') -1.0*dsecondary3(1)
            !   WRITE(sd_saponite,'(F25.10)') -1.0*dsecondary3(2)
            !   WRITE(sd_celadonite,'(F25.10)') -1.0*dsecondary3(3)
            !   WRITE(sd_clinoptilolite,'(F25.10)') -1.0*dsecondary3(4)
            !   WRITE(sd_pyrite,'(F25.10)') -1.0*dsecondary3(5)
            !   WRITE(sd_mont_na,'(F25.10)') -1.0*dsecondary3(6)
            !   WRITE(sd_goethite,'(F25.10)') -1.0*dsecondary3(7)
            !   WRITE(sd_smectite,'(F25.10)') -1.0*dsecondary3(8)
            !   WRITE(sd_calcite,'(F25.10)') -1.0*dsecondary3(9)
            !   WRITE(sd_kspar,'(F25.10)') -1.0*dsecondary3(10)
            !   WRITE(sd_saponite_na,'(F25.10)') -1.0*dsecondary3(11) !!!!
            !   WRITE(sd_nont_na,'(F25.10)') -1.0*dsecondary3(12)
            !   WRITE(sd_nont_mg,'(F25.10)') -1.0*dsecondary3(13)
            !   WRITE(sd_fe_celadonite,'(F25.10)') -1.0*dsecondary3(14)
            !   WRITE(sd_nont_ca,'(F25.10)') -1.0*dsecondary3(15)
            !   WRITE(sd_mesolite,'(F25.10)') -1.0*dsecondary3(16)
            !   WRITE(sd_hematite,'(F25.10)') -1.0*dsecondary3(17)
            !   WRITE(sd_mont_ca,'(F25.10)') -1.0*dsecondary3(18)
            !   WRITE(sd_verm_ca,'(F25.10)') -1.0*dsecondary3(19)
            !   WRITE(sd_analcime,'(F25.10)') -1.0*dsecondary3(20)
            !   WRITE(sd_phillipsite,'(F25.10)') -1.0*dsecondary3(21)
            !   WRITE(sd_mont_mg,'(F25.10)') -1.0*dsecondary3(22)
            !   WRITE(sd_gismondine,'(F25.10)') -1.0*dsecondary3(23)
            !   WRITE(sd_verm_mg,'(F25.10)') -1.0*dsecondary3(24)
            !   WRITE(sd_natrolite,'(F25.10)') -1.0*dsecondary3(25)
            !   WRITE(sd_talc,'(F25.10)') -1.0*dsecondary3(26) !!!!!!!!!
            !   WRITE(sd_smectite_low,'(F25.10)') -1.0*dsecondary3(27)
            !   WRITE(sd_prehnite,'(F25.10)') -1.0*dsecondary3(28)
            !   WRITE(sd_chlorite,'(F25.10)') -1.0*dsecondary3(29) !!!!!!!!
            !   WRITE(sd_scolecite,'(F25.10)') -1.0*dsecondary3(30)
            !   WRITE(sd_clinochlore14a,'(F25.10)') -1.0*dsecondary3(31)
            !   WRITE(sd_clinochlore7a,'(F25.10)') -1.0*dsecondary3(32)
            !   WRITE(sd_saponite_ca,'(F25.10)') -1.0*dsecondary3(33)
            !   WRITE(sd_verm_na,'(F25.10)') -1.0*dsecondary3(34)
            !   WRITE(sd_pyrrhotite,'(F25.10)') -1.0*dsecondary3(35)
            !   WRITE(sd_fe_saponite_ca,'(F25.10)') -1.0*dsecondary3(36) !!!
            !   WRITE(sd_fe_saponite_mg,'(F25.10)') -1.0*dsecondary3(37) !!!
            !   WRITE(sd_daphnite_7a,'(F25.10)') -1.0*dsecondary3(38) !!!
            !   WRITE(sd_daphnite_14a,'(F25.10)') -1.0*dsecondary3(39) !!!
            !   WRITE(sd_epidote,'(F25.10)') -1.0*dsecondary3(40) !!!


            ! WRITE(sd_kaolinite,'(F25.10)') -0.9*dsecondary3(1)
            ! WRITE(sd_saponite,'(F25.10)') -0.9*dsecondary3(2)
            ! WRITE(sd_celadonite,'(F25.10)') -0.9*dsecondary3(3)
            ! WRITE(sd_clinoptilolite,'(F25.10)') -0.9*dsecondary3(4)
            ! WRITE(sd_pyrite,'(F25.10)') -0.9*dsecondary3(5)
            ! WRITE(sd_mont_na,'(F25.10)') -0.9*dsecondary3(6)
            ! WRITE(sd_goethite,'(F25.10)') -0.9*dsecondary3(7)
            ! WRITE(sd_smectite,'(F25.10)') -0.9*dsecondary3(8)
            ! WRITE(sd_calcite,'(F25.10)') -0.9*dsecondary3(9)
            ! WRITE(sd_kspar,'(F25.10)') -0.9*dsecondary3(10)
            ! WRITE(sd_saponite_na,'(F25.10)') -0.9*dsecondary3(11) !!!!
            ! WRITE(sd_nont_na,'(F25.10)') -0.9*dsecondary3(12)
            ! WRITE(sd_nont_mg,'(F25.10)') -0.9*dsecondary3(13)
            ! WRITE(sd_fe_celadonite,'(F25.10)') -0.9*dsecondary3(14)
            ! WRITE(sd_nont_ca,'(F25.10)') -0.9*dsecondary3(15)
            ! WRITE(sd_mesolite,'(F25.10)') -0.9*dsecondary3(16)
            ! WRITE(sd_hematite,'(F25.10)') -0.9*dsecondary3(17)
            ! WRITE(sd_mont_ca,'(F25.10)') -0.9*dsecondary3(18)
            ! WRITE(sd_verm_ca,'(F25.10)') -0.9*dsecondary3(19)
            ! WRITE(sd_analcime,'(F25.10)') -0.9*dsecondary3(20)
            ! WRITE(sd_phillipsite,'(F25.10)') -0.9*dsecondary3(21)
            ! WRITE(sd_mont_mg,'(F25.10)') -0.9*dsecondary3(22)
            ! WRITE(sd_gismondine,'(F25.10)') -0.9*dsecondary3(23)
            ! WRITE(sd_verm_mg,'(F25.10)') -0.9*dsecondary3(24)
            ! WRITE(sd_natrolite,'(F25.10)') -0.9*dsecondary3(25)
            ! WRITE(sd_talc,'(F25.10)') -0.9*dsecondary3(26) !!!!!!!!!
            ! WRITE(sd_smectite_low,'(F25.10)') -0.9*dsecondary3(27)
            ! WRITE(sd_prehnite,'(F25.10)') -0.9*dsecondary3(28)
            ! WRITE(sd_chlorite,'(F25.10)') -0.9*dsecondary3(29) !!!!!!!!
            ! WRITE(sd_scolecite,'(F25.10)') -0.9*dsecondary3(30)
            ! WRITE(sd_clinochlore14a,'(F25.10)') -0.9*dsecondary3(31)
            ! WRITE(sd_clinochlore7a,'(F25.10)') -0.9*dsecondary3(32)
            ! WRITE(sd_saponite_ca,'(F25.10)') -0.9*dsecondary3(33)
            ! WRITE(sd_verm_na,'(F25.10)') -0.9*dsecondary3(34)
            ! WRITE(sd_pyrrhotite,'(F25.10)') -0.9*dsecondary3(35)
            ! WRITE(sd_fe_saponite_ca,'(F25.10)') -0.9*dsecondary3(36) !!!
            ! WRITE(sd_fe_saponite_mg,'(F25.10)') -0.9*dsecondary3(37) !!!
            ! WRITE(sd_daphnite_7a,'(F25.10)') -0.9*dsecondary3(38) !!!
            ! WRITE(sd_daphnite_14a,'(F25.10)') -0.9*dsecondary3(39) !!!
            ! WRITE(sd_epidote,'(F25.10)') -0.9*dsecondary3(40) !!!

            ! WRITE(sd_kaolinite,'(F25.10)') -0.1*dsecondary3(1)
            ! WRITE(sd_saponite,'(F25.10)') -0.1*dsecondary3(2)
            ! WRITE(sd_celadonite,'(F25.10)') -0.1*dsecondary3(3)
            ! WRITE(sd_clinoptilolite,'(F25.10)') -0.1*dsecondary3(4)
            ! WRITE(sd_pyrite,'(F25.10)') -0.1*dsecondary3(5)
            ! WRITE(sd_mont_na,'(F25.10)') -0.1*dsecondary3(6)
            ! WRITE(sd_goethite,'(F25.10)') -0.1*dsecondary3(7)
            ! WRITE(sd_smectite,'(F25.10)') -0.1*dsecondary3(8)
            ! WRITE(sd_calcite,'(F25.10)') -0.1*dsecondary3(9)
            ! WRITE(sd_kspar,'(F25.10)') -0.1*dsecondary3(10)
            ! WRITE(sd_saponite_na,'(F25.10)') -0.1*dsecondary3(11) !!!!
            ! WRITE(sd_nont_na,'(F25.10)') -0.1*dsecondary3(12)
            ! WRITE(sd_nont_mg,'(F25.10)') -0.1*dsecondary3(13)
            ! WRITE(sd_fe_celadonite,'(F25.10)') -0.1*dsecondary3(14)
            ! WRITE(sd_nont_ca,'(F25.10)') -0.1*dsecondary3(15)
            ! WRITE(sd_mesolite,'(F25.10)') -0.1*dsecondary3(16)
            ! WRITE(sd_hematite,'(F25.10)') -0.1*dsecondary3(17)
            ! WRITE(sd_mont_ca,'(F25.10)') -0.1*dsecondary3(18)
            ! WRITE(sd_verm_ca,'(F25.10)') -0.1*dsecondary3(19)
            ! WRITE(sd_analcime,'(F25.10)') -0.1*dsecondary3(20)
            ! WRITE(sd_phillipsite,'(F25.10)') -0.1*dsecondary3(21)
            ! WRITE(sd_mont_mg,'(F25.10)') -0.1*dsecondary3(22)
            ! WRITE(sd_gismondine,'(F25.10)') -0.1*dsecondary3(23)
            ! WRITE(sd_verm_mg,'(F25.10)') -0.1*dsecondary3(24)
            ! WRITE(sd_natrolite,'(F25.10)') -0.1*dsecondary3(25)
            ! WRITE(sd_talc,'(F25.10)') -0.1*dsecondary3(26) !!!!!!!!!
            ! WRITE(sd_smectite_low,'(F25.10)') -0.1*dsecondary3(27)
            ! WRITE(sd_prehnite,'(F25.10)') -0.1*dsecondary3(28)
            ! WRITE(sd_chlorite,'(F25.10)') -0.1*dsecondary3(29) !!!!!!!!
            ! WRITE(sd_scolecite,'(F25.10)') -0.1*dsecondary3(30)
            ! WRITE(sd_clinochlore14a,'(F25.10)') -0.1*dsecondary3(31)
            ! WRITE(sd_clinochlore7a,'(F25.10)') -0.1*dsecondary3(32)
            ! WRITE(sd_saponite_ca,'(F25.10)') -0.1*dsecondary3(33)
            ! WRITE(sd_verm_na,'(F25.10)') -0.1*dsecondary3(34)
            ! WRITE(sd_pyrrhotite,'(F25.10)') -0.1*dsecondary3(35)
            ! WRITE(sd_fe_saponite_ca,'(F25.10)') -0.1*dsecondary3(36) !!!
            ! WRITE(sd_fe_saponite_mg,'(F25.10)') -0.1*dsecondary3(37) !!!
            ! WRITE(sd_daphnite_7a,'(F25.10)') -0.1*dsecondary3(38) !!!
            ! WRITE(sd_daphnite_14a,'(F25.10)') -0.1*dsecondary3(39) !!!
            ! WRITE(sd_epidote,'(F25.10)') -0.1*dsecondary3(40) !!!

            ! WRITE(sd_kaolinite,'(F25.10)') -0.999*dsecondary3(1)
            ! WRITE(sd_saponite,'(F25.10)') -0.999*dsecondary3(2)
            ! WRITE(sd_celadonite,'(F25.10)') -0.999*dsecondary3(3)
            ! WRITE(sd_clinoptilolite,'(F25.10)') -0.999*dsecondary3(4)
            ! WRITE(sd_pyrite,'(F25.10)') -0.999*dsecondary3(5)
            ! WRITE(sd_mont_na,'(F25.10)') -0.999*dsecondary3(6)
            ! WRITE(sd_goethite,'(F25.10)') -0.999*dsecondary3(7)
            ! WRITE(sd_smectite,'(F25.10)') -0.999*dsecondary3(8)
            ! WRITE(sd_calcite,'(F25.10)') -0.999*dsecondary3(9)
            ! WRITE(sd_kspar,'(F25.10)') -0.999*dsecondary3(10)
            ! WRITE(sd_saponite_na,'(F25.10)') -0.999*dsecondary3(11) !!!!
            ! WRITE(sd_nont_na,'(F25.10)') -0.999*dsecondary3(12)
            ! WRITE(sd_nont_mg,'(F25.10)') -0.999*dsecondary3(13)
            ! WRITE(sd_fe_celadonite,'(F25.10)') -0.999*dsecondary3(14)
            ! WRITE(sd_nont_ca,'(F25.10)') -0.999*dsecondary3(15)
            ! WRITE(sd_mesolite,'(F25.10)') -0.999*dsecondary3(16)
            ! WRITE(sd_hematite,'(F25.10)') -0.999*dsecondary3(17)
            ! WRITE(sd_mont_ca,'(F25.10)') -0.999*dsecondary3(18)
            ! WRITE(sd_verm_ca,'(F25.10)') -0.999*dsecondary3(19)
            ! WRITE(sd_analcime,'(F25.10)') -0.999*dsecondary3(20)
            ! WRITE(sd_phillipsite,'(F25.10)') -0.999*dsecondary3(21)
            ! WRITE(sd_mont_mg,'(F25.10)') -0.999*dsecondary3(22)
            ! WRITE(sd_gismondine,'(F25.10)') -0.999*dsecondary3(23)
            ! WRITE(sd_verm_mg,'(F25.10)') -0.999*dsecondary3(24)
            ! WRITE(sd_natrolite,'(F25.10)') -0.999*dsecondary3(25)
            ! WRITE(sd_talc,'(F25.10)') -0.999*dsecondary3(26) !!!!!!!!!
            ! WRITE(sd_smectite_low,'(F25.10)') -0.999*dsecondary3(27)
            ! WRITE(sd_prehnite,'(F25.10)') -0.999*dsecondary3(28)
            ! WRITE(sd_chlorite,'(F25.10)') -0.999*dsecondary3(29) !!!!!!!!
            ! WRITE(sd_scolecite,'(F25.10)') -0.999*dsecondary3(30)
            ! WRITE(sd_clinochlore14a,'(F25.10)') -0.999*dsecondary3(31)
            ! WRITE(sd_clinochlore7a,'(F25.10)') -0.999*dsecondary3(32)
            ! WRITE(sd_saponite_ca,'(F25.10)') -0.999*dsecondary3(33)
            ! WRITE(sd_verm_na,'(F25.10)') -0.999*dsecondary3(34)
            ! WRITE(sd_pyrrhotite,'(F25.10)') -0.999*dsecondary3(35)
            ! WRITE(sd_fe_saponite_ca,'(F25.10)') -0.999*dsecondary3(36) !!!
            ! WRITE(sd_fe_saponite_mg,'(F25.10)') -0.999*dsecondary3(37) !!!
            ! WRITE(sd_daphnite_7a,'(F25.10)') -0.999*dsecondary3(38) !!!
            ! WRITE(sd_daphnite_14a,'(F25.10)') -0.999*dsecondary3(39) !!!
            ! WRITE(sd_epidote,'(F25.10)') -0.999*dsecondary3(40) !!!


            !IF ((se_toggle .EQ. 1) .OR. (MOD(j_root-(4*mstep),mstep*se_factor) .EQ. 0) .OR. (MOD(j_root-(5*mstep),mstep*se_factor) .EQ. 0)) THEN

                WRITE(sd_kaolinite,'(F25.10)') -0.99*dsecondary3(1)
                WRITE(sd_saponite,'(F25.10)') -0.99*dsecondary3(2)
                WRITE(sd_celadonite,'(F25.10)') -0.99*dsecondary3(3)
                WRITE(sd_clinoptilolite,'(F25.10)') -0.99*dsecondary3(4)
                WRITE(sd_pyrite,'(F25.10)') -0.99*dsecondary3(5)
                WRITE(sd_mont_na,'(F25.10)') -0.99*dsecondary3(6)
                WRITE(sd_goethite,'(F25.10)') -0.99*dsecondary3(7)
                WRITE(sd_smectite,'(F25.10)') -0.99*dsecondary3(8)
                WRITE(sd_calcite,'(F25.10)') -0.99*dsecondary3(9)
                WRITE(sd_kspar,'(F25.10)') -0.99*dsecondary3(10)
                WRITE(sd_saponite_na,'(F25.10)') -0.99*dsecondary3(11) !!!!
                WRITE(sd_nont_na,'(F25.10)') -0.99*dsecondary3(12)
                WRITE(sd_nont_mg,'(F25.10)') -0.99*dsecondary3(13)
                WRITE(sd_fe_celadonite,'(F25.10)') -0.99*dsecondary3(14)
                WRITE(sd_nont_ca,'(F25.10)') -0.99*dsecondary3(15)
                WRITE(sd_mesolite,'(F25.10)') -0.99*dsecondary3(16)
                WRITE(sd_hematite,'(F25.10)') -0.99*dsecondary3(17)
                WRITE(sd_mont_ca,'(F25.10)') -0.99*dsecondary3(18)
                WRITE(sd_verm_ca,'(F25.10)') -0.99*dsecondary3(19)
                WRITE(sd_analcime,'(F25.10)') -0.99*dsecondary3(20)
                WRITE(sd_phillipsite,'(F25.10)') -0.99*dsecondary3(21)
                WRITE(sd_mont_mg,'(F25.10)') -0.99*dsecondary3(22)
                WRITE(sd_gismondine,'(F25.10)') -0.99*dsecondary3(23)
                WRITE(sd_verm_mg,'(F25.10)') -0.99*dsecondary3(24)
                WRITE(sd_natrolite,'(F25.10)') -0.99*dsecondary3(25)
                WRITE(sd_talc,'(F25.10)') -0.99*dsecondary3(26) !!!!!!!!!
                WRITE(sd_smectite_low,'(F25.10)') -0.99*dsecondary3(27)
                WRITE(sd_prehnite,'(F25.10)') -0.99*dsecondary3(28)
                WRITE(sd_chlorite,'(F25.10)') -0.99*dsecondary3(29) !!!!!!!!
                WRITE(sd_scolecite,'(F25.10)') -0.99*dsecondary3(30)
                WRITE(sd_clinochlore14a,'(F25.10)') -0.99*dsecondary3(31)
                WRITE(sd_clinochlore7a,'(F25.10)') -0.99*dsecondary3(32)
                WRITE(sd_saponite_ca,'(F25.10)') -0.99*dsecondary3(33)
                WRITE(sd_verm_na,'(F25.10)') -0.99*dsecondary3(34)
                WRITE(sd_pyrrhotite,'(F25.10)') -0.99*dsecondary3(35)
                WRITE(sd_fe_saponite_ca,'(F25.10)') -0.99*dsecondary3(36) !!!
                WRITE(sd_fe_saponite_mg,'(F25.10)') -0.99*dsecondary3(37) !!!
                WRITE(sd_daphnite_7a,'(F25.10)') -0.99*dsecondary3(38) !!!
                WRITE(sd_daphnite_14a,'(F25.10)') -0.99*dsecondary3(39) !!!
                WRITE(sd_epidote,'(F25.10)') -0.99*dsecondary3(40) !!!

            !END IF

            ! WRITE(sd_kaolinite,'(F25.10)') -0.999*dsecondary3(1)
            ! WRITE(sd_saponite,'(F25.10)') -0.999*dsecondary3(2)
            ! WRITE(sd_celadonite,'(F25.10)') -0.999*dsecondary3(3)
            ! WRITE(sd_clinoptilolite,'(F25.10)') -0.999*dsecondary3(4)
            ! WRITE(sd_pyrite,'(F25.10)') -0.999*dsecondary3(5)
            ! WRITE(sd_mont_na,'(F25.10)') -0.999*dsecondary3(6)
            ! WRITE(sd_goethite,'(F25.10)') -0.999*dsecondary3(7)
            ! WRITE(sd_smectite,'(F25.10)') -0.999*dsecondary3(8)
            ! WRITE(sd_calcite,'(F25.10)') -0.999*dsecondary3(9)
            ! WRITE(sd_kspar,'(F25.10)') -0.999*dsecondary3(10)
            ! WRITE(sd_saponite_na,'(F25.10)') -0.999*dsecondary3(11) !!!!
            ! WRITE(sd_nont_na,'(F25.10)') -0.999*dsecondary3(12)
            ! WRITE(sd_nont_mg,'(F25.10)') -0.999*dsecondary3(13)
            ! WRITE(sd_fe_celadonite,'(F25.10)') -0.999*dsecondary3(14)
            ! WRITE(sd_nont_ca,'(F25.10)') -0.999*dsecondary3(15)
            ! WRITE(sd_mesolite,'(F25.10)') -0.999*dsecondary3(16)
            ! WRITE(sd_hematite,'(F25.10)') -0.999*dsecondary3(17)
            ! WRITE(sd_mont_ca,'(F25.10)') -0.999*dsecondary3(18)
            ! WRITE(sd_verm_ca,'(F25.10)') -0.999*dsecondary3(19)
            ! WRITE(sd_analcime,'(F25.10)') -0.999*dsecondary3(20)
            ! WRITE(sd_phillipsite,'(F25.10)') -0.999*dsecondary3(21)
            ! WRITE(sd_mont_mg,'(F25.10)') -0.999*dsecondary3(22)
            ! WRITE(sd_gismondine,'(F25.10)') -0.999*dsecondary3(23)
            ! WRITE(sd_verm_mg,'(F25.10)') -0.999*dsecondary3(24)
            ! WRITE(sd_natrolite,'(F25.10)') -0.999*dsecondary3(25)
            ! WRITE(sd_talc,'(F25.10)') -0.999*dsecondary3(26) !!!!!!!!!
            ! WRITE(sd_smectite_low,'(F25.10)') -0.999*dsecondary3(27)
            ! WRITE(sd_prehnite,'(F25.10)') -0.999*dsecondary3(28)
            ! WRITE(sd_chlorite,'(F25.10)') -0.999*dsecondary3(29) !!!!!!!!
            ! WRITE(sd_scolecite,'(F25.10)') -0.999*dsecondary3(30)
            ! WRITE(sd_clinochlore14a,'(F25.10)') -0.999*dsecondary3(31)
            ! WRITE(sd_clinochlore7a,'(F25.10)') -0.999*dsecondary3(32)
            ! WRITE(sd_saponite_ca,'(F25.10)') -0.999*dsecondary3(33)
            ! WRITE(sd_verm_na,'(F25.10)') -0.999*dsecondary3(34)
            ! WRITE(sd_pyrrhotite,'(F25.10)') -0.999*dsecondary3(35)
            ! WRITE(sd_fe_saponite_ca,'(F25.10)') -0.999*dsecondary3(36) !!!
            ! WRITE(sd_fe_saponite_mg,'(F25.10)') -0.999*dsecondary3(37) !!!
            ! WRITE(sd_daphnite_7a,'(F25.10)') -0.999*dsecondary3(38) !!!
            ! WRITE(sd_daphnite_14a,'(F25.10)') -0.999*dsecondary3(39) !!!
            ! WRITE(sd_epidote,'(F25.10)') -0.999*dsecondary3(40) !!!

            !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
            !       write(*,*) "sd_kaolinite" , sd_kaolinite
            !       write(*,*) "sd_saponite" , sd_saponite
            !       write(*,*) "sd_celadonite" , sd_celadonite
            !       write(*,*) "sd_clinoptilolite" , sd_clinoptilolite
            !       write(*,*) "sd_pyrite" , sd_pyrite
            !   end if

            !   WRITE(sd_mont_na,'(F25.10)') dsecondary3(6)
            !   WRITE(sd_goethite,'(F25.10)') dsecondary3(7)
            !   WRITE(sd_smectite,'(F25.10)') dsecondary3(8)
            !   WRITE(sd_calcite,'(F25.10)') dsecondary3(9)
            !   WRITE(sd_kspar,'(F25.10)') dsecondary3(10)
            !   WRITE(sd_saponite_na,'(F25.10)') dsecondary3(11) !!!!
            !   WRITE(sd_nont_na,'(F25.10)') dsecondary3(12)
            !   WRITE(sd_nont_mg,'(F25.10)') dsecondary3(13)
            !   WRITE(sd_fe_celadonite,'(F25.10)') dsecondary3(14)
            !   WRITE(sd_nont_ca,'(F25.10)') dsecondary3(15)
            !   WRITE(sd_mesolite,'(F25.10)') dsecondary3(16)
            !   WRITE(sd_hematite,'(F25.10)') dsecondary3(17)
            !   WRITE(sd_mont_ca,'(F25.10)') dsecondary3(18)
            !   WRITE(sd_verm_ca,'(F25.10)') dsecondary3(19)
            !   WRITE(sd_analcime,'(F25.10)') dsecondary3(20)
            !   WRITE(sd_phillipsite,'(F25.10)') dsecondary3(21)
            !   WRITE(sd_mont_mg,'(F25.10)') dsecondary3(22)
            !   WRITE(sd_gismondine,'(F25.10)') dsecondary3(23)
            !   WRITE(sd_verm_mg,'(F25.10)') dsecondary3(24)
            !   WRITE(sd_natrolite,'(F25.10)') dsecondary3(25)
            !   WRITE(sd_talc,'(F25.10)') dsecondary3(26) !!!!!!!!!
            !   WRITE(sd_smectite_low,'(F25.10)') dsecondary3(27)
            !   WRITE(sd_prehnite,'(F25.10)') dsecondary3(28)
            !   WRITE(sd_chlorite,'(F25.10)') dsecondary3(29) !!!!!!!!
            !   WRITE(sd_scolecite,'(F25.10)') dsecondary3(30)
            !   WRITE(sd_clinochlore14a,'(F25.10)') dsecondary3(31)
            !   WRITE(sd_clinochlore7a,'(F25.10)') dsecondary3(32)
            !   WRITE(sd_saponite_ca,'(F25.10)') dsecondary3(33)
            !   WRITE(sd_verm_na,'(F25.10)') dsecondary3(34)
            !   WRITE(sd_pyrrhotite,'(F25.10)') dsecondary3(35)
            !   WRITE(sd_fe_saponite_ca,'(F25.10)') dsecondary3(36) !!!
            !   WRITE(sd_fe_saponite_mg,'(F25.10)') dsecondary3(37) !!!
            !   WRITE(sd_daphnite_7a,'(F25.10)') dsecondary3(38) !!!
            !   WRITE(sd_daphnite_14a,'(F25.10)') dsecondary3(39) !!!
            !   WRITE(sd_epidote,'(F25.10)') dsecondary3(40) !!!






              ! OTHER INFORMATION TO STRINGS
              WRITE(s_temp,'(F25.10)') temp3
              WRITE(s_precip,'(F25.10)') medium3(2)
              WRITE(s_reactive,'(F25.10)') medium3(4)


            !-equilibrium input file
            IF (se_toggle .EQ. 1) THEN

              ! EQ solution
              inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
                   &"    units   mol/kgw" //NEW_LINE('')// &
                   &"    temp" // TRIM(s_temp) //NEW_LINE('')// &
                   &"    Ca " // TRIM(s_ca) //NEW_LINE('')// &
                   &"    Mg " // TRIM(s_mg) //NEW_LINE('')// &
                   &"    Na " // TRIM(s_na) //NEW_LINE('')// &
                   &"    K " // TRIM(s_k) //NEW_LINE('')// &
                   &"    Fe " // TRIM(s_fe) //NEW_LINE('')// &
                   &"    S "// TRIM(s_s)  //NEW_LINE('')// &
                   &"    Si " // TRIM(s_si) //NEW_LINE('')// &
                   &"    Cl " // TRIM(s_cl) //NEW_LINE('')// &
                   &"    Al " // TRIM(s_al) //NEW_LINE('')// &
                   &"    C " // TRIM(s_co2) //NEW_LINE('')// &
                   &"    Alkalinity " // TRIM(s_alk) //NEW_LINE('')// &
                   &"    -water "// TRIM(s_water) // " # kg" //NEW_LINE('')// &
                   &" "  //NEW_LINE('')

              ! EQ equilibrium phases
              if (medium3(2) .eq. precip_th) then

              inputz0 = TRIM(inputz0) // "EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
                   &"    Goethite " // TRIM(s_precip) // TRIM(s_goethite) // kinetics //NEW_LINE('')// &
                   &"    Celadonite " // TRIM(s_precip) // TRIM(s_celadonite) // kinetics //NEW_LINE('')// & ! mica
                   &"    Saponite-Mg " // TRIM(s_precip) // TRIM(s_saponite) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Pyrite " // TRIM(s_precip) // TRIM(s_pyrite) // kinetics //NEW_LINE('')// &
                   &"    Saponite-Na " // TRIM(s_precip) // TRIM(s_saponite_na) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Nontronite-Na " // TRIM(s_precip) // TRIM(s_nont_na) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Nontronite-Mg " // TRIM(s_precip) // TRIM(s_nont_mg) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Fe-Celadonite " // TRIM(s_precip) // TRIM(s_fe_celadonite) // kinetics //NEW_LINE('')// & ! mica
                   &"    Nontronite-Ca " // TRIM(s_precip) // TRIM(s_nont_ca) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Analcime " // TRIM(s_precip) // TRIM(s_analcime) // kinetics //NEW_LINE('')// & ! zeolite
                   &"    Phillipsite " // TRIM(s_precip) // TRIM(s_phillipsite) // kinetics //NEW_LINE('')// & ! zeolite
                   &"    Natrolite " // TRIM(s_precip) // TRIM(s_natrolite) // kinetics //NEW_LINE('')// & ! zeolite
                   &"    Talc " // TRIM(s_precip) // TRIM(s_talc) // kinetics //NEW_LINE('')// &
                   &"    Chlorite(14A) " // TRIM(s_precip) // TRIM(s_chlorite) // kinetics //NEW_LINE('')// & ! chlorite
                   &"    Clinochlore-14A " // TRIM(s_precip) // TRIM(s_clinochlore14a) // kinetics //NEW_LINE('')// & ! chlorite
                   &"    Clinochlore-7A " // TRIM(s_precip) // TRIM(s_clinochlore7a) // kinetics //NEW_LINE('')// & ! chlorite
                   &"    Saponite-Ca " // TRIM(s_precip) // TRIM(s_saponite_ca) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Pyrrhotite " // TRIM(s_precip) // TRIM(s_pyrrhotite) // kinetics //NEW_LINE('')//& ! sulfide
                   &"    Fe-Saponite-Ca " // TRIM(s_precip) // TRIM(s_fe_saponite_ca) // kinetics //NEW_LINE('')// & ! sap smec
                   &"    Fe-Saponite-Mg " // TRIM(s_precip) // TRIM(s_fe_saponite_mg) // kinetics //NEW_LINE('')// &! sap smec
                   ! 		!&"    Calcite " // trim(s_precip) // trim(s_calcite) // kinetics //NEW_LINE('')// & ! .135
                   &"    Montmor-Na " // TRIM(s_precip) // TRIM(s_mont_na) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Montmor-Mg " // TRIM(s_precip) // TRIM(s_mont_mg) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Montmor-Ca " // TRIM(s_precip) // TRIM(s_mont_ca) // kinetics //NEW_LINE('')// & ! smectite
                   &"    Smectite-high-Fe-Mg " // trim(s_precip) // trim(s_smectite) // kinetics //NEW_LINE('')// & ! smectite
 	 	           &"    Vermiculite-Na " // TRIM(s_precip) // TRIM(s_verm_na) // kinetics //NEW_LINE('')// & ! clay
                   &"    Vermiculite-Ca " // TRIM(s_precip) // TRIM(s_verm_ca) // kinetics //NEW_LINE('')// & ! clay
                   &"    Vermiculite-Mg " // TRIM(s_precip) // TRIM(s_verm_mg) // kinetics //NEW_LINE('')//& ! clay
                   &"    Hematite " // TRIM(s_precip) // TRIM(s_hematite) // kinetics //NEW_LINE('')//& ! iron oxide
                   &"    Epidote  " // trim(s_precip) // trim(s_epidote) // kinetics //NEW_LINE('')// &
                    		&"    Smectite-low-Fe-Mg 0.0 " // trim(s_smectite_low) // kinetics //NEW_LINE('')// & ! smectite
                   &"   Daphnite-7a " // trim(s_precip) // trim(s_daphnite_7a) // kinetics //NEW_LINE('')// & ! chlorite
              	   &"   Daphnite-14a " // trim(s_precip) // trim(s_daphnite_14a) // kinetics //NEW_LINE('')// &! chlorite
                   !&"    Kaolinite " // trim(s_precip) // trim(s_kaolinite) // kinetics //NEW_LINE('')// & ! clay
                   &"    Clinoptilolite-Ca " // trim(s_precip) // trim(s_clinoptilolite) // kinetics //NEW_LINE('')// & ! zeolite
                   !&"    K-Feldspar " // trim(s_precip) // trim(s_kspar) // kinetics //NEW_LINE('')// &
                   !&"    Mesolite " // trim(s_precip) // trim(s_mesolite) // kinetics //NEW_LINE('')// & ! zeolite
                   &"    Prehnite " // trim(s_precip) // trim(s_prehnite) // kinetics //NEW_LINE('')// &
                   "    Scolecite " // trim(s_precip) // trim(s_scolecite) // kinetics //NEW_LINE('')// & ! zeolite
                   !&"    Gismondine " // trim(s_precip) // trim(s_gismondine) // kinetics //NEW_LINE('')// & ! zeolite
                   &" "  //NEW_LINE('')

               end if

              ! EQ rates
              inputz0 = TRIM(inputz0) // "RATES" //NEW_LINE('')// &

               &"BGlass" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"	 10 base0 = 1e-10" //NEW_LINE('')// &
               &"	 20 if (ACT('Al+3') > 1e-10) then base0 = ACT('Al+3')" //NEW_LINE('')// &
               &"    30 rate0=M*110.0*(1.52e-5)*" // TRIM(param_exp_string) // "*(1.0e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))*(((ACT('H+')^3)/(ACT('Al+3')))^.33333)" //NEW_LINE('')// &
               &"    40 save rate0 * time" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &

               ! olivine
               &"Basalt1" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"    10 rate0=M*140.7*(1.52e-5)*" // TRIM(exp_ol) //"*(" //TRIM(ol_k1)//"*(ACT('H+')^"//TRIM(ol_n1)//")*exp(-("//TRIM(ol_e1)//"/.008314)*((1.0/TK) - (1.0/298.0))) + "//TRIM(ol_k2)//"*exp(-("//TRIM(ol_e2)//"/.008314)*((1.0/TK) - (1.0/298.0))))" //NEW_LINE('')// &
               &"    20 save rate0 * time" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &

               ! pyroxene
               &"Basalt2" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"    10 rate0=M*250.0*(1.52e-5)*" // TRIM(exp_pyr) //"*(" //TRIM(pyr_k1)//"*(ACT('H+')^"//TRIM(pyr_n1)//")*exp(-("//TRIM(pyr_e1)//"/.008314)*((1.0/TK) - (1.0/298.0))) + "//TRIM(pyr_k2)//"*exp(-("//TRIM(pyr_e2)//"/.008314)*((1.0/TK) - (1.0/298.0))))" //NEW_LINE('')// &
               &"    20 save rate0 * time" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &

               ! plagioclase
               &"Basalt3" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"    10 rate0=M*270.0*(1.52e-5)*" // TRIM(exp_plag) //"*" //TRIM(plag_k1)//"*(ACT('H+')^"//TRIM(plag_n1)//")*exp(-("//TRIM(plag_e1)//"/.008314)*((1.0/TK) - (1.0/298.0)))" //NEW_LINE('')// &
               &"    20 save rate0 * time" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &



               ! EQ kinetics
               &"KINETICS 1" //NEW_LINE('')// &

               &"BGlass" //NEW_LINE('')// &
               &"-f CaO .1997 SiO2 .847 Al2O3 .138 " //&
               & "Fe2O3 .149 MgO .1744 K2O .002 " //&
               & "Na2O .043" //NEW_LINE('')// &
               &"-m0 " // TRIM(s_glass) //NEW_LINE('')// &

               &"Basalt1 " //NEW_LINE('')// &
               & TRIM(param_ol_string) //NEW_LINE('')// &
               &"-m0 " // TRIM(s_basalt1) //NEW_LINE('')// &

               &"Basalt2 " //NEW_LINE('')// &
               & TRIM(param_pyr_string) //NEW_LINE('')// &
               &"-m0 " // TRIM(s_basalt2) //NEW_LINE('')// &

               &"Basalt3 " //NEW_LINE('')// &
               & TRIM(param_plag_string) //NEW_LINE('')// &
               &"-m0 " // TRIM(s_basalt3) //NEW_LINE('')// &
               &"    -step " // TRIM(s_timestep) // " in 1" //NEW_LINE('')// &

               &"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &



               &"CALCULATE_VALUES" //NEW_LINE('')// &

               &"R(sum)" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"10 sum = 1.0" //&
               &"" //NEW_LINE('')// &
               &"100 SAVE sum" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &

               &"R(s_sp)" //NEW_LINE('')// &
               &"-start" //NEW_LINE('')// &
               &"10 s_sp = 1.53e-5" //&
               &"" //NEW_LINE('')// &
               &"100 SAVE s_sp" //NEW_LINE('')// &
               &"-end" //NEW_LINE('')// &


               &"SELECTED_OUTPUT" //NEW_LINE('')// &
               &"    -reset false" //NEW_LINE('')// &
               &"    -high_precision true" //NEW_LINE('')// &
               &"    -k basalt3 basalt2 basalt1 bglass" //NEW_LINE('')// &
               &"    -ph" //NEW_LINE('')// &
               &"    -pe false" //NEW_LINE('')// &
               &"    -totals C" //NEW_LINE('')// &
               &"    -totals Ca Mg Na K Fe S Si Cl Al " //NEW_LINE('')// &
               &"    -molalities HCO3-" //NEW_LINE('')// &
               &"    -water true" //NEW_LINE('')// &
               &"    -alkalinity" //NEW_LINE('')// &
               &"    -p Kaolinite Saponite-Mg Celadonite Clinoptilolite-Ca Pyrite Montmor-Na Goethite" //NEW_LINE('')// & ! 7
               &"    -p Smectite-high-Fe-Mg Calcite K-Feldspar Saponite-Na Nontronite-Na Nontronite-Mg" //NEW_LINE('')// & ! 6
               &"    -p Fe-Celadonite Nontronite-Ca Mesolite Hematite Montmor-Ca Vermiculite-Ca Analcime" //NEW_LINE('')// & ! 7
               &"    -p Phillipsite Montmor-Mg Gismondine Vermiculite-Mg Natrolite Talc Smectite-low-Fe-Mg " //NEW_LINE('')// & ! 7
               &"    -p Prehnite Chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A Saponite-Ca" //NEW_LINE('')// & ! 6
               &"    -p Vermiculite-Na Pyrrhotite Fe-Saponite-Ca Fe-Saponite-Mg" //NEW_LINE('')// & ! 4

               &"    -p Daphnite-7a Daphnite-14a Epidote Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// & ! 6
               &"    -p prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// & ! 6
               &"    -p prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// & ! 6
               &"    -s kaolinite" //NEW_LINE('')// &		! 1
               ! 		&"    -s kaolinite saponite-mg celadonite Clinoptilolite-Ca pyrite montmor-na goethite" //NEW_LINE('')// &
               ! 		&"    -s Smectite-high-Fe-Mg calcite k-feldspar saponite-na nontronite-na nontronite-mg" //NEW_LINE('')// &
               ! 		&"    -s fe-celadonite nontronite-ca mesolite hematite montmor-ca vermiculite-ca analcime" //NEW_LINE('')// &
               ! 		&"    -s phillipsite diopside gismondine vermiculite-mg natrolite talc Smectite-low-Fe-Mg " //NEW_LINE('')// &
               ! 		&"    -s prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// &
               ! 		&"    -s vermiculite-na pyrrhotite Fe-Saponite-Ca Fe-Saponite-Mg" //NEW_LINE('')// &
               &"    -calculate_values R(sum) R(s_sp)" //NEW_LINE('')// &
               &"    -time" //NEW_LINE('')// &
               &"END"



              !-phreeq equilibrium
              id = CreateIPhreeqc()


              ! IF (SetErrorFileOn(id,.TRUE.).NE.IPQ_OK) THEN
              ! 	CALL OutputErrorString(id)
              ! 	!STOP
              ! END IF

              ! IF (seterrorstringon(id, .true.).NE.ipq_ok) THEN
              !      CALL outputerrorstring(id)
              !      stop
              !   END IF


              !-slave writes to variables

              ! IF (id.LT.0) THEN
              ! 	write(*,*) "weird stop?"
              ! 	STOP
              ! END IF

              IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
                 CALL OutputErrorString(id)
                 WRITE(*,*) "primary"
                 WRITE(*,*) primary3
                 WRITE(*,*) "secondary"
                 WRITE(*,*) secondary3
                 WRITE(*,*) "solute"
                 WRITE(*,*) solute3
                 WRITE(*,*) "medium"
                 WRITE(*,*) medium3
                 WRITE(*,*) "temp"
                 WRITE(*,*) temp3
                 !STOP
              END IF

              IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
                 CALL OutputErrorString(id)
                 !STOP
              END IF

              !IF (LoadDatabase(id, 'l5.dat').NE.0) THEN
              IF (LoadDatabaseString(id, TRIM(L5)).NE.0) THEN
                 CALL OutputErrorString(id)
                 WRITE(*,*) "primary"
                 WRITE(*,*) primary3
                 WRITE(*,*) "secondary"
                 WRITE(*,*) secondary3
                 WRITE(*,*) "solute"
                 WRITE(*,*) solute3
                 WRITE(*,*) "medium"
                 WRITE(*,*) medium3
                 WRITE(*,*) "temp"
                 WRITE(*,*) temp3
                 !STOP
              END IF

              ! RUN INPUT
              IF (RunString(id, TRIM(inputz0)).NE.0) THEN
                 WRITE(*,*) "issue is:" , RunString(id, TRIM(inputz0))
                 CALL OutputErrorString(id)
                 WRITE(*,*) "primary"
                 WRITE(*,*) primary3
                 WRITE(*,*) "secondary"
                 WRITE(*,*) secondary3
                 WRITE(*,*) "solute"
                 WRITE(*,*) solute3
                 WRITE(*,*) "medium"
                 WRITE(*,*) medium3
                 WRITE(*,*) "temp"
                 WRITE(*,*) temp3
                 IF (RunString(id, TRIM(inputz0)).NE.0) THEN
                    WRITE(*,*) "another chance 2"
                    CALL OutputErrorString(id)
                 END IF
                 !STOP
              END IF

              ! WRITE AWAY
              DO i=1,GetSelectedOutputStringLineCount(id)
                 CALL GetSelectedOutputStringLine(id, i, line)
                 ! HEADER BITS YOU MAY WANT
                !  	if (i .eq. 1) then
                !   	   write(12,*) trim(line)
                !  ! 	   !if ((medium3(6) .gt. 24000.0) .and. (medium3(7) .gt. -100.0)) then
                !   	   write(*,*) trim(line) ! PRINT LABELS FOR EVERY FIELD (USEFUL)
                !  ! 	   !end if
                !  	end if

                 ! MEAT
                 IF (i .GT. 1) THEN
                    READ(line,*) outmat(i,:)
!!!!write(12,*) outmat(i,:) ! this writes to file, which i don't need (USEFUL)
                    ! 		if ((medium3(6) .gt. 23000.0) .and. (medium3(7) .gt. -200.0)) then
                    ! 		write(*,*) i
                    ! 		write(*,*) trim(line) ! PRINT EVERY GOD DAMN LINE
                    ! 		write(*,*) ""
                    ! ! 		! write(*,*) solute3
                    ! ! ! 		write(*,*) ""
                    ! ! 		write(*,*) ""
                    ! 		end if
                 END IF
              END DO


              ! OUTPUT TO THE MAIN MASSACR METHOD
              alt0(1,:) = outmat(3,:)
              alt_mat(m,:) = alt0(1,:)

              !write(*,*) "an output alt0: ", alt0

              IF (GetSelectedOutputStringLineCount(id) .NE. 3) THEN
                 alt0(1,:) = 0.0
                 WRITE(*,*) "not 3 lines error"
              END IF


              ! DESTROY INSTANCE
              IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
                 CALL OutputErrorString(id)
                 WRITE(*,*) "cannot be destroyed error"
                 STOP
              END IF


              if (alt0(1,2) .GT. 1.0) then

              solLocal(m,:) = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
                   alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
                   alt0(1,13), alt0(1,14), alt0(1,15), 0.0/)

              priLocal(m,:) = (/ 0.0*alt0(1,136), alt0(1,127), alt0(1,129), alt0(1,131), alt0(1,133)/)

              dpriLocal(m,:) = (/ 0.0*alt0(1,136), alt0(1,128), alt0(1,130), alt0(1,132), alt0(1,134)/)
            !   if (my_id .EQ. 10) then
            !       write(*,*) "dpriLocal:" , dpriLocal(m,:)
            !   end if

              end if

              IF (alt0(1,2) .LT. 1.0) THEN
                 !medLocal(m,5) = 0.0
                 solLocal(m,:) = (/ solute3(1), solute3(2), solute3(3), solute3(4), solute3(5), &
                      solute3(6), solute3(7), solute3(8), solute3(9), solute3(10), solute3(11), &
                      solute3(12), solute3(13), solute3(14), 0.0/)
              END IF


           END IF ! end if-cell-is-on loop (medLocl 5 == 1)


           if (alt0(1,2) .GT. 1.0) then

               if (medium3(2) .eq. precip_th) then
               !secLocal = 0.0
                   DO ii=1,g_sec/2
                      secLocal(m,ii) = alt_mat(m,2*ii+14)
                      dsecLocal(m,ii) = secLocal(m,ii) - secondary3(ii)
                   END DO
               end if

           end if

            medLocal(:,1) = 0.0
            DO ii=1,g_sec/2
               medLocal(:,1) = medLocal(:,1) + secLocal(:,ii)*sec_molar(ii)/sec_density(ii)
            END DO

            DO ii=1,g_pri
               medLocal(:,1) = medLocal(:,1) + priLocal(:,ii)*pri_molar(ii)/pri_density(ii)
            END DO

            IF ((medLocal(m,1) + solLocal(m,3)) .GT. 0.0) THEN
               medLocal(m,1) = medLocal(m,1)/(medLocal(m,1) + 1000.0*solLocal(m,3))
            END IF




        END IF ! if se_toggle == 1





        !-kinetic input file
        IF (se_toggle .EQ. 0) THEN
            ! !# FE FIX
            ! if ((my_id .EQ. 34) .or. (my_id .EQ. 35)) then
            !     s_fe = "0.000009"
            ! end if



            inputz0 = "SOLUTION" //NEW_LINE('')// &
                 &"    units   mol/kgw" //NEW_LINE('')// &
                 &"    temp" // TRIM(s_temp) //NEW_LINE('')// &
                 &"    Ca " // TRIM(s_ca) //NEW_LINE('')// &
                 &"    Mg " // TRIM(s_mg) //NEW_LINE('')// &
                 &"    Na " // TRIM(s_na) //NEW_LINE('')// &
                 &"    K " // TRIM(s_k) //NEW_LINE('')// &
                 &"    Fe " // TRIM(s_fe) //NEW_LINE('')// &
                 &"    S "// TRIM(s_s)  //NEW_LINE('')// &
                 &"    Si " // TRIM(s_si) //NEW_LINE('')// &
                 &"    Cl " // TRIM(s_cl) //NEW_LINE('')// &
                 &"    Al " // TRIM(s_al) //NEW_LINE('')// &
                 &"    C " // TRIM(s_co2) //NEW_LINE('')// &
                 &"    Alkalinity " // TRIM(s_alk) //NEW_LINE('')// &
                 &"    -water "// TRIM(s_water) // " # kg" //NEW_LINE('')// &

                 &" "  //NEW_LINE('')


                 inputz0 = TRIM(inputz0) // "RATES" //NEW_LINE('')// &

                  &"BGlass" //NEW_LINE('')// &
                  &"-start" //NEW_LINE('')// &
                  &"    10 save " // TRIM(sd_dglass) //NEW_LINE('')// &
                  &"-end" //NEW_LINE('')// &

                  ! olivine
                  &"Basalt1" //NEW_LINE('')// &
                  &"-start" //NEW_LINE('')// &
                  &"    10 save " // TRIM(sd_dbasalt1) //NEW_LINE('')// &
                  &"-end" //NEW_LINE('')// &

                  ! pyroxene
                  &"Basalt2" //NEW_LINE('')// &
                  &"-start" //NEW_LINE('')// &
                  &"    10 save " // TRIM(sd_dbasalt2) //NEW_LINE('')// &
                  &"-end" //NEW_LINE('')// &

                  ! plagioclase
                  &"Basalt3" //NEW_LINE('')// &
                  &"-start" //NEW_LINE('')// &
                  &"    10 save " // TRIM(sd_dbasalt3) //NEW_LINE('')// &
                  &"-end" //NEW_LINE('')// &



                ! kaolinite
                &"Kaolinite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_kaolinite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! saponite-mg
                &"Saponite-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_saponite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! celadonite
                &"Celadonite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_celadonite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! clinoptilolite
                &"Clinoptilolite-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_clinoptilolite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! pyrite
                &"Pyrite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_pyrite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &


                ! mont-na
                &"Montmor-Na" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_mont_na) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! goethite
                &"Goethite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_goethite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! smectite
                &"Smectite-high-Fe-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_smectite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! calcite
                &"Calcite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_calcite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! kspar
                &"K-Feldspar" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_kspar) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! saponite-na
                &"Saponite-Na" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_saponite_na) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! nont na
                &"Nontronite-Na" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_nont_na) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! nont mg
                &"Nontronite-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_nont_mg) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! fe celad
                &"Fe-Celadonite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_fe_celadonite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! nont ca
                &"Nontronite-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_nont_ca) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! mesolite
                &"Mesolite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_mesolite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! hematite
                &"Hematite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_hematite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! mont ca
                &"Montmor-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_mont_ca) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! verm ca
                &"Vermiculite-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_verm_ca) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! analcime
                &"Analcime" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_analcime) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! phillipsite
                &"Phillipsite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_phillipsite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! mont mg
                &"Montmor-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_mont_mg) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! gismondine
                &"Gismondine" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_gismondine) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! verm mg
                &"Vermiculite-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_verm_mg) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! natrolite
                &"Natrolite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_natrolite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! talc
                &"Talc" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_talc) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! smectite low
                &"Smectite-low-Fe-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_smectite_low) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! prehnite
                &"Prehnite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_prehnite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! chlorite
                &"Chlorite(14a)" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_chlorite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! scolecite
                &"Scolecite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_scolecite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! clinochlorte
                &"Clinochlore-14A" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_clinochlore14a) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! Clinochlore-7A
                &"Clinochlore-7A" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_clinochlore7a) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! saponite_ca
                &"Saponite-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_saponite_ca) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! verm-na
                &"Vermiculite-Na" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_verm_na) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! pyrrhotite
                &"Pyrrhotite" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_pyrrhotite) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! fe-sap-ca
                &"Fe-Saponite-Ca" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_fe_saponite_ca) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! fe sap mg
                &"Fe-Saponite-Mg" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_fe_saponite_mg) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! daphnite-7a
                &"Daphnite-7a" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_daphnite_7a) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! daphnite-14a
                &"Daphnite-14a" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_daphnite_14a) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &

                ! epidote
                &"Epidote" //NEW_LINE('')// &
                &"-start" //NEW_LINE('')// &
                &"    10 save " // TRIM(sd_epidote) //NEW_LINE('')// &
                &"-end" //NEW_LINE('')// &




                  ! EQ kinetics
                  &"KINETICS" //NEW_LINE('')// &

                  &"BGlass" //NEW_LINE('')// &
                  &"-f CaO .1997 SiO2 .847 Al2O3 .138 " //&
                  & "Fe2O3 .149 MgO .1744 K2O .002 " //&
                  & "Na2O .043" //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_glass) //NEW_LINE('')// &

                  &"Basalt1 " //NEW_LINE('')// &
                  & TRIM(param_ol_string) //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_basalt1) //NEW_LINE('')// &

                  &"Basalt2 " //NEW_LINE('')// &
                  & TRIM(param_pyr_string) //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_basalt2) //NEW_LINE('')// &

                  &"Basalt3 " //NEW_LINE('')// &
                  & TRIM(param_plag_string) //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_basalt3) //NEW_LINE('')// &


                  &"Kaolinite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_kaolinite) //NEW_LINE('')// &

                  &"Saponite-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_saponite) //NEW_LINE('')// &

                  &"Celadonite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_celadonite) //NEW_LINE('')// &

                  &"Clinoptilolite-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_clinoptilolite) //NEW_LINE('')// &

                  &"Pyrite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_pyrite) //NEW_LINE('')// &


                  &"Montmor-Na " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_mont_na) //NEW_LINE('')// &

                  &"Goethite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_goethite) //NEW_LINE('')// &

                  &"Smectite-high-Fe-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_smectite) //NEW_LINE('')// &

                  &"Calcite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_calcite) //NEW_LINE('')// &

                  &"K-Feldspar " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_kspar) //NEW_LINE('')// &

                  &"Saponite-Na " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_saponite_na) //NEW_LINE('')// &

                  &"Nontronite-Na " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_nont_na) //NEW_LINE('')// &

                  &"Nontronite-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_nont_mg) //NEW_LINE('')// &

                  &"Fe-Celadonite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_fe_celadonite) //NEW_LINE('')// &

                  &"Nontronite-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_nont_ca) //NEW_LINE('')// &

                  &"Mesolite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_mesolite) //NEW_LINE('')// &

                  &"Hematite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_hematite) //NEW_LINE('')// &

                  &"Montmor-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_mont_ca) //NEW_LINE('')// &

                  &"Vermiculite-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_verm_ca) //NEW_LINE('')// &

                  &"Analcime " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_analcime) //NEW_LINE('')// &

                  &"Phillipsite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_phillipsite) //NEW_LINE('')// &

                  &"Montmor-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_mont_mg) //NEW_LINE('')// &

                  &"Gismondine " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_gismondine) //NEW_LINE('')// &

                  &"Vermiculite-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_verm_mg) //NEW_LINE('')// &

                  &"Natrolite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_natrolite) //NEW_LINE('')// &

                  &"Talc " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_talc) //NEW_LINE('')// &

                  &"Smectite-low-Fe-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_smectite_low) //NEW_LINE('')// &

                  &"Prehnite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_prehnite) //NEW_LINE('')// &

                  &"Chlorite(14a) " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_chlorite) //NEW_LINE('')// &

                  &"Scolecite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_scolecite) //NEW_LINE('')// &

                  &"Clinochlore-14a " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_clinochlore14a) //NEW_LINE('')// &

                  &"Clinochlore-7a " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_clinochlore7a) //NEW_LINE('')// &

                  &"Saponite-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_saponite_ca) //NEW_LINE('')// &

                  &"Vermiculite-Na " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_verm_na) //NEW_LINE('')// &

                  &"Pyrrhotite " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_pyrrhotite) //NEW_LINE('')// &

                  &"Fe-Saponite-Ca " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_fe_saponite_ca) //NEW_LINE('')// &

                  &"Fe-Saponite-Mg " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_fe_saponite_mg) //NEW_LINE('')// &

                  &"Daphnite-7a " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_daphnite_7a) //NEW_LINE('')// &

                  &"Daphnite-14a " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_daphnite_14a) //NEW_LINE('')// &

                  &"Epidote " //NEW_LINE('')// &
                  &"-m0 " // TRIM(s_epidote) //NEW_LINE('')// &

                  &"    -step "//TRIM(s_timestep)//" in 1" //NEW_LINE('')// &

                  &"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &



                  &"SELECTED_OUTPUT" //NEW_LINE('')// &
                  &"    -reset false" //NEW_LINE('')// &
                  &"    -high_precision true" //NEW_LINE('')// &
                   &"    -k basalt3 basalt2 basalt1 bglass" //NEW_LINE('')// &
                  &"    -ph" //NEW_LINE('')// &
                  &"    -pe false" //NEW_LINE('')// &
                  &"    -totals C" //NEW_LINE('')// &
                  &"    -totals Ca Mg Na K Fe S Si Cl Al " //NEW_LINE('')// &
                  &"    -molalities HCO3-" //NEW_LINE('')// &
                  &"    -water true" //NEW_LINE('')// &
                  &"    -alkalinity" //NEW_LINE('')// &
                  &"    -k Kaolinite Saponite-Mg Celadonite Clinoptilolite-Ca Pyrite Montmor-Na Goethite" //NEW_LINE('')// & ! 7
                  &"    -k Smectite-high-Fe-Mg Calcite K-Feldspar Saponite-Na Nontronite-Na Nontronite-Mg" //NEW_LINE('')// & ! 6
                  &"    -k Fe-Celadonite Nontronite-Ca Mesolite Hematite Montmor-Ca Vermiculite-Ca Analcime" //NEW_LINE('')// & ! 7
                  &"    -k Phillipsite Montmor-Mg Gismondine Vermiculite-Mg Natrolite Talc Smectite-low-Fe-Mg " //NEW_LINE('')// & ! 7
                  &"    -k Prehnite Chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A Saponite-Ca" //NEW_LINE('')// & ! 6
                  &"    -k Vermiculite-Na Pyrrhotite Fe-Saponite-Ca Fe-Saponite-Mg" //NEW_LINE('')// & ! 4
                   &"    -k Daphnite-7a Daphnite-14a Epidote" //NEW_LINE('')// & ! 3



                !   &"    -p prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// & ! 6
                !   &"    -p prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// & ! 6
                !   &"    -s kaolinite" //NEW_LINE('')// &		! 1
                  ! 		&"    -s kaolinite saponite-mg celadonite Clinoptilolite-Ca pyrite montmor-na goethite" //NEW_LINE('')// &
                  ! 		&"    -s Smectite-high-Fe-Mg calcite k-feldspar saponite-na nontronite-na nontronite-mg" //NEW_LINE('')// &
                  ! 		&"    -s fe-celadonite nontronite-ca mesolite hematite montmor-ca vermiculite-ca analcime" //NEW_LINE('')// &
                  ! 		&"    -s phillipsite diopside gismondine vermiculite-mg natrolite talc Smectite-low-Fe-Mg " //NEW_LINE('')// &
                  ! 		&"    -s prehnite chlorite(14a) scolecite Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// &
                  ! 		&"    -s vermiculite-na pyrrhotite Fe-Saponite-Ca Fe-Saponite-Mg" //NEW_LINE('')// &
                  &"    -time" //NEW_LINE('')// &
                  &"END"

                !-kinetic print statements


                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "plag pyr ol glass"
                !       write(*,*) s_basalt3, s_basalt2, s_basalt1, s_glass
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "sd_ plag pyr ol glass"
                !       write(*,*) sd_dbasalt3, sd_dbasalt2, sd_dbasalt1, sd_dglass
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "amount sap" , s_saponite
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "sd_saponite" , sd_saponite
                !   end if




                  !-kinetic PHREEQ
                  id = CreateIPhreeqc()



                  IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
                     CALL OutputErrorString(id)
                     WRITE(*,*) "primary"
                     WRITE(*,*) primary3
                     WRITE(*,*) "secondary"
                     WRITE(*,*) secondary3
                     WRITE(*,*) "solute"
                     WRITE(*,*) solute3
                     WRITE(*,*) "medium"
                     WRITE(*,*) medium3
                     WRITE(*,*) "temp"
                     WRITE(*,*) temp3
                     !STOP
                  END IF


                  IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
                     CALL OutputErrorString(id)
                     !STOP
                  END IF




                  if (my_id .EQ. 15) then
                      call system_clock(counti, count_rate, count_max)
                  end if


                  !IF (LoadDatabase(id, 'l5.dat').NE.0) THEN
                  IF (LoadDatabaseString(id, TRIM(L5)).NE.0) THEN
                     CALL OutputErrorString(id)
                     WRITE(*,*) "primary"
                     WRITE(*,*) primary3
                     WRITE(*,*) "secondary"
                     WRITE(*,*) secondary3
                     WRITE(*,*) "solute"
                     WRITE(*,*) solute3
                     WRITE(*,*) "medium"
                     WRITE(*,*) medium3
                     WRITE(*,*) "temp"
                     WRITE(*,*) temp3
                     !STOP
                  END IF

                  if (my_id .EQ. 15) then
                      call system_clock(countf, count_rate, count_max)
                      write(*,*) "      END LOAD DB" , countf - counti
                  end if

                !   if (my_id .EQ. 34) then
                !       write(*,*) "34b jjj" , jjj , slave_vector(jjj) , medium3(6) , medium3(7)
                !   end if
                  !
                !   if (my_id .EQ. 35) then
                !       write(*,*) "35b jjj" , jjj , slave_vector(jjj) , medium3(6) , medium3(7)
                !   end if


                if (my_id .EQ. 15) then
                    call system_clock(counti, count_rate, count_max)
                end if

                  ! RUN INPUT
                  IF (RunString(id, TRIM(inputz0)).NE.0) THEN
                     WRITE(*,*) "issue is:" , RunString(id, TRIM(inputz0))
                     CALL OutputErrorString(id)
                     WRITE(*,*) "primary"
                     WRITE(*,*) primary3
                     WRITE(*,*) "secondary"
                     WRITE(*,*) secondary3
                     WRITE(*,*) "solute"
                     WRITE(*,*) solute3
                     WRITE(*,*) "medium"
                     WRITE(*,*) medium3
                     WRITE(*,*) "temp"
                     WRITE(*,*) temp3
                     IF (RunString(id, TRIM(inputz0)).NE.0) THEN
                        WRITE(*,*) "another chance 2"
                        CALL OutputErrorString(id)
                     END IF
                     !STOP
                  END IF

                  if (my_id .EQ. 15) then
                      call system_clock(countf, count_rate, count_max)
                      write(*,*) "      END RUNSTRING" , countf - counti
                  end if

                !   if (my_id .EQ. 34) then
                !       write(*,*) "34 after jjj" , jjj , slave_vector(jjj)
                !   end if
                  !
                !   if (my_id .EQ. 35) then
                !       write(*,*) "35 after jjj" , jjj , slave_vector(jjj)
                !   end if


                  ! WRITE AWAY
                  DO i=1,GetSelectedOutputStringLineCount(id)
                     CALL GetSelectedOutputStringLine(id, i, line)
                     ! HEADER BITS YOU MAY WANT
                    !  if (my_id .EQ. 10) then
                    !       	if (i .eq. 1) then
                    !        	   write(12,*) trim(line)
                    !        	   write(*,*) trim(line) ! PRINT LABELS FOR EVERY FIELD (USEFUL)
                    !       	end if
                    ! end if

                     ! MEAT
                     IF (i .GT. 1) THEN
                        READ(line,*) outmat(i,1:103)
    !!!!write(12,*) outmat(i,:) ! this writes to file, which i don't need (USEFUL)
                        ! 		if ((medium3(6) .gt. 23000.0) .and. (medium3(7) .gt. -200.0)) then
                        ! 		write(*,*) i
                        ! 		write(*,*) trim(line) ! PRINT EVERY GOD DAMN LINE
                        ! 		write(*,*) ""
                        ! ! 		! write(*,*) solute3
                        ! ! ! 		write(*,*) ""
                        ! ! 		write(*,*) ""
                        ! 		end if
                     END IF
                  END DO





                  ! OUTPUT TO THE MAIN MASSACR METHOD
                  alt0(1,1:103) = outmat(3,1:103)
                  alt_mat(m,1:103) = alt0(1,1:103)

                !   if (my_id .EQ. 10) then
                !       write(*,*) "all of outmat : " , alt_mat(m,1:103)
                !   end if

                  !write(*,*) "an output alt0: ", alt0

                  IF (GetSelectedOutputStringLineCount(id) .NE. 3) THEN
                     alt0(1,:) = 0.0
                     WRITE(*,*) "not 3 lines error"
                  END IF


                  ! DESTROY INSTANCE
                  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
                     CALL OutputErrorString(id)
                     WRITE(*,*) "cannot be destroyed error"
                     STOP
                  END IF

                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "BEFORE priLocal:" , priLocal(m,2:)
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "BEFORE dpriLocal:" , dpriLocal(m,2:)
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "BEFORE secLocal:" , secLocal(m,2:5)
                !   end if
                  !
                !   if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !       write(*,*) "BEFORE dsecLocal:" , dsecLocal(m,2:5)
                !   end if


                  if (alt0(1,2) .GT. 1.0) then

                  solLocal(m,:) = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
                       alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
                       alt0(1,13), alt0(1,14), alt0(1,15), 0.0/)

                !   priLocal(m,:) = (/ 0.0*alt0(1,136), alt0(1,127), alt0(1,129), alt0(1,131), alt0(1,133)/)
                priLocal(m,:) = (/ 0.0*alt0(1,16), alt0(1,16), alt0(1,18), alt0(1,20), alt0(1,22)/)
                !priLocal(m,:) = (/0.0, primary3(2) + dprimary3(2), primary3(3) + dprimary3(3), primary3(4) + dprimary3(4), primary3(5) + dprimary3(5)/)

                ! if (my_id .EQ. 10) then
                !     write(*,*) "priLocal:" , priLocal(m,:)
                ! end if

                !   dpriLocal(m,:) = (/ 0.0*alt0(1,136), alt0(1,128), alt0(1,130), alt0(1,132), alt0(1,134)/)
                dpriLocal(m,:) = (/ 0.0*alt0(1,16), alt0(1,17), alt0(1,19), alt0(1,21), alt0(1,23)/)
                !dpriLocal(m,:) = (/0.0, dprimary3(2), dprimary3(3), dprimary3(4), dprimary3(5)/)


                  end if

                  IF (alt0(1,2) .LT. 1.0) THEN
                     !medLocal(m,5) = 0.0
                     solLocal(m,:) = (/ solute3(1), solute3(2), solute3(3), solute3(4), solute3(5), &
                          solute3(6), solute3(7), solute3(8), solute3(9), solute3(10), solute3(11), &
                          solute3(12), solute3(13), solute3(14), 0.0/)
                  END IF


               if (alt0(1,2) .GT. 1.0) then

                   !if (medium3(2) .eq. precip_th) then
                   !secLocal = 0.0
                       DO ii=1,g_sec/2
                          secLocal(m,ii) = alt_mat(m,2*ii+22)
                          !secLocal(m,ii) = secondary3(ii) + dsecondary3(ii)
                          dsecLocal(m,ii) = alt_mat(m,2*ii+23)!secLocal(m,ii) - secondary3(ii)
                          !dsecLocal(m,ii) = dsecondary3(ii)
                       END DO
                   !end if

               end if

                medLocal(:,1) = 0.0
                DO ii=1,g_sec/2
                   medLocal(:,1) = medLocal(:,1) + secLocal(:,ii)*sec_molar(ii)/sec_density(ii)
                END DO

                DO ii=1,g_pri
                   medLocal(:,1) = medLocal(:,1) + priLocal(:,ii)*pri_molar(ii)/pri_density(ii)
                END DO

                IF ((medLocal(m,1) + solLocal(m,3)) .GT. 0.0) THEN
                   medLocal(m,1) = medLocal(m,1)/(medLocal(m,1) + 1000.0*solLocal(m,3))
                END IF


                ! if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !     write(*,*) "AFTER priLocal:" , priLocal(m,2:)
                ! end if
                !
                ! if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !     write(*,*) "AFTER dpriLocal:" , dpriLocal(m,2:)
                ! end if
                !
                ! if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !     write(*,*) "AFTER secLocal:" , secLocal(m,2:5)
                ! end if
                !
                ! if ((my_id .EQ. 10) .AND. (jjj .EQ. 1)) then
                !     write(*,*) "AFTER dsecLocal:" , dsecLocal(m,2:5)
                ! end if



                  !priLocal(m,:) = primary3
                  !secLocal(m,:) = secondary3
                  !solLocal(m,:) = solute3



        END IF ! if se_toggle == 0



    END DO ! end jjj = 1,end_loop



    call system_clock(countf, count_rate, count_max)
    write(*,*) "my_id" , my_id , "jjj loop time:" , countf - counti

    !#GEOCHEM: slave sends to master

    CALL MPI_SEND( end_loop, 1, MPI_INTEGER, root_process, return_data_tag, MPI_COMM_WORLD, ierr)

    CALL MPI_SEND( slave_vector, end_loop, MPI_INTEGER, root_process, return_data_tag, MPI_COMM_WORLD, ierr)

    ! send primary array chunk back to root process
    DO ii = 1,g_pri
       CALL MPI_SEND( priLocal(slave_vector(1:end_loop),ii), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    END DO

    ! send secondary array chunk back to root process
    DO ii = 1,g_sec/2
       CALL MPI_SEND( secLocal(slave_vector(1:end_loop),ii), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    END DO

    ! ! send solute array chunk back to root process
    ! CALL MPI_SEND( solLocal(slave_vector(1:end_loop),3), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    !
    ! CALL MPI_SEND( solLocal(slave_vector(1:end_loop),2), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    !
    ! CALL MPI_SEND( solLocal(slave_vector(1:end_loop),1), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)

    DO ii = 1,g_sol
       CALL MPI_SEND( solLocal(slave_vector(1:end_loop),ii), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    END DO

    ! send medium array chunk back to root process
    DO ii = 1,g_med
       CALL MPI_SEND( medLongBitFull(slave_vector(1:end_loop),ii), end_loop, MPI_REAL4, root_process, return_data_tag, MPI_COMM_WORLD, ierr)
    END DO



    END DO ! end do jj = tn/mstep





  END IF ! end loop through processors




  ! close up shop
  CALL MPI_FINALIZE ( ierr )


END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: computes the 2D temperature profile for the current timestep
!
! INPUTS: h(xn,yn) : temperature profile of previous timestep
!         psi(xn,yn) : 2D streamfunction array
!         rho_in(xn,yn) : 2D density array
!         flux(xn,2) : top and bottom heat flux boundary conditions
!
! RETURNS: h_next(xn,yn) : temperature profile of current timestep
!
! ----------------------------------------------------------------------------------%%

FUNCTION h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! interface
  !
  ! 	function velocities(psi)
  ! 		use globals
  ! 		use initialize
  ! 		implicit none
  ! 		real(4) :: velocities(xn,2*yn), psi(xn,yn)
  ! 		real(4) :: u0(xn,yn), v0(xn,yn)
  ! 	end function velocities
  !
  ! end interface

  ! declare errthing

  ! integers
  INTEGER :: i, j, n, ii, m=3
  ! inputs
  REAL(4) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2), phi_in(xn,yn)
  ! velocity stuff
  REAL(4) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
  REAL(4) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
  REAL(4) ::  velocities0(xn,2*yn)
  ! matrix stuff
  REAL(4) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
  ! real(4) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  REAL(4) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
  ! real(4) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  REAL(4) :: h0(xn,yn), hMid(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
  REAL(4) :: kMatLong((xn-2)*(yn-2))
  REAL(4) :: mn(xn,yn)
  REAL(4) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
  REAL(4) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
  REAL(4) :: frac6_in(yn,2), temp6_in(yn,2), dt_in


  ! ! calculate velocities from streamfunction values
  ! velocities0 = velocities(psi)
  ! uf = velocities0(1:xn,1:yn)
  ! vf = velocities0(1:xn,yn+1:2*yn)
  u = -1.0*u_in
  v = -1.0*v_in
  ! uf = -1.0*uf!*rho_in
  ! vf = -1.0*vf!*rho_in


  ! do ii=2,yn-1
  ! 	do i=1,xn
  ! 		if ((maskP(i,ii) .eq. 2.5)) then
  ! 			u(i,ii) = 0.0
  ! 			v(i,ii) = 0.0
  ! 		end if
  ! 		if ((maskP(i,ii) .eq. 7.5)) then
  ! 			u(i,ii) = 0.0
  ! 			v(i,ii) = 0.0
  ! 		end if
  ! 	end do
  ! end do

  ! do ii=2,yn
  ! 	do i=1,xn
  ! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i,ii-1) .ne. 0.0)) then
  ! 			u(i,ii) = (uf(i,ii) + uf(i,ii-1))/(2.0)
  ! 		end if
  ! 	end do
  ! end do
  !
  ! do ii=1,yn
  ! 	do i=2,xn-2
  ! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i-1,ii) .ne. 0.0)) then
  ! 			v(i,ii) = (vf(i,ii) + vf(i-1,ii))/(2.0)
  ! 		end if
  ! 	end do
  ! end do
  !
  !
  !
  ! do ii=1,yn
  ! 	do i=2,xn-2
  ! 		if ((maskP(i,ii) .eq. 9.0)) then
  ! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i+1,ii) - psi(i,ii))/dx
  ! 		end if
  ! 		if ((maskP(i,ii) .eq. 3.0)) then
  ! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i,ii) - psi(i-1,ii))/dx
  ! 		end if
  ! 	end do
  ! end do
  !
  ! do ii=1,yn
  ! 	do i=2,xn-2
  ! 		if ((maskP(i,ii) .eq. 6.0)) then
  ! 			v(i,ii) = (frac6_in(ii,2) - frac6_in(ii,1))/(2.0*param_f_dx)
  ! 		end if
  ! 	end do
  ! end do

  ! do ii=1,yn
  ! 	do i=2,xn-2
  ! 		if ((maskP(i,ii) .eq. 6.0)) then
  ! 			v(i,ii) = -1.0*(frac6_in(ii,2) - frac6_in(ii,1))/(rho_in(i,ii)*2.0*param_f_dx)
  ! 		end if
  ! 	end do
  ! end do

  ! ! fracture vertical velocity
  ! v(xn-1,:) = v(xn-1,:)*param_f_dx/dx






  uLong = -1.0*RESHAPE(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  vLong = -1.0*RESHAPE(TRANSPOSE(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))


  mn = h



  h0 = h
  hMid = h

  qx = dt/(dx)
  qy = dt/(dy)
  sx = (2.0*dt_in*lambdaMat(1,1))/(dx*dx*rho_fluid*cp)
  sy = (2.0*dt_in*lambdaMat(1,1))/(dy*dy*rho_fluid*cp)

  ! qxMat = 1.0*dt*4179.0*phi_in/(dx*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
  ! qyMat = 1.0*dt*4179.0*phi_in/(dy*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
  ! ! sxMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
  ! ! syMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
  ! sxMat = 2.0*dt*(lambdaMat)/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
  ! syMat = 2.0*dt*(lambdaMat)/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)


  qxMat = dt_in*4179.0/(dx*cp)
  qyMat = dt_in*4179.0/(dy*cp)
  sxMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*rho_in*cp)
  syMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*rho_in*cp)


  ! ! print stability conditions at each timestep
  ! write(*,*) " "
  ! write(*,*) "velocity check"
  ! write(*,"(F10.5)") maxval(abs(u))*maxval(abs(qxMat))
  ! write(*,"(F10.5)") maxval(abs(v))*maxval(abs(qyMat))
  ! write(*,*) "conduction check"
  ! write(*,"(F10.5)") maxval(abs(syMat))
  ! write(*,*) " "
  !
  ! write(*,*) "velocity check"
  ! write(*,*) maxval(abs(u))
  ! write(*,*) maxval(abs(v))

  ! vertical boundary conditions

  ! stretch
  stretch = h0
  stretch(2,:) = stretch(1,:)
  stretch(xn-1,:) = stretch(xn,:)
  DO ii=2,yn-2
     DO i=2,xn-1
        IF ((mask(i,ii) .EQ. 5.0) .OR. (mask(i,ii) .EQ. 12.5) ) THEN
           stretch(i,ii) = stretch(i+1,ii)
        END IF

        IF ((mask(i,ii) .EQ. 3.0) .OR. (mask(i,ii) .EQ. 3.5) .OR. (mask(i,ii) .EQ. 3.1) .OR. (mask(i,ii) .EQ. 3.05)) THEN
           stretch(i,ii) = temp6_in(ii,1)
        END IF

        IF ((mask(i,ii) .EQ. 10.0) .OR. (mask(i,ii) .EQ. 17.5)) THEN
           stretch(i,ii) = stretch(i-1,ii)
        END IF

        IF ((mask(i,ii) .EQ. 6.0) .OR. (mask(i,ii) .EQ. 6.5) .OR. (mask(i,ii) .EQ. 6.1) .OR. (mask(i,ii) .EQ. 6.05)) THEN
           stretch(i,ii) = temp6_in(ii,2)
        END IF
     END DO
  END DO

  stretchLong = RESHAPE(stretch(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

  DO ii=2,yn-1
     DO i=2,xn-1
	! right outcrop (left boundary)
	IF ((mask(i,ii) .EQ. 17.5)) THEN
           h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	END IF
	IF ((mask(i,ii) .EQ. 10.0)) THEN
           h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	END IF
	IF ((mask(i,ii) .EQ. 6.0) .OR. (mask(i,ii) .EQ. 6.5) .OR. (mask(i,ii) .EQ. 6.1) .OR. (mask(i,ii) .EQ. 6.05)) THEN
           h(i,ii) = h(i,ii) + temp6_in(ii,2)*sxMat(i,ii)/2.0
	END IF

 ! left outcrop (right boundary)
	IF ((mask(i,ii) .EQ. 12.5)) THEN
           h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	END IF
	IF ((mask(i,ii) .EQ. 5.0)) THEN
           h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	END IF
	IF ((mask(i,ii) .EQ. 3.0) .OR. (mask(i,ii) .EQ. 3.5) .OR. (mask(i,ii) .EQ. 3.1) .OR. (mask(i,ii) .EQ. 3.05)) THEN
           h(i,ii) = h(i,ii) + temp6_in(ii,1)*sxMat(i,ii)/2.0
	END IF



     END DO

     IF (mask(1,ii) .EQ. 1.0) THEN
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! left
     END IF
     IF (mask(1,ii) .EQ. 25.0) THEN
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! top left corner
     END IF

     IF (mask(xn,ii) .EQ. 1.0) THEN
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! right
     END IF
     IF (mask(xn,ii) .EQ. 50.0) THEN
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! top right corner
     END IF
  END DO



  !h(2,2) = h(2,2) + h0(1,2)*sxMat(2,2)/2.0  ! bottom left corner
  !h(xn-1,2) = h(xn-1,2) + h0(xn,2)*sxMat(xn-1,2)/2.0  ! bottom right corner

  uVec = RESHAPE(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  h0Long = RESHAPE(h0(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  sxLong = RESHAPE(sxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  syLong = RESHAPE(syMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  qxLong = RESHAPE(qxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
  qyLong = RESHAPE(qyMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

  ! make the band
  aBand = 0.0
  aBand(:,1) = -sxLong/2.0 - uLong*qxLong/2.0
  aBand(:,2) = 1.0+sxLong
  aBand(:,3) = -sxLong/2.0 + uLong*qxLong/2.0

  ! bottom left corner, 2 and 3
  aBand(1,1) =  0.0
  aBand(1,2) = 1.0 + sxLong(1)/1.0 - uLong(1)*qxLong(1)!/2.0
  aBand(1,3) = -sxLong(1)/2.0 + uLong(1)*qxLong(1)!/2.0

  ! top right corner, 1 and 2
  aBand((xn-2)*(yn-2),1) = -sxLong((xn-2)*(yn-2))/2.0 - uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
  aBand((xn-2)*(yn-2),2) = 1.0 + sxLong((xn-2)*(yn-2))/1.0 + uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
  aBand((xn-2)*(yn-2),3) =  0.0



  DO i = 2,(xn-2)*(yn-2)-1

     ! flow left anywhere, 2 and 3
     IF (uLong(i) .LT. 0.0) THEN

        aBand(i,1) = -sxLong(i)/2.0
        aBand(i,2) = 1.0+sxLong(i) - uLong(i)*qxLong(i)
        aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)

     END IF


     ! flow right anywhere, 1 and 2
     IF (uLong(i) .GT. 0.0) THEN

        aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
        aBand(i,2) = 1.0+sxLong(i) + uLong(i)*qxLong(i)
        aBand(i,3) = -sxLong(i)/2.0

     END IF



     ! left edges, default 2 and 3
     IF (((MOD(i-1,xn-2).EQ.0)) .OR. (maskLong(i).EQ.10.0) .OR. (maskLong(i).EQ.17.5) .OR. (maskLong(i).EQ.6.0) .OR. (maskLong(i).EQ.6.5) .OR. (maskLong(i).EQ.6.1) .OR. (maskLong(i).EQ.6.05)) THEN
        aBand(i,1) =  0.0
        aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
        aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)
     END IF



     ! left edge flowing to the right, 1 & 2
     IF (uLong(i) .GT. 0.0) THEN

        ! left edge of right outcrop
        IF ((maskLong(i).EQ.10.0) .OR. (maskLong(i).EQ.17.5) .OR. (maskLong(i).EQ.6.0) .OR. (maskLong(i).EQ.6.5) .OR. (maskLong(i).EQ.6.1) .OR. (maskLong(i).EQ.6.05)) THEN
           aBand(i,1) =  0.0
           aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
           aBand(i,3) = -sxLong(i)/2.0
           uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
        END IF

        ! left edge but not uppper left corner
        IF ((MOD(i-1,xn-2).EQ.0) .AND. (maskLong(i) .NE. 25.0)) THEN
           aBand(i,1) =  0.0
           aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
           aBand(i,3) = -sxLong(i)/2.0
           uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
        END IF


     END IF


     ! right edge, 1 and 2 by default
     IF ((MOD(i,xn-2) .EQ. 0) .OR. (maskLong(i).EQ.5.0) .OR. (maskLong(i).EQ.12.5) .OR. (maskLong(i).EQ.3.0) .OR. (maskLong(i).EQ.3.5) .OR. (maskLong(i).EQ.3.1) .OR. (maskLong(i).EQ.3.05)) THEN
        aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
        aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
        aBand(i,3) =  0.0
     END IF


     ! right edge flowing to the left, 2 and 3
     IF (uLong(i) .LT. 0.0) THEN

        ! right edge of left outcrop
        IF ((maskLong(i).EQ.5.0) .OR. (maskLong(i).EQ.12.5) .OR. (maskLong(i).EQ.3.0) .OR. (maskLong(i).EQ.3.5) .OR. (maskLong(i).EQ.3.1) .OR. (maskLong(i).EQ.3.05)) THEN
           aBand(i,1) = -sxLong(i)/2.0
           aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
           aBand(i,3) =  0.0
           uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
        END IF

        ! right edge but not upper right corner
        IF ((MOD(i,xn-2) .EQ. 0) .AND. (maskLong(i) .NE. 25.0)) THEN
           aBand(i,1) = -sxLong(i)/2.0
           aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
           aBand(i,3) =  0.0
           uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
        END IF


     END IF


     ! upper left corner, default 2 and 3
     IF ((MOD(i-1,xn-2).EQ.0) .AND. (maskLong(i) .EQ. 25.0)) THEN
        aBand(i,1) =  0.0
        aBand(i,2) = 1.0 + sxLong(i)/1.0 - uLong(i)*qxLong(i)!/2.0
        aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)!/2.0
     END IF

     ! upper right corner, default 1 and 2
     IF ((MOD(i,xn-2).EQ.0) .AND. (maskLong(i) .EQ. 25.0)) THEN
        aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
        aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
        aBand(i,3) =  0.0
     END IF

     ! bottom right corner, 1 and 2
     IF (i.EQ.xn-2) THEN
        aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
        aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
        aBand(i,3) =  0.0
     END IF






     !

  END DO

  DO i=1,(xn-2)*(yn-2)
     ! mask
     IF ((maskLong(i) .EQ. 0.0) .OR. (maskLong(i) .EQ. 600.0)) THEN
        aBand(i,2) = 1.0
        aBand(i,1) = 0.0
        aBand(i,3) = 0.0
     END IF

  END DO

  ! make sure solver doesn't go out of bounds
  DO i=1,((yn-2)-1)
     ii = i*(xn-2)
     aBand(ii,3) = 0.0
     aBand(ii+1,1) = 0.0
  END DO

!!!!!!!!!!!! THIS !!!!!!!!!!!
  h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
  h(2:xn-1,2:yn-1) = RESHAPE(h_nextRow, (/xn-2, yn-2/))




  stretch = h0
  stretch(:,2) = stretch(:,1)
  DO ii=2,yn-1
     DO i=2,xn-1
        IF ((mask(i,ii) .EQ. 25.0) .OR. (mask(i,ii) .EQ. 50.0)  .OR. (mask(i,ii) .EQ. 17.5)  .OR. (mask(i,ii) .EQ. 12.5)  .OR. (mask(i,ii) .EQ. 3.5)  .OR. (mask(i,ii) .EQ. 6.5)) THEN
           stretch(i,ii) = stretch(i,ii+1)
        END IF

     END DO
  END DO

  stretchLongT = RESHAPE(TRANSPOSE(stretch(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))




  !h0LongT = reshape(transpose(h0(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
  !h0T = transpose(h0)
  sxLong = RESHAPE(TRANSPOSE(sxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
  syLong = RESHAPE(TRANSPOSE(syMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
  qxLong = RESHAPE(TRANSPOSE(qxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
  qyLong = RESHAPE(TRANSPOSE(qyMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

  ! horizontal boundary conditions
  h(3:xn-2,2) = h(3:xn-2,2) + h0(3:xn-2,1)*syMat(3:xn-2,2)/2.0 ! bottom


  ! top of sediment
  DO ii=2,yn-1
     DO i=2,xn-1
	! top of sediment and any short outcrops
	IF ((mask(i,ii) .EQ. 50.0) .OR. (mask(i,ii) .EQ. 2.0) .OR. (mask(i,ii) .EQ. 7.0)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	END IF

	IF ((mask(i,ii) .EQ. 25.0) .AND. (i .GE. 3) .AND. (i .LE. xn-2)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	END IF
	IF ((mask(i,ii) .EQ. 25.0) .AND. (i .EQ. 2)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top left corner
	END IF
	IF ((mask(i,ii) .EQ. 25.0) .AND. (i .EQ. xn-1)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top right corner
	END IF

	IF ((mask(i,ii) .EQ. 1.0) .AND. (ii .EQ. 2) .AND. (i .EQ. 2)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom left corner
	END IF
	IF ((mask(i,ii) .EQ. 1.0) .AND. (ii .EQ. 2) .AND. (i .EQ. xn-1)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom right corner
	END IF


	IF ((mask(i,ii) .EQ. 12.5) .OR. (mask(i,ii) .EQ. 17.5) .OR. (mask(i,ii) .EQ. 3.5) .OR. (mask(i,ii) .EQ. 6.5)) THEN
           h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	END IF


     END DO
  END DO

  !genTrans = transpose(h)
  !h_nextRow = reshape(genTrans(2:yn-1,2:xn-1), (/(xn-2)*(yn-2)/))
  h_nextRow = RESHAPE(TRANSPOSE(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

  ! make the band
  bBand = 0.0
  bBand(:,1) = -syLong/2.0 - vLong(:)*qyLong/2.0
  bBand(:,2) = 1.0+syLong
  bBand(:,3) = -syLong/2.0 + vLong*qyLong/2.0

  ! bottom left corner
  bBand(1,1) =  0.0
  bBand(1,2) = 1.0 + syLong(1)/1.0 - vLong(1)*qyLong(1)!/2.0
  bBand(1,3) = -syLong(1)/2.0 + vLong(1)*qyLong(1)!/2.0

  ! top right corner
  bBand((xn-2)*(yn-2),1) = -syLong((xn-2)*(yn-2))/2.0 - vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
  bBand((xn-2)*(yn-2),2) = 1.0 + syLong((xn-2)*(yn-2))/1.0 + vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
  bBand((xn-2)*(yn-2),3) =  0.0
  DO i = 2,(xn-2)*(yn-2)-1


     ! flow going down anywhere, 2 and 3
     IF (vLong(i) .LT. 0.0) THEN

        bBand(i,1) = -syLong(i)/2.0
        bBand(i,2) = 1.0+syLong(i) - vLong(i)*qyLong(i)
        bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)

     END IF

     ! flow coming up anywhere, 1 and 2
     IF (vLong(i) .GT. 0.0) THEN

        bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
        bBand(i,2) = 1.0+syLong(i) + vLong(i)*qyLong(i)
        bBand(i,3) = -syLong(i)/2.0

     END IF


!!!!! TOP EDGES


     ! bottom rrow, default 2 and 3 !!
     IF (MOD(i-1,yn-2) .EQ. 0) THEN
        bBand(i,1) =  0.0
        bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
        bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)
     END IF

     ! bottom rrow, if flow is coming up (1 and 2)
     IF (vLong(i) .GT. 0.0) THEN

        ! bottom row but not bottom right corner
        IF ((MOD(i-1,yn-2) .EQ. 0) .AND. (i .NE. (xn-2)*(yn-2)-(yn-2)+1)) THEN
           bBand(i,1) =  0.0
           bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
           bBand(i,3) = -syLong(i)/2.0
           h_nextRow(i) = h_nextRow(i) + vLong(i)*qyLong(i)*stretchLongT(i)
        END IF

     END IF


     ! last/top edge, default 1 and 2
     IF ((maskLongT(i).EQ.25.0) .OR. (maskLongT(i).EQ.50.0) .OR. (maskLongT(i).EQ.2.0) .OR. (maskLongT(i).EQ.7.0)) THEN
        bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
        bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
        bBand(i,3) =  0.0

     END IF

     IF ((maskLongT(i).EQ.12.5).OR.(maskLongT(i).EQ.17.5).OR.(maskLongT(i).EQ.3.5).OR.(maskLongT(i).EQ.6.5)) THEN
        bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
        bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
        bBand(i,3) =  0.0
     END IF


     ! top edges if flow coming down (2, 3)
     IF (vLong(i) .LT. 0.0) THEN

        IF ((maskLongT(i).EQ.50.0) .OR. (maskLongT(i).EQ.2.0) .OR. (maskLongT(i).EQ.7.0)) THEN
           bBand(i,1) = -syLong(i)/2.0
           bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
           bBand(i,3) =  0.0
           h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)
        END IF

        ! top but not top left or top right corner
        IF ((maskLongT(i).EQ.25.0) .AND. (i .GT. yn-2) .AND. (i .LE. (xn-2)*(yn-2)-(yn-2))) THEN
           bBand(i,1) = -syLong(i)/2.0
           bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
           bBand(i,3) =  0.0
           h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)

        END IF

        IF ((maskLongT(i).EQ.12.5).OR.(maskLongT(i).EQ.17.5).OR.(maskLongT(i).EQ.3.5).OR.(maskLongT(i).EQ.6.5)) THEN
           bBand(i,1) = -syLong(i)/2.0
           bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
           bBand(i,3) =  0.0
           h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)!/2.0
        END IF

     END IF




     ! upper left corner
     IF ((i .LE. (yn-2)) .AND. (maskLongT(i) .EQ. 25.0)) THEN
        bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
        bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
        bBand(i,3) =  0.0
     END IF

     ! upper right corner
     IF ((i .GT. (yn-2)*(xn-2)-(yn-2)) .AND. (maskLongT(i) .EQ. 25.0)) THEN
        bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
        bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
        bBand(i,3) =  0.0

     END IF



     ! bottom right corner
     IF (i .EQ. (yn-2)*(xn-2)-(yn-2)+1) THEN
        bBand(i,1) =  0.0
        bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
        bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)!/2.0
     END IF




  END DO

  DO i=1,(xn-2)*(yn-2)
     ! mask
     IF ((maskLongT(i) .EQ. 0.0) .OR. (maskLongT(i) .EQ. 600.0)) THEN
        bBand(i,2) = 1.0
        bBand(i,1) = 0.0
        bBand(i,3) = 0.0
     END IF

  END DO

  ! make sure solver doesn't go out of bounds
  DO i=1,((xn-2)-1)
     ii = i*(yn-2)
     bBand(ii,3) = 0.0
     bBand(ii+1,1) = 0.0
  END DO

  h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))
  h_next(2:xn-1,2:yn-1) = TRANSPOSE(RESHAPE(h_nextRow, (/yn-2, xn-2/)))
  hMid = h_next
  !
  !
  !

  ! ! something awful right here
  ! ! possible not-doing-upwind correction?
  ! do ii=2,yn-1
  ! do i = 2,xn-1
  ! 	if ( (hMid(i,ii).gt.h0(i-1,ii)).and.(hMid(i,ii).gt.h0(i+1,ii)).and. &
  ! 	(hMid(i,ii).gt.h0(i,ii-1)).and.(hMid(i,ii).gt.h0(i,ii+1))) then
  ! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
  ! 	end if
  !
  ! 	if ( (hMid(i,ii).lt.h0(i-1,ii)).and.(hMid(i,ii).lt.h0(i+1,ii)).and. &
  ! 	(hMid(i,ii).lt.h0(i,ii-1)).and.(hMid(i,ii).lt.h0(i,ii+1))) then
  ! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
  ! 	end if
  !
  !
  ! end do
  ! end do




  ! check out how this equation is converging to steady-state
  !write(*,*) "deltaT"
  !write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-h_next(2:xn-1,2:yn-1))/h_next(2:xn-1,2:yn-1)))

  RETURN

END FUNCTION h_next


! ----------------------------------------------------------------------------------%%
!
! PSI_NEXT
!
! SUMMARY: computes the 2D streamfunction array of the current timestep
!
! INPUTS: h(xn,yn) : temperature profile
!         rhs0(xn,yn) : right hand side of streamfunction-vorticity equation
!         psi(xn,yn) : 2D streamfunction array of previous timestep
!         top_in(xn,1) : permeable upper boundary
!         rho_in(xn,yn) : 2D density array
!
! RETURNS: psi_next(xn,yn): 2D streamfunction array for current timestep
!
! ----------------------------------------------------------------------------------%%

FUNCTION psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage, frac6_in)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! interface
  !
  ! 	function partial(array,rows,cols,d1,d2,dim)
  ! 		use globals
  ! 		use initialize
  ! 		implicit none
  ! 		integer :: rows, cols, dim, i, j, ii, jj
  ! 		real(4) :: array(rows,cols), d1, d2, d
  ! 		real(4) :: partial(rows,cols)
  ! 	end function partial
  !
  ! end interface

  ! declare errthing

  ! integers
  INTEGER :: i, j, ii, n, m, stage
  ! inputs
  REAL(4) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
  REAL(4) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
  ! matrix stuff
  REAL(4) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
  REAL(4) :: psi_next(xn,yn)
  REAL(4) :: mn(xn,yn)
  ! back to band
  REAL(4) :: aBand0(longP,2*((yn/2)-2) + 1), band_in(longP,2*((yn/2)-2) + 1)
  REAL(4) :: rhoLong(longP)
  REAL(4) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
  REAL(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
  REAL(4) :: psi_f

  psi_next = 0.0

  mn = psi

  m = 2*((yn/2)-2) + 1

  rhs1 = rhs0

  rho_in = rho_fluid

  !phi_in = 1.0

  permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  DO ii=2,yn-1
     DO i=2,xn-1
	permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
	permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)

	IF (maskP(i,ii) .EQ. 5.0) THEN
           permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
	END IF

	permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
	permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)

	IF ((maskP(i,ii) .EQ. 6.0) .OR. (maskP(i,ii) .EQ. 6.5) .OR. (maskP(i,ii) .EQ. 6.1) .OR. (maskP(i,ii) .EQ. 6.05) .OR. (maskP(i,ii) .EQ. 6.01)) THEN
           permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	END IF
	IF ((maskP(i,ii) .EQ. 3.0) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 3.1) .OR. (maskP(i,ii) .EQ. 3.05) .OR. (maskP(i,ii) .EQ. 3.01)) THEN
           permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	END IF
     END DO
  END DO

  DO ii=2,yn-1
     DO i=2,xn-1
	IF ((maskP(i,ii) .EQ. 50.0) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 6.5)) THEN
           permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
	END IF
     END DO
  END DO


  ! do ii=yn/2,yn-1
  ! do i=2,xn-1
  ! 	if ((y(ii) .ge. sed(i)-dy) .and. (any(maskP(i,:) .eq. 50.0))) then
  ! 		permy_top(i,ii) = phi_in(i,ii) / ((1e-16 + 1e-16) / 2.0)
  ! 	end if
  ! end do
  ! end do

  DO ii=yn/2,yn
     IF (maskP(1,ii) .EQ. 1.0) THEN
        ! left
        rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
     END IF

     IF (maskP(xn,ii) .EQ. 1.0) THEN
	! right
        rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
     END IF

     ! bottom left corner
     IF (maskP(1,ii) .EQ. 100.0) THEN
        ! left
        rhs1(2,ii) = rhs1(2,ii) + psi(2,ii-1)*permy_bottom(2,ii)/(dy*dy)
        rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
     END IF

     ! bottom right corner
     IF (maskP(xn,ii) .EQ. 100.0) THEN
        ! left
        rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn-1,ii-1)*permy_bottom(xn-1,ii)/(dy*dy)
        rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
     END IF


  END DO




  ! mask boundary conditions
  ! vertical boundary conditions
  DO ii=yn/2,yn-1
     DO i=2,xn-1

	IF ((maskP(i,ii) .EQ. 17.5)) THEN
           rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 10.0)) THEN
    !rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
           rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i+1,ii))*permx_left(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 6.0) .OR. (maskP(i,ii) .EQ. 6.5) .OR. (maskP(i,ii) .EQ. 6.1) .OR. (maskP(i,ii) .EQ. 6.05) .OR. (maskP(i,ii) .EQ. 6.01)) THEN
           rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,2)*permx_left(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 12.5)) THEN
           rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 5.0)) THEN
    !rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
           rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i-1,ii))*permx_right(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 3.05) .OR. (maskP(i,ii) .EQ. 3.01)) THEN
           rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
	END IF

	IF ((maskP(i,ii) .EQ. 3.0) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 3.1)) THEN
           rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
           !rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*phi_in(i,ii)*2.0*viscosity/(grav*rho_fluid*(rho_fluid*perm_in(i,ii) + (rho_fluid*param_f_dx*param_f_dx/12.0) ))
           !rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*(phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/12.0)*rho_fluid) / 2.0))
	END IF


     END DO
  END DO




  !
  ! ! at the bottom of everything
  DO ii=yn/2,yn-1
     DO i=3,xn-2
	IF ((maskP(i,ii) .EQ. 100.0) .OR. (maskP(i,ii) .EQ. 3.01) .OR. (maskP(i,ii) .EQ. 6.01)) THEN
           rhs1(i,ii) = rhs1(i,ii) + psi(i,ii-1)*permy_bottom(i,ii)/(dy*dy)
	END IF
     END DO
  END DO

  ! do ii=yn/2,yn
  ! do i=2,xn-1
  !
  ! if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 25.0) .or. (maskP(i,ii) .eq. 12.5) .or. (maskP(i,ii) .eq. 17.5)) then
  ! 	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
  ! end if
  !
  ! end do
  ! end do

  DO ii=yn/2,yn
     DO i=2,xn-1

        IF ((maskP(i,ii) .EQ. 25.0) .OR. (maskP(i,ii) .EQ. 12.5) .OR. (maskP(i,ii) .EQ. 17.5) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 6.5)) THEN
           rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
        END IF

     END DO
  END DO

  DO ii=yn/2,yn
     DO i=2,xn-1

        IF ((maskP(i,ii) .EQ. 50.0)) THEN
           rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
        END IF

     END DO
  END DO

  ! ! CORNER 2
  ! do ii=yn/2,yn
  ! do i=2,xn-1
  !
  ! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i-1,ii) .ne. 50.0)) then
  ! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
  ! end if
  !
  ! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i+1,ii) .ne. 50.0)) then
  ! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
  ! end if
  !
  ! end do
  ! end do


  DO ii=1,yn
     DO i=1,xn
	IF ((maskP(i,ii) .EQ. 0.0) .OR. (maskP(i,ii) .EQ. 600.0)) THEN
           rhs1(i,ii) = 0.0
	END IF
     END DO
  END DO


  uVec = RESHAPE(TRANSPOSE(rhs1( 2:xn-1 , (yn/2)+2:yn-1 )),(/longP/))

  psi_next = 0.0





  ! THIS IS WHERE THE BAND IS MADE
  aband0 = band_in

  ! use the banded solver here
  psi_nextRow = solve(aBand0,uVec,2*((yn/2)-2) + 1,longP)
  psi_next(2:xn-1,(yn/2)+2:yn-1) = TRANSPOSE(RESHAPE(psi_nextRow, (/(yn/2)-2, xn-2/)))




  !write(*,*) "deltaPSI"
  !write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))


  RETURN

END FUNCTION psi_next




! ----------------------------------------------------------------------------------%%
!
! MAKE BAND
!
! SUMMARY: only make the big matrix every once in a while
!
!
! ----------------------------------------------------------------------------------%%





FUNCTION make_band(perm_in,phi_in,permx,permy,rho_in)


  USE globals
  USE initialize
  IMPLICIT NONE


  INTERFACE

     FUNCTION partial(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial(rows,cols)
     END FUNCTION partial

     FUNCTION partial_edge(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge(rows,cols)
     END FUNCTION partial_edge

     FUNCTION partial_edge_p(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge_p(rows,cols)
     END FUNCTION partial_edge_p

  END INTERFACE


  INTEGER :: i, j, ii, n, m
  REAL(4) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
  REAL(4) :: permx(xn,yn), permy(xn,yn), permLong(longP)
  REAL(4) :: permxLong(longP), permyLong(longP)
  REAL(4) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
  REAL(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
  REAL(4) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
  REAL(4) :: perm_long(longP)

  !	phi_in = 1.0

  rho_in = rho_fluid

  permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
  DO ii=2,yn-1
     DO i=2,xn-1
        permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
        permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)

        IF (maskP(i,ii) .EQ. 5.0) THEN
           permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
        END IF

        permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
        permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)

        IF ((maskP(i,ii) .EQ. 6.0) .OR. (maskP(i,ii) .EQ. 6.5) .OR. (maskP(i,ii) .EQ. 6.1) .OR. (maskP(i,ii) .EQ. 6.05) .OR. (maskP(i,ii) .EQ. 6.01)) THEN
           permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
        END IF
        IF ((maskP(i,ii) .EQ. 3.0) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 3.1) .OR. (maskP(i,ii) .EQ. 3.05) .OR. (maskP(i,ii) .EQ. 3.01)) THEN
           permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
        END IF
     END DO
  END DO

  DO ii=2,yn-1
     DO i=2,xn-1
        IF ((maskP(i,ii) .EQ. 50.0) .OR. (maskP(i,ii) .EQ. 3.5) .OR. (maskP(i,ii) .EQ. 6.5)) THEN
           permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
        END IF
     END DO
  END DO


  permx_left_long = RESHAPE(TRANSPOSE(permx_left(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
  permx_right_long = RESHAPE(TRANSPOSE(permx_right(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
  permy_bottom_long = RESHAPE(TRANSPOSE(permy_bottom(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
  permy_top_long = RESHAPE(TRANSPOSE(permy_top(2:xn-1,(yn/2)+2:yn-1)),(/longP/))


  ! 	do ii=2,yn-1
  ! 	do i=2,xn-1
  ! 		if (maskP(i,ii) .eq. 50.0) then
  ! 			permy_top(:,ii) = perm_in(:,ii)! phi_in(i,ii) / ((perm_in(i,ii) + perm_in(i,ii+1)) / 2.0)
  ! 		end if
  ! 	end do
  ! 	end do

  ! make the band
  innerBand = 0.0
  make_band = 0.0
  m = 2*((yn/2)-2) + 1


  ! trial (7 days) inner band
  innerBand = 0.0
  DO i = 1,longP
     ! diagonal
     innerBand(i,(m+1)/2) = (permx_right_long(i) + permx_left_long(i))/(dx*dx) + (permy_top_long(i) + permy_bottom_long(i))/(dy*dy)


     ! off-diagonals
     ! all but left edge
     IF ((i .GT. ((yn/2)-2)) .AND. (maskPLongT(i) .NE. 10.0) .AND. (maskPLongT(i) .NE. 17.5) .AND. (maskPLongT(i) .NE. 6.0) .AND. (maskPLongT(i) .NE. 6.5) .AND. (maskPLongT(i) .NE. 6.1) .AND. (maskPLongT(i) .NE. 6.05) .AND. (maskPLongT(i) .NE. 6.01)) THEN
        innerBand(i,1) = -permx_left_long(i)/(dx*dx)
     END IF
     ! all but right edge
     IF ((i .LE. (longP)-((yn/2)-2)) .AND. (maskPLongT(i) .NE. 5.0) .AND. (maskPLongT(i) .NE. 12.5) .AND. (maskPLongT(i) .NE. 3.0) .AND. (maskPLongT(i) .NE. 3.5) .AND. (maskPLongT(i) .NE. 3.1) .AND. (maskPLongT(i) .NE. 3.05) .AND. (maskPLongT(i) .NE. 3.01)) THEN
        innerBand(i,m) = -permx_right_long(i)/(dx*dx)
     END IF


     ! all but bottom
     IF ((maskPLongT(i) .NE. 100.0) .OR. (maskPLongT(i) .NE. 3.01) .OR. (maskPLongT(i) .NE. 6.01)) THEN
        innerBand(i,(m+1)/2-1) = -permy_bottom_long(i)/(dy*dy)
     END IF

     ! all but top
     IF ((maskPLongT(i) .NE. 50.0) .AND. (maskPLongT(i) .NE. 25.0) .AND. (maskPLongT(i) .NE. 12.5) .AND. (maskPLongT(i) .NE. 17.5) .AND. (maskPLongT(i) .NE. 3.5) .AND. (maskPLongT(i) .NE. 6.5)) THEN
        innerBand(i,(m+1)/2+1) = -permy_top_long(i)/(dy*dy)
     END IF

     ! 			if (maskPLongT(i) .eq. 2.0) then
     ! 				innerBand(i,m) = 0.0 ! skip right
     ! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
     ! 			end if
     !
     ! 			if (maskPLongT(i) .eq. 7.0) then
     ! 				innerBand(i,1) = 0.0 ! skip left
     ! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
     ! 			end if

  END DO

  DO i = 1,longP

     ! mask
     IF ((maskPLongT(i) .EQ. 0.0) .OR. (maskPLongT(i) .EQ. 600.0)) THEN
        innerBand(i,:) = 0.0
        innerBand(i,(m+1)/2) = 1.0!(2.0)/(permLong(i)*dx*dx) + (2.0)/(permLong(i)*dy*dy)
     END IF

  END DO

  !write(*,*) innerBand(1,:)

  make_band = innerBand

  RETURN
END FUNCTION make_band








! ----------------------------------------------------------------------------------%%
!
! PARTICLES NEXT
!
! SUMMARY: transports particles by euler or 4th order runge kutta
!
! INPUTS:
!
!
! RETURNS:
!
! ----------------------------------------------------------------------------------%%

FUNCTION particles_next (trace_in, uTransport, vTransport, inval, num, num_sat)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing

  ! integers
  INTEGER :: i, j, ii, n, m, mm, nn, num, num_sat
  ! inputs
  REAL(4) :: trace_in(5,num), particles_next(5,num)
  REAL(4) :: uTransport(xn,yn), vTransport(xn,yn)
  REAL(4) :: u_wt, v_wt, inval
  REAL(4) :: rando
  REAL(4) :: r_i, r_ii


  DO mm=1,cstep

     ! delete particles out of bound

     ! do n = 1,num
     ! 	if ((trace_in(3,n) .ne. 0.0)) then
     ! 		i = floor(trace_in(1,n)/dx) + 1
     ! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
     ! 		if ((i .lt. 1) .or. (i .ge. xn) .or.  (ii .lt. 1) .or. (ii .ge. yn) )then
     ! 			trace_in(:,n) = 0.0
     ! 		end if
     ! 	end if
     ! end do
     !
     ! do n = 1,num
     ! 	if ((trace_in(3,n) .ne. 0.0)) then
     ! 		i = floor(trace_in(1,n)/dx) + 1
     ! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
     ! 		if (maskP(i,ii) .eq. 0.0) then
     ! 			trace_in(:,n) = 0.0
     ! 		end if
     !
     ! 		if ((maskP(i,ii+1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii+1)) .lt. dy/2.0)) then
     ! 			trace_in(:,n) = 0.0
     ! 		end if
     !
     ! 		if ((maskP(i,ii-1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii-1)) .lt. dy/2.0)) then
     ! 			trace_in(:,n) = 0.0
     ! 		end if
     !
     ! 	end if
     ! end do


     ! put in new particles in inflow cells

     DO nn = 1, xn

        IF (vTransport(twentyfives(1,nn),twentyfives(2,nn)) .LT. 0.0) THEN

           ! generate a new particle so the cell is at saturation
           DO j=1,num_sat
              m = 0
              n = 0
              ! find an empty slot to put the particle in
              DO WHILE (m .EQ. 0)
                 n = n +1
                 IF ((trace_in(3,n) .EQ. 0.0)) THEN
                    CALL RANDOM_NUMBER(rando)
                    trace_in(1,n) = dx*rando + (x(twentyfives(1,nn)) - dx/2)
                    trace_in(2,n) = dy*rando + (y(twentyfives(2,nn)) - dy/2)
                    trace_in(3,n) = inval

                    trace_in(4,n) = twentyfives(1,nn)
                    trace_in(5,n) = twentyfives(2,nn)
                    m = 1

                 END IF
              END DO
           END DO

        END IF

     END DO


     DO nn = yn/2, yn

        IF (uTransport(fives(1,nn),fives(2,nn)) .LT. 0.0) THEN

           ! generate a new particle so the cell is at saturation
           DO j=1,num_sat
              m = 0
              n = 0
              ! find an empty slot to put the particle in
              DO WHILE (m .EQ. 0)
                 n = n +1
                 IF ((trace_in(3,n) .EQ. 0.0)) THEN
                    CALL RANDOM_NUMBER(rando)
                    trace_in(1,n) = dx*rando + (x(fives(1,nn)) - dx/2)
                    trace_in(2,n) = dy*rando + (y(fives(2,nn)) - dy/2)
                    trace_in(3,n) = inval

                    trace_in(4,n) = fives(1,nn)
                    trace_in(5,n) = fives(2,nn)
                    m = 1
                 END IF
              END DO
           END DO

        END IF


        IF (uTransport(tens(1,nn),tens(2,nn)) .GT. 0.0) THEN

           ! generate a new particle so the cell is at saturation
           DO j=1,num_sat
              m = 0
              n = 0
              ! find an empty slot to put the particle in
              DO WHILE (m .EQ. 0)
                 n = n +1
                 IF ((trace_in(3,n) .EQ. 0.0)) THEN
                    CALL RANDOM_NUMBER(rando)
                    trace_in(1,n) = dx*rando + (x(tens(1,nn)) - dx/2)
                    trace_in(2,n) = dy*rando + (y(tens(2,nn)) - dy/2)
                    trace_in(3,n) = inval

                    trace_in(4,n) = tens(1,nn)
                    trace_in(5,n) = tens(2,nn)
                    m = 1
                 END IF
              END DO
           END DO

        END IF

     END DO


     ! do nn = 1, xn
     !
     ! 		if (vTransport(fifties(1,nn),fifties(2,nn)) .lt. 0.0) then
     !
     ! 		! generate a new particle so the cell is at saturation
     ! 		do j=1,num_sat
     ! 			m = 0
     ! 			n = 0
     ! 			! find an empty slot to put the particle in
     ! 			do while (m .eq. 0)
     ! 				n = n +1
     ! 				if ((trace_in(3,n) .eq. 0.0)) then
     ! 					call random_number(rando)
     ! 					trace_in(1,n) = dx*rando + (x(fifties(1,nn)) - dx/2)
     ! 					trace_in(2,n) = dy*rando + (y(fifties(2,nn)) - dy/2)
     ! 					trace_in(3,n) = inval
     !
     ! 					trace_in(4,n) = fifties(1,nn)
     ! 					trace_in(5,n) = fifties(2,nn)
     ! 					m = 1
     ! 				end if
     ! 			end do
     ! 		end do
     !
     ! 		end if
     !
     ! end do



     ! advect all of the particles

     DO n = 1,num

	IF ((trace_in(3,n) .NE. 0.0)) THEN
    ! 		i = floor(trace_in(1,n)/dx) + 1
    ! 		i = min(i,xn-1)
    ! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
    ! 		ii = min(ii,yn-1)
           i = trace_in(4,n)
           ii = trace_in(5,n)


           u_wt = ( uTransport(i,ii) * (y(ii+1) - trace_in(2,n)) * (x(i+1) - trace_in(1,n)) + &
                uTransport(i+1,ii) * (y(ii+1) - trace_in(2,n)) * (trace_in(1,n) - x(i)) + &
                uTransport(i+1,ii+1) * (trace_in(2,n) - y(ii)) * (trace_in(1,n) - x(i)) + &
                uTransport(i,ii+1) * (trace_in(2,n) - y(ii)) * (x(i+1) - trace_in(1,n)) ) / (dx*dy)

           v_wt = ( vTransport(i,ii) * (y(ii+1) - trace_in(2,n)) * (x(i+1) - trace_in(1,n)) + &
                vTransport(i+1,ii) * (y(ii+1) - trace_in(2,n)) * (trace_in(1,n) - x(i)) + &
                vTransport(i+1,ii+1) * (trace_in(2,n) - y(ii)) * (trace_in(1,n) - x(i)) + &
                vTransport(i,ii+1) * (trace_in(2,n) - y(ii)) * (x(i+1) - trace_in(1,n)) ) / (dx*dy)

           ! 					u_wt = uTransport(i,ii)
           ! 					v_wt = vTransport(i,ii)


           ! advect the thing
           trace_in(1,n) = trace_in(1,n) + u_wt*dt*mstep/(cstep)!/1000.0
           trace_in(2,n) = trace_in(2,n) + v_wt*dt*mstep/(cstep)!/1000.0

           IF ((trace_in(1,n).GT.x_max) .OR. (trace_in(1,n).LT.x_min) .OR. (trace_in(2,n).GT.y_max) .OR. (trace_in(2,n).LT.y_min)) THEN
              trace_in(:,n) = 0.0
           END IF

           ! 		i = 1.0*floor(trace_in(1,n)/dx) + 1
           ! 		i = 1.0*min(i,xn-1)
           ! 		trace_in(4,n) = i
           !
           ! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
           ! 		ii = 1.0*min(ii,yn-1)
           ! 		trace_in(5,n) = ii

           trace_in(4,n) = MIN(FLOOR(trace_in(1,n)/dx) + 1,xn-1)
           trace_in(5,n) = MIN(yn - FLOOR(-1.0*trace_in(2,n)/dy) - 1,yn-1)
           i = trace_in(4,n)
           ii = trace_in(5,n)

           IF (maskP(i,ii) .EQ. 0.0) THEN
              trace_in(:,n) = 0.0
           END IF

           IF ((maskP(i,ii+1) .EQ. 0.0) .AND. (ABS(trace_in(2,n) - y(ii+1)) .LT. dy/2.0)) THEN
              trace_in(:,n) = 0.0
           END IF

           IF ((maskP(i,ii-1) .EQ. 0.0) .AND. (ABS(trace_in(2,n) - y(ii-1)) .LT. dy/2.0)) THEN
              trace_in(:,n) = 0.0
           END IF
           !
           ! 		if (permeability(i,ii) .le. 1e-13) then
           ! 			trace_in(:,n) = 0.0
           ! 		end if
           !
	END IF
     END DO


  END DO

  particles_next = trace_in

  RETURN


END FUNCTION particles_next


! ----------------------------------------------------------------------------------%%
!
! SOLUTE_NEXT
!
! SUMMARY: transports solutes on coarse mesh
!
! INPUTS: sol(xn/cell,yn/cell) : 2D array of initial solute concentrations
!         uTransport(xn/cell,yn/cell) : lateral velocities (coarse mesh)
!         vTransport(xn/cell,yn/cell) : vertical velocities (coarse mesh)
!
! RETURNS: solute_next(xn/cell,yn/cell): 2D array of solute concentrations
!
! ----------------------------------------------------------------------------------%%

FUNCTION solute_next (sol, uTransport, vTransport, seaw)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing

  ! integers
  INTEGER :: i, j, ii, n, m
  ! inputs
  REAL(4) :: sol(xn,yn), sol0(xn,yn)
  REAL(4) :: uTransport(xn,yn), vTransport(xn,yn)
  ! solver stuff
  REAL(4) :: uLong(((xn)-2)*((yn)-0)), vLong(((xn)-0)*((yn)-2))
  REAL(4) :: aBand(((xn)-2)*((yn)-0),5), bBand(((xn)-0)*((yn)-2),5)
  REAL(4) :: qx, qy, solute_next(xn,yn), vec(((xn)-2)*((yn)-0))
  REAL(4) :: sol_nextRow(((xn)-2)*((yn)-0)), sol_nextRowB(((xn)-0)*((yn)-2))
  REAL(4) :: seaw
  REAL(4) :: bm1(xn,yn), b0(xn,yn), bp1(xn,yn), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
  REAL(4) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6

  ! call init_mini()


  ! do i = 1,xn
  ! 	do j = 1,yn
  ! 		if ((maskP(i,j) .eq. 0.0)) then
  ! 			sol(i,j) = seaw
  ! 		end if
  ! 	end do
  ! end do

  sol(1,:) = seaw!(4.0/3.0)*sol(2,:) - (1.0/3.0)*sol(3,:)
  sol(xn,:) = (4.0/3.0)*sol(xn-1,:) - (1.0/3.0)*sol(xn-2,:)




  qx = dt*mstep/(cstep*dx)
  qy = dt*mstep/(cstep*dy)

  !write(*,*) "qx, qy" , qx , qy

  ! uLong = reshape(uTransport(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
  ! !! transpose coarse needed!
  ! vLong = reshape(transpose(vTransport(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))

  uTransport(1,:) = 0.0
  !uTransport(xn,:) = 0.0

  ! write(*,*) qx*maxval(abs(uTransport))
  ! write(*,*) qy*maxval(abs(vTransport))



  DO i = 1,xn

     ! 	if (i .eq. f_index1-1) then
     ! 		sol(i,:) = (4.0/3.0)*sol(i-1,:) - (1.0/3.0)*sol(i-2,:)
     ! 		do j = yn/2,yn
     ! 			sol(i+1:,j) = sol(i,j)
     ! 		end do
     !
     ! 	end if

     DO j = yn/2,yn

        ! 		if ((maskP(i,j) .eq. 0.0)) then
        ! 			sol(i,j) = seaw
        ! 		end if

        IF ((maskP(i,j) .EQ. 5.0)) THEN

           IF (uTransport(i,j) .LT. 0.0) THEN
              sol(i+1,j) = seaw
              uTransport(i+1,j) = uTransport(i,j)
           ELSE
              sol(i+1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i-1,j)
           END IF

        END IF


        ! 		if ((maskP(i,j) .eq. 10.0) .or. (maskP(i,j) .eq. 170.5)) then
        ! 			if (uTransport(i,j) .gt. 0.0) then
        ! 				sol(i-1,j) = seaw
        ! 				uTransport(i-1,j) = uTransport(i,j)
        ! 			else
        ! 				sol(i-1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i+1,j)
        ! 			end if
        ! 		end if
     END DO
  END DO

  sol0 = sol
  solute_next = sol

  DO j = yn/2,yn-1
     ! do i = 2,xn-1
     solute_next(2,j) = sol0(2,j) - qx*uTransport(2,j)*( sol0(2,j) - sol0(1,j) )
     DO i = 3,xn-2
        IF (uTransport(i,j) .GT. 1e-10) THEN
           !do i = 3,f_index1-2



           ! if (uTransport(i,j) .gt. 0.0) then
           !if (uTransport(i,j) .gt. 0.0) then
           ! upwind including LHS value
           solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i,j) - sol0(i-1,j) )


           ! correction loop: sort of a mess
           !if (i .gt. 2) then
           !if (maskP(i-2,j) .ne. 0.0) then
           !sigma1 = 0.0
           sigma1a = (sol0(i+1,j) - sol0(i,j))/dx
           sigma1b = 2.0*(sol0(i,j) - sol0(i-1,j))/dx

           !if (sigma1a*sigma1b .gt. 0.0) then
           !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
           ! 							if (sigma1 .eq. abs(sigma1a)) then
           ! 								sigma1 = sign(1.0,sigma1a)*sigma1
           ! 							end if
           ! 							if (sigma1 .eq. abs(sigma1b)) then
           ! 								sigma1 = sign(1.0,sigma1b)*sigma1
           ! 							end if
           ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
           ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
           ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
           sigma1 = ((MINLOC((/ABS(sigma1a), ABS(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((MINLOC((/ABS(sigma1b), ABS(sigma1a)/),DIM=1)-1.0)*sigma1a)
           !end if

           !sigma3 = 0.0
           sigma3a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           sigma3b = (sol0(i,j) - sol0(i-1,j))/dx

           !if (sigma3a*sigma3b .gt. 0.0) then
           !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
           ! 						if (sigma3 .eq. abs(sigma3a)) then
           ! 							sigma3 = sign(1.0,sigma3a)*sigma3
           ! 						end if
           ! 						if (sigma3 .eq. abs(sigma3b)) then
           ! 							sigma3 = sign(1.0,sigma3b)*sigma3
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
           ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
           ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
           sigma3 = ((MINLOC((/ABS(sigma3a), ABS(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((MINLOC((/ABS(sigma3b), ABS(sigma3a)/),DIM=1)-1.0)*sigma3a)
           !end if

           ! choosing sigma5
           sigma5 = 0.0
           IF (sigma1*sigma3 .GT. 0.0) THEN
              sigma5 = SIGN(1.0,sigma1)*MAXVAL((/ABS(sigma1), ABS(sigma3)/))
           END IF

           !sigma2 = 0.0
           sigma2a = (sol0(i,j) - sol0(i-1,j))/dx
           sigma2b = 2.0*(sol0(i-1,j) - sol0(i-2,j))/dx

           !if (sigma2a*sigma2b .gt. 0.0) then
           !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
           ! 						if (sigma2 .eq. abs(sigma2a)) then
           ! 							sigma2 = sign(1.0,sigma2a)*sigma2
           ! 						end if
           ! 						if (sigma2 .eq. abs(sigma2b)) then
           ! 							sigma2 = sign(1.0,sigma2b)*sigma2
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
           ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
           ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
           sigma2 = ((MINLOC((/ABS(sigma2a), ABS(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((MINLOC((/ABS(sigma2b), ABS(sigma2a)/),DIM=1)-1.0)*sigma2a)
           !end if

           !sigma4 = 0.0
           sigma4a = 2.0*(sol0(i,j) - sol0(i-1,j))/dx
           sigma4b = (sol0(i-1,j) - sol0(i-2,j))/dx

           !if (sigma4a*sigma4b .gt. 0.0) then
           !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
           ! 						if (sigma4 .eq. abs(sigma4a)) then
           ! 							sigma4 = sign(1.0,sigma4a)*sigma4
           ! 						end if
           ! 						if (sigma4 .eq. abs(sigma4b)) then
           ! 							sigma4 = sign(1.0,sigma4b)*sigma4
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
           ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
           ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
           sigma4 = ((MINLOC((/ABS(sigma4a), ABS(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((MINLOC((/ABS(sigma4b), ABS(sigma4a)/),DIM=1)-1.0)*sigma4a)
           !end if


           ! choosing sigma6
           sigma6 = 0.0
           IF (sigma2*sigma4 .GT. 0.0) THEN
              sigma6 = SIGN(1.0,sigma2)*MAXVAL((/ABS(sigma2), ABS(sigma4)/))
           END IF

           ! 						write(*,*) "sigma5"
           ! 						write(*,*) sigma5
           ! 						write(*,*) "sigma6"
           ! 						write(*,*) sigma6
           correction = (uTransport(i,j)*qx*0.5) * (sigma5 - sigma6) * (dx - uTransport(i,j)*qx*dx)
           solute_next(i,j) = solute_next(i,j) - correction

           !end if ! end if maskP i-2,j .eq. 0
           !end if ! end if i .gt. 2
           ! end correction loop

           !end if ! end if u .gt. 0.0

           ! 			if (uTransport(i,j) .lt. 0.0) then
           ! 				! upwind including RHS value
           ! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i+1,j) - sol0(i,j) )
           !
           !
           ! 				! correction loop: sort of a mess
           ! 				!if (i .lt. xn-1) then
           ! 					if (maskP(i+2,j) .ne. 0.0) then
           ! 						sigma1 = 0.0
           ! 						sigma1a = (sol0(i+2,j) - sol0(i+1,j))/dx
           ! 						sigma1b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma1a*sigma1b .gt. 0.0) then
           ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
           ! ! 							if (sigma1 .eq. abs(sigma1a)) then
           ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
           ! ! 							end if
           ! ! 							if (sigma1 .eq. abs(sigma1b)) then
           ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
           ! ! 							end if
           ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
           ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
           ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
           ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
           ! 							!end if
           !
           ! 						sigma3 = 0.0
           ! 						sigma3a = 2.0*(sol0(i+2,j) - sol0(i+1,j))/dx
           ! 						sigma3b = (sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma3a*sigma3b .gt. 0.0) then
           ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
           ! ! 						if (sigma3 .eq. abs(sigma3a)) then
           ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
           ! ! 						end if
           ! ! 						if (sigma3 .eq. abs(sigma3b)) then
           ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
           ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
           ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
           ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
           ! 							!end if
           !
           !
           ! 						! choosing sigma5
           ! 						sigma5 = 0.0
           ! 						if (sigma1*sigma3 .gt. 0.0) then ! not sure about these conditionals...?
           ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
           ! 						end if
           !
           !
           !
           !
           !
           ! 						sigma2 = 0.0
           ! 						sigma2a = (sol0(i+1,j) - sol0(i,j))/dx
           ! 						sigma2b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma2a*sigma2b .gt. 0.0) then
           ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
           ! ! 						if (sigma2 .eq. abs(sigma2a)) then
           ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
           ! ! 						end if
           ! ! 						if (sigma2 .eq. abs(sigma2b)) then
           ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
           ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
           ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
           ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
           ! 							!end if
           !
           ! 						sigma4 = 0.0
           ! 						sigma4a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           ! 						sigma4b = (sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma4a*sigma4b .gt. 0.0) then
           ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
           ! ! 						if (sigma4 .eq. abs(sigma4a)) then
           ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
           ! ! 						end if
           ! ! 						if (sigma4 .eq. abs(sigma4b)) then
           ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
           ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
           ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
           ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
           ! 							!end if
           !
           !
           ! 						! choosing sigma6
           ! 						sigma6 = 0.0
           ! 						if (sigma2*sigma4 .gt. 0.0) then
           ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
           ! 						end if
           !
           ! ! 						write(*,*) "sigma5"
           ! ! 						write(*,*) sigma5
           ! ! 						write(*,*) "sigma6"
           ! ! 						write(*,*) sigma6
           ! 						correction = (uTransport(i,j)*qx*0.5) * (sigma6 - sigma5) * (dx - uTransport(i,j)*qx*dx)
           ! 						solute_next(i,j) = solute_next(i,j) - correction
           ! 					end if
           ! 				!end if ! end if i .lt. xn-1
           ! 				! end correction loop
           !
           ! 			end if
           !
        END IF ! end mod thing

        ! if ((solute_next(i,j) .gt. sol0(i,j)*1.1) .and. (solute_next(i,j) .gt. sol0(i-1,j)*1.1)) then
        !     solute_next(i,j) = sol0(i,j)
        ! end if

     END DO
     !solute_next(xn-1,j) = sol0(xn-1,j) - qx*uTransport(xn-1,j)*( sol0(xn-1,j) - sol0(xn-2,j) )
  END DO


  ! sol = solute_next
  !
  !
  !
  ! do j = yn/2,yn
  ! 	do i = 1,xn
  ! 		if ((maskP(i,j) .eq. 25.0)) then
  !
  ! 			if (vTransport(i,j+1) .lt. 0.0) then
  ! 				sol(i,j+1) = seaw
  ! 			else
  ! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
  ! 			end if
  !
  ! 		end if
  !
  ! 		if ((maskP(i,j) .eq. 50.0)) then
  !
  ! ! 			if (vTransport(i,j+1) .lt. 0.0) then
  ! 				sol(i,j+1) = seaw
  ! ! 			else
  ! ! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
  ! ! 			end if
  !
  ! 		end if
  !
  !
  !
  ! 		if (maskP(i,j) .eq. 100.0) then
  !
  ! 				sol(i,j-1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j+1)
  !
  ! 		end if
  !
  !
  ! 	end do
  ! end do
  !
  ! sol0 = sol
  !
  !
  !
  ! do j = yn/2,yn-1
  ! 	! do i = 1,xn
  ! 	do i = 1,f_index1-1
  !
  ! 			if (vTransport(i,j) .gt. 1.0e-9) then
  ! 				! upwind including bottom value
  ! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j) - sol0(i,j-1) )
  ! !
  ! 			! correction loop: sort of a mess
  ! 				!if (maskP(i,j-2) .ne. 0.0) then
  ! 						!sigma1 = 0.0
  ! 						sigma1a = (sol0(i,j+1) - sol0(i,j))/dy
  ! 						sigma1b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
  !
  ! 						!if (sigma1a*sigma1b .gt. 0.0) then
  ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
  ! ! 							if (sigma1 .eq. abs(sigma1a)) then
  ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
  ! ! 							end if
  ! ! 							if (sigma1 .eq. abs(sigma1b)) then
  ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
  ! ! 							end if
  ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
  ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
  ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
  ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
  ! 							!end if
  !
  ! 						!sigma3 = 0.0
  ! 						sigma3a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
  ! 						sigma3b = (sol0(i,j) - sol0(i,j-1))/dy
  !
  ! 						!if (sigma3a*sigma3b .gt. 0.0) then
  ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
  ! ! 						if (sigma3 .eq. abs(sigma3a)) then
  ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
  ! ! 						end if
  ! ! 						if (sigma3 .eq. abs(sigma3b)) then
  ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
  ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
  ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
  ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma5
  ! 						sigma5 = 0.0
  ! 						if (sigma1*sigma3 .gt. 0.0) then
  ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
  ! 						end if
  !
  !
  !
  !
  !
  ! 						!sigma2 = 0.0
  ! 						sigma2a = (sol0(i,j) - sol0(i,j-1))/dy
  ! 						sigma2b = 2.0*(sol0(i,j-1) - sol0(i,j-2))/dy
  !
  ! 						!if (sigma2a*sigma2b .gt. 0.0) then
  ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
  ! ! 						if (sigma2 .eq. abs(sigma2a)) then
  ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
  ! ! 						end if
  ! ! 						if (sigma2 .eq. abs(sigma2b)) then
  ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
  ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
  ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
  ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
  ! 							!end if
  !
  ! 						!sigma4 = 0.0
  ! 						sigma4a = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
  ! 						sigma4b = (sol0(i,j-1) - sol0(i,j-2))/dy
  !
  ! 						!if (sigma4a*sigma4b .gt. 0.0) then
  ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
  ! ! 						if (sigma4 .eq. abs(sigma4a)) then
  ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
  ! ! 						end if
  ! ! 						if (sigma4 .eq. abs(sigma4b)) then
  ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
  ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
  ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
  ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma6
  ! 						sigma6 = 0.0
  ! 						if (sigma2*sigma4 .gt. 0.0) then
  ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
  ! 						end if
  !
  ! ! 						write(*,*) "sigma5"
  ! ! 						write(*,*) sigma5
  ! ! 						write(*,*) "sigma6"
  ! ! 						write(*,*) sigma6
  ! 						correction = (vTransport(i,j)*qy*0.5) * (sigma5 - sigma6) * (dy - vTransport(i,j)*qy*dy)
  ! 						solute_next(i,j) = solute_next(i,j) - correction
  ! 				!end if ! end if maskP i,j-2 .ne. 0.0
  ! 				! end correction loop
  !
  !
  ! 			end if
  !
  ! 			if (vTransport(i,j) .lt. -1.0e-9) then
  ! 				! upwind including top value
  ! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j+1) - sol0(i,j) )
  !
  ! 				! correction loop: sort of a mess
  ! 				!if (j .lt. yn-1) then
  ! 					if (maskP(i,j+2) .ne. 0.0) then
  ! 						!sigma1 = 0.0
  ! 						sigma1a = (sol0(i,j+2) - sol0(i,j+1))/dy
  ! 						sigma1b = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
  !
  ! 						!if (sigma1a*sigma1b .gt. 0.0) then
  ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
  ! ! 							if (sigma1 .eq. abs(sigma1a)) then
  ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
  ! ! 							end if
  ! ! 							if (sigma1 .eq. abs(sigma1b)) then
  ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
  ! ! 							end if
  ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
  ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
  ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
  ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
  ! 							!end if
  !
  ! 						!sigma3 = 0.0
  ! 						sigma3a = 2.0*(sol0(i,j+2) - sol0(i,j+1))/dy
  ! 						sigma3b = (sol0(i,j+1) - sol0(i,j))/dy
  !
  ! 						!if (sigma3a*sigma3b .gt. 0.0) then
  ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
  ! ! 						if (sigma3 .eq. abs(sigma3a)) then
  ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
  ! ! 						end if
  ! ! 						if (sigma3 .eq. abs(sigma3b)) then
  ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
  ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
  ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
  ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma5
  ! 						sigma5 = 0.0
  ! 						if (sigma1*sigma3 .gt. 0.0) then
  ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
  ! 						end if
  !
  !
  !
  !
  !
  ! 						!sigma2 = 0.0
  ! 						sigma2a = (sol0(i,j+1) - sol0(i,j))/dy
  ! 						sigma2b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
  !
  ! 						!if (sigma2a*sigma2b .gt. 0.0) then
  ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
  ! ! 						if (sigma2 .eq. abs(sigma2a)) then
  ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
  ! ! 						end if
  ! ! 						if (sigma2 .eq. abs(sigma2b)) then
  ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
  ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
  ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
  ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
  ! 							!end if
  !
  ! 						!sigma4 = 0.0
  ! 						sigma4a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
  ! 						sigma4b = (sol0(i,j) - sol0(i,j-1))/dy
  !
  ! 						!if (sigma4a*sigma4b .gt. 0.0) then
  ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
  ! ! 						if (sigma4 .eq. abs(sigma4a)) then
  ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
  ! ! 						end if
  ! ! 						if (sigma4 .eq. abs(sigma4b)) then
  ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
  ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
  ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
  ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma6
  ! 						sigma6 = 0.0
  ! 						if (sigma2*sigma4 .gt. 0.0) then
  ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
  ! 						end if
  !
  ! ! 						write(*,*) "sigma5"
  ! ! 						write(*,*) sigma5
  ! ! 						write(*,*) "sigma6"
  ! ! 						write(*,*) sigma6
  ! 						correction = (vTransport(i,j)*qy*0.5) * (sigma6 - sigma5) * (dy - vTransport(i,j)*qy*dy)
  ! 						solute_next(i,j) = solute_next(i,j) - correction
  !
  ! 					end if
  ! 				!end if ! end if j .lt. yn-1
  ! 				! end correction loop
  !
  !
  ! 			end if
  !
  !
  ! 	end do
  ! 	solute_next(i,yn-1) = sol0(i,yn-1) - qy*vTransport(i,yn-1)*( sol0(i,yn) - sol0(i,yn-1) )
  ! end do
  !



  RETURN

END FUNCTION solute_next












! ----------------------------------------------------------------------------------%%
!
! SOLUTE_NEXT_COARSE
!
! SUMMARY: transports solutes on coarse mesh
!
! INPUTS: sol(xn/cell,yn/cell) : 2D array of initial solute concentrations
!         uTransport(xn/cell,yn/cell) : lateral velocities (coarse mesh)
!         vTransport(xn/cell,yn/cell) : vertical velocities (coarse mesh)
!
! RETURNS: solute_next(xn/cell,yn/cell): 2D array of solute concentrations
!
! ----------------------------------------------------------------------------------%%

FUNCTION solute_next_coarse (sol, uTransport, vTransport, phiTransport, seaw)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing

  ! integers
  INTEGER :: i, j, ii, n, m
  ! inputs
  REAL(4) :: sol((xn-1)/cellx,yn/(2*celly)), sol0((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: uTransport((xn-1)/cellx,yn/(2*celly)), vTransport((xn-1)/cellx,yn/(2*celly))
  real(4) :: phiTransport((xn-1)/cellx,yn/(2*celly))
  ! solver stuff
  REAL(4) :: uLong((((xn-1)/cellx)-2)*((yn/(2*celly))-0)), vLong((((xn-1)/cellx)-0)*((yn/(2*celly))-2))
  REAL(4) :: aBand((((xn-1)/cellx)-2)*((yn/(2*celly))-0),5), bBand((((xn-1)/cellx)-0)*((yn/(2*celly))-2),5)
  REAL(4) :: qx, qy, solute_next_coarse((xn-1)/cellx,yn/(2*celly)), vec((((xn-1)/cellx)-2)*((yn/(2*celly))-0))
  REAL(4) :: sol_nextRow((((xn-1)/cellx)-2)*((yn/(2*celly))-0)), sol_nextRowB((((xn-1)/cellx)-0)*((yn/(2*celly))-2))
  REAL(4) :: seaw
  REAL(4) :: bm1((xn-1)/cellx,yn/(2*celly)), b0((xn-1)/cellx,yn/(2*celly)), bp1((xn-1)/cellx,yn/(2*celly)), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
  REAL(4) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6

  ! call init_mini()



  sol(1,:) = seaw!(4.0/3.0)*sol(2,:) - (1.0/3.0)*sol(3,:)
  !sol((xn-1)/cellx,:) = sol((xn-1)/cellx-1,:)
  !sol((xn-1)/cellx,:) = (4.0/3.0)*sol((xn-1)/cellx-1,:) - (1.0/3.0)*sol((xn-1)/cellx-2,:)




  qx = dt*mstep/(cstep*dx*cellx)
  qy = dt*mstep/(cstep*dy*celly)


  uTransport(1,:) = 0.0



  !sol((xn-1)/cellx,:) = (4.0/3.0)*sol((xn-1)/cellx-1,:) - (1.0/3.0)*sol((xn-1)/cellx-2,:)

  sol0 = sol
  solute_next_coarse = sol

  DO j = 1,yn/(2*celly)
     ! do i = 2,xn-1
      !solute_next_coarse(2,j) = sol0(2,j)-(qx*uTransport(2,j)/phiTransport(2,j))*(sol0(2,j)-sol0(1,j))! - qx*uTransport(2,j)*sol0(2,j)*((1.0/phiTransport(2,j))-(1.0/phiTransport(1,j)))

    ! right before i looked at more equations (never mind no real changes?)
    ! solute_next_coarse(2,j) = sol0(2,j) - (qx*uTransport(2,j))*(sol0(2,j)-sol0(1,j)) - qx*(uTransport(2,j)/phiTransport(2,j))*sol0(2,j)*(phiTransport(2,j)-phiTransport(1,j))
    !
    ! solute_next_coarse((xn-1)/cellx,j) = sol0((xn-1)/cellx,j) - (qx*uTransport((xn-1)/cellx,j))*(sol0((xn-1)/cellx,j)-sol0((xn-1)/cellx-1,j)) - qx*(uTransport((xn-1)/cellx,j)/phiTransport((xn-1)/cellx,j))*sol0((xn-1)/cellx,j)*(phiTransport((xn-1)/cellx,j)-phiTransport((xn-1)/cellx-1,j))
    !
    ! solute_next_coarse((xn-1)/cellx-1,j) = sol0((xn-1)/cellx-1,j) - (qx*uTransport((xn-1)/cellx-1,j))*(sol0((xn-1)/cellx-1,j)-sol0((xn-1)/cellx-2,j)) - qx*(uTransport((xn-1)/cellx-1,j)/phiTransport((xn-1)/cellx-1,j))*sol0((xn-1)/cellx-1,j)*(phiTransport((xn-1)/cellx-1,j)-phiTransport((xn-1)/cellx-2,j))

    !solute_next_coarse(2,j) = sol0(2,j) - (qx*0.12866E-06)*(sol0(2,j)-sol0(1,j)) - qx*(0.12866E-06/phiTransport(2,j))*sol0(2,j)*(phiTransport(2,j)-phiTransport(1,j))















    solute_next_coarse(2,j) = sol0(2,j) - (qx*uTransport(2,j))*(sol0(2,j)-sol0(1,j)) - qx*(uTransport(2,j)/phiTransport(2,j))*sol0(2,j)*(phiTransport(2,j)-phiTransport(1,j))

     solute_next_coarse((xn-1)/cellx,j) = sol0((xn-1)/cellx,j) - (qx*uTransport((xn-1)/cellx,j))*(sol0((xn-1)/cellx,j)-sol0((xn-1)/cellx-1,j)) - qx*(uTransport((xn-1)/cellx,j)/phiTransport((xn-1)/cellx,j))*sol0((xn-1)/cellx,j)*(phiTransport((xn-1)/cellx,j)-phiTransport((xn-1)/cellx-1,j))

    !solute_next_coarse((xn-1)/cellx,j) = (4.0/3.0)*sol0((xn-1)/cellx-1,j) - (1.0/3.0)*sol0((xn-1)/cellx-2,j)

    ! solute_next_coarse((xn-1)/cellx-1,j) = sol0((xn-1)/cellx-1,j) - (qx*uTransport((xn-1)/cellx-1,j))*(sol0((xn-1)/cellx-1,j)-sol0((xn-1)/cellx-2,j)) - qx*(uTransport((xn-1)/cellx-1,j)/phiTransport((xn-1)/cellx-1,j))*sol0((xn-1)/cellx-1,j)*(phiTransport((xn-1)/cellx-1,j)-phiTransport((xn-1)/cellx-2,j))

     DO i = 3,(xn-1)/cellx-1
        IF (uTransport(i,j) .GT. 1e-9) THEN
           !do i = 3,f_index1-2



           ! if (uTransport(i,j) .gt. 0.0) then
           !if (uTransport(i,j) .gt. 0.0) then
           ! upwind including LHS value
            !solute_next_coarse(i,j) = sol0(i,j)-(qx*uTransport(i,j)/phiTransport(i,j))*(sol0(i,j)-sol0(i-1,j))! - qx*uTransport(i,j)*sol0(i,j)*((1.0/phiTransport(i,j))-(1.0/phiTransport(i-1,j)))

        ! new advection scheme
        solute_next_coarse(i,j) = sol0(i,j) - (qx*uTransport(i,j))*(sol0(i,j)-sol0(i-1,j)) - qx*(uTransport(i,j)/phiTransport(i,j))*sol0(i,j)*(phiTransport(i,j)-phiTransport(i-1,j))

        ! solute_next_coarse(i,j) = sol0(i,j) - (qx*uTransport(i,j))*(sol0(i,j)-sol0(i-1,j)) - qx*(uTransport(i,j)/((phiTransport(i,j)+phiTransport(i-1,j))/2.0))*sol0(i,j)*(phiTransport(i,j)-phiTransport(i-1,j))

        !solute_next_coarse(i,j) = sol0(i,j) - (qx*0.12866E-06)*(sol0(i,j)-sol0(i-1,j)) - qx*(0.12866E-06/phiTransport(i,j))*sol0(i,j)*(phiTransport(i,j)-phiTransport(i-1,j))

           !
           ! correction loop: sort of a mess
           !if (i .gt. 2) then
           !if (maskP(i-2,j) .ne. 0.0) then
           !sigma1 = 0.0
           sigma1a = (sol0(i+1,j) - sol0(i,j))/(dx*cellx)
           sigma1b = 2.0*(sol0(i,j) - sol0(i-1,j))/(dx*cellx)

           !if (sigma1a*sigma1b .gt. 0.0) then
           !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
           ! 							if (sigma1 .eq. abs(sigma1a)) then
           ! 								sigma1 = sign(1.0,sigma1a)*sigma1
           ! 							end if
           ! 							if (sigma1 .eq. abs(sigma1b)) then
           ! 								sigma1 = sign(1.0,sigma1b)*sigma1
           ! 							end if
           ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
           ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
           ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
           sigma1 = ((MINLOC((/ABS(sigma1a), ABS(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((MINLOC((/ABS(sigma1b), ABS(sigma1a)/),DIM=1)-1.0)*sigma1a)
           !end if

           !sigma3 = 0.0
           sigma3a = 2.0*(sol0(i+1,j) - sol0(i,j))/(dx*cellx)
           sigma3b = (sol0(i,j) - sol0(i-1,j))/(dx*cellx)

           !if (sigma3a*sigma3b .gt. 0.0) then
           !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
           ! 						if (sigma3 .eq. abs(sigma3a)) then
           ! 							sigma3 = sign(1.0,sigma3a)*sigma3
           ! 						end if
           ! 						if (sigma3 .eq. abs(sigma3b)) then
           ! 							sigma3 = sign(1.0,sigma3b)*sigma3
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
           ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
           ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
           sigma3 = ((MINLOC((/ABS(sigma3a), ABS(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((MINLOC((/ABS(sigma3b), ABS(sigma3a)/),DIM=1)-1.0)*sigma3a)
           !end if

           ! choosing sigma5
           sigma5 = 0.0
           IF (sigma1*sigma3 .GT. 0.0) THEN
              sigma5 = SIGN(1.0,sigma1)*MAXVAL((/ABS(sigma1), ABS(sigma3)/))
           END IF

           !sigma2 = 0.0
           sigma2a = (sol0(i,j) - sol0(i-1,j))/(dx*cellx)
           sigma2b = 2.0*(sol0(i-1,j) - sol0(i-2,j))/(dx*cellx)

           !if (sigma2a*sigma2b .gt. 0.0) then
           !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
           ! 						if (sigma2 .eq. abs(sigma2a)) then
           ! 							sigma2 = sign(1.0,sigma2a)*sigma2
           ! 						end if
           ! 						if (sigma2 .eq. abs(sigma2b)) then
           ! 							sigma2 = sign(1.0,sigma2b)*sigma2
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
           ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
           ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
           sigma2 = ((MINLOC((/ABS(sigma2a), ABS(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((MINLOC((/ABS(sigma2b), ABS(sigma2a)/),DIM=1)-1.0)*sigma2a)
           !end if

           !sigma4 = 0.0
           sigma4a = 2.0*(sol0(i,j) - sol0(i-1,j))/(dx*cellx)
           sigma4b = (sol0(i-1,j) - sol0(i-2,j))/(dx*cellx)

           !if (sigma4a*sigma4b .gt. 0.0) then
           !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
           ! 						if (sigma4 .eq. abs(sigma4a)) then
           ! 							sigma4 = sign(1.0,sigma4a)*sigma4
           ! 						end if
           ! 						if (sigma4 .eq. abs(sigma4b)) then
           ! 							sigma4 = sign(1.0,sigma4b)*sigma4
           ! 						end if
           ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
           ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
           ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
           sigma4 = ((MINLOC((/ABS(sigma4a), ABS(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((MINLOC((/ABS(sigma4b), ABS(sigma4a)/),DIM=1)-1.0)*sigma4a)
           !end if


           ! choosing sigma6
           sigma6 = 0.0
           IF (sigma2*sigma4 .GT. 0.0) THEN
              sigma6 = SIGN(1.0,sigma2)*MAXVAL((/ABS(sigma2), ABS(sigma4)/))
           END IF

           ! 						write(*,*) "sigma5"
           ! 						write(*,*) sigma5
           ! 						write(*,*) "sigma6"
           ! 						write(*,*) sigma6
           correction = (uTransport(i,j)*qx*0.5) * (sigma5 - sigma6) * (dx*cellx - uTransport(i,j)*qx*dx*cellx)
           solute_next_coarse(i,j) = solute_next_coarse(i,j) - correction











           !
        !    ! second term correction?
        !    !if (i .gt. 2) then
        !    !if (maskP(i-2,j) .ne. 0.0) then
        !    !sigma1 = 0.0
        !    sigma1a = (phiTransport(i+1,j) - phiTransport(i,j))/(dx*cellx)
        !    sigma1b = 2.0*(phiTransport(i,j) - phiTransport(i-1,j))/(dx*cellx)
           !
           !
        !    sigma1 = ((MINLOC((/ABS(sigma1a), ABS(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((MINLOC((/ABS(sigma1b), ABS(sigma1a)/),DIM=1)-1.0)*sigma1a)
        !    !end if
           !
        !    !sigma3 = 0.0
        !    sigma3a = 2.0*(phiTransport(i+1,j) - phiTransport(i,j))/(dx*cellx)
        !    sigma3b = (phiTransport(i,j) - phiTransport(i-1,j))/(dx*cellx)
           !
           !
        !    sigma3 = ((MINLOC((/ABS(sigma3a), ABS(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((MINLOC((/ABS(sigma3b), ABS(sigma3a)/),DIM=1)-1.0)*sigma3a)
        !    !end if
           !
        !    ! choosing sigma5
        !    sigma5 = 0.0
        !    IF (sigma1*sigma3 .GT. 0.0) THEN
        !       sigma5 = SIGN(1.0,sigma1)*MAXVAL((/ABS(sigma1), ABS(sigma3)/))
        !    END IF
           !
        !    !sigma2 = 0.0
        !    sigma2a = (phiTransport(i,j) - phiTransport(i-1,j))/(dx*cellx)
        !    sigma2b = 2.0*(phiTransport(i-1,j) - phiTransport(i-2,j))/(dx*cellx)
           !
           !
        !    sigma2 = ((MINLOC((/ABS(sigma2a), ABS(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((MINLOC((/ABS(sigma2b), ABS(sigma2a)/),DIM=1)-1.0)*sigma2a)
        !    !end if
           !
        !    !sigma4 = 0.0
        !    sigma4a = 2.0*(phiTransport(i,j) - phiTransport(i-1,j))/(dx*cellx)
        !    sigma4b = (phiTransport(i-1,j) - phiTransport(i-2,j))/(dx*cellx)
           !
           !
        !    sigma4 = ((MINLOC((/ABS(sigma4a), ABS(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((MINLOC((/ABS(sigma4b), ABS(sigma4a)/),DIM=1)-1.0)*sigma4a)
        !    !end if
           !
           !
        !    ! choosing sigma6
        !    sigma6 = 0.0
        !    IF (sigma2*sigma4 .GT. 0.0) THEN
        !       sigma6 = SIGN(1.0,sigma2)*MAXVAL((/ABS(sigma2), ABS(sigma4)/))
        !    END IF
           !
        !    ! 						write(*,*) "sigma5"
        !    ! 						write(*,*) sigma5
        !    ! 						write(*,*) "sigma6"
        !    ! 						write(*,*) sigma6
        !    correction = ((uTransport(i,j)/phiTransport(i,j))*sol0(i,j)*qx*0.5) * (sigma5 - sigma6) * (dx*cellx - (uTransport(i,j)/phiTransport(i,j))*sol0(i,j)*qx*dx*cellx)
        !    solute_next_coarse(i,j) = solute_next_coarse(i,j) - correction
           !
           !
           !






















           !end if ! end if maskP i-2,j .eq. 0
           !end if ! end if i .gt. 2
           ! end correction loop

           !end if ! end if u .gt. 0.0

           ! 			if (uTransport(i,j) .lt. 0.0) then
           ! 				! upwind including RHS value
           ! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i+1,j) - sol0(i,j) )
           !
           !
           ! 				! correction loop: sort of a mess
           ! 				!if (i .lt. xn-1) then
           ! 					if (maskP(i+2,j) .ne. 0.0) then
           ! 						sigma1 = 0.0
           ! 						sigma1a = (sol0(i+2,j) - sol0(i+1,j))/dx
           ! 						sigma1b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma1a*sigma1b .gt. 0.0) then
           ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
           ! ! 							if (sigma1 .eq. abs(sigma1a)) then
           ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
           ! ! 							end if
           ! ! 							if (sigma1 .eq. abs(sigma1b)) then
           ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
           ! ! 							end if
           ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
           ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
           ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
           ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
           ! 							!end if
           !
           ! 						sigma3 = 0.0
           ! 						sigma3a = 2.0*(sol0(i+2,j) - sol0(i+1,j))/dx
           ! 						sigma3b = (sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma3a*sigma3b .gt. 0.0) then
           ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
           ! ! 						if (sigma3 .eq. abs(sigma3a)) then
           ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
           ! ! 						end if
           ! ! 						if (sigma3 .eq. abs(sigma3b)) then
           ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
           ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
           ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
           ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
           ! 							!end if
           !
           !
           ! 						! choosing sigma5
           ! 						sigma5 = 0.0
           ! 						if (sigma1*sigma3 .gt. 0.0) then ! not sure about these conditionals...?
           ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
           ! 						end if
           !
           !
           !
           !
           !
           ! 						sigma2 = 0.0
           ! 						sigma2a = (sol0(i+1,j) - sol0(i,j))/dx
           ! 						sigma2b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma2a*sigma2b .gt. 0.0) then
           ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
           ! ! 						if (sigma2 .eq. abs(sigma2a)) then
           ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
           ! ! 						end if
           ! ! 						if (sigma2 .eq. abs(sigma2b)) then
           ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
           ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
           ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
           ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
           ! 							!end if
           !
           ! 						sigma4 = 0.0
           ! 						sigma4a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
           ! 						sigma4b = (sol0(i+1,j) - sol0(i,j))/dx
           !
           ! 						!if (sigma4a*sigma4b .gt. 0.0) then
           ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
           ! ! 						if (sigma4 .eq. abs(sigma4a)) then
           ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
           ! ! 						end if
           ! ! 						if (sigma4 .eq. abs(sigma4b)) then
           ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
           ! ! 						end if
           ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
           ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
           ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
           ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
           ! 							!end if
           !
           !
           ! 						! choosing sigma6
           ! 						sigma6 = 0.0
           ! 						if (sigma2*sigma4 .gt. 0.0) then
           ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
           ! 						end if
           !
           ! ! 						write(*,*) "sigma5"
           ! ! 						write(*,*) sigma5
           ! ! 						write(*,*) "sigma6"
           ! ! 						write(*,*) sigma6
           ! 						correction = (uTransport(i,j)*qx*0.5) * (sigma6 - sigma5) * (dx - uTransport(i,j)*qx*dx)
           ! 						solute_next(i,j) = solute_next(i,j) - correction
           ! 					end if
           ! 				!end if ! end if i .lt. xn-1
           ! 				! end correction loop
           !
           ! 			end if
           !
        END IF ! end mod thing
     END DO
     !solute_next(xn-1,j) = sol0(xn-1,j) - qx*uTransport(xn-1,j)*( sol0(xn-1,j) - sol0(xn-2,j) )
  END DO


  !sol = solute_next_coarse


  !
  ! do j = 1,yn/(2*celly)
  ! 	do i = 1,(xn-1)/cellx
  ! ! 		if ((coarse_mask(i,j) .eq. 0.0)) then
  ! !
  ! ! 			if (vTransport(i,j) .lt. 0.0) then
  ! ! 				sol(i,j+1) = seaw
  ! ! 			else
  ! ! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
  ! ! 			end if
  ! !
  ! ! 		end if
  !
  !
  !
  !
  !
  !
  !
  ! ! 		if ((maskP(i,j) .eq. 50.0)) then
  ! !
  ! ! ! 			if (vTransport(i,j+1) .lt. 0.0) then
  ! ! 				sol(i,j+1) = seaw
  ! ! ! 			else
  ! ! ! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
  ! ! ! 			end if
  ! !
  ! ! 		end if
  !
  !
  ! !
  ! ! 		if (maskP(i,j) .eq. 100.0) then
  ! !
  ! ! 				sol(i,j-1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j+1)
  ! !
  ! ! 		end if
  !
  !
  ! 	end do
  ! end do
  !
  ! sol0 = sol
  !
  !
  !
  ! do j = 1,(yn/(2*celly))-1
  ! 	! do i = 1,xn
  ! 	do i = 1,(xn-1)/cellx
  !
  ! 			if (vTransport(i,j) .gt. 1.0e-9) then
  ! 				! upwind including bottom value
  ! 				solute_next_coarse(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j) - sol0(i,j-1) )
  ! !
  ! 			! correction loop: sort of a mess
  ! 				!if (maskP(i,j-2) .ne. 0.0) then
  ! 						!sigma1 = 0.0
  ! 						sigma1a = (sol0(i,j+1) - sol0(i,j))/(dy*celly)
  ! 						sigma1b = 2.0*(sol0(i,j) - sol0(i,j-1))/(dy*celly)
  !
  ! 						!if (sigma1a*sigma1b .gt. 0.0) then
  ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
  ! ! 							if (sigma1 .eq. abs(sigma1a)) then
  ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
  ! ! 							end if
  ! ! 							if (sigma1 .eq. abs(sigma1b)) then
  ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
  ! ! 							end if
  ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
  ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
  ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
  ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
  ! 							!end if
  !
  ! 						!sigma3 = 0.0
  ! 						sigma3a = 2.0*(sol0(i,j+1) - sol0(i,j))/(dy*celly)
  ! 						sigma3b = (sol0(i,j) - sol0(i,j-1))/(dy*celly)
  !
  ! 						!if (sigma3a*sigma3b .gt. 0.0) then
  ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
  ! ! 						if (sigma3 .eq. abs(sigma3a)) then
  ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
  ! ! 						end if
  ! ! 						if (sigma3 .eq. abs(sigma3b)) then
  ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
  ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
  ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
  ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma5
  ! 						sigma5 = 0.0
  ! 						if (sigma1*sigma3 .gt. 0.0) then
  ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
  ! 						end if
  !
  !
  !
  !
  !
  ! 						!sigma2 = 0.0
  ! 						sigma2a = (sol0(i,j) - sol0(i,j-1))/(dy*celly)
  ! 						sigma2b = 2.0*(sol0(i,j-1) - sol0(i,j-2))/(dy*celly)
  !
  ! 						!if (sigma2a*sigma2b .gt. 0.0) then
  ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
  ! ! 						if (sigma2 .eq. abs(sigma2a)) then
  ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
  ! ! 						end if
  ! ! 						if (sigma2 .eq. abs(sigma2b)) then
  ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
  ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
  ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
  ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
  ! 							!end if
  !
  ! 						!sigma4 = 0.0
  ! 						sigma4a = 2.0*(sol0(i,j) - sol0(i,j-1))/(dy*celly)
  ! 						sigma4b = (sol0(i,j-1) - sol0(i,j-2))/(dy*celly)
  !
  ! 						!if (sigma4a*sigma4b .gt. 0.0) then
  ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
  ! ! 						if (sigma4 .eq. abs(sigma4a)) then
  ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
  ! ! 						end if
  ! ! 						if (sigma4 .eq. abs(sigma4b)) then
  ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
  ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
  ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
  ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma6
  ! 						sigma6 = 0.0
  ! 						if (sigma2*sigma4 .gt. 0.0) then
  ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
  ! 						end if
  !
  ! ! 						write(*,*) "sigma5"
  ! ! 						write(*,*) sigma5
  ! ! 						write(*,*) "sigma6"
  ! ! 						write(*,*) sigma6
  ! 						correction = (vTransport(i,j)*qy*0.5) * (sigma5 - sigma6) * ((dy*celly) - vTransport(i,j)*qy*(dy*celly))
  ! 						solute_next_coarse(i,j) = solute_next_coarse(i,j) - correction
  ! 				!end if ! end if maskP i,j-2 .ne. 0.0
  ! 				! end correction loop
  !
  !
  ! 			end if
  !
  ! 			if (vTransport(i,j) .lt. -1.0e-9) then
  ! 				! upwind including top value
  ! 				solute_next_coarse(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j+1) - sol0(i,j) )
  !
  ! 				! correction loop: sort of a mess
  ! 				!if (j .lt. yn-1) then
  ! 					!if (maskP(i,j+2) .ne. 0.0) then
  ! 						!sigma1 = 0.0
  ! 						sigma1a = (sol0(i,j+2) - sol0(i,j+1))/(dy*celly)
  ! 						sigma1b = 2.0*(sol0(i,j+1) - sol0(i,j))/(dy*celly)
  !
  ! 						!if (sigma1a*sigma1b .gt. 0.0) then
  ! !							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
  ! ! 							if (sigma1 .eq. abs(sigma1a)) then
  ! ! 								sigma1 = sign(1.0,sigma1a)*sigma1
  ! ! 							end if
  ! ! 							if (sigma1 .eq. abs(sigma1b)) then
  ! ! 								sigma1 = sign(1.0,sigma1b)*sigma1
  ! ! 							end if
  ! ! 							sig_bool_a = th_bool(sigma1 .eq. abs(sigma1a))
  ! ! 							sig_bool_b = th_bool(sigma1 .eq. abs(sigma1b))
  ! ! 							sigma1 = sig_bool_a*sign(1.0,sigma1a)*sigma1 + sig_bool_b*sign(1.0,sigma1b)*sigma1
  ! 							sigma1 = ((minloc((/abs(sigma1a), abs(sigma1b)/),DIM=1)-1.0)*sigma1b) + ((minloc((/abs(sigma1b), abs(sigma1a)/),DIM=1)-1.0)*sigma1a)
  ! 							!end if
  !
  ! 						!sigma3 = 0.0
  ! 						sigma3a = 2.0*(sol0(i,j+2) - sol0(i,j+1))/(dy*celly)
  ! 						sigma3b = (sol0(i,j+1) - sol0(i,j))/(dy*celly)
  !
  ! 						!if (sigma3a*sigma3b .gt. 0.0) then
  ! !						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
  ! ! 						if (sigma3 .eq. abs(sigma3a)) then
  ! ! 							sigma3 = sign(1.0,sigma3a)*sigma3
  ! ! 						end if
  ! ! 						if (sigma3 .eq. abs(sigma3b)) then
  ! ! 							sigma3 = sign(1.0,sigma3b)*sigma3
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma3 .eq. abs(sigma3a))
  ! ! 						sig_bool_b = th_bool(sigma3 .eq. abs(sigma3b))
  ! ! 						sigma3 = sig_bool_a*sign(1.0,sigma3a)*sigma3 + sig_bool_b*sign(1.0,sigma3b)*sigma3
  ! 						sigma3 = ((minloc((/abs(sigma3a), abs(sigma3b)/),DIM=1)-1.0)*sigma3b) + ((minloc((/abs(sigma3b), abs(sigma3a)/),DIM=1)-1.0)*sigma3a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma5
  ! 						sigma5 = 0.0
  ! 						if (sigma1*sigma3 .gt. 0.0) then
  ! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
  ! 						end if
  !
  !
  !
  !
  !
  ! 						!sigma2 = 0.0
  ! 						sigma2a = (sol0(i,j+1) - sol0(i,j))/(dy*celly)
  ! 						sigma2b = 2.0*(sol0(i,j) - sol0(i,j-1))/(dy*celly)
  !
  ! 						!if (sigma2a*sigma2b .gt. 0.0) then
  ! !						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
  ! ! 						if (sigma2 .eq. abs(sigma2a)) then
  ! ! 							sigma2 = sign(1.0,sigma2a)*sigma2
  ! ! 						end if
  ! ! 						if (sigma2 .eq. abs(sigma2b)) then
  ! ! 							sigma2 = sign(1.0,sigma2b)*sigma2
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma2 .eq. abs(sigma2a))
  ! ! 						sig_bool_b = th_bool(sigma2 .eq. abs(sigma2b))
  ! ! 						sigma2 = sig_bool_a*sign(1.0,sigma2a)*sigma2 + sig_bool_b*sign(1.0,sigma2b)*sigma2
  ! 						sigma2 = ((minloc((/abs(sigma2a), abs(sigma2b)/),DIM=1)-1.0)*sigma2b) + ((minloc((/abs(sigma2b), abs(sigma2a)/),DIM=1)-1.0)*sigma2a)
  ! 							!end if
  !
  ! 						!sigma4 = 0.0
  ! 						sigma4a = 2.0*(sol0(i,j+1) - sol0(i,j))/(dy*celly)
  ! 						sigma4b = (sol0(i,j) - sol0(i,j-1))/(dy*celly)
  !
  ! 						!if (sigma4a*sigma4b .gt. 0.0) then
  ! !						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
  ! ! 						if (sigma4 .eq. abs(sigma4a)) then
  ! ! 							sigma4 = sign(1.0,sigma4a)*sigma4
  ! ! 						end if
  ! ! 						if (sigma4 .eq. abs(sigma4b)) then
  ! ! 							sigma4 = sign(1.0,sigma4b)*sigma4
  ! ! 						end if
  ! ! 						sig_bool_a = th_bool(sigma4 .eq. abs(sigma4a))
  ! ! 						sig_bool_b = th_bool(sigma4 .eq. abs(sigma4b))
  ! ! 						sigma4 = sig_bool_a*sign(1.0,sigma4a)*sigma4 + sig_bool_b*sign(1.0,sigma4b)*sigma4
  ! 						sigma4 = ((minloc((/abs(sigma4a), abs(sigma4b)/),DIM=1)-1.0)*sigma4b) + ((minloc((/abs(sigma4b), abs(sigma4a)/),DIM=1)-1.0)*sigma4a)
  ! 							!end if
  !
  !
  ! 						! choosing sigma6
  ! 						sigma6 = 0.0
  ! 						if (sigma2*sigma4 .gt. 0.0) then
  ! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
  ! 						end if
  !
  ! ! 						write(*,*) "sigma5"
  ! ! 						write(*,*) sigma5
  ! ! 						write(*,*) "sigma6"
  ! ! 						write(*,*) sigma6
  ! 						correction = (vTransport(i,j)*qy*0.5) * (sigma6 - sigma5) * ((dy*celly) - vTransport(i,j)*qy*(dy*celly))
  ! 						solute_next_coarse(i,j) = solute_next_coarse(i,j) - correction
  !
  ! 					!end if
  ! 				!end if ! end if j .lt. yn-1
  ! 				! end correction loop
  !
  !
  ! 			end if
  !
  !
  ! 	end do
  ! 	solute_next_coarse(i,yn/(2*celly)-1) = sol0(i,yn/(2*celly)-1) - qy*vTransport(i,yn/(2*celly)-1)*( sol0(i,yn/(2*celly)) - sol0(i,yn/(2*celly)-1) )
  ! end do




  RETURN

END FUNCTION solute_next_coarse









! function solute_next_coarse (sol_coarse, uTransport, vTransport, seaw)
!
! use globals
! use initialize
! implicit none
!
! ! declare errthing
!
! ! integers
! integer :: i, j, ii, n, m
! ! inputs
! real(4) :: sol((xn-1)/cellx,yn/celly), sol0((xn-1)/cellx,yn/celly)
! real(4) :: uTransport((xn-1)/cellx,yn/celly), vTransport((xn-1)/cellx,yn/celly)
! ! solver stuff
! real(4) :: uLong((((xn-1)/cellx)-2)*((yn/celly)-0)), vLong((((xn-1)/cellx)-0)*((yn/celly)-2))
! real(4) :: aBand((((xn-1)/cellx)-2)*((yn/celly)-0),5), bBand((((xn-1)/cellx)-0)*((yn/celly)-2),5)
! real(4) :: qx, qy, solute_next_coarse((xn-1)/cellx,yn/celly), vec((((xn-1)/cellx)-2)*((yn/celly)-0))
! real(4) :: sol_nextRow((((xn-1)/cellx)-2)*((yn/celly)-0)), sol_nextRowB((((xn-1)/cellx)-0)*((yn/celly)-2))
! real(4) :: seaw
! real(4) :: bm1((xn-1)/cellx,yn/celly), b0((xn-1)/cellx,yn/celly), bp1((xn-1)/cellx,yn/celly), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
! real(4) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6
!
!
!
! do i = 1,(xn-1)/cellx
! 	do j = 1,yn/celly
! 		if ((coarse_mask(i,j) .eq. 0.0)) then
! 			sol(i,j) = seaw
! 		end if
! 	end do
! end do
!
! sol(1,:) = (4.0/3.0)*sol(2,:) - (1.0/3.0)*sol(3,:)
! sol((xn-1)/cellx,:) = (4.0/3.0)*sol((xn-1)/cellx-1,:) - (1.0/3.0)*sol((xn-1)/cellx-2,:)
!
! sol0 = sol
! solute_next_coarse = sol
!
! qx = dt*mstep/(cstep*dx*cellx)
! qy = dt*mstep/(cstep*dy*celly)
!
! ! uLong = reshape(uTransport(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
! ! !! transpose coarse needed!
! ! vLong = reshape(transpose(vTransport(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))
!
! uTransport(1,:) = 0.0
! uTransport(xn/cell,:) = 0.0
!
! ! write(*,*) qx*maxval(abs(uTransport))
! ! write(*,*) qy*maxval(abs(vTransport))
!
!
!
! do i = 2,xn-1
! 	do j = yn/2,yn
! 		if ((coarse_mask(i,j) .eq. 5.0) .or. (maskP(i,j) .eq. 120.5)) then
!
! 			if (uTransport(i,j) .lt. 0.0) then
! 				sol(i+1,j) = seaw
! 				uTransport(i+1,j) = uTransport(i,j)
! 			else
! 				sol(i+1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i-1,j)
! 			end if
!
! 		end if
!
! 		if ((maskP(i,j) .eq. 10.0) .or. (maskP(i,j) .eq. 170.5)) then
! 			if (uTransport(i,j) .gt. 0.0) then
! 				sol(i-1,j) = seaw
! 				uTransport(i-1,j) = uTransport(i,j)
! 			else
! 				sol(i-1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i+1,j)
! 			end if
! 		end if
! 	end do
! end do
!
!
!
!
! do i = 2,xn-1
! 	do j = yn/2,yn-1
!
!
!
! 			if (uTransport(i,j) .gt. 0.0) then
! 				! upwind including LHS value
! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i,j) - sol0(i-1,j) )
!
!
! 				! correction loop: sort of a mess
! 				if (i .gt. 2) then
! 					if (maskP(i-2,j) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i+1,j) - sol0(i,j))/dx
! 						sigma1b = 2.0*(sol0(i,j) - sol0(i-1,j))/dx
!
! 						!if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
! 							if (sigma1 .eq. abs(sigma1a)) then
! 								sigma1 = sign(1.0,sigma1a)*sigma1
! 							end if
! 							if (sigma1 .eq. abs(sigma1b)) then
! 								sigma1 = sign(1.0,sigma1b)*sigma1
! 							end if
! 							!end if
!
! 						sigma3 = 0.0
! 						sigma3a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
! 						sigma3b = (sol0(i,j) - sol0(i-1,j))/dx
!
! 						!if (sigma3a*sigma3b .gt. 0.0) then
! 						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
! 						if (sigma3 .eq. abs(sigma3a)) then
! 							sigma3 = sign(1.0,sigma3a)*sigma3
! 						end if
! 						if (sigma3 .eq. abs(sigma3b)) then
! 							sigma3 = sign(1.0,sigma3b)*sigma3
! 						end if
! 							!end if
!
!
! 						! choosing sigma5
! 						sigma5 = 0.0
! 						!if (sigma1*sigma3 .gt. 0.0) then
! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
! 							!end if
!
!
!
!
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i,j) - sol0(i-1,j))/dx
! 						sigma2b = 2.0*(sol0(i-1,j) - sol0(i-2,j))/dx
!
! 						!if (sigma2a*sigma2b .gt. 0.0) then
! 						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
! 						if (sigma2 .eq. abs(sigma2a)) then
! 							sigma2 = sign(1.0,sigma2a)*sigma2
! 						end if
! 						if (sigma2 .eq. abs(sigma2b)) then
! 							sigma2 = sign(1.0,sigma2b)*sigma2
! 						end if
! 							!end if
!
! 						sigma4 = 0.0
! 						sigma4a = 2.0*(sol0(i,j) - sol0(i-1,j))/dx
! 						sigma4b = (sol0(i-1,j) - sol0(i-2,j))/dx
!
! 						!if (sigma4a*sigma4b .gt. 0.0) then
! 						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
! 						if (sigma4 .eq. abs(sigma4a)) then
! 							sigma4 = sign(1.0,sigma4a)*sigma4
! 						end if
! 						if (sigma4 .eq. abs(sigma4b)) then
! 							sigma4 = sign(1.0,sigma4b)*sigma4
! 						end if
! 							!end if
!
!
! 						! choosing sigma6
! 						sigma6 = 0.0
! 						!if (sigma2*sigma4 .gt. 0.0) then
! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
! 							!end if
!
! ! 						write(*,*) "sigma5"
! ! 						write(*,*) sigma5
! ! 						write(*,*) "sigma6"
! ! 						write(*,*) sigma6
! 						correction = (uTransport(i,j)*qx*0.5) * (sigma5 - sigma6) * (dx - uTransport(i,j)*qx*dx)
! 						solute_next(i,j) = solute_next(i,j) - correction
!
! 					end if
! 				end if
! 				! end correction loop
!
! 			end if
!
! 			if (uTransport(i,j) .lt. 0.0) then
! 				! upwind including RHS value
! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i+1,j) - sol0(i,j) )
!
!
! 				! correction loop: sort of a mess
! 				if (i .lt. xn-1) then
! 					if (maskP(i+2,j) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i+2,j) - sol0(i+1,j))/dx
! 						sigma1b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
!
! 						!if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
! 							if (sigma1 .eq. abs(sigma1a)) then
! 								sigma1 = sign(1.0,sigma1a)*sigma1
! 							end if
! 							if (sigma1 .eq. abs(sigma1b)) then
! 								sigma1 = sign(1.0,sigma1b)*sigma1
! 							end if
! 							!end if
!
! 						sigma3 = 0.0
! 						sigma3a = 2.0*(sol0(i+2,j) - sol0(i+1,j))/dx
! 						sigma3b = (sol0(i+1,j) - sol0(i,j))/dx
!
! 						!if (sigma3a*sigma3b .gt. 0.0) then
! 						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
! 						if (sigma3 .eq. abs(sigma3a)) then
! 							sigma3 = sign(1.0,sigma3a)*sigma3
! 						end if
! 						if (sigma3 .eq. abs(sigma3b)) then
! 							sigma3 = sign(1.0,sigma3b)*sigma3
! 						end if
! 							!end if
!
!
! 						! choosing sigma5
! 						sigma5 = 0.0
! 						if (sigma1*sigma3 .gt. 0.0) then
! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
! 						end if
!
!
!
!
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i+1,j) - sol0(i,j))/dx
! 						sigma2b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
!
! 						!if (sigma2a*sigma2b .gt. 0.0) then
! 						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
! 						if (sigma2 .eq. abs(sigma2a)) then
! 							sigma2 = sign(1.0,sigma2a)*sigma2
! 						end if
! 						if (sigma2 .eq. abs(sigma2b)) then
! 							sigma2 = sign(1.0,sigma2b)*sigma2
! 						end if
! 							!end if
!
! 						sigma4 = 0.0
! 						sigma4a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
! 						sigma4b = (sol0(i+1,j) - sol0(i,j))/dx
!
! 						!if (sigma4a*sigma4b .gt. 0.0) then
! 						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
! 						if (sigma4 .eq. abs(sigma4a)) then
! 							sigma4 = sign(1.0,sigma4a)*sigma4
! 						end if
! 						if (sigma4 .eq. abs(sigma4b)) then
! 							sigma4 = sign(1.0,sigma4b)*sigma4
! 						end if
! 							!end if
!
!
! 						! choosing sigma6
! 						sigma6 = 0.0
! 						if (sigma2*sigma4 .gt. 0.0) then
! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
! 						end if
!
! ! 						write(*,*) "sigma5"
! ! 						write(*,*) sigma5
! ! 						write(*,*) "sigma6"
! ! 						write(*,*) sigma6
! 						correction = (uTransport(i,j)*qx*0.5) * (sigma6 - sigma5) * (dx - uTransport(i,j)*qx*dx)
! 						solute_next(i,j) = solute_next(i,j) - correction
! 					end if
! 				end if
! 				! end correction loop
!
! 			end if
!
!
! 	end do
! end do
!
!
! sol = solute_next
!
!
!
!
! do i = 1,xn
! 	do j = yn/2,yn
! 		if ((maskP(i,j) .eq. 25.0) .or. (maskP(i,j) .eq. 51.0) .or. (maskP(i,j) .eq. 120.5) .or. (maskP(i,j) .eq. 170.5) ) then
!
! 			if (vTransport(i,j+1) .lt. 0.0) then
! 				sol(i,j+1) = seaw
! 			else
! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
! 			end if
!
! 		end if
!
! 		if ((maskP(i,j) .eq. 50.0)) then
!
! ! 			if (vTransport(i,j+1) .lt. 0.0) then
! 				sol(i,j+1) = seaw
! ! 			else
! ! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
! ! 			end if
!
! 		end if
!
!
!
! 		if (maskP(i,j) .eq. 100.0) then
!
! 				sol(i,j-1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j+1)
!
! 		end if
!
!
! 	end do
! end do
!
! sol0 = sol
!
!
! do i = 1,xn
! 	do j = yn/2,yn-1
!
!
!
! 			if (vTransport(i,j) .gt. 0.0) then
! 				! upwind including bottom value
! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j) - sol0(i,j-1) )
! !
! 			! correction loop: sort of a mess
! 				if (maskP(i,j-2) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i,j+1) - sol0(i,j))/dy
! 						sigma1b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
!
! 						!if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
! 							if (sigma1 .eq. abs(sigma1a)) then
! 								sigma1 = sign(1.0,sigma1a)*sigma1
! 							end if
! 							if (sigma1 .eq. abs(sigma1b)) then
! 								sigma1 = sign(1.0,sigma1b)*sigma1
! 							end if
! 							!end if
!
! 						sigma3 = 0.0
! 						sigma3a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
! 						sigma3b = (sol0(i,j) - sol0(i,j-1))/dy
!
! 						!if (sigma3a*sigma3b .gt. 0.0) then
! 						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
! 						if (sigma3 .eq. abs(sigma3a)) then
! 							sigma3 = sign(1.0,sigma3a)*sigma3
! 						end if
! 						if (sigma3 .eq. abs(sigma3b)) then
! 							sigma3 = sign(1.0,sigma3b)*sigma3
! 						end if
! 							!end if
!
!
! 						! choosing sigma5
! 						sigma5 = 0.0
! 						if (sigma1*sigma3 .gt. 0.0) then
! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
! 						end if
!
!
!
!
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i,j) - sol0(i,j-1))/dy
! 						sigma2b = 2.0*(sol0(i,j-1) - sol0(i,j-2))/dy
!
! 						!if (sigma2a*sigma2b .gt. 0.0) then
! 						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
! 						if (sigma2 .eq. abs(sigma2a)) then
! 							sigma2 = sign(1.0,sigma2a)*sigma2
! 						end if
! 						if (sigma2 .eq. abs(sigma2b)) then
! 							sigma2 = sign(1.0,sigma2b)*sigma2
! 						end if
! 							!end if
!
! 						sigma4 = 0.0
! 						sigma4a = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
! 						sigma4b = (sol0(i,j-1) - sol0(i,j-2))/dy
!
! 						!if (sigma4a*sigma4b .gt. 0.0) then
! 						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
! 						if (sigma4 .eq. abs(sigma4a)) then
! 							sigma4 = sign(1.0,sigma4a)*sigma4
! 						end if
! 						if (sigma4 .eq. abs(sigma4b)) then
! 							sigma4 = sign(1.0,sigma4b)*sigma4
! 						end if
! 							!end if
!
!
! 						! choosing sigma6
! 						sigma6 = 0.0
! 						if (sigma2*sigma4 .gt. 0.0) then
! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
! 						end if
!
! ! 						write(*,*) "sigma5"
! ! 						write(*,*) sigma5
! ! 						write(*,*) "sigma6"
! ! 						write(*,*) sigma6
! 						correction = (vTransport(i,j)*qy*0.5) * (sigma5 - sigma6) * (dy - vTransport(i,j)*qy*dy)
! 						solute_next(i,j) = solute_next(i,j) - correction
! 				end if
! 				! end correction loop
!
!
! 			end if
!
! 			if (vTransport(i,j) .lt. 0.0) then
! 				! upwind including top value
! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j+1) - sol0(i,j) )
!
! 				! correction loop: sort of a mess
! 				if (j .lt. yn-1) then
! 					if (maskP(i,j-2) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i,j+2) - sol0(i,j+1))/dy
! 						sigma1b = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
!
! 						!if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
! 							if (sigma1 .eq. abs(sigma1a)) then
! 								sigma1 = sign(1.0,sigma1a)*sigma1
! 							end if
! 							if (sigma1 .eq. abs(sigma1b)) then
! 								sigma1 = sign(1.0,sigma1b)*sigma1
! 							end if
! 							!end if
!
! 						sigma3 = 0.0
! 						sigma3a = 2.0*(sol0(i,j+2) - sol0(i,j+1))/dy
! 						sigma3b = (sol0(i,j+1) - sol0(i,j))/dy
!
! 						!if (sigma3a*sigma3b .gt. 0.0) then
! 						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
! 						if (sigma3 .eq. abs(sigma3a)) then
! 							sigma3 = sign(1.0,sigma3a)*sigma3
! 						end if
! 						if (sigma3 .eq. abs(sigma3b)) then
! 							sigma3 = sign(1.0,sigma3b)*sigma3
! 						end if
! 							!end if
!
!
! 						! choosing sigma5
! 						sigma5 = 0.0
! 						if (sigma1*sigma3 .gt. 0.0) then
! 							sigma5 = sign(1.0,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
! 						end if
!
!
!
!
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i,j+1) - sol0(i,j))/dy
! 						sigma2b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
!
! 						!if (sigma2a*sigma2b .gt. 0.0) then
! 						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
! 						if (sigma2 .eq. abs(sigma2a)) then
! 							sigma2 = sign(1.0,sigma2a)*sigma2
! 						end if
! 						if (sigma2 .eq. abs(sigma2b)) then
! 							sigma2 = sign(1.0,sigma2b)*sigma2
! 						end if
! 							!end if
!
! 						sigma4 = 0.0
! 						sigma4a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
! 						sigma4b = (sol0(i,j) - sol0(i,j-1))/dy
!
! 						!if (sigma4a*sigma4b .gt. 0.0) then
! 						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
! 						if (sigma4 .eq. abs(sigma4a)) then
! 							sigma4 = sign(1.0,sigma4a)*sigma4
! 						end if
! 						if (sigma4 .eq. abs(sigma4b)) then
! 							sigma4 = sign(1.0,sigma4b)*sigma4
! 						end if
! 							!end if
!
!
! 						! choosing sigma6
! 						sigma6 = 0.0
! 						if (sigma2*sigma4 .gt. 0.0) then
! 							sigma6 = sign(1.0,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
! 						end if
!
! ! 						write(*,*) "sigma5"
! ! 						write(*,*) sigma5
! ! 						write(*,*) "sigma6"
! ! 						write(*,*) sigma6
! 						correction = (vTransport(i,j)*qy*0.5) * (sigma6 - sigma5) * (dy - vTransport(i,j)*qy*dy)
! 						solute_next(i,j) = solute_next(i,j) - correction
!
! 					end if
! 				end if
! 				! end correction loop
!
!
! 			end if
!
!
! 	end do
! end do
!
!
!
!
! return
!
! end function solute_next_coarse



! ! ----------------------------------------------------------------------------------%%
! !
! ! ALT_NEXT
! !
! ! SUMMARY: solves for equilibrium at a single grid cell using PHREEQC
! !
! ! INPUTS: temp : temperature of grid cell
! !         timestep : time elapsed
! !         primaryList(5) : amounts of primary minerals
! !         secondaryList(16) : amounts of secondary minerals
! !         soluteList(11) : concentrations of solutes
! !
! ! RETURNS: alt_next(1,altnum): returns everything from PHREEQC in a big pile
! !          and it gets parsed in the main method's geochem loop
! !
! ! ----------------------------------------------------------------------------------%%
!
! function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
! use globals
! use initialize
! use alteration
! implicit none
!
! interface
!
! end interface
!
! ! declare errthing
!
! integer :: order
! real(4) :: temp, timestep
! real(4) :: alt_next(1,167)
! real(4) :: alter0(1,167)
! real(4) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)
!
! ! use the alteration module
! alter0 = alter(temp-272.9, timestep, primaryList, secondaryList, soluteList, mediumList)
!
! ! rename it for a reason that i now forget
! alt_next = alter0
!
! end function alt_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

FUNCTION rho_next(h_in)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: i,j
  REAL(4) :: h_in(xn,yn), rho_next(xn,yn)


  DO i=1,xn
     DO j = 1,yn
        rho_next(i,j) = rho_fluid*(1.0 - alpha*(h_in(i,j)-273.0))
     END DO
  END DO

  RETURN

END FUNCTION rho_next



! ----------------------------------------------------------------------------------%%
!
! RHO_ONE
!
!
! ----------------------------------------------------------------------------------%%

FUNCTION rho_one(h_in)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: i,j
  REAL(4) :: h_in, rho_one


  rho_one = rho_fluid*(1.0 - alpha*((h_in-273.0)))


  RETURN

END FUNCTION rho_one


! ----------------------------------------------------------------------------------%%
!
! VISC_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

FUNCTION visc_next(h_in)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: i,j
  REAL(4) :: h_in(xn,yn), visc_next(xn,yn)


  DO i=1,xn
     DO j = 1,yn
        !visc_next(i,j) = (2.44e-5)*10.0**(247.8/(h_in(i,j)-140.0))
        visc_next(i,j) = .001!.001!.0002 !.00350*exp(-(h_in(i,j)-273.16)/35.0) + .0002
        ! from nist.gov
     END DO
  END DO

  RETURN

END FUNCTION visc_next




! ----------------------------------------------------------------------------------%%
!
! VELOCITIES
!
! SUMMARY : computes the darcy velocity (specific discharge) from the streamfunction
!           using finite difference partial derivatives
!
! INPUTS : psi(xn,yn) : 2D streamfunction array of current timestep
!
! RETURNS : velocities(xn,2*yn) : both u and v velocities in one matrix
!
! ----------------------------------------------------------------------------------%%


FUNCTION velocities(psi)

  USE globals
  USE initialize
  IMPLICIT NONE

  INTERFACE

     FUNCTION partial_edge(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge(rows,cols)
     END FUNCTION partial_edge

     FUNCTION partial_edge_p(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_edge_p(rows,cols)
     END FUNCTION partial_edge_p

  END INTERFACE

  ! declare errthing
  INTEGER :: i,ii
  REAL(4) :: velocities(xn,2*yn), psi(xn,yn)
  REAL(4) :: u0(xn,yn), v0(xn,yn)

  u0 = partial_edge_p(psi,xn,yn,dx,dy,2)
  v0 = -partial_edge_p(psi,xn,yn,dx,dy,1)

  velocities(1:xn,1:yn) = u0
  velocities(1:xn,yn+1:2*yn) = v0

  RETURN
END FUNCTION velocities




! ----------------------------------------------------------------------------------%%
!
! VELOCITIES_COARSE
!
! ----------------------------------------------------------------------------------%%



! calculates velocities from COARSE streamfunction values
FUNCTION velocities_coarse(psi_coarse)
  USE globals
  USE initialize
  IMPLICIT NONE
  INTEGER :: i,ii
  REAL(4) :: velocities_coarse((xn-1)/cellx,yn/celly), psi_coarse((xn-1)/cellx,yn/(2*celly))
  REAL(4) :: u0((xn-1)/cellx,yn/(2*celly)), v0((xn-1)/cellx,yn/(2*celly))

  INTERFACE

     FUNCTION partial_coarse(array,rows,cols,d1,d2,dim)
       USE globals
       USE initialize
       IMPLICIT NONE
       INTEGER :: rows, cols, dim, i, j, ii, jj
       REAL(4) :: array(rows,cols), d1, d2, d
       REAL(4) :: partial_coarse(rows,cols)
     END FUNCTION partial_coarse

  END INTERFACE

  u0 = partial_coarse(psi_coarse,(xn-1)/cellx,yn/(2*celly),dx*cellx,dy*celly,2)
  v0 = -partial_coarse(psi_coarse,(xn-1)/cellx,yn/(2*celly),dx*cellx,dy*celly,1)

  ! do i =1,xn
  ! 	do ii = 1,yn
  ! 		if (mask(i,ii) .eq. 50.0) then
  ! 			u0(i,ii+1) = 0.0
  ! 		end if
  ! 	end do
  ! end do

  velocities_coarse(1:(xn-1)/cellx,1:yn/(2*celly)) = u0
  velocities_coarse(1:(xn-1)/cellx,yn/(2*celly)+1:yn/celly) = v0

  RETURN
END FUNCTION velocities_coarse



! ----------------------------------------------------------------------------------%%
!
! PARTIAL
!
! ----------------------------------------------------------------------------------%%

FUNCTION partial(array,rows,cols,d1,d2,dim)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: rows, cols, dim, i, j, ii, jj
  REAL(4) :: array(rows,cols), d1, d2, d
  REAL(4) :: partial(rows,cols)

  ! write(*,*) "dim"
  ! write(*,*) rows
  ! write(*,*) cols

  ! figure out which direction derivative goes (dx or dy)

  partial = 0.0

  IF (dim .EQ. 1) THEN
     ii = 1
     jj = 0
     d = d1

     ! compute edges beforehand
     partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
     partial(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

     DO i = 2,rows-1
        DO j = 1,cols
           partial(i,j) = (array(i+1,j) - array(i-1,j))/(2.0*d)
        END DO
     END DO

     DO i = 2,rows-1
        DO j = 1,cols
           IF ((maskP(i,j) .EQ. 3.0) .OR. (maskP(i,j) .EQ. 3.5) .OR. (maskP(i,j) .EQ. 3.1) .OR. (maskP(i,j) .EQ. 3.05) .OR. (maskP(i,j) .EQ. 3.01) .OR. (mask(i,j) .EQ. 3.05)) THEN
              partial(i,j) = (array(i,j) - array(i-1,j))/d
              !partial(i,j) = (3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j))/(2.0*d)
           END IF
           IF ((maskP(i,j) .EQ. 6.0) .OR. (maskP(i,j) .EQ. 6.5) .OR. (maskP(i,j) .EQ. 6.1) .OR. (maskP(i,j) .EQ. 6.05) .OR. (maskP(i,j) .EQ. 6.01) .OR. (mask(i,j) .EQ. 6.05)) THEN
              partial(i,j) = (array(i+1,j) - array(i,j))/d
              !partial(i,j) = (-3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j))/(2.0*d)
           END IF
        END DO
     END DO

  END IF


  IF (dim .EQ. 2) THEN
     ii = 0
     jj = 1
     d = d2

     ! compute edges beforehand
     partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
     partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
  END IF

  ! ! use central difference method ignoring edges (already done)
  ! do i = 2-jj,rows-1+jj
  !     do j = 2-ii,cols-1+ii
  !     	partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
  ! 	end do
  ! end do


  ! if (dim .eq. 1) then
  ! 	do i=2,rows-1
  ! 		do j=2,cols-1
  !
  ! 			if (maskP(i,j) .eq. 5.0) then
  ! 				partial(i,j) = (array(i,j) - array(i-1,j))/d
  ! 				partial(i+1,j) = (array(i+1,j) - array(i,j))/d
  ! 			end if
  !
  ! 			if (maskP(i,j) .eq. 10.0) then
  ! 				partial(i,j) = (array(i+1,j) - array(i,j))/d
  ! 				partial(i-1,j) = (array(i,j) - array(i-1,j))/d
  ! 			end if
  !
  !
  ! 		end do
  ! 	end do
  !
  !
  ! end if

  ! ! if jj=1, dim=2, cols=yn, rows=xn

  ! ! compute edges beforehand
  ! partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
  ! partial(:,yn) = ( 3.0*array(:,yn) - 4.0*array(:,yn-1) + array(:,yn-2) ) / (2.0*d)


  ! do i = 1,xn
  !     do j = 2,yn-1
  !     	partial(i,j) = (array(i,j+1)-array(i,j-1))/(2.0*d)
  ! 	end do
  ! end do






  ! ! if ii=1, dim=1, cols=yn, rows=xn

  ! ! compute edges beforehand
  ! partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
  ! partial(xn,:) = ( 3.0*array(xn,:) - 4.0*array(xn-1,:) + array(xn-2,:) ) / (2.0*d)


  ! do i = 2,xn-1
  !     do j = 1,yn
  !     	partial(i,j) = (array(i+1,j)-array(i-1,j))/(2.0*d)
  ! 	end do
  ! end do

  RETURN
END FUNCTION partial










! ----------------------------------------------------------------------------------%%
!
! PARTIAL_COARSE
!
! ----------------------------------------------------------------------------------%%

FUNCTION partial_coarse(array,rows,cols,d1,d2,dim)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: rows, cols, dim, i, j, ii, jj
  REAL(4) :: array(rows,cols), d1, d2, d
  REAL(4) :: partial_coarse(rows,cols)


  partial_coarse = 0.0

  IF (dim .EQ. 1) THEN
     ii = 1
     jj = 0
     d = d1

     ! compute edges beforehand
     partial_coarse(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
     partial_coarse(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

     DO i = 2,rows-1
        DO j = 1,cols
           partial_coarse(i,j) = (array(i+1,j) - array(i-1,j))/(2.0*d)
        END DO
     END DO


  END IF




  IF (dim .EQ. 2) THEN
     ii = 0
     jj = 1
     d = d2

     ! compute edges beforehand
     partial_coarse(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
     partial_coarse(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)

     DO i = 1,rows
        DO j = 2,cols-1
           IF ((coarse_mask(i,j) .NE. 0.0) .AND. (coarse_mask(i,j+1) .NE. 0.0) .AND. (coarse_mask(i,j-1) .NE. 0.0)) THEN
              partial_coarse(i,j) = (array(i,j+1) - array(i,j-1))/(2.0*d)
           END IF

           IF ((coarse_mask(i,j) .NE. 0.0) .AND. (coarse_mask(i,j+1) .NE. 0.0) .AND. (coarse_mask(i,j-1) .EQ. 0.0)) THEN
              partial_coarse(i,j) = (array(i,j+1) - array(i,j))/(1.0*d)
           END IF

           IF ((coarse_mask(i,j) .NE. 0.0) .AND. (coarse_mask(i,j+1) .EQ. 0.0) .AND. (coarse_mask(i,j-1) .NE. 0.0)) THEN
              partial_coarse(i,j) = (array(i,j) - array(i,j-1))/(1.0*d)
           END IF

        END DO
     END DO

  END IF




  RETURN
END FUNCTION partial_coarse










! ----------------------------------------------------------------------------------%%
!
! PARTIAL_EDGE
!
! ----------------------------------------------------------------------------------%%

FUNCTION partial_edge(array,rows,cols,d1,d2,dim)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: rows, cols, dim, i, j, ii, jj
  REAL(4) :: array(rows,cols), d1, d2, d
  REAL(4) :: partial_edge(rows,cols)

  ! write(*,*) "dim"
  ! write(*,*) rows
  ! write(*,*) cols


  ! figure out which direction derivative goes (dx or dy)
  IF (dim .EQ. 1) THEN
     ii = 1
     jj = 0
     d = d1
  END IF

  IF (dim .EQ. 2) THEN
     ii = 0
     jj = 1
     d = d2
  END IF
  partial_edge = 0.0
  ! use central difference method ignoring edges (already done)
  DO i = 2-jj,rows-1+jj
     DO j = 2-ii,cols-1+ii
        IF ((mask(i,j) .NE. 0.0)) THEN
           partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
        END IF
     END DO
  END DO

  ! do i = 2,rows-1
  !     do j = 2,cols-1
  ! 		if ((mask(i,j) .ne. 0.0) .or. (mask(i+jj,j) .ne. 0.0) .or. (mask(i-jj,j-ii) .ne. 0.0)) then
  !     		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
  ! 		end if
  ! 	end do
  ! end do

  DO i = 2,rows-1
     DO j = 2,cols-1
        IF ((mask(i,j) .NE. 0.0) .OR. (mask(i-1,j) .NE. 0.0) .OR. (mask(i+1,j) .NE. 0.0) .OR. (mask(i,j-1) .NE. 0.0) .OR. (mask(i,j+1) .NE. 0.0) ) THEN
           partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
        END IF
     END DO
  END DO


  ! figure out which direction derivative goes (dx or dy)
  IF (dim .EQ. 1) THEN
     ! compute edges beforehand
     partial_edge(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
     partial_edge(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

     ! inner edges
     DO j=2,cols-1
        DO i=2,rows-1
           IF ((mask(i-1,j) .EQ. 12.5)) THEN
              partial_edge(i,j) =( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
           END IF
           IF ((mask(i-1,j) .EQ. 5.0) .AND. (mask(i-1,j) .EQ. 5.0)) THEN
              partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
           END IF

           IF ((mask(i+1,j) .EQ. 17.5)) THEN
              partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)
           END IF
           IF ((mask(i+1,j) .EQ. 10.0) .AND. (mask(i+1,j) .EQ. 10.0)) THEN
              partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d)
           END IF


           IF (mask(i,j-1) .EQ. 12.5) THEN
              partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
           END IF

           IF (mask(i,j-1) .EQ. 17.5) THEN
              partial_edge(i,j) =( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d)
           END IF


        END DO
     END DO


     DO j=1,cols
        DO i=2,rows-1

           IF (mask(i-1,j-1) .EQ. 12.5) THEN
              partial_edge(i,j) = partial_edge(i-1,j)
           END IF

           IF (mask(i+1,j-1) .EQ. 17.5) THEN
              partial_edge(i,j) = partial_edge(i+1,j)
           END IF

        END DO
     END DO


  END IF




  IF (dim .EQ. 2) THEN
     ! compute edges beforehand
     partial_edge(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
     partial_edge(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)


     ! inner edges
     DO j=2,cols
        DO i=1,rows
           IF ((mask(i,j-1).EQ.25.0) .OR. (mask(i,j-1).EQ.12.5)  .OR. (mask(i,j-1).EQ.17.5)) THEN
              partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
           END IF
        END DO
     END DO

     DO j=2,cols-1
        DO i=2,rows-1
           IF ((mask(i,j-1).EQ.50.0) .AND. (mask(i,j-1).EQ.50.0) .AND. (mask(i,j-1).EQ.50.0)) THEN
              partial_edge(i,j) =  ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)!
           END IF
        END DO
     END DO



     ! JANUARY DELETION...
     DO j=2,cols
        DO i=2,rows-1

           IF (mask(i-1,j) .EQ. 12.5) THEN
              partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
           END IF

           IF (mask(i+1,j) .EQ. 17.5) THEN
              partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
           END IF

        END DO
     END DO

     DO j=2,cols
        DO i=2,rows-1

           IF ((mask(i-1,j-1) .EQ. 12.5)) THEN
              partial_edge(i,j) = partial_edge(i,j-1)
           END IF

           IF ((mask(i+1,j-1) .EQ. 17.5)) THEN
              partial_edge(i,j) = partial_edge(i,j-1)
           END IF


        END DO
     END DO



  END IF





  RETURN
END FUNCTION partial_edge










! ----------------------------------------------------------------------------------%%
!
! PARTIAL_EDGE_P
!
! ----------------------------------------------------------------------------------%%


FUNCTION partial_edge_p(array,rows,cols,d1,d2,dim)

  USE globals
  USE initialize
  IMPLICIT NONE

  ! declare errthing
  INTEGER :: rows, cols, dim, i, j, ii, jj
  REAL(4) :: array(rows,cols), d1, d2, d
  REAL(4) :: partial_edge_p(rows,cols)

  ! write(*,*) "dim"
  ! write(*,*) rows
  ! write(*,*) cols

  partial_edge_p = 0.0

  ! figure out which direction derivative goes (dx or dy)
  IF (dim .EQ. 1) THEN
     ii = 1
     jj = 0
     d = d1

     DO i=2,rows-1
        DO j=1,cols
           IF (maskP(i,j) .NE. 0.0) THEN
              partial_edge_p(i,j) = (array(i+1,j) - array(i-1,j)) / (2.0*d)
           END IF
        END DO
     END DO

     DO i=2,rows-1
        DO j=1,cols
           IF ((maskP(i,j) .EQ. 5.0) .OR. (maskP(i,j) .EQ. 2.5)) THEN
              partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / (1.0*d)
           END IF
        END DO
     END DO

     DO i=1,rows-1
        DO j=1,cols
           ! 			if (maskP(i,j) .eq. 6.0) then
           ! 				partial_edge_p(i,j) = 0.0
           ! 			end if
           IF ((maskP(i,j) .EQ. 6.0) .OR. (maskP(i,j) .EQ. 6.5) .OR. (maskP(i,j) .EQ. 6.1) .OR. (maskP(i,j) .EQ. 6.05) .OR. (maskP(i,j) .EQ. 6.01)) THEN
              partial_edge_p(i,j) = (array(i+1,j) - array(i,j)) / d
           END IF
           IF ((maskP(i,j) .EQ. 3.0) .OR. (maskP(i,j) .EQ. 3.5) .OR. (maskP(i,j) .EQ. 3.1) .OR. (maskP(i,j) .EQ. 3.05) .OR. (maskP(i,j) .EQ. 3.01)) THEN
              partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / d
           END IF
        END DO
     END DO
     !
     ! fracture vertical velocity
     !partial_edge_p(xn-1,:) = (array(xn,:) - array(xn-1,:)) / param_f_dx


  END IF

  IF (dim .EQ. 2) THEN
     ii = 0
     jj = 1
     d = d2

     DO i=1,rows
        DO j=2,cols-1
           IF ((maskP(i,j) .NE. 0.0)) THEN
              partial_edge_p(i,j) = (array(i,j+1) - array(i,j-1)) / (2.0*d)
           END IF
        END DO
     END DO

     DO i=1,rows
        DO j=2,cols-1
           IF ((maskP(i,j) .EQ. 50.0) .OR. (maskP(i,j) .EQ. 3.5) .OR. (maskP(i,j) .EQ. 6.5)) THEN
              partial_edge_p(i,j) = (array(i,j) - array(i,j-1)) / (1.0*d)
           END IF
        END DO
     END DO



  END IF




  RETURN
END FUNCTION partial_edge_p




! ----------------------------------------------------------------------------------%%
!
! WRITE_VEC
!
! ----------------------------------------------------------------------------------%%

FUNCTION write_vec ( n, vector, filename )
  USE globals
  IMPLICIT NONE
  INTEGER :: n, j, output_status, unit0
  CHARACTER ( len = * ) filename
  REAL(4)  :: vector(n), write_vec



  unit0 = get_unit ()
  OPEN ( unit = unit0, file = filename, status = 'replace', iostat = output_status )
  IF ( output_status /= 0 ) THEN
     WRITE ( *, '(a,i8)' ) 'COULD NOT OPEN OUTPUT FILE "' // &
          TRIM ( filename ) // '" USING UNIT ', unit0
     unit0 = -1
     STOP
  END IF


  IF ( 0 < n ) THEN
     DO j = 1, n
        WRITE ( unit0, '(2x,g24.16)' ) vector(j)
     END DO

  END IF


  CLOSE ( unit = unit0 )
  write_vec = 1.0
  RETURN
END FUNCTION write_vec




! ----------------------------------------------------------------------------------%%
!
! WRITE_MATRIX
!
! ----------------------------------------------------------------------------------%%

FUNCTION write_matrix ( m, n, table, filename )
  USE globals
  IMPLICIT NONE
  INTEGER :: m, n, j, output_status, unit0, reclen
  CHARACTER ( len = * ) filename
  CHARACTER ( len = 30 ) string
  REAL(4)  :: table(m,n) , write_matrix



  INQUIRE(iolength=reclen)table
  unit0 = get_unit ()
  OPEN ( unit = unit0, file = filename, &
       status = 'replace', iostat = output_status, buffered='YES', buffercount=500)

  IF ( output_status /= 0 ) THEN
     WRITE ( *, '(a)' ) ' '
     WRITE ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
     WRITE ( *, '(a,i8)' ) 'Could not open the output file "' // &
          TRIM ( filename ) // '" on unit ', unit0
     unit0 = -1
     STOP
  END IF



  !
  WRITE ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
  ! 	!write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
  !
  !     do j = 1, n
  !       write ( unit0, string) table(1:m,j)
  !     end do

  DO j = 1, n
     WRITE ( unit0, 400) table(1:m,j)
  END DO
400 FORMAT(<m>g24.16)


  CLOSE ( unit = unit0 )
  write_matrix = 2.0
  RETURN
END FUNCTION write_matrix
