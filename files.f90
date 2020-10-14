MODULE files

  USE mpi, ONLY: mype
  USE parameters, ONLY: slab_mype

  IMPLICIT NONE

  !********************!                                                                                                                    
  !** Plotting stuff **!                                                                                                                                             
  !********************!                                                                          
                                                                         
  !unit no of files for output                                                                                                                            
  integer, parameter :: unit_slices  = 1
  integer, parameter :: unit_slicesv = 2
  integer, parameter :: unit_slices2 = 3
  integer, parameter :: unit_slices2v= 4
  integer, parameter :: unit_slices3 = 5
  integer, parameter :: unit_slices3v= 6
  integer, parameter :: unit_energy =50
!  integer, parameter :: unit_etrop  =51
!  integer, parameter :: unit_estra  =52
  integer, parameter :: unit_energyr=53
  integer, parameter :: unit_we=5441235
  integer, parameter :: unit_ce=5441234
  integer, parameter :: unit_conv=5451234
  integer, parameter :: unit_gamma=5461234

!  integer, parameter :: unit_hbot   =61
!  integer, parameter :: unit_hmid   =62
!  integer, parameter :: unit_htop   =63
  integer, parameter :: unit_ez     =64
  integer, parameter :: unit_wz     =65
  integer, parameter :: unit_wvave  =66

  integer, parameter :: unit_rco    =67
  integer, parameter :: unit_bs     =68
  integer, parameter :: unit_bs_dim =69
  integer, parameter :: unit_ell    =70

!  integer, parameter :: unit_hbotg  =71
!  integer, parameter :: unit_hmidg  =72
!  integer, parameter :: unit_htopg  =73

!  integer, parameter :: unit_ro_bot =76
!  integer, parameter :: unit_ro_mid =77
!  integer, parameter :: unit_ro_top =78

  integer, parameter :: unit_rotz   =81
  integer, parameter :: unit_pv     =82
  integer, parameter :: unit_pvz    =83
  integer, parameter :: unit_ens    =84
  integer, parameter :: unit_ensz   =85
  integer, parameter :: unit_specz  =86

  integer, parameter :: unit_h0     =100
  integer, parameter :: unit_h1     =101
  integer, parameter :: unit_h2     =102
  integer, parameter :: unit_h3     =103
  integer, parameter :: unit_h4     =104
  integer, parameter :: unit_h5     =105
  integer, parameter :: unit_h6     =106
  integer, parameter :: unit_h7     =107
  integer, parameter :: unit_h8     =108
  integer, parameter :: unit_h9     =109

  integer, parameter :: unit_h0w     =1001
  integer, parameter :: unit_h1w     =1011
  integer, parameter :: unit_h2w     =1012
  integer, parameter :: unit_h3w     =1013
  integer, parameter :: unit_h4w     =1014
  integer, parameter :: unit_h5w     =1015
  integer, parameter :: unit_h6w     =1016
  integer, parameter :: unit_h7w     =1017
  integer, parameter :: unit_h8w     =1018
  integer, parameter :: unit_h9w     =1019

  integer, parameter :: unit_hg0    =110
  integer, parameter :: unit_hg1    =111
  integer, parameter :: unit_hg2    =112
  integer, parameter :: unit_hg3    =113
  integer, parameter :: unit_hg4    =114
  integer, parameter :: unit_hg5    =115
  integer, parameter :: unit_hg6    =116
  integer, parameter :: unit_hg7    =117
  integer, parameter :: unit_hg8    =118
  integer, parameter :: unit_hg9    =119
  
  integer, parameter :: unit_ro0     =120
  integer, parameter :: unit_ro1     =121
  integer, parameter :: unit_ro2     =122
  integer, parameter :: unit_ro3     =123
  integer, parameter :: unit_ro4     =124
  integer, parameter :: unit_ro5     =125
  integer, parameter :: unit_ro6     =126
  integer, parameter :: unit_ro7     =127
  integer, parameter :: unit_ro8     =128
  integer, parameter :: unit_ro9     =129

  integer, parameter :: unit_frac    =130
  integer, parameter :: unit_cond1   =131
  integer, parameter :: unit_cond2   =132
  integer, parameter :: unit_cond3   =133
  integer, parameter :: unit_cond1j   =1311
  integer, parameter :: unit_cond2j   =1322
  integer, parameter :: unit_cond3j   =1333
  integer, parameter :: unit_condwz   =1344
  integer, parameter :: unit_grow    =134
  integer, parameter :: unit_hb      =135
  integer, parameter :: unit_hbt     =136
  integer, parameter :: unit_hb2    =1352
  integer, parameter :: unit_hbt2   =1362

  integer, parameter :: unit_psi     =137
  integer, parameter :: unit_psiz    =138


  integer, parameter :: unit_run    = 141

  integer, parameter :: unit_integrate    = 151

  integer, parameter :: unit_divtot =152
  integer, parameter :: unit_dwhere =153

  integer, parameter :: unit_slab = 160

  integer, parameter :: unit_tspec =161


  CONTAINS

    SUBROUTINE init_files

      !Energy
      if(mype==0) open (unit=unit_energy   ,file="energy.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_ez       ,file="ez.dat"       ,action="write",status="replace")
      if(mype==0) open (unit=unit_wz       ,file="wz.dat"       ,action="write",status="replace")
      if(mype==0) open (unit=unit_energyr  ,file="erot.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_we       ,file="we.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_ce       ,file="ce.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_conv     ,file="conv.dat"  ,action="write",status="replace")
      if(mype==0) open (unit=unit_gamma     ,file="gamma.dat"  ,action="write",status="replace")

!      if(mype==0) open (unit=unit_hbot     ,file="hbot.dat"     ,action="write",status="replace")
!      if(mype==0) open (unit=unit_hmid     ,file="hmid.dat"     ,action="write",status="replace")
!      if(mype==0) open (unit=unit_htop     ,file="htop.dat"     ,action="write",status="replace")
!      if(mype==0) open (unit=unit_ro_bot   ,file="ro_bot.dat"   ,action="write",status="replace")
!      if(mype==0) open (unit=unit_ro_mid   ,file="ro_mid.dat"   ,action="write",status="replace")
!      if(mype==0) open (unit=unit_ro_top   ,file="ro_top.dat"   ,action="write",status="replace")
!      if(mype==0) open (unit=unit_hbotg    ,file="hbotg.dat"    ,action="write",status="replace")
!      if(mype==0) open (unit=unit_hmidg    ,file="hmidg.dat"    ,action="write",status="replace")
!      if(mype==0) open (unit=unit_htopg    ,file="htopg.dat"    ,action="write",status="replace")

      


      !Base-state structure
      if(mype==0) open (unit=unit_rco      ,file="rcoeff.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_bs       ,file="bs.dat"       ,action="write",status="replace")
      if(mype==0) open (unit=unit_bs_dim   ,file="bs_dim.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_ell      ,file="ell.dat"      ,action="write",status="replace")

      if(mype==0) open (unit=unit_specz    ,file="specz.dat"    ,action="write",status="replace")

      if(mype==0) open (unit=unit_rotz     ,file="rotz.dat"     ,action="write",status="replace")  !Rotational energy as a function of z... 
      if(mype==0) open (unit=unit_pv       ,file="pv.dat"       ,action="write",status="replace")  !Unused
      if(mype==0) open (unit=unit_pvz      ,file="pvz.dat"      ,action="write",status="replace")  !Unused
      if(mype==0) open (unit=unit_ens      ,file="ens.dat"      ,action="write",status="replace")
      if(mype==0) open (unit=unit_ensz     ,file="ensz.dat"     ,action="write",status="replace")

      !Energy spectrum at various levels
      open (unit=unit_h0       ,file="h0.dat"        ,action="write",status="replace")
      open (unit=unit_h1       ,file="h1.dat"        ,action="write",status="replace")
      open (unit=unit_h2       ,file="h2.dat"        ,action="write",status="replace")
      open (unit=unit_h3       ,file="h3.dat"        ,action="write",status="replace")
      open (unit=unit_h4       ,file="h4.dat"        ,action="write",status="replace")
      open (unit=unit_h5       ,file="h5.dat"        ,action="write",status="replace")
      open (unit=unit_h6       ,file="h6.dat"        ,action="write",status="replace")
      open (unit=unit_h7       ,file="h7.dat"        ,action="write",status="replace")
      open (unit=unit_h8       ,file="h8.dat"        ,action="write",status="replace")
      open (unit=unit_h9       ,file="h9.dat"        ,action="write",status="replace")

      open (unit=unit_h0w       ,file="h0w.dat"        ,action="write",status="replace")
      open (unit=unit_h1w       ,file="h1w.dat"        ,action="write",status="replace")
      open (unit=unit_h2w       ,file="h2w.dat"        ,action="write",status="replace")
      open (unit=unit_h3w       ,file="h3w.dat"        ,action="write",status="replace")
      open (unit=unit_h4w       ,file="h4w.dat"        ,action="write",status="replace")
      open (unit=unit_h5w       ,file="h5w.dat"        ,action="write",status="replace")
      open (unit=unit_h6w       ,file="h6w.dat"        ,action="write",status="replace")
      open (unit=unit_h7w       ,file="h7w.dat"        ,action="write",status="replace")
      open (unit=unit_h8w       ,file="h8w.dat"        ,action="write",status="replace")
      open (unit=unit_h9w       ,file="h9w.dat"        ,action="write",status="replace")

      open (unit=unit_tspec    ,file="tspec.dat"     ,action="write",status="replace")

      !Geostrophic energy spectrum at various levels
      open (unit=unit_hg0      ,file="hg0.dat"      ,action="write",status="replace")
      open (unit=unit_hg1      ,file="hg1.dat"      ,action="write",status="replace")
      open (unit=unit_hg2      ,file="hg2.dat"      ,action="write",status="replace")
      open (unit=unit_hg3      ,file="hg3.dat"      ,action="write",status="replace")
      open (unit=unit_hg4      ,file="hg4.dat"      ,action="write",status="replace")
      open (unit=unit_hg5      ,file="hg5.dat"      ,action="write",status="replace")
      open (unit=unit_hg6      ,file="hg6.dat"      ,action="write",status="replace")
      open (unit=unit_hg7      ,file="hg7.dat"      ,action="write",status="replace")
      open (unit=unit_hg8      ,file="hg8.dat"      ,action="write",status="replace")
      open (unit=unit_hg9      ,file="hg9.dat"      ,action="write",status="replace")

      !Rossby number at various levels
      open (unit=unit_ro0      ,file="ro0.dat"      ,action="write",status="replace")
      open (unit=unit_ro1      ,file="ro1.dat"      ,action="write",status="replace")
      open (unit=unit_ro2      ,file="ro2.dat"      ,action="write",status="replace")
      open (unit=unit_ro3      ,file="ro3.dat"      ,action="write",status="replace")
      open (unit=unit_ro4      ,file="ro4.dat"      ,action="write",status="replace")
      open (unit=unit_ro5      ,file="ro5.dat"      ,action="write",status="replace")
      open (unit=unit_ro6      ,file="ro6.dat"      ,action="write",status="replace")
      open (unit=unit_ro7      ,file="ro7.dat"      ,action="write",status="replace")
      open (unit=unit_ro8      ,file="ro8.dat"      ,action="write",status="replace")
      open (unit=unit_ro9      ,file="ro9.dat"      ,action="write",status="replace")

      !Conditions of integrability
      if(mype==0) open (unit=unit_frac     ,file="frac.dat"     ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond1    ,file="cond1.dat"    ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond2    ,file="cond2.dat"    ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond3    ,file="cond3.dat"    ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond1j   ,file="cond1j.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond2j   ,file="cond2j.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_cond3j   ,file="cond3j.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_condwz   ,file="condwz.dat"   ,action="write",status="replace")
      if(mype==0) open (unit=unit_grow     ,file="growth.dat"   ,action="write",status="replace")

      !Vertical scale of buoyancy and vertical velocity
      if(mype==0) open (unit=unit_hb       ,file="hb.dat"       ,action="write",status="replace")
      if(mype==0) open (unit=unit_hbt      ,file="hbt.dat"      ,action="write",status="replace")
      if(mype==0) open (unit=unit_hb2      ,file="hb2.dat"       ,action="write",status="replace")
      if(mype==0) open (unit=unit_hbt2     ,file="hbt2.dat"      ,action="write",status="replace")

      !Psi illustration
      if(mype==0) open (unit=unit_psi      ,file="psi.dat"      ,action="write",status="replace")
      if(mype==0) open (unit=unit_psiz     ,file="psiz.dat"     ,action="write",status="replace")

      !Specifications and parameters of the run
      if(mype==0) open (unit=unit_run    ,file="run_specs.dat"    ,action="write",status="replace")

      if(mype==0) open (unit=unit_integrate ,file="sum.dat"    ,action="write",status="replace")
      if(mype==0) open (unit=unit_divtot ,file="div.dat"    ,action="write",status="replace")
      if(mype==0) open (unit=unit_dwhere ,file="dwhere.dat"    ,action="write",status="replace")

      !Slab for slow/fast spectrum
      if(mype==slab_mype) open (unit=unit_slab  ,file="slab.dat"   ,action="write",status="replace")

    END SUBROUTINE init_files


    SUBROUTINE kill_files

      if(mype==0) close(unit=unit_energy)
      if(mype==0) close(unit=unit_energyr)
      if(mype==0) close(unit=unit_we)

!      if(mype==0) close(unit=unit_hbot)
!      if(mype==0) close(unit=unit_hmid)
!      if(mype==0) close(unit=unit_htop)
!      if(mype==0) close(unit=unit_ro_bot)
!      if(mype==0) close(unit=unit_ro_mid)
!      if(mype==0) close(unit=unit_ro_top)
      if(mype==0) close(unit=unit_ez)
!      if(mype==0) close(unit=unit_hbotg)
!      if(mype==0) close(unit=unit_hmidg)
!      if(mype==0) close(unit=unit_htopg)
      if(mype==0) close(unit=unit_rco)
      if(mype==0) close(unit=unit_bs)
      if(mype==0) close(unit=unit_bs_dim)
      if(mype==0) close(unit=unit_ell)
      if(mype==0) close(unit=unit_rotz)
      if(mype==0) close(unit=unit_pv)
      if(mype==0) close(unit=unit_pvz)
      if(mype==0) close(unit=unit_ens)
      if(mype==0) close(unit=unit_ensz)
      if(mype==0) close(unit=unit_specz)

     if(mype==0) close(unit=unit_tspec)

     if(mype==0) close(unit=unit_h0)
     if(mype==0) close(unit=unit_h1)
     if(mype==0) close(unit=unit_h2)
     if(mype==0) close(unit=unit_h3)
     if(mype==0) close(unit=unit_h4)
     if(mype==0) close(unit=unit_h5)
     if(mype==0) close(unit=unit_h6)
     if(mype==0) close(unit=unit_h7)
     if(mype==0) close(unit=unit_h8)
     if(mype==0) close(unit=unit_h9)
     
     close(unit=unit_hg0)
     close(unit=unit_hg1)
     close(unit=unit_hg2)
     close(unit=unit_hg3)
     close(unit=unit_hg4)
     close(unit=unit_hg5)
     close(unit=unit_hg6)
     close(unit=unit_hg7)
     close(unit=unit_hg8)
     close(unit=unit_hg9)

     close(unit=unit_ro0)
     close(unit=unit_ro1)
     close(unit=unit_ro2)
     close(unit=unit_ro3)
     close(unit=unit_ro4)
     close(unit=unit_ro5)
     close(unit=unit_ro6)
     close(unit=unit_ro7)
     close(unit=unit_ro8)
     close(unit=unit_ro9)

     if(mype==0) close(unit=unit_frac)
     if(mype==0) close(unit=unit_cond1)
     if(mype==0) close(unit=unit_cond2)
     if(mype==0) close(unit=unit_cond3)
     if(mype==0) close(unit=unit_condwz)
     if(mype==0) close(unit=unit_cond1j)
     if(mype==0) close(unit=unit_cond2j)
     if(mype==0) close(unit=unit_cond3j)
     if(mype==0) close(unit=unit_grow)
     if(mype==0) close(unit=unit_hb)
     if(mype==0) close(unit=unit_hbt)
     if(mype==0) close(unit=unit_hb2)
     if(mype==0) close(unit=unit_hbt2)

     if(mype==0) close(unit=unit_psi)
     if(mype==0) close(unit=unit_psiz)

     if(mype==0) close(unit=unit_run)
     if(mype==0) close(unit=unit_integrate)
     if(mype==0) close(unit=unit_divtot)
     if(mype==0) close(unit=unit_dwhere)

     if(mype==slab_mype) close(unit=unit_slab)




    END SUBROUTINE kill_files


END MODULE files

