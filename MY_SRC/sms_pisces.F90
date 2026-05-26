MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables are declared and allocated
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   CHARACTER(:), ALLOCATABLE ::   numnatp_ref   !! Character buffer for reference namelist pisces
   CHARACTER(:), ALLOCATABLE ::   numnatp_cfg   !! Character buffer for configuration namelist pisces
   INTEGER ::   numonp      = -1                !! Logical unit for namelist pisces output

   !!* Model used
   LOGICAL  ::  ln_p2z            !: Flag to use LOBSTER model
   LOGICAL  ::  ln_p4z            !: Flag to use PISCES  model
   LOGICAL  ::  ln_p5z            !: Flag to use PISCES  quota model
   LOGICAL  ::  ln_ligand         !: Flag to enable organic ligands
   LOGICAL  ::  ln_sediment       !: Flag to enable sediment module
   LOGICAL  ::  ln_bait           !: Flag to use BAIT Fe module 
   LOGICAL  ::  ln_felith         !: Flag to use BAIT lithogenic Fe module
   LOGICAL  ::  ln_feauth         !: Flag to use BAIT authigenic Fe module
   LOGICAL  ::  ln_p6z            !: Flag to use PISCES  quota explicit diazotrophy
   LOGICAL  ::  ln_tricho         !: switch explicit diazotroph between Trichodesmium (TRUE) 
                                  !: and Crocosphaera (FALSE) parameterisation
   LOGICAL  ::  ln_facul          !: Flag to enable facultative diazotrophy
   LOGICAL  ::  ln_evolve         !: Flag to allow dynamic evolution of croco

   !!*  Time variables
   INTEGER  ::   nrdttrc          !: ???
   REAL(wp) ::   rfact , rfactr   !: time step duration (in seconds)
   REAL(wp) ::   rfact2, rfact2r  !: time step duration (in seconds) when timesplitting is activated for PISCES
   REAL(wp) ::   xstep            !: Time step duration for biology
   REAL(wp) ::   ryyss            !: number of seconds per year 
   REAL(wp) ::   r1_ryyss         !: inverse number of seconds per year 

   !!*  Biological parameters 
   REAL(wp) ::   rno3             !: C/N stoichiometric ratio
   REAL(wp) ::   o2ut             !: O2/N stoichiometric ratio for ammonification
   REAL(wp) ::   po4r             !: C/P stoichiometric ratio
   REAL(wp) ::   rdenit           !: C/N ratio for denitrification
   REAL(wp) ::   rdenita          !: C/N ratio for denitrification
   REAL(wp) ::   o2nit            !: O2/N ratio for nitrification
   REAL(wp) ::   wsbio, wsbio2    !: Sinking speeds of particles
   REAL(wp) ::   wsbio2max        !: Maximum sinking speed of the largest particles
   REAL(wp) ::   wsbio2scale      !: Length scale for the variations of wsbio2
   REAL(wp) ::   oxymin           !:  half saturation constant for anoxia
   REAL(wp) ::   xkmort           !: Mortality half-saturation constant
   REAL(wp) ::   feratz           !: Fe/C in microzooplankton
   REAL(wp) ::   feratm           !: Fe/C in mesozooplankton
   REAL(wp) ::   ldocp            !: Ligand production ratio during PP
   REAL(wp) ::   ldocz            !: Ligand production ratio by grazing
   REAL(wp) ::   lthet            !: Uptake of ligand by phytoplankton
   REAL(wp) ::   no3rat3          !: C/N ratio of zooplankton
   REAL(wp) ::   po4rat3          !: C/P ratio of zooplankton
   REAL(wp) ::   ragglfe          !: Aggregation rate of LFe
   REAL(wp) ::   aggrat           !: Aggregation scaling
   REAL(wp) ::   rdislfea         !: Diaggregation rate of LFea
   REAL(wp) ::   wslfe            !: lfe sinking rate
   REAL(wp) ::   wslfea           !: lfea sinking rate
   REAL(wp) ::   wslfeamax        !: lfea max sinking rate
   REAL(wp) ::   cfeagg        !: aggregation of small authigenic Fe
   REAL(wp) ::   cfeagg2        !: aggregation of small authigenic Fe mediated by big afe
   REAL(wp) ::   wsafes
   REAL(wp) ::   wsafeb
   REAL(wp) ::   wsafebmax

   !!*  diagnostic parameters 
   REAL(wp) ::  tpp               !: total primary production
   REAL(wp) ::  t_oce_co2_exp     !: total carbon export
   REAL(wp) ::  t_oce_co2_flx     !: Total ocean carbon flux
   REAL(wp) ::  t_oce_co2_flx_cum !: Cumulative Total ocean carbon flux
   REAL(wp) ::  t_atm_co2_flx     !: global mean of atmospheric pco2

   !!* Mass conservation
   LOGICAL  ::  ln_check_mass     !: Flag to check mass conservation
   LOGICAL, PUBLIC ::   ln_ironice   !: boolean for Fe input from sea ice

   !!* Diurnal cycle in PISCES
   LOGICAL  ::  ln_p4z_dcyc       !: Flag to activate diurnal cycle in PISCES

   !!*  Biological fluxes for light : variables shared by pisces & lobster
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::  strn  !: Day duration in hours
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enanom, ediatm !: mean PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epico          !: PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epicom         !: mean PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoy, etotm    !: averaged PAR in the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  Half-saturation con,stant for diatoms
   ! Explicit Diazotroph
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  ediaz
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  ediazm  !: PAR for diazotrophs
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  vnfmax_dz   !: Maximum nutrient uptake rate diazotroph

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    ::   xksimax    !: Maximum half-saturation constant over the year (Si)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   biron      !: bioavailable fraction of iron
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   plig       !: proportion of iron organically complexed

   !!*  Sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio5   !: POC sinking
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio6   !: GOC sinking
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio7   !: POC sinking
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio8   !: GOC sinking

   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: Fraction of nanophytoplankton that are calcifying organisms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: OMZ 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac2   !: N depleted indice
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   orem       !: oxic remineralisation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: Shear rate
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodpoc    !: POC production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   conspoc    !: POC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodgoc    !: GOC production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consgoc    !: GOC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consfe3    !: GOC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   blim       !: bacterial production factor
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizen, logsizen    !: size of nanophyto
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizep, logsizep    !: size of picophyto
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sized, logsized    !: size of diatoms   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizena     !: size of nanophytoplankton, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizepa     !: size of picophyto, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizeda     !: size of diatomss, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfecolagg  !: Refractory diagnostic concentration of ligands
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xcoagfe    !: Coagulation rate of colloidal Fe/ligands
   ! Explicit Diazotroph
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizedz, logsizedz     !: size of diazotrophs
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizedza    !: size of diazotrophs, after

   !!* Evolutionary clock for diazos
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   evoclock !: Forwards Evolutionary clock
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   bckevoclock !: Backwards Evolutionary clock

   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: Carbonate chemistry constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: Carbonate chemistry constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: Solubility product of CaCO3
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: Proton concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: CO3 saturation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: 


   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   remintpoc  ! Mean remineralisation rate of POC
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   remintgoc  ! Mean remineralisation rate of GOC   

   LOGICAL, SAVE :: lk_sed

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: sms_pisces.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      INTEGER ::   ierr(18)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for pisces & lobster
      ALLOCATE( xksi(jpi,jpj), strn(jpi,jpj),  STAT=ierr(1) )

      IF( ln_p4z .OR. ln_p5z .OR. ln_p6z ) THEN

         !* Optics
         ALLOCATE(  enano(jpi,jpj,jpk) , ediat(jpi,jpj,jpk) ,   &
           &        enanom(jpi,jpj,jpk), ediatm(jpi,jpj,jpk),   &
           &        emoy(jpi,jpj,jpk)  , etotm(jpi,jpj,jpk),   STAT=ierr(2) )

         !* Biological SMS
         ALLOCATE( xksimax(jpi,jpj)  , biron(jpi,jpj,jpk)      ,  STAT=ierr(3) )

         ! Biological SMS
         ALLOCATE( xfracal  (jpi,jpj,jpk), orem    (jpi,jpj,jpk),  &
            &      nitrfac  (jpi,jpj,jpk), nitrfac2(jpi,jpj,jpk),  &
            &      prodcal  (jpi,jpj,jpk), xdiss   (jpi,jpj,jpk),  &
            &      prodpoc  (jpi,jpj,jpk), conspoc (jpi,jpj,jpk),  &
            &      prodgoc  (jpi,jpj,jpk), consgoc (jpi,jpj,jpk),  &
            &      blim     (jpi,jpj,jpk), consfe3 (jpi,jpj,jpk),  &
            &      xfecolagg(jpi,jpj,jpk), xcoagfe (jpi,jpj,jpk), STAT=ierr(4) )

         !* Carbonate chemistry
         ALLOCATE( ak13  (jpi,jpj,jpk)  ,                          &
            &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,    &
            &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,    &
            &      aphscale(jpi,jpj,jpk),                         STAT=ierr(5) )
         !
         !* Temperature dependency of SMS terms
         ALLOCATE( tgfunc (jpi,jpj,jpk) , tgfunc2(jpi,jpj,jpk),   STAT=ierr(6) )
         !
         !* Sinking speed
         ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4 (jpi,jpj,jpk),   STAT=ierr(7) )

         !*  Size of phytoplankton cells
         ALLOCATE( sizen (jpi,jpj,jpk), sized (jpi,jpj,jpk),      &
           &       sizena(jpi,jpj,jpk), sizeda(jpi,jpj,jpk),      &
           &       logsizen(jpi,jpj,jpk), logsized(jpi,jpj,jpk), STAT=ierr(8) )
         ! 
         ALLOCATE( plig(jpi,jpj,jpk)  ,                           STAT=ierr(9) )
         !
         ALLOCATE( remintpoc(jpi,jpj,jpk), remintgoc(jpi,jpj,jpk), STAT=ierr(10) )
         !      
      ENDIF
      !
      IF( ln_p5z .OR. ln_p6z) THEN
         ! PISCES-QUOTA specific part      
         ALLOCATE( epico(jpi,jpj,jpk)   , epicom(jpi,jpj,jpk) ,   STAT=ierr(11) ) 

         !*  Size of phytoplankton cells
         ALLOCATE( sizep(jpi,jpj,jpk), sizepa(jpi,jpj,jpk), logsizep(jpi,jpj,jpk), STAT=ierr(12) )
      ENDIF
      IF( ln_felith ) THEN
         ALLOCATE( wsbio5 (jpi,jpj,jpk) , wsbio6 (jpi,jpj,jpk),   STAT=ierr(13) )
      ENDIF
      IF( ln_feauth ) THEN
         ALLOCATE( wsbio7 (jpi,jpj,jpk) , wsbio8 (jpi,jpj,jpk),   STAT=ierr(14) )
      ENDIF
      IF( ln_p6z ) THEN
         ! PISCES QUOTA Explicit Diazotroph
         ALLOCATE(ediaz(jpi,jpj,jpk)   , ediazm(jpi,jpj,jpk),  STAT=ierr(15) )
         !* size of diazotroph cells
         ALLOCATE(  sizedz(jpi,jpj,jpk), sizedza(jpi,jpj,jpk), logsizedz(jpi,jpj,jpk), STAT=ierr(16) )
         !* Maximum nutrient uptake by diazotroph
         ALLOCATE( vnfmax_dz(jpi,jpj,jpk),  STAT=ierr(17) )
         !* Evolution Clock for Croco
         ALLOCATE( evoclock(jpi,jpj,jpk), bckevoclock(jpi,jpj,jpk),   STAT=ierr(18) )
      ENDIF
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sms_pisces_alloc: failed to allocate arrays' ) 
      !
   END FUNCTION sms_pisces_alloc

   !!======================================================================   
END MODULE sms_pisces    
