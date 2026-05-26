MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :   PISCES Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p4z_fechem       : Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  : Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc : Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zche          ! chemical model
   USE p4zbc           ! Boundary conditions from sediments
   USE prtctl          ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem        ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init   ! called in trcsms_pisces.F90

   LOGICAL          ::   ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::   xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::   xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::   ligand       !: ligand concentration in the ocean 
   REAL(wp), PUBLIC ::   kfep         !: rate constant for nanoparticle formation
   REAL(wp), PUBLIC ::   scaveff      !: Fraction of scavenged iron that is considered as being subject to solubilization
   LOGICAL          ::   ln_fixlogk   !: boolean for fixed logK for ligands
   REAL(wp), PUBLIC ::   logk        !: fixed logk
   REAL(wp), PUBLIC ::   xlamdust1     !: scavenging rate of Iron by lfe
   REAL(wp), PUBLIC ::   xlamdust2     !: scavenging rate of Iron by lfe aggregates
   REAL(wp), PUBLIC ::   xlamafe1     !: scavenging rate of Iron by lfe
   REAL(wp), PUBLIC ::   xlamafe2     !: scavenging rate of Iron by lfe
   REAL(wp), PUBLIC ::   mincolfe     !:
   REAL(wp), PUBLIC ::   collf        !:
   REAL(wp), PUBLIC ::   kcfe
   REAL(wp), PUBLIC ::   kbcfe        !:
   REAL(wp), PUBLIC ::   coagf        !:
   REAL(wp), PUBLIC ::   siscav      !:
   REAL(wp), PUBLIC ::   calscav     !:
   REAL(wp), PUBLIC ::   scaveff2    !: Fraction of scavenged iron that is considered as being subject to solubilizatio from si and cal

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zfechem.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   A simple chemistry model of iron from Aumont and Bopp (2006)
      !!                based on one ligand and one inorganic form
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   ji, jj, jk, jic, jn
      REAL(wp) ::   zlam1a, zlam1b
      REAL(wp) ::   zkeq, zfesatur, zligco !fe3sol, zligco
      REAL(wp) ::   zscave, zaggdfea, zaggdfeb, ztrc, zdust, zklight
      REAL(wp) ::   ztfe, zhplus, zxlam, zaggliga, zaggligb
      REAL(wp) ::   zprecip, zprecipno3,  zconsfe, za1
      REAL(wp) ::   zrfact2
      REAL(wp) ::   zscaveb, zscave3, zscave4, zscave1, zscave2, zdust2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zcfe
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zTL1, zFe3, ztotlig, zfeprecip, zFeL1, zfecoll
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zcoll3d, zscav3d, zlcoll3d
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zscava3d, zscavl3d, fe3sol
      REAL(wp) :: zaggdfec, zaggdfed, zlam1c, zlam1d, zaggdfee, tbio, biof
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_fechem')
      !
      zFe3     (:,:,jpk) = 0.
      zFeL1    (:,:,jpk) = 0.
      zTL1     (:,:,jpk) = 0.
      zfeprecip(:,:,jpk) = 0.
      zcoll3d  (:,:,jpk) = 0.
      zscav3d  (:,:,jpk) = 0.
      zlcoll3d (:,:,jpk) = 0.
      zfecoll  (:,:,jpk) = 0.
      xfecolagg(:,:,jpk) = 0.
      xcoagfe  (:,:,jpk) = 0.
      !
      ! Derive colloidal Fe
      ! -------------------------------------------------
      IF ( ln_bait ) THEN
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol(ji,jj,jk)  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
               &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
               &         + fesol(ji,jj,jk,5) / zhplus )
               zcfe(ji,jj,jk) = MIN( MAX( ( tr(ji,jj,jk,jpfer,Kbb) - fe3sol(ji,jj,jk) ) , mincolfe  &
               &         * tr(ji,jj,jk,jpfer,Kbb) ) , ( 0.9 * tr(ji,jj,jk,jpfer,Kbb) ) )
      END_3D
      ELSE
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
               zcfe(ji,jj,jk) = 0.5 * tr(ji,jj,jk,jpfer,Kbb)
      END_3D
      ENDIF
      !
      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Pham and Ito (2018)
      ! -------------------------------------------------
      xfecolagg(:,:,:) = ligand * 1E9 + MAX(0., chemo2(:,:,:) - tr(:,:,:,jpoxy,Kbb) ) / 400.E-6
      !
      IF( ln_ligand ) THEN  
         ztotlig(:,:,:) = tr(:,:,:,jplgw,Kbb) * 1E9
      ELSE 
         IF( ln_ligvar ) THEN
            ztotlig(:,:,:) =  0.09 * 0.667 * tr(:,:,:,jpdoc,Kbb) * 1E6 + xfecolagg(:,:,:)
            ztotlig(:,:,:) =  MIN( ztotlig(:,:,:), 10. )
         ELSE
            ztotlig(:,:,:) = ligand * 1E9 
         ENDIF
      ENDIF

      ! ------------------------------------------------------------
      !  from Aumont and Bopp (2006)
      ! This model is based on one ligand, Fe2+ and Fe3+ 
      ! Chemistry is supposed to be fast enough to be at equilibrium
      ! ------------------------------------------------------------
        IF( ln_bait ) THEN
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
          zTL1(ji,jj,jk)  = ztotlig(ji,jj,jk)
          zkeq            = fekeq(ji,jj,jk)
          IF( ln_fixlogk )  fekeq(ji,jj,jk) = 10**logk
          zklight         = 4.77E-7 * etot(ji,jj,jk) * 0.5 / ( 10**(-6.3) )
          zconsfe         = consfe3(ji,jj,jk) / ( 10**(-6.3) )
          zfesatur        = zTL1(ji,jj,jk) * 1E-9
          ztfe            = max(0., tr(ji,jj,jk,jpfer,Kbb) - zcfe(ji,jj,jk) )
         ! Fe' is the root of a 2nd order polynom
          za1 =  1. + zfesatur * zkeq + zklight +  zconsfe - zkeq * tr(ji,jj,jk,jpfer,Kbb)
          zFe3 (ji,jj,jk) = ( -1 * za1 + SQRT( za1**2 + 4. * ztfe * zkeq) ) / ( 2. * zkeq + rtrn )
      END_3D
        ELSE
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
          zTL1(ji,jj,jk)  = ztotlig(ji,jj,jk)
          zkeq            = fekeq(ji,jj,jk)
          zklight         = 4.77E-7 * etot(ji,jj,jk) * 0.5 / ( 10**(-6.3) )
          zconsfe         = consfe3(ji,jj,jk) / ( 10**(-6.3) )
          zfesatur        = zTL1(ji,jj,jk) * 1E-9
          ztfe            = (1.0 + zklight) * tr(ji,jj,jk,jpfer,Kbb) 
          ! Fe' is the root of a 2nd order polynom
          za1 =  1. + zfesatur * zkeq + zklight +  zconsfe - zkeq * tr(ji,jj,jk,jpfer,Kbb)
          zFe3 (ji,jj,jk) = ( -1 * za1 + SQRT( za1**2 + 4. * ztfe * zkeq) ) / ( 2. * zkeq + rtrn )
          zFeL1(ji,jj,jk) = MAX( 0., tr(ji,jj,jk,jpfer,Kbb) - zFe3(ji,jj,jk) )
      END_3D
        ENDIF
      !
        IF( ln_bait ) THEN
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
         zFeL1(ji,jj,jk) = MAX( 0., tr(ji,jj,jk,jpfer,Kbb) - zcfe(ji,jj,jk)  - zFe3(ji,jj,jk) )
      END_3D
        ELSE
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
         zFeL1(ji,jj,jk) = MAX( 0., tr(ji,jj,jk,jpfer,Kbb) - zFe3(ji,jj,jk) )
      END_3D
        ENDIF
      !
      plig(:,:,:) =  MAX( 0., ( zFeL1(:,:,:) / ( tr(:,:,:,jpfer,Kbb) + rtrn ) ) )
      !
      zdust = 0.         ! if no dust available

      ! Computation of the colloidal fraction that is subjecto to coagulation
      ! The assumption is that 50% of complexed iron is colloidal. Furthermore
      ! The refractory part is supposed to be non sticky. The refractory
      ! fraction is supposed to equal to the background concentration + 
      ! the fraction that accumulates in the deep ocean. AOU is taken as a 
      ! proxy of that accumulation following numerous studies showing 
      ! some relationship between weak ligands and AOU.
      ! An issue with that parameterization is that when ligands are not
      ! prognostic or non variable, all the colloidal fraction is supposed
      ! to coagulate
      ! ----------------------------------------------------------------------
      IF (ln_bait) THEN
      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
           zfecoll(ji,jj,jk) = zcfe(ji,jj,jk)
          END_3D
      ELSE
      IF( ln_ligand ) THEN
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
            zfecoll(ji,jj,jk) = 0.5 * zFeL1(ji,jj,jk) * MAX(0., ztotlig(ji,jj,jk) - xfecolagg(ji,jj,jk) ) &
                  &              / ( ztotlig(ji,jj,jk) + rtrn )
         END_3D
      ELSE
         IF( ln_ligvar ) THEN
            DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
               zfecoll(ji,jj,jk) = 0.5 * zFeL1(ji,jj,jk) * MAX(0., ztotlig(ji,jj,jk) - xfecolagg(ji,jj,jk) ) &
                  &              / ( ztotlig(ji,jj,jk) + rtrn )
            END_3D
         ELSE
            DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
               zfecoll(ji,jj,jk) = 0.5 * zFeL1(ji,jj,jk) * MAX(0., 0.09 * 0.667 * tr(ji,jj,jk,jpdoc,Kbb) * 1E6 ) &
                  &             / ( ztotlig(ji,jj,jk) + rtrn )
            END_3D
         ENDIF
      ENDIF
      ENDIF ! ln_bait

      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1)
         ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
         ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
         ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
         ! --------------------------------------------------------------------------------------
         zhplus  = max( rtrn, hi(ji,jj,jk) )
         fe3sol(ji,jj,jk)  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
         &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
         &         + fesol(ji,jj,jk,5) / zhplus )
         !
         ! precipitation of Fe3+, creation of nanoparticles
         zprecip = MAX( 0., ( zFe3(ji,jj,jk) - fe3sol(ji,jj,jk) ) ) * kfep * xstep * ( 1.0 - nitrfac(ji,jj,jk) ) 
         ! Precipitation of Fe2+ due to oxidation by NO3 (Croot et al., 2019)
         ! This occurs in anoxic waters only
         zprecipno3 = 2.0 * 130.0 * tr(ji,jj,jk,jpno3,Kbb) * nitrfac(ji,jj,jk) * xstep * zFe3(ji,jj,jk)
         !
         zfeprecip(ji,jj,jk) = zprecip + zprecipno3
         !!
!         ztrc   = ( tr(ji,jj,jk,jppoc,Kbb) + tr(ji,jj,jk,jpgoc,Kbb) + tr(ji,jj,jk,jpcal,Kbb) + tr(ji,jj,jk,jpgsi,Kbb) ) * 1.e6 
         ztrc   = ( tr(ji,jj,jk,jppoc,Kbb) + tr(ji,jj,jk,jpgoc,Kbb) + tr(ji,jj,jk,jpcal,Kbb) + tr(ji,jj,jk,jpgsi,Kbb) )
         IF ( ln_bait ) THEN
!            ztrc = ( tr(ji,jj,jk,jppoc,Kbb) + tr(ji,jj,jk,jpgoc,Kbb) &
!         &   + (calscav * tr(ji,jj,jk,jpcal,Kbb)) + (siscav * tr(ji,jj,jk,jpgsi,Kbb) ) ) * 1.e6
            ztrc = ( tr(ji,jj,jk,jppoc,Kbb) + tr(ji,jj,jk,jpgoc,Kbb) &
         &   + (calscav * tr(ji,jj,jk,jpcal,Kbb)) + (siscav * tr(ji,jj,jk,jpgsi,Kbb) ) )
         ENDIF 
         ztrc = MAX( rtrn, ztrc )
         !
         !  Compute the coagulation of colloidal iron. This parameterization
         !  could be thought as an equivalent of colloidal pumping.
         !  It requires certainly some more work as it is very poorly
         !  constrained.
         !  ----------------------------------------------------------------
         IF (ln_bait ) THEN
               ! amount of living biomass
                IF( ln_p5z ) THEN
                 tbio = tr(ji,jj,jk,jpphy,Kbb) + tr(ji,jj,jk,jpdia,Kbb) + tr(ji,jj,jk,jppic,Kbb) + rtrn
                ELSE
                 tbio = tr(ji,jj,jk,jpphy,Kbb) + tr(ji,jj,jk,jpdia,Kbb) + rtrn
                ENDIF
               ! coagulation with DOC incr. by coagf so background biof should
               ! be = 1/coagf
               biof = max( 1/coagf , tbio / ( tbio + kbcfe ) )
               !! colloidal shunt associated with small particles and DOC
               zlam1a   = ( ((12.0*coagf)*biof)  * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) &
                   &   + 9.05  * tr(ji,jj,jk,jppoc,Kbb) ) * xdiss(ji,jj,jk)    &
                   &    + ( 2.49  * tr(ji,jj,jk,jppoc,Kbb) )     &
                   &    + ( ((127.8*coagf)*biof) * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) &
                   &    + 725.7 * tr(ji,jj,jk,jppoc,Kbb) )
               zaggdfea = zlam1a * xstep * zcfe(ji,jj,jk)
               ! autocatalytic removal of cFe produces authigenic Fe, lower in
               ! the light due to photochem dissolution of small authigenic
               ! particles
               zaggdfee = ( (kfep*collf) * ( zcfe(ji,jj,jk)**4 / ( zcfe(ji,jj,jk)**4 + kcfe**4 ) ) ) &
               &          * (1 - (etot(ji,jj,jk)**2 / ( etot(ji,jj,jk)**2 + 10**2) ) ) &
               &          * zcfe(ji,jj,jk) * xstep * xdiss(ji,jj,jk)
         ELSE
         zlam1a   = ( 12.0  * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) + 9.05  * tr(ji,jj,jk,jppoc,Kbb) ) * xdiss(ji,jj,jk)    &
             &    + ( 2.49  * tr(ji,jj,jk,jppoc,Kbb) )     &
             &    + ( 127.8 * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) + 725.7 * tr(ji,jj,jk,jppoc,Kbb) )
         zaggdfea = zlam1a * xstep * zfecoll(ji,jj,jk)
         zaggdfee = 0.
         ENDIF
         ! define the dust concentration
         IF( ll_dust )  zdust  = dust(ji,jj) / ( wdust / rday ) * tmask(ji,jj,jk)
         IF( ln_felith ) THEN
                   zdust  = tr(ji,jj,jk,jplfe,Kbb) / 0.035 * 55.85 ! converted to g/L
                   zdust2 = tr(ji,jj,jk,jplfa,Kbb) / 0.035 * 55.85 ! converted to g/L
         ENDIF
         ! role of oxygen as a proxy for Fe2/Fe3 redox
         zxlam  = MAX( 1.E-3, (1. - EXP(-2 * tr(ji,jj,jk,jpoxy,Kbb) / 100.E-6 ) ))
         IF ( .NOT. ln_bait ) THEN
         zlam1b = 3.e-5 + ( xlamdust * zdust + xlam1 * ztrc ) * zxlam
         zscave = zFe3(ji,jj,jk) * zlam1b * xstep
         ! additional scavenging terms set to zero
         zscave1 = 0. ; zscave2 = 0. ; zscave3 = 0. ; zscave4 = 0. ; zscaveb = 0.
         ELSE
          ! scavenging by organics
          zxlam  = MAX( 1.E-3, (1. - EXP(-2 * tr(ji,jj,jk,jpoxy,Kbb) / 100.E-6 ) ) )
!          zlam1b = 3.e-5 + ( xlam1 * ztrc ) * zxlam ! organic particles
          zlam1b = ( xlam1 * ztrc ) * zxlam ! organic particles
          zscave = zFe3(ji,jj,jk) * zlam1b * xstep
          zscaveb = zFe3(ji,jj,jk) * 3.e-5 * xstep ! base scavenging rate
               ! lithogenic particles
               IF (ln_felith) THEN
               zscave1 = zFe3(ji,jj,jk) * ( xlamdust1 * zdust ) * zxlam * xstep
               zscave2 = zFe3(ji,jj,jk) * ( xlamdust2 * zdust2 ) * zxlam * xstep
               ELSE
               zscave1 = zFe3(ji,jj,jk) * (xlamdust * zdust ) * zxlam * xstep
               zscave2 = 0.
               ENDIF
               ! authgenic particles
               IF (ln_feauth ) then
               zscave3 = zFe3(ji,jj,jk) * ( xlamafe1 * tr(ji,jj,jk,jpafs,Kbb) ) * zxlam * xstep
               zscave4 = zFe3(ji,jj,jk) * ( xlamafe2 * tr(ji,jj,jk,jpafb,Kbb) ) * zxlam * xstep
               ELSE
               zscave3 = 0.
               zscave4 = 0.
               ENDIF
         ENDIF
            ! this bit the same for bait and Stnd
            zlam1b   = ( 1.94 * xdiss(ji,jj,jk) + 1.37 ) * tr(ji,jj,jk,jpgoc,Kbb)
            zaggdfeb = zlam1b * xstep * zfecoll(ji,jj,jk)
           ! correct even if bait closures are used:
           ! this is used in p4zligand
            xcoagfe(ji,jj,jk) =  zlam1a + zlam1b

          ! Scavenged iron is supposed to be released back to seawater
          ! when POM is solubilized. This is highly uncertain as probably
          ! a significant part of it may be rescavenged back onto 
          ! the particles. An efficiency factor is applied that is read
          ! in the namelist. 
          ! See for instance Tagliabue et al. (2019).
          ! Aggregated FeL is considered as biogenic Fe as it 
          ! probably remains  complexed when the particle is solubilized.
          ! -------------------------------------------------------------
          tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zscave * scaveff * tr(ji,jj,jk,jppoc,Kbb) / ztrc
          tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zscave * scaveff * tr(ji,jj,jk,jpgoc,Kbb) / ztrc
            IF ( ln_bait ) THEN
            tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zscave * scaveff2 * ( tr(ji,jj,jk,jpcal,Kbb) * calscav ) / ztrc
            tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zscave * scaveff2 * ( tr(ji,jj,jk,jpgsi,Kbb) * siscav ) / ztrc
            !tr(ji,jj,jk,jpgsf,Krhs) = tr(ji,jj,jk,jpgsf,Krhs) + zscave * scaveff2 * ( tr(ji,jj,jk,jpgsi,Kbb) * siscav ) / ztrc
            ELSE
          tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zscave * scaveff2 * tr(ji,jj,jk,jpcal,Kbb) / ztrc
          tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zscave * scaveff2 * tr(ji,jj,jk,jpgsi,Kbb)  / ztrc
            ENDIF
          !
          zscav3d(ji,jj,jk)   = zscave + zscaveb
         ! 
            ! STANDARD ACROSS ALL FE MODELS:
            ! fate of colloidal coag Fe (if not ln_bait then zscave3, zscave4,
            ! zaggdfee = 0 )
            IF (ln_feauth) THEN
               tr(ji,jj,jk,jpafb,Krhs) = tr(ji,jj,jk,jpafb,Krhs) + zaggdfeb
               tr(ji,jj,jk,jpafs,Krhs) = tr(ji,jj,jk,jpafs,Krhs) + zaggdfea + zaggdfee
            ELSE ! if we don't resolve the authigenics then add to sFe and bFe
               tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zaggdfea + zaggdfee
               tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zaggdfeb
            ENDIF
         ! now STANDARD ACROSS ALL FE MODELS:
         ! terms in bait model 0 by default if not used
          tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) - ( zscave + zscaveb) - ( zaggdfea + zaggdfeb ) &
              &          - ( zprecip + zprecipno3 ) - ( zscave1 + zscave2 ) - ( zscave3 + zscave4 ) &
              &          - zaggdfee
               IF (ln_felith) THEN
                 tr(ji,jj,jk,jplfe,Krhs) = tr(ji,jj,jk,jplfe,Krhs) + zscave1
                 tr(ji,jj,jk,jplfa,Krhs) = tr(ji,jj,jk,jplfa,Krhs) + zscave2
               ENDIF
               IF (ln_feauth) THEN
                 tr(ji,jj,jk,jpafs,Krhs) = tr(ji,jj,jk,jpafs,Krhs) + zscave3
                 tr(ji,jj,jk,jpafb,Krhs) = tr(ji,jj,jk,jpafb,Krhs) + zscave4
               ENDIF
          !
            zcoll3d(ji,jj,jk)   = zaggdfea + zaggdfeb + zaggdfee
            zscavl3d(ji,jj,jk)   = zscave1 + zscave2
            zscava3d(ji,jj,jk)   = zscave3 + zscave4
          !
      END_3D
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      biron(:,:,:) = tr(:,:,:,jpfer,Kbb) 
      !
      !  Output of some diagnostics variables
      !     ---------------------------------
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
         IF( iom_use("Fe3")    )  CALL iom_put("Fe3"    , zFe3   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
         IF( iom_use("FeL1")   )  CALL iom_put("FeL1"   , zFeL1  (:,:,:)       * tmask(:,:,:) )   ! FeL1
         IF( iom_use("TL1")    )  CALL iom_put("TL1"    , zTL1   (:,:,:)       * tmask(:,:,:) )   ! TL1
         IF( iom_use("Totlig") )  CALL iom_put("Totlig" , ztotlig(:,:,:)       * tmask(:,:,:) )   ! TL
         IF( iom_use("Biron")  )  CALL iom_put("Biron"  , biron  (:,:,:)  * 1e9 * tmask(:,:,:) )   ! biron
         IF( iom_use("FESCAV") )  CALL iom_put("FESCAV" , zscav3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FECOLL") )  CALL iom_put("FECOLL" , zcoll3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FEPREC") )  CALL iom_put("FEPREC" , zfeprecip(:,:,:) *1e9*tmask(:,:,:)*zrfact2 )
         IF( iom_use("FESCAVA") )  CALL iom_put("FESCAVA" , zscava3d(:,:,:) *1e9*tmask(:,:,:)*zrfact2 )
         IF( iom_use("FESCAVL") )  CALL iom_put("FESCAVL" , zscavl3d(:,:,:)*1e9*tmask(:,:,:)*zrfact2 )
         IF( iom_use("CFE") )      CALL iom_put("CFE" , zcfe(:,:,:) *1e6 * tmask(:,:,:) )
         IF( iom_use("FE3SOL") )   CALL iom_put("FE3SOL", fe3sol(:,:,:) *1e6 * tmask(:,:,:) )
      ENDIF

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_fechem')
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !!
      NAMELIST/nampisfer/ ln_ligvar, xlam1, xlamdust, ligand, kfep, scaveff,  &
      &      ln_fixlogk, logk, xlamdust1, xlamdust2, xlamafe1, xlamafe2, mincolfe, collf, kcfe &
      &      , kbcfe, coagf, siscav, calscav, scaveff2
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of iron chemistry parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnatp_ref, nampisfer, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisfer in reference namelist' )
      READ  ( numnatp_cfg, nampisfer, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisfer in configuration namelist' )
      IF(lwm) WRITE( numonp, nampisfer )

      IF(lwp) THEN                     ! control print
         WRITE(numout,*) '   Namelist : nampisfer'
         WRITE(numout,*) '      variable concentration of ligand          ln_ligvar    =', ln_ligvar
         WRITE(numout,*) '      scavenging rate of Iron                   xlam1        =', xlam1
         WRITE(numout,*) '      scavenging rate of Iron by dust           xlamdust     =', xlamdust
         WRITE(numout,*) '      ligand concentration in the ocean         ligand       =', ligand
         WRITE(numout,*) '      rate constant for nanoparticle formation  kfep         =', kfep
         WRITE(numout,*) '      Scavenged iron that is added to POFe      scaveff      =', scaveff
      IF( ln_bait) THEN
         WRITE(numout,*) '      fixed logK                                ln_fixlogk   =', ln_fixlogk
         WRITE(numout,*) '      logk                                      logk         =', logk
         WRITE(numout,*) '      scavenging rate of Iron by lfe            xlamdust1    =', xlamdust1
         WRITE(numout,*) '      scavenging rate of Iron by lfe aggregates = xlamdust2  =', xlamdust2
         WRITE(numout,*) ' minimum colloidal Fe fraction                  mincolfe     =',mincolfe
         WRITE(numout,*) ' factor enhancement of colloidal aggregation rate collf      =', collf
         WRITE(numout,*) '  shape function for cfe self aggregation       kcfe         =',kcfe
         WRITE(numout,*) ' shape function for the control of cfe agg with doc by biology kbcfe=', kbcfe
         WRITE(numout,*) ' factor modulation of cFe coagulation           coagf        =',coagf
         WRITE(numout,*) ' factor modulation of scavenging by biogenic silica siscav   =', siscav
         WRITE(numout,*) ' factor modulation of scavenging by calcite     calscav      =', calscav
         WRITE(numout,*) '      Scavenged iron that solubilised from BSi and CaCO3  scaveff2      =', scaveff2

       ENDIF
      ENDIF
      !
   END SUBROUTINE p4z_fechem_init
   
   !!======================================================================
END MODULE p4zfechem
