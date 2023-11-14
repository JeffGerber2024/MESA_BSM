! The 'extras_finish_step' and 'low_mass_wind_scheme' subroutines were taken from the MIST run_star_extras.f file. These methods include new winds definition, a key for turning up Blocker scaling, a key to start diffusion, and a way to save a model immediately following the AGB. The rest of the subroutines are taken from the 1M_pre_ms_to_wd run_star_extras.f90 file with some slight adjustments (they are basically MESAs default methods for the test suite).
! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def
      use auto_diff
      
      implicit none

      real(dp) :: original_diffusion_dt_limit
      real(dp) :: burn_check = 0.0
      real(dp) :: postAGB_check = 0.0
      Logical :: is_TPAGB = .false.
      real(dp) :: rot_set_check = 0.0
      logical :: wd_diffusion = .false.
      real(dp) :: X_C_init, X_N_init

      include "test_suite_extras_def.inc"
      include 'xtra_coeff_os/xtra_coeff_os_def.inc'

      !NMM Parameters
      real(dp) :: mu_12 = 0.316

        !ALP Parameters
        real(dp) :: axion_g10
        real(dp) :: axion_g
        real(dp) :: m_axion !in eV

        !Global arrays used for interpolation
        real(dp) :: mary(1001), kary(61), finit(1001*61), params(2)
        real(dp), pointer :: f1(:)
        real(dp), target :: farray(4*1001*61)
        integer :: MD

      !HP Parameters
      real(dp) :: chi = 1d-11 ! eV/m kinetic mixing parameter constrained by solar neutrino measurements
      real(dp) :: m = 1d-1  !eV, energy of hidden photon

      ! these routines are called by the standard run_star check_model
      
      contains

      include "test_suite_extras.inc"
      include 'xtra_coeff_os/xtra_coeff_os.inc'
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         original_diffusion_dt_limit = s% diffusion_dt_limit
         include 'xtra_coeff_os/xtra_coeff_os_controls.inc'
	 s% other_neu => other_neu_HP !other_neu_axions, other_neu_HP, other_neu_NMM
         s% other_wind => low_mass_wind_scheme    
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls

      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
        integer :: j, cid
        real(dp) :: frac, rot_full_off, rot_full_on, frac2, vct30, vct100
        character(len=256) :: photosphere_summary, tau100_summary
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      
!     set OPACITIES: Zbase for Type 2 Opacities automatically to the Z for the star
      s% kap_rq% Zbase = s% initial_z
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) 'Zbase for Type 2 Opacities: ', s% kap_rq% Zbase
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) 'Initial y: ', s% initial_y
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


!     set ROTATION: extra param are set in inlist: star_job
      rot_full_off = s% job% extras_rpar(1) !1.2
      rot_full_on = s% job% extras_rpar(2) !1.8

      
      if (s% job% extras_rpar(3) > 0.0) then
         if (s% star_mass < rot_full_off) then
            frac2 = 0
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'no rotation'
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         else if (s% star_mass >= rot_full_off .and. s% star_mass <= rot_full_on) then
            frac2 = (s% star_mass - rot_full_off) / &
            (rot_full_on - rot_full_off)
            frac2 = 0.5d0*(1 - cos(pi*frac2))
            s% job% set_near_zams_omega_div_omega_crit_steps = 10
            s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         else
            frac2 = 1.0
            s% job% set_near_zams_omega_div_omega_crit_steps = 10
            s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2 !nominally 0.4
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         end if
      else
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'no rotation'
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      end if
      
      
!     set VARCONTROL: for massive stars, turn up varcontrol gradually to help them evolve
      vct30 = 1e-4
      vct100 = 3e-3
      
      if (s% initial_mass > 30.0) then
         frac = (s% initial_mass-30.0)/(100.0-30.0)
         frac = 0.5d0*(1 - cos(pi*frac))
         s% varcontrol_target = vct30 + (vct100-vct30)*frac
         
         if (s% initial_mass > 100.0) then
            s% varcontrol_target = vct100
         end if
         
         !CONVERGENCE TEST CHANGING C
         s% varcontrol_target = s% varcontrol_target * 1.0 

         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'varcontrol_target is set to ', s% varcontrol_target
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      end if
      end subroutine extras_startup

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going 
      !     define STOPPING CRITERION: stopping criterion for C burning, massive stars.
      if ((s% center_h1 < 1d-4) .and. (s% center_he4 < 1d-4)) then
         if ((s% center_c12 < 1d-4) .and. (s% initial_mass >= 10.0)) then
            termination_code_str(t_xtra1) = 'central C12 mass fraction below 1e-4'
            s% termination_code = t_xtra1
            extras_check_model = terminate
         else if ((s% center_c12 < 1d-2) .and. (s% initial_mass < 10.0)) then
            termination_code_str(t_xtra2) = 'central C12 mass fraction below 1e-2'
            s% termination_code = t_xtra2
            extras_check_model = terminate
	 else if ((PostAGB_check == 1) .and. (s% initial_mass < 10.0)) then
            termination_code_str(t_xtra3) = 'Reached Post AGB Phase, terminating'
            s% termination_code = t_xtra3
            extras_check_model = terminate
	 else if ((s% initial_mass > 6.5) .and. (s% initial_mass < 8.5)) then
            termination_code_str(t_xtra4) = 'Reached End of Helium Burning Phase, terminating'
            s% termination_code = t_xtra4
            extras_check_model = terminate
         end if
      end if
	  end function extras_check_model
      

      ! returns either keep_going, retry, or terminate.
integer function extras_finish_step(id)
	use chem_def
	integer, intent(in) :: id
	integer :: ierr
	real(dp) :: vct100, vct30, envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
            critmass, feh, frac2, mass_difference
        real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
        real(dp), parameter :: new_varcontrol_target = 1d-3
	logical :: diff_test1, diff_test2, diff_test3
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_finish_step = keep_going
! postAGB: save a model
        envelope_mass_fraction = 1d0 - max(s% he_core_mass, s% co_core_mass, s% one_core_mass)/s% star_mass
        if ((s% initial_mass < 10) .and. (envelope_mass_fraction < 0.1) .and. (s% center_h1 < 1d-4) .and. (s% center_he4 < 1d-4) &
        .and. (s% L_phot > 3.0) .and. (s% Teff > 7000.0)) then
	    	if (postAGB_check == 0.0) then !only print the first time
	    		write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	    		write(*,*) 'now at post AGB phase, saving a model'
	    		write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	    		
	    		!save a model
	    		call star_write_model(id, s% job% save_model_filename, ierr)
	    		postAGB_check = 1.0
	    	end if		  
        end if


! check DIFFUSION: to determine whether or not diffusion should happen
! no diffusion for fully convective, post-MS, and mega-old models
! do diffusion during the WD phase
	    min_center_h1_for_diff = 1d-10
	    diff_test1 = abs(s% mass_conv_core - s% star_mass) < 1d-2 !fully convective
	    diff_test2 = s% star_age > 5d10 !really old
	    diff_test3 = s% center_h1 < min_center_h1_for_diff !past the main sequence
	    if( diff_test1 .or. diff_test2 .or. diff_test3 )then
            s% diffusion_dt_limit = huge_dt_limit
        else
            s% diffusion_dt_limit = original_diffusion_dt_limit
	    end if
        
        if (wd_diffusion) then
            s% diffusion_dt_limit = original_diffusion_dt_limit
        end if

      
      end function extras_finish_step
      
 	subroutine low_mass_wind_scheme(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
        use star_def
        use chem_def, only: ih1, ihe4
        integer, intent(in) :: id
        real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
!       NOTE: surface is outermost cell. not necessarily at photosphere.
        real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
        integer, intent(out) :: ierr
        integer :: h1, he4
        real(dp) :: plain_reimers, reimers_w, blocker_w, vink_w, center_h1, center_he4
        real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc, Zsurf
        real(dp), parameter :: Zsolar_V = 0.019d0 ! for Vink et al formula (Look into this)
	    type (star_info), pointer :: s
	    ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
	    
        h1 = s% net_iso(ih1)
        he4 = s% net_iso(ihe4)
        center_h1 = s% xa(h1,s% nz)
        center_he4 = s% xa(he4,s% nz)
        Zsurf = 1.0 - (s% surface_h1 + s% surface_he3 + s% surface_he4)
        
        !reimers
        plain_reimers = 4d-13*(Lsurf*Rsurf/Msurf)/(Lsun*Rsun/Msun)
      
        reimers_w = plain_reimers * s% Reimers_scaling_factor
        
        !blocker
        blocker_w = plain_reimers * s% Blocker_scaling_factor * &
            4.83d-9 * pow(Msurf/Msun,-2.1d0) * pow(Lsurf/Lsun,2.7d0)
            
        !vink
        ! alfa = 1 for hot side, = 0 for cool side
        if (s% Teff > 27500d0) then
            alfa = 1
        else if (s% Teff < 22500d0) then
            alfa = 0
        else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
            Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Zsurf/Zsolar_V)))
            dT = 100d0
            if (s% Teff > Teff_jump + dT) then
                alfa = 1
            else if (s% Teff < Teff_jump - dT) then
                alfa = 0
            else
                alfa = (s% Teff - (Teff_jump - dT)) / (2*dT)
            end if
        end if
        
        if (alfa > 0) then ! eval hot side wind (eqn 24)
            vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
            vinf_div_vesc = vinf_div_vesc*pow(Zsurf/Zsolar_V,0.13d0) ! corrected for Z
            logMdot = &
               - 6.697d0 &
               + 2.194d0*log10(s% photosphere_L/1d5) &
               - 1.313d0*log10(s% photosphere_m/30) &
               - 1.226d0*log10(vinf_div_vesc/2d0) &
               + 0.933d0*log10(s% Teff/4d4) &
               - 10.92d0*pow2(log10(s% Teff/4d4)) &
               + 0.85d0*log10(Zsurf/Zsolar_V)
            w1 = exp10(logMdot)
        else
                w1 = 0
        end if
        
        if (alfa < 1) then ! eval cool side wind (eqn 25)
            vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
            vinf_div_vesc = vinf_div_vesc*pow(Zsurf/Zsolar_V,0.13d0) ! corrected for Z
            logMdot = &
               - 6.688d0 &
               + 2.210d0*log10(s% photosphere_L/1d5) &
               - 1.339d0*log10(s% photosphere_m/30) &
               - 1.601d0*log10(vinf_div_vesc/2d0) &
               + 1.07d0*log10(s% Teff/2d4) &
               + 0.85d0*log10(Zsurf/Zsolar_V)
            w2 = exp10(logMdot)
        else
            w2 = 0
        end if
        
        vink_w = alfa*w1 + (1 - alfa)*w2
        
        !for hot high mass MS stars, use V, then transition to R/B post-MS.
        !V is for 12.5k - 50k
        !for lower mass MS stars, use R or B.
        
        if (s% Teff > 12500d0) then
            w = vink_w
        else if ((s% Teff > 11500d0) .and. (s% Teff <= 12500d0)) then
            !transition from V to cool
            w = ((12500d0 - s% Teff)/(12500d0 - 11500d0)) * reimers_w &
            + (1.0 - (12500d0 - s% Teff)/(12500d0 - 11500d0)) * vink_w
        else
            !for below 11500
            !don't use B prior to AGB
            if (center_h1 < 0.01d0 .and. center_he4 < s% RGB_to_AGB_wind_switch) then
                w = max(reimers_w, blocker_w)
            else
                w = reimers_w
            end if
        end if
	
	end subroutine low_mass_wind_scheme 

      ! This function interpolates between the pre-processed grid of values stored in the work directory using MESA's interp2d function.
    subroutine other_neu_NMM(  &
            id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
         use neu_lib, only: neu_get
         use neu_def
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: log10_Tlim 
         logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: ierr
         
         real(dp) :: ye,ratio_mu_plas,mu_plas,mu_pair,omega_pl
         type (star_info), pointer :: s
         
         include 'formats'

         ierr = 0         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call neu_get(  &
            T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
         if (ierr /= 0) return

         !..   Add neutrino magnetic moment (Added by K. Mori on Jan. 2020)
         
         ye = zbar/abar
		 
	 omega_pl=28.7*sqrt(ye*Rho)/((1d0+(1.019d-6*ye*rho)**(2.0/3.0))**0.25)
	 ratio_mu_plas=0.318*(10d3/omega_pl)**2d0*mu_12**2d0
	 mu_plas=sources(plas_neu_type,ineu)*ratio_mu_plas	!Eq. (B1) in Heger et al. (2009)
		 
	 mu_pair=1.6d11*(mu_12*0.01)**2d0/(Rho/1d4)*exp(-118.5/(T/1d8)) !Eq. (B18) in Heger et al. (2009)
         
         loss(ineu) = loss(ineu) + mu_plas+mu_pair       
            
      end subroutine other_neu_NMM

      ! This function interpolates between the pre-processed grid of values stored in the work directory using MESA's interp2d function.
      ! It takes as input the value of log10(ma/kb*T) and log10(k_s/kb*T) specified in the other_neu subroutine below
      subroutine faxion_get(m_a, k_a, result, ierr, xary, yary, ND, f1)
          use interp_2d_lib_db, only: interp_rgbi3p_db, interp_mkbicub_db, interp_evbicub_db
          real(dp), intent(in) :: m_a, k_a, xary(1001), yary(61) ! Input: m_a (log10(ma/kb*T)), (k_a log10(k_s/kb*T))
                                                                 ! xary, yary - arrays of x and y values in interpolation range
          integer, intent(in) :: ND !=1 for new data, = 2 for old data
          real(dp), intent(out) :: result(6) ! output array - see interp_evbicub_db documentation for more information
          integer, intent(out) :: ierr
          real(dp), intent(inout), pointer :: f1(:) !Internal workspace for interp_evbicub_db
          integer :: i, j
          real(dp) :: xval(1), yval(1), zval(1), xmax(61), xmin(61), ymax(1001), ymin(1001)
          integer :: xsize, ysize, ilinx, iliny, ict(6)
          ict = [1,0,0,0,0,0] ! This array specifies what we would like interp_evbicub_db to calculate. See its documentation for more
                            ! information
          ! Boundary values, used for extrapolation if necessary
          ymin(:) = 0.0
          !ymax(:) = zary(:, 61)
          !xmin(:) = zary(1,:)
          xmax(:) = 0.0
          ! The point for which we want to know the z value. Has to be in an array for interp_evbicub_db to work.
          xval(1) = m_a
          yval(1) = k_a
          ! The size of each array
          xsize = SIZE(xary)
          ysize = SIZE(yary)
          ! Tells interp_mkbicub_db that xary and yary are linearly spaced.
          ilinx = 1
          iliny = 1
          if (ND == 1) then ! If ND==1, it implies that this is the first time this data is being used. interp_mkbicub_db needs to be run
                            ! to set up the internal workspace.
              call interp_mkbicub_db(xary, 1001, yary, 61, f1, 1001, 0, xmin, 0, xmax, 0, ymin, 0, ymax, ilinx, iliny, ierr)
          end if
          ! Call the interp_evbicub_db routine to determine the value of our function at m_a, k_a
          call interp_evbicub_db(m_a, k_a, xary, 1001, yary, 61, ilinx, iliny, f1, 1001, ict, result, ierr)
      end subroutine faxion_get


      subroutine other_neu_axions(  &
                  id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, loss, sources, ierr)
          use neu_lib, only: neu_get
          use neu_def
          integer, intent(in) :: id ! id for star
          integer, intent(in) :: k ! cell number or 0 if not for a particular cell
          real(dp), intent(in) :: T ! temperature
          real(dp), intent(in) :: log10_T ! log10 of temperature
          real(dp), intent(in) :: Rho ! density
          real(dp), intent(in) :: log10_Rho ! log10 of density
          real(dp), intent(in) :: abar ! mean atomic weight
          real(dp), intent(in) :: zbar ! mean charge
          real(dp), intent(in) :: log10_Tlim
          logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
          real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
          real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
          integer, intent(out) :: ierr
          ! Parameters important for calculating energy-loss to ALP production
          real(dp) :: ye, axionz2ye, axioncsi, faxioncsi, d_faxioncsi_dT, &
             sprimakoff, d_sprimakoff_dT, d_sprimakoff_dRho, m_T, result, k_a, res(6), log10mT, log10ka
          integer :: i,j
          type (star_info), pointer :: s
          include 'formats'
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          ! Call the standard neu_get routine. We add the contribution of ALPs to this.
          call neu_get(  &
             T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
             loss, sources, ierr)
          if (ierr /= 0) return
          ! Step One: importing data and writing it to arrays
          ! Perform this once the first time the routine is called
          if (MD/=2) then
              ! Define arrays of log10(mu) and log10(xi) values as defined in the paper
              do i=1, 61
                  kary(i) = -4 + 0.1*(i-1)
              end do
              do i=1, 1001
                  mary(i) = -5 + 0.01*(i-1)
              end do
              ! Open the data file - this contains the array of (F(mu, xi)+G(mu)) values from the paper
              ! NB file path should be absolute if run on a cluster
              open(3, file = "ALP_Data.txt", &
                          status = 'old', access = 'sequential', form = 'formatted', action = 'read')
              read(3,*, end = 102) finit !Write the contents into the array finit
102           close(3)
              ! Rewrite the contents of finit into the required form for the internal workspace of interp_mkbicub_db
              do i=1, 1001*61
                 farray(4*i-3) = finit(i)
                 farray(4*i-2) = 0.0
    		     farray(4*i-1) = 0.0
    		     farray(4*i) = 0.0
    	      end do
    	      f1 => farray
              ! The file bsm.txt contains the relevant ALP parameters (m_a and g10) - this was useful when running on a cluster
              ! This can be incorporated into an inlist using extras_rpar
              open(4, file = "bsm.txt", status = 'old', access = 'sequential', &
                             form = 'formatted', action = 'read')
    	      read(4,*) params
    	      close(4)
    		    axion_g10=params(1)
    		    m_axion=params(2)
    	      MD = 1
    	  end if
          ! ****** End Writing to Arrays *********
          ! Step Two: Calculate the value of epsilon_a
          ! Calculate log10(xi) as per our paper
          ! This section is modified from 1210.1271
          ye = s% ye(k)
          axionz2ye=zbar*zbar/abar+ye
          axioncsi=  1.658d20*axionz2ye*Rho/(T*T*T)
          !.. coefficient = 4*pi*alpha/(4*(1 Kelvin)^3) * N_A * cm^(-3)
          !..    pi*(1/137.0.35)/(1/11604.5)^3 * (6.022*10^23) * (197.326*10^6*10^(-13))^3
          m_T = m_axion/(8.61733d-5*T)
          log10mT = log10(m_T)
          k_a = axioncsi**(0.5)
          if (log10(k_a) < 2) then
              log10ka = log10(k_a)
          else
              log10ka = 1.99999 !Specify upper value for log10ka (fusion is dominant for large mass)
          end if
          ! Call routine if within the appropriate range - (don't want NANs in areas that aren't physically interesting)
          if (log10mT .LE. 5 .AND. log10mT .GE. -5) then
              call faxion_get(log10mT, log10ka, res, ierr, mary, kary, MD, f1)
              MD = 2 ! Arrays have been loaded and mkbicub run.
                     ! Change the value of global integer MD so this section will not be called again.
              faxioncsi = res(1)
              sprimakoff = 8.87736d-51*T**7*axion_g10**2*faxioncsi/(Rho)
              loss(ineu) = loss(ineu) + sprimakoff
              !loss(idneu_dT) = loss(idneu_dT) + d_sprimakoff_dT
              !loss(idneu_dRho) = loss(idneu_dRho) + d_sprimakoff_dRho
          end if
       end subroutine other_neu_axions   

    subroutine other_neu_HP(  &
            id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
         use neu_lib, only: neu_get
         use neu_def
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: log10_Tlim
         logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: ierr
         
         real(dp) :: ye, hprate, res, y0, xm, zm
         type (star_info), pointer :: s
         
         include 'formats'

         ierr = 0         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call neu_get(  &
            T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
         if (ierr /= 0) return

         !..   Add hidden photon production and energy loss.
         !..   Fitting formula by Ayala...
         
         ! ym = 0.d0
         ye = zbar/abar
         zm = zbar*zbar/abar

         ! omega=3.35d5*(Rho*ye)**0.5/((1.0d0+(1.019d-6*Rho*ye)**0.66667)**0.25)

         y0=3.33D5*(Rho*ye)**0.5/(T*(1.d0+(1.019d-6*Rho*ye)**0.666667)**0.25)
         xm=m/T
         xm=xm*1.16009d4
         res=y0**2-xm**2
         ! res=sqrt(res)

         if(res<0.0d0) then
         hprate=0.0d0
         else if(res>0.0d0) then
         hprate=1.62d4*(chi**2*m**2*y0**2*T**3)*sqrt(res)
         hprate=hprate/(exp(y0)-1.0d0)
         hprate=hprate/Rho
         ! print *, "hprate =", hprate
         ! print(hprate)
         end if

         loss(ineu) = loss(ineu) + hprate
         !loss(idneu_dT) = loss(idneu_dT) + d_hprate_dT
         !loss(idneu_dRho) = loss(idneu_dRho) + d_hprate_dRho

      end subroutine other_neu_HP
      end module run_star_extras