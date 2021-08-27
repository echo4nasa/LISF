!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_qc_snowobs
! \label{noahmp401_qc_snowobs}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  30 Jan 2015: Yuqiong Liu; added additional QC
!  03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
!  14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SNOW
!
! !INTERFACE:
subroutine noahmp401_qc_snowobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the snow observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  ground is covered with vegatation (more than 50%).  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_snow_field

  real, pointer            :: snowobs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))
  real                     :: fveg_obs(LIS_rc%obs_ngrid(k))
  real                     :: tv_obs(LIS_rc%obs_ngrid(k))
  real                     :: snowliq(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: snowliq_obs(LIS_rc%obs_ngrid(k))
  real                     :: isnow_obs(LIS_rc%obs_ngrid(k))
!  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))

print *, 'Reading /da_snow/noahmp401_qc_snowobs'

  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noahmp401_qc_snowobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=snowobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in noahmp401_qc_snowobs")
  
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     !stc1(t) = noahmp401_struc(n)%noahmp401(t)%sstc(1) ! get snow/veg temp.
     stc1(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1) ! get snow/veg temp.
     snowliq(t) = noahmp401_struc(n)%noahmp401(t)%snowliq(1) ! get snow layer liquid water [mm] 
     vegt(t) = LIS_surface(n,1)%tile(t)%vegt
  enddo

  call LIS_convertPatchSpaceToObsSpace(n,k,LIS_rc%lsm_index, & 
       noahmp401_struc(n)%noahmp401(:)%tv,tv_obs) !tv: vegetation temperature. unit: K 
  call LIS_convertPatchSpaceToObsSpace(n,k,LIS_rc%lsm_index, &    !fveg: green vegetation fraction. unit: - 
       noahmp401_struc(n)%noahmp401(:)%fveg,fveg_obs)
!  call LIS_convertPatchSpaceToObsSpace(n,k,LIS_rc%lsm_index, &    !isnow: no. of snow layer.  unit: - 
!       noahmp401_struc(n)%noahmp401(:)%isnow,isnow_obs) 
!! add rainf flag by Eunsang Cho 08/05/2021
!  call LIS_convertPatchSpaceToObsSpace(n,k,LIS_rc%lsm_index, & 
!       noahmp401_struc(n)%noahmp401(:)%prcp,prcp_obs)  ! prcp: total precip [in], rainf: rainfall rate [out] 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,snowliq,snowliq_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,stc1,stc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,vegt,vegt_obs)

!! Eunsang: check for rain -----
!  do t = 1,LIS_rc%obs_ngrid(k)
!        if(rainf_obs(t).gt.3E-6) then  
!           snowobs(t) = LIS_rc%udef
!           print*, 'rainf ',gid,t,noahmp401_struc(n)%noahmp401(t)%prcp
!     endif
!  enddo

  do t=1,LIS_rc%obs_ngrid(k)
     if(snowobs(t).ne.LIS_rc%udef) then
       if(fveg_obs(t).gt.0.7) then ! green vegetation fraction
           snowobs(t) = LIS_rc%udef
         elseif(vegt_obs(t).le.4) then ! forest types
           snowobs(t) = LIS_rc%udef
!assume that snow will not form at 3 deg. celcius or higher soil (stc1) or vegetation (tv) temp. 
         elseif(tv_obs(t).ge.276.15) then
           snowobs(t) = LIS_rc%udef
         elseif(stc1_obs(t).ge.276.15) then
           snowobs(t) = LIS_rc%udef
        endif
! if snow layer liquid water [mm] > 0, mask by Eunsang Cho 8/24
       if(snowliq_obs(t).gt.0) then
           snowobs(t) = LIS_rc%udef
        endif
     endif
  enddo

!! ---------------------------
end subroutine noahmp401_qc_snowobs


