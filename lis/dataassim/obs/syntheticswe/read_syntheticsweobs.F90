!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"  !es
#include "LIS_NetCDF_inc.h" !es 
!BOP
! !ROUTINE: read_syntheticsweobs
!  \label{read_syntheticsweobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_syntheticsweobs(n, k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_historyMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod

  use LIS_mpiMod
  use LIS_dataAssimMod
  use LIS_pluginIndices
  use syntheticsweobs_module
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf !es
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic SWE observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField
  integer             :: fnd,p  !EC 08/06/2021
  integer             :: ftn !es 
  integer             :: c,r !es
  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)
 
  logical             :: data_upd_flag_local !EC 08/06/2021
  logical             :: data_upd_flag(LIS_npes) !EC 08/06/2021
  logical             :: data_upd !EC 08/06/2021
  integer             :: snodid !es
  real                :: swe_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)) !EC 08/06/2021
  real                :: snodobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)) !es
  character*100       :: sweobsdir
  logical             :: data_update
  logical             :: file_exists
  character*80        :: name

  logical             :: readflag
  integer             :: status

  integer             :: t

print *, 'read_syntheticswe_code' 
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call synswe_filename(name,sweobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
print *, readflag , name
  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     write(LIS_logunit,*)  'Reading syn data ',name
     
     call ESMF_StateGet(OBS_State,"Observation01",sweField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
 !    open(90, file=trim(name),form='unformatted')
 !    do t=1,1
 !       if(t==1) then 
 !          call LIS_readvar_gridded(90,n,obsl)
 !       else 
 !          call LIS_readvar_gridded(90,n,dummy)
 !       endif
 !    end do
 !   close(90)
!! es start
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     write(LIS_logunit,*)  '[INFO] Reading syn data ',name

     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SWE_tavg',snodid),&
          'Error nf90_inq_varid: SWE_tavg')

     call LIS_verify(nf90_get_var(ftn,snodid,snodobs,&
          start=(/LIS_ews_obs_halo_ind(n,LIS_localPet+1),&
          LIS_nss_obs_halo_ind(n,LIS_localPet+1)/),&
          count = (/LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)/)),&
          'Error in nf90_get_var')
     call LIS_verify(nf90_close(ftn))

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   snodobs(c,r)
           end if
        end do
     end do

#endif

!-------------------------------------------------------------------------
!  Apply LSM-based QC of observations
!-------------------------------------------------------------------------     
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_synsweId)//char(0),n, k, OBS_state)

     call LIS_checkForValidObs(n, k,obsl,fnd,swe_current)

     if(fnd.eq.0) then
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif

!#if (defined SPMD)
!     call MPI_ALLGATHER(data_upd_flag_local,1, &
!          MPI_LOGICAL, data_upd_flag(:),&
!          1, MPI_LOGICAL, LIS_mpi_comm, status)
!#endif
     data_upd = .false.
     do p=1,LIS_npes
        data_upd = data_upd.or.data_upd_flag(p)
     enddo

    readflag = .false.
    
    if(data_upd) then     
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
           fnd = 1 
        else
           assimflag(t) = 0 
        endif
     enddo

print *, 'assimflag', assimflag(t)

     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)
    
    if(LIS_rc%obs_ngrid(k).gt.0) then !es

     call ESMF_AttributeSet(sweField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
   
     call ESMF_AttributeSet(sweField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
   endif
  !   deallocate(dummy)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     endif
     write(LIS_logunit,*)  '[INFO] Finished reading syn data ',name
  else 
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)   
  return
  end if
  
  do t=1,LIS_rc%obs_ngrid(k)
     if(obsl(t).ne.-9999.0) then 
        if(obsl(t).lt.0.0) obsl(t) = 0.0
!        if(obsl(t).gt.200.0 ) obsl(t) = 200.0
     endif
  enddo


end subroutine read_syntheticsweobs

subroutine synswe_filename(name, ndir, yr, mo,da,hr,mn)
  
  implicit none
  character*80      :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  
  
!  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
!       trim(fmn)//'.d01.gs4r'

  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.nc'

end subroutine synswe_filename

!BOP
!
! !ROUTINE: readobsvar_1dgridded
! \label{readobsvar_1dgridded}
!
! !INTERFACE:
subroutine readobsvar_1dgridded_swe(ftn,n,k,var)
! !USES:
  use LIS_coreMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS:
  integer              :: ftn
  integer              :: n
  integer              :: k
  real                 :: var(LIS_rc%obs_ngrid(k))
!
! !DESCRIPTION:
!  This routine reads the observation data and subsets to the
!  local processor's domain, in a 1-d vector formulation.
!
!EOP

  real,  allocatable   :: gobs(:,:)
  integer              :: nc,c1,r1,c,r,gid

  allocate(gobs(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
  read(ftn) gobs

  nc = (LIS_ewe_obs_halo_ind(k,LIS_localPet+1)-&
       LIS_ews_obs_halo_ind(k,LIS_localPet+1))+1

  do r=LIS_nss_obs_halo_ind(k,LIS_localPet+1),&
       LIS_nse_obs_halo_ind(k,LIS_localPet+1)
     do c=LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
          LIS_ewe_obs_halo_ind(k,LIS_localPet+1)
        c1 = c-LIS_ews_obs_halo_ind(k,LIS_localPet+1)+1
        r1 = r-LIS_nss_obs_halo_ind(k,LIS_localPet+1)+1
        gid = LIS_obs_domain(n,k)%gindex(c1,r1)
        if(gid.ne.-1) then
           var(gid) = gobs(c,r)
        endif
     enddo
  enddo
  deallocate(gobs)

end subroutine readobsvar_1dgridded_swe


