!-----------------------------------------------------------------------------------!
!             Group on Data Assimilation Development - GDAD/CPTEC/INPE              !
!                                                                                   !
!                                                                                   !
!       AUTHORS: Arletis Roque Carrasco                                             !
!              	 Luiz Fernando Sapucci                                              !
!										    !
!       Adapted from the work of: Arletis Roque Carrasco                            !
!				  Maibys Sierra Lorenzo                             !
!				  Israel Borrajero Montejo                          !
!				  Camilo Rodríguez Geno                             !
!                               						    ! 
!-----------------------------------------------------------------------------------!

! MODULE: m_mode.f90
!
! DESCRIPTON:
! Module for calculating the Method for Object-based Diagnostic Evaluation (MODE)
! The MODE is based on the objects identification by a convolution procedure
! whereby the fields are first smoothed over space and then thresholded by applying 
! an intensity threshold to the field. 

! Once objects are identified (contiguous nonzero pixels), they are merged and 
! matched by an algorithm of Fuzzy Logic utilizing information about various  
! attributes (centroid position, total area, area overlap, intensity distribution,  
! orientation angle, and boundary separation). 
!


MODULE mode

   USE scamtec_module                 		! module where the structure scantec is defined
   USE SCAM_dataMOD, only : scamdata  		! SCANTEC data matrix
   USE SCAM_Utils
   USE time_module, only: jul2cal, cal2jul 	! Time operations
   USE m_string                      		! string manipulations
   USE m_die                         		! Error Messages
   USE m_mode_objects		      		! module where objects are identified
   USE m_mode_singleAttrib	      		! module where single objects attributes are calculated
   USE m_mode_pairAttrib	      		! module where pair objects attributes are calculated

   IMPLICIT NONE
   PRIVATE

   integer, allocatable :: Idx(:)     ! Array to save undefined points index
   integer              :: nidx       ! Total undefined points
 
   real, pointer     :: prefield(:,:) ! Pointer to save precipitation observation data
   real, pointer     :: expfield(:,:) ! Pointer to save precipitation experiment data 
   
   real, allocatable :: precOriginalField(:,:)  ! Matrix to save precipitation observation data 
   real, allocatable :: expOriginalField(:,:)   ! Matrix to save precipitation experiment data 

   ! Statistical indices
   integer			:: hits, false_alarms, misses
   real				:: CSI, POD, FAR, BIAS

   !character(len=512) :: filename, fmt
   character(len=512) :: FNameOut = '%iy4%im2%id2%ih2%fy4%fm2%fd2%fh2'
   integer            :: FUnitOut = 30

   type statistic
     integer(I4B) 	:: atime
     integer, allocatable  :: fcst_time(:)
     integer, allocatable  :: misses(:), falseAlarms(:), hits(:)
     real, allocatable	:: csi(:)
     real, allocatable	:: pod(:)
     real, allocatable	:: far(:)
     real, allocatable	:: vies(:)
   end type statistic

   type(statistic), allocatable :: indices(:,:)

   integer, public, Parameter :: NumIndices = 8
   character(len=8), public, parameter ::   IndicesName(1:NumIndices) = (/ &
                                           'Forecast',& ! Virtual Temperature @ 925 hPa [K]
                                           'Misses',& ! Virtual Temperature @ 850 hPa [K]
                                           'FAlarms',& ! Virtual Temperature @ 500 hPa [K]                                           
                                           'Hits',& ! Absolute Temperature @ 850 hPa [K]
                                           'CSI',& ! Absolute Temperature @ 500 hPa [K]
                                           'POD',& ! Absolute Temperature @ 250 hPa [K]
                                           'FAR',& ! Pressure reduced to MSL [hPa]
					   'BIAS' & ! CONVECTIVE PRECIPITATION @ 1000 hPa [kg/m2/day]
                                          /)

   public :: mode_init
   public :: mode_ObjectIdentf
   public :: mode_run



   Contains

    !**************************************************************************************************************************************
      Subroutine mode_init(nexp)  
      ! Subroutine where variables and structures memory is allocated.     

         Implicit None

	 ! Input Parameter         
	 integer, intent(in) :: nexp      ! experiment number

         integer            :: i, npts      

    	 npts = scamtec%nxpt*scamtec%nypt                 ! scamtec -> variable type definida em scamtec_module.f90

	 Allocate(precOriginalField(scamtec%nypt,scamtec%nxpt))
	 Allocate(expOriginalField(scamtec%nypt,scamtec%nxpt))          
        
         ! transferindo dados de precipitacion 
	 prefield => scamdata(1)%prefield             
         expfield => scamdata(nexp)%expfield           

	 ! Convirtiendo el vector de los datos de precipitacion para una matriz 
	 precOriginalField = RESHAPE(prefield(:,21), (/scamtec%nypt,scamtec%nxpt/))
	 expOriginalField = RESHAPE(expfield(:,hist%tipo_precip), (/scamtec%nypt,scamtec%nxpt/))      

         !open(46,file=trim(scamtec%output_dir)//'/'//'EXP_precip'//'.bin',form='unformatted',status='unknown',access = 'sequential')         
         !write(46)expOriginalField 

         !print*,'mode_init' 

	 !print*,'Min/Max PREFIELD_MODE: ',minval(precOriginalField(:,:)),maxval(precOriginalField(:,:))  

         !print*,'Min/Max EXPFIELD_MODE: ',minval(expOriginalField(:,:)),maxval(expOriginalField(:,:))  
         !stop  

         if(.NOT.Allocated(indices))Allocate(indices(scamtec%nexp,scamtec%ntime_forecast))

         !print*, 'scamtec%ntime_forecast', scamtec%ntime_forecast
 	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%fcst_time))Allocate(indices(nexp,scamtec%ftime_count(1))%fcst_time(scamtec%ntime_forecast-1))
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%misses))Allocate(indices(nexp,scamtec%ftime_count(1))%misses(scamtec%ntime_forecast-1))
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%falseAlarms))Allocate(indices(nexp,scamtec%ftime_count(1))%falseAlarms(scamtec%ntime_forecast-1))
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%hits))Allocate(indices(nexp,scamtec%ftime_count(1))%hits(scamtec%ntime_forecast-1))	
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%csi))Allocate(indices(nexp,scamtec%ftime_count(1))%csi(scamtec%ntime_forecast-1))
         if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%pod))Allocate(indices(nexp,scamtec%ftime_count(1))%pod(scamtec%ntime_forecast-1))
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%far))Allocate(indices(nexp,scamtec%ftime_count(1))%far(scamtec%ntime_forecast-1))
	 if(.NOT.Allocated(indices(nexp,scamtec%ftime_count(1))%vies))Allocate(indices(nexp,scamtec%ftime_count(1))%vies(scamtec%ntime_forecast-1))

      End Subroutine mode_init
    !**************************************************************************************************************************************



    !************************************************************************************************************************************** 
      Subroutine mode_finalize(nexp)
      ! Subroutine where variables and structures memory is released.

	 Implicit None

         ! INPUT PARAMETER:
         integer, intent(in) :: nexp ! experiment number

         ! Desassociando ponteiros
         if (associated(expfield)) nullify(expfield)
         if (associated(prefield)) nullify(prefield)

         ! Desalocando variaveis
         !DeAllocate(Idx)
         DeAllocate(precOriginalField)
         DeAllocate(expOriginalField)

      End Subroutine mode_finalize 
    !**************************************************************************************************************************************



    !**************************************************************************************************************************************
      Subroutine mode_run(nexp)
         Implicit None
         integer, intent(in) :: nexp ! experiment number

			   ! dimensions of the fields to compare and loop variables
         integer 	:: rowSize, colSize, i, j, t, f, idField  

	 		   ! Weights of the attributes used in the fuzzy logic
	 real		:: min_boundary_dist_weight, dif_centroid_weight, area_ratio_weight, perimeter_ratio_weight, dif_angle_weight, aspect_ratio_weight, complexity_ratio_weight, int_area_ratio_weight, weight(8), total_interest_tresh, grid_res  

	 real, allocatable		:: obsConvField(:,:)    ! Field resulting of convolution process
	 real, allocatable		:: expConvField(:,:)    ! Field resulting of convolution process         
	 real, allocatable		:: obsRestoreField(:,:) ! Matrix to save original field values where the mask is 1
	 real, allocatable		:: expRestoreField(:,:) ! Matrix to save original field values where the mask is 1
       

					  ! Mask to count objects -> resulting of Object Identification Algorithm
         integer, allocatable   	:: mask(:,:), obs_maskObj(:,:), exp_maskObj(:,:) 

					  ! Identified objects total (Observation and Forecast)
         integer 			:: prec_nobj, exp_nobj

					  ! List to save objects and attributes (Observation and Forecast)
         type(attrs), pointer		:: prec_objects(:), exp_objects(:)
         type(atrib_pair), pointer	:: atrib_matched(:)

					  ! Masks resulting of Matching Algorithm (Objects pairs have the same id in each field)
	 integer, allocatable   	:: obsMatch_mask(:,:), expMatch_mask(:,:)
					  ! Auxiliary variables used in the Matching Algorithm
	 integer			:: fcst_id, obs_id, num, x, y, cont 
         type(attrs)			:: objaux

!					  ! Statistical indices
!	 integer			:: hits, false_alarms, misses
!	 real				:: CSI, POD, FAR, BIAS

         call mode_init(nexp)

         rowSize = scamtec%nypt
         colSize = scamtec%nxpt          

         grid_res = 0.400  ! Este valor creo q es dom(I)%x  

         ! Attributes weight used in Merging and Matching process
         min_boundary_dist_weight = 0.0
	 dif_centroid_weight = 4.0
	 area_ratio_weight = 2.0
	 perimeter_ratio_weight = 0.0
	 dif_angle_weight = 1.0
	 aspect_ratio_weight = 0.0
	 complexity_ratio_weight = 0.0
	 int_area_ratio_weight = 2.0
	 total_interest_tresh = 0.5

         weight(1) = min_boundary_dist_weight
	 weight(2) = dif_centroid_weight 
	 weight(3) = area_ratio_weight
	 weight(4) = perimeter_ratio_weight
	 weight(5) = dif_angle_weight
	 weight(6) = aspect_ratio_weight
	 weight(7) = complexity_ratio_weight
	 weight(8) = int_area_ratio_weight	 
	        
         f = scamtec%ftime_idx 
         
	      ! Subroutine defined in m_mode_objects where convolution, tresholding, identification of objects, attributes calculation and 
	      ! merging algorithms are made
	 ! Observation
         call mode_ObjectIdentf(rowSize, colSize, precOriginalField, obsConvField, obsRestoreField, weight, total_interest_tresh, grid_res, mask, obs_maskObj, prec_nobj, prec_objects)
            !print*
	    !print*, 'prec_nobj', prec_nobj
            ! imprimiendo campo objeto.
	    !DO i=1, rowSize            
               !write(*,*) (obsConvField(i,j), j=1, colSize)         
            !ENDDO

	    !print*
	    !DO i=1, rowSize            
               !write(*,*) (obs_maskObj(i,j), j=1, colSize)         
            !ENDDO
	 idField = 1  
	 call mode_writeFields(nexp, idField, precOriginalField, obsConvField, mask, obsRestoreField) 
	 !stop                     
         
	 ! Forecast
         call mode_ObjectIdentf(rowSize, colSize, expOriginalField, expConvField, expRestoreField, weight, total_interest_tresh, grid_res, mask, exp_maskObj, exp_nobj, exp_objects)
            !print*
	    !print*, 'exp_nobj', exp_nobj
	    !DO i=1, rowSize            
               !write(*,*) (expOriginalField(i,j), j=1, colSize)         
            !ENDDO

	    !print*
	    !DO i=1, rowSize            
               !write(*,*) (exp_maskObj(i,j), j=1, colSize)         
            !ENDDO
	 idField = 0	    
         call mode_writeFields(nexp, idField, expOriginalField, expConvField, mask, expRestoreField) 
	 stop                     
	 
         If  (f .GT. 1) then
	     ! Subroutine defined in m_mode_pairAttrib where observation objects attributes and forecast objects attributes are compared to select pair objects         
           call object_matching(prec_nobj, prec_objects, exp_nobj, exp_objects, weight, grid_res, total_interest_tresh, atrib_matched, cont)
         !stop
           allocate(obsMatch_mask(rowSize,colSize))
	   allocate(expMatch_mask(rowSize,colSize))
	 
!          obsMatch_mask = obs_maskObj
!	   expMatch_mask = exp_maskObj

           obsMatch_mask = 0
	   expMatch_mask = 0	 
         
           num = 1
	   do i=1, cont
	     fcst_id = atrib_matched(i)%id1
	     obs_id = atrib_matched(i)%id2           

	     objaux = prec_objects(obs_id)
	     !print*, 'prec_objects(obs_id)total_pts', prec_objects(obs_id)%total_pts
             do j=1, objaux%area
	       x = objaux%total_pts(j)%x
	       y = objaux%total_pts(j)%y
	     
	       obsMatch_mask(y,x) = num
	     enddo

	     objaux = exp_objects(fcst_id)

             do j=1, objaux%area
	       x = objaux%total_pts(j)%x
	       y = objaux%total_pts(j)%y
	     
	       expMatch_mask(y,x) = num
	     enddo
             num = num + 1	   
           enddo	   

	   !print*
           !DO i=1, rowSize            
             !write(*,*) (obsMatch_mask(i,j), j=1, colSize)         
           !ENDDO

	   !print*
	   !DO i=1, rowSize            
             !write(*,*) (expMatch_mask(i,j), j=1, colSize)         
           !ENDDO	 

	   call mode_finalize(nexp)         

	 !*********** Object-based Statistical Índices *********************************

	 misses = abs(prec_nobj - cont)
	 false_alarms = abs(exp_nobj - cont)
	 hits = cont

	 print*
	 print*, 'misses', misses, 'false_alarms', false_alarms, 'hits', hits

	 CSI = REAL(hits) / (REAL(hits) + REAL(misses) + REAL(false_alarms))
	 print*
	 print*,  'CSI', CSI

	 POD = REAL(hits) / (REAL(hits) + REAL(misses))
	 print*
	 print*,  'POD', POD

	 FAR = REAL(false_alarms) / (REAL(hits) + REAL(false_alarms))
	 print*
	 print*,  'FAR', FAR

	 BIAS = (REAL(hits) + REAL(false_alarms)) / (REAL(hits) + REAL(misses))
	 print*
	 print*, 'BIAS', BIAS        

	indices(nexp,scamtec%ftime_count(1))%atime = scamtec%atime
	indices(nexp,scamtec%ftime_count(1))%fcst_time(scamtec%ftime_idx-1) = int(abs(cal2jul(scamtec%atime)-cal2jul(scamtec%ftime))*24)
	indices(nexp,scamtec%ftime_count(1))%misses(scamtec%ftime_idx-1) = misses
	indices(nexp,scamtec%ftime_count(1))%falseAlarms(scamtec%ftime_idx-1) = false_alarms
	indices(nexp,scamtec%ftime_count(1))%hits(scamtec%ftime_idx-1) = hits
        indices(nexp,scamtec%ftime_count(1))%csi(scamtec%ftime_idx-1) = CSI
        indices(nexp,scamtec%ftime_count(1))%pod(scamtec%ftime_idx-1) = POD
	indices(nexp,scamtec%ftime_count(1))%far(scamtec%ftime_idx-1) = FAR
	indices(nexp,scamtec%ftime_count(1))%vies(scamtec%ftime_idx-1) = BIAS

	if (scamtec%ftime_idx .EQ. scamtec%ntime_forecast) then
	   !nprint = 3
	   call mode_write(nexp)	   
        endif
      Endif
        

      End Subroutine mode_run
    !**************************************************************************************************************************************



    !**************************************************************************************************************************************
      Subroutine mode_writeFields(nexp, id, Original, Convolution, Mask, Restored)
        Implicit None
        integer, intent(in) :: nexp,  id ! experiment number
	real, allocatable, intent(in) 		:: Original(:,:)
	real, allocatable, intent(in)		:: Convolution(:,:)    ! Field resulting of convolution process
	real, allocatable, intent(in)		:: Restored(:,:)
        integer, allocatable, intent(in)   	:: Mask(:,:) ! Mask to count objects -> resulting of Object Identification Algorithm
	 
	integer            :: iret, i,j, ier
	character(len=512) :: filename, fname, fmt
	integer            :: nymd, nhms
    	integer            :: fymd, fhms	

        nymd = scamtec%atime/100
        nhms = MOD(scamtec%atime,100) * 10000
        fymd = scamtec%ftime/100
        fhms = MOD(scamtec%ftime,100) * 10000

	If (id .EQ. 1) then

	  fname = 'PrecipField'
	  inquire(unit=FUnitOut, opened=iret)
	  if(.not.iret) then 
	    filename = trim(fname)//'_'//trim(FNameOut)
            call str_template(filename, nymd, nhms, fymd, fhms, label=num2str(nexp,'(I2.2)'))

	    open(unit   = FUnitOut+0,	&
	         File   = trim(scamtec%output_dir)//'/Original'//Trim(filename)//'.bin',   &
	         status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+1,	&
	         File   = trim(scamtec%output_dir)//'/Convolution'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+2,	&
	         File   = trim(scamtec%output_dir)//'/Mask'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+3,	&
	         File   = trim(scamtec%output_dir)//'/Restored'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)

	    write(FUnitOut+0)Original	        
	    write(FUnitOut+1)Convolution
	    write(FUnitOut+2)Mask
	    write(FUnitOut+3)Restored

	    Close(FUnitOut+0)
            Close(FUnitOut+1)
    	    Close(FUnitOut+2)
    	    Close(FUnitOut+3)
	  endif

	Else

	  fname = 'ExpField'
	  inquire(unit=FUnitOut, opened=iret)
	  if(.not.iret) then 
	    filename = trim(fname)//'_'//trim(FNameOut)
            call str_template(filename, nymd, nhms, fymd, fhms, label=num2str(nexp,'(I2.2)'))

	    open(unit   = FUnitOut+4,	&
	         File   = trim(scamtec%output_dir)//'/Original'//Trim(filename)//'.bin',   &
	         status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+5,	&
	         File   = trim(scamtec%output_dir)//'/Convolution'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+6,	&
	         File   = trim(scamtec%output_dir)//'/Mask'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    open(unit   = FUnitOut+7,	&
	         File   = trim(scamtec%output_dir)//'/Restored'//Trim(filename)//'.bin',   &
                 status='unknown', &
                 form   = 'unformatted', &
                 access = 'sequential'  &                 
		)
	    write(FUnitOut+4)Original	        
	    write(FUnitOut+5)Convolution
	    write(FUnitOut+6)Mask
	    write(FUnitOut+7)Restored

	    Close(FUnitOut+4)
            Close(FUnitOut+5)
    	    Close(FUnitOut+6)
    	    Close(FUnitOut+7)
	  endif

	Endif
      End Subroutine mode_writeFields
    !**************************************************************************************************************************************
      
  

    !**************************************************************************************************************************************
      Subroutine mode_write(nexp)
	Implicit None
        integer, intent(in) :: nexp ! experiment number
	integer            :: iret, nparameters, i
	character(len=512) :: filename, fname, fmt
	integer            :: nymd, nhms
    	integer            :: fymd, fhms

	integer(I4B) 	:: a
     	integer	:: b,c,d,e     	
     	real	:: f,g,h,j

        nymd = scamtec%atime/100
        nhms = MOD(scamtec%atime,100) * 10000
        fymd = scamtec%ftime/100
        fhms = MOD(scamtec%ftime,100) * 10000	

        !If (nprint .EQ. 3) then
          fname = 'StatisticIndices'
	  nparameters = 9

	  inquire(unit=FUnitOut+8, opened=iret)
	  if(.not.iret) then 
	    filename = trim(fname)//'_'//trim(FNameOut)
            call str_template(filename, nymd, nhms, fymd, fhms, label=num2str(nexp,'(I2.2)'))

	    open(unit   = FUnitOut+8,	&
	         File   = trim(scamtec%output_dir)//Trim(filename)//'T.scam',   &
                 access = 'sequential',  &
                 Form   = 'formatted', &
                 Status = 'replace'      &
                )

	    write(FUnitOut+5,'(A)')'%Analysis    Forecast    Misses   FAlarms  Hits      CSI         POD       FAR       BIAS'

	    Do i=1, scamtec%ftime_idx-1
	      a = indices(nexp,scamtec%ftime_count(1))%atime
	      b = indices(nexp,scamtec%ftime_count(1))%fcst_time(i)
	      c = indices(nexp,scamtec%ftime_count(1))%misses(i)
	      d = indices(nexp,scamtec%ftime_count(1))%falseAlarms(i)
	      e = indices(nexp,scamtec%ftime_count(1))%hits(i)
	      f = indices(nexp,scamtec%ftime_count(1))%csi(i)
	      g = indices(nexp,scamtec%ftime_count(1))%pod(i)
	      h = indices(nexp,scamtec%ftime_count(1))%far(i)
	      j = indices(nexp,scamtec%ftime_count(1))%vies(i)
	    
	      write(FUnitOut+8,95)a,b,c,d,e,f,g,h,j
95            FORMAT(I10,6X,I2,3X,3(7X,I1),2X,4(3X,F8.6))
	    Enddo
  	    close(FUnitOut+8)	  
          endif	
	!Endif

      End Subroutine mode_write
    !**************************************************************************************************************************************
  




END MODULE mode

     

