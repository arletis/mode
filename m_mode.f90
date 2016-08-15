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

   USE scamtec_module                 ! module where the structure samtec is defined
   USE SCAM_dataMOD, only : scamdata  ! SCANTEC data matrix
   USE m_mode_objects
   USE m_mode_singleAttrib

   IMPLICIT NONE
   PRIVATE

   integer, allocatable :: Idx(:)     ! Array to save undefined points index
   integer              :: nidx       ! Total undefined points
 
   real, pointer     :: prefield(:,:) ! Pointer to save precipitation observation data
   real, pointer     :: expfield(:,:) ! Pointer to save precipitation experiment data 
   
   real, allocatable :: precOriginalField(:,:)  ! Matrix to save precipitation observation data 
   real, allocatable :: expOriginalField(:,:)   ! Matrix to save precipitation experiment data 


   public :: mode_init
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
	 precOriginalField = RESHAPE(prefield, (/scamtec%nypt,scamtec%nxpt/))
	 expOriginalField = RESHAPE(expfield, (/scamtec%nypt,scamtec%nxpt/))         

         print*,'mode_init'

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
      SUBROUTINE mode_run(nexp)
         Implicit None
         integer, intent(in) :: nexp ! experiment number

	 real, allocatable    		:: cfilter(:,:)      ! Matrix to save circular filter's values
	 real, allocatable    		:: convField(:,:)    ! Field resulting of convolution process
	 integer, allocatable 		:: maskField(:,:)    ! Binary field -mask- resulting of tresholding process 
	 real, allocatable    		:: restoreField(:,:) ! Matrix to save original field values where the mask is 1
	 real    			:: treshold	     ! Precipitation threshold defined by the user

         ! Radio and threshold values should be defined in scamtec.conf
	 integer 			:: i, j, rowSize, colSize, radio, binary_treshold=1, is_valid, totalObj, cont	 

	 integer, allocatable   	:: mask(:,:), maskObj(:,:) ! Mask to count objects -> resulting of Object Identification Algorithm 	 
	 ! objects attributes 
	 integer			:: perimeter, area, xcent, ycent
	 real				:: angle, aspect_ratio, complexity, area_hull 
	 type(mark), pointer		:: per_linked, total_linked 

	 real				:: area_tresh  ! Area threshold defined according to the area of interest to the user (scamtec.conf)

	 type(attrs), pointer		:: objects(:)
	 type(attrs_linked), pointer	:: objects_linked => NULL()
	 type(attrs_linked), pointer	:: attrs_aux => NULL()
	 type(attrs_linked), pointer	:: dir_aux => NULL()
	 type(mark), pointer		:: dir

         real, allocatable    		:: TESTE(:,:) ! EXEMPLOS PARA COMPROVAR OS ALGORITMOS
	 
	 call mode_init(nexp)


	 ! EXEMPLOS PARA COMPROVAR OS ALGORITMOS	 
	 rowSize=15
	 colSize=5
	 ! os valores do radio e o limiar devem ser definidos no scamtec.conf
	 radio=1
         treshold=0.5
	 area_tresh=2

	 allocate(TESTE(rowSize,colSize))
	 
	 allocate( convField(rowSize,colSize), maskField(rowSize,colSize), restoreField(rowSize,colSize) )	 	 
	 allocate( mask(rowSize,colSize), maskObj(rowSize,colSize) )


	 	 
 
	 !**** TESTE  *************

	 TESTE = 0
	 mask = 0
         maskObj = 0

	 TESTE(1,1)=60.0
	 TESTE(1,2)=60.0
	 TESTE(2,1)=60.0
	 TESTE(2,2)=59.9
	 TESTE(2,3)=59
	 TESTE(2,4)=60.0
	 TESTE(3,2)=59.0
	 TESTE(3,3)=50
	 TESTE(5,1)=52
	 TESTE(5,4)=52
	 TESTE(5,5)=30
	 TESTE(6,1)=30.01
	 TESTE(6,5)=30.01
	 TESTE(7,1)=32.1
	 TESTE(7,2)=32.1
	 TESTE(7,3)=30.01
	 TESTE(7,4)=32.1
	 TESTE(7,5)=32.2
	 TESTE(8,2)=32.1
	 TESTE(8,3)=32.1
	 TESTE(8,5)=32.25 
	 TESTE(9,1)=32.25 
	 TESTE(9,2)=32.25
         TESTE(9,5)=42
	 TESTE(10,5)=32.25	
	 TESTE(11,5)=42
	 TESTE(12,4)=42
	 TESTE(12,5)=32.25
	 TESTE(13,4)=32.25
	 TESTE(13,5)=32.251
	 TESTE(14,3)=32.25
	 TESTE(14,5)=32.25

	 print*
         print*, ' **** TESTE *****'
	 Do i=1, rowSize
	    write(*,*) (TESTE(i,j), j=1, colSize)
         ENDDO
 	
	 print*
         print*,' **** CIRCULAR FILTER ***** '

	 ! Calculate Circular filter
	 call circular_filter(cfilter, radio)
	 
            DO i=1, 3
               write(*,*) (cfilter(i,j), j=1, 3)
            ENDDO

	 print*	 
	 print*,' **** CONVOLUTION **** '

	 ! Convolution process
         call convolution(convField, TESTE, cfilter, rowSize, colSize)
         
         DO i=1, rowSize            
               write(*,*) (convField(i,j), j=1, colSize)            
         ENDDO

	 print*
	 print*,' **** tresholding ****'

	 ! Tresholding process
         call tresholding(TESTE, convField, maskField, restoreField, rowSize, colSize, treshold)

	 DO i=1, rowSize            
               write(*,*) (maskField(i,j), j=1, colSize)             
         ENDDO
	 print*
	 DO i=1, rowSize            
               write(*,*) (restoreField(i,j), j=1, colSize)         
         ENDDO

	 ! Loop to identify objects (rain areas) and attributes of each object
         print*
	 print*,' **** Object Identification ****'
         totalObj=0
	 DO j=1, colSize
	    DO i=1, rowSize            
               call valid(restoreField, rowSize, colSize, i, j, treshold, is_valid)		
	       if (is_valid .and. (mask(i,j) .EQ. 0) ) then		  
                  call singleObj_Ident_Attrib(i, j, restoreField, rowSize, colSize, treshold, mask, totalObj, maskObj, perimeter, area, xcent, ycent, angle, aspect_ratio, area_hull, complexity, per_linked, total_linked)		  
	
		  if (area .GE. area_tresh) then ! Creating temporal list to save objects attributes
		     allocate(attrs_aux)
		     allocate(attrs_aux%next)
		     attrs_aux%next => objects_linked
		     objects_linked => attrs_aux

		     objects_linked%xcent = xcent
		     objects_linked%ycent = ycent
		     objects_linked%area = area
		     objects_linked%perimeter = perimeter
		     objects_linked%angle = angle
		     objects_linked%aspect_ratio = aspect_ratio
		     objects_linked%complexity = complexity
		     objects_linked%area_hull = area_hull

		     allocate(objects_linked%pts_per(perimeter))
		     dir => per_linked
	    	     cont = 0 
		     do while (associated(dir))
			objects_linked%pts_per(cont+1)%x = dir%x 
			objects_linked%pts_per(cont+1)%y = dir%y
			dir => dir%next
			cont = cont+1
		     enddo

		     allocate(objects_linked%total_pts(area))
		     dir => total_linked
	    	     cont = 0 
		     do while (associated(dir))
			objects_linked%total_pts(cont+1)%x = dir%x 
			objects_linked%total_pts(cont+1)%y = dir%y
			dir => dir%next
			cont = cont+1
		     enddo

		     totalObj = totalObj + 1		     		     
		  endif          	       
	       endif	       
            ENDDO
         ENDDO

	 ! Save the attributes into an array
	 allocate(objects(totalObj))
	 dir_aux => objects_linked
	 i=0
	 Do while ()

	 Enddo

	 DO i=1, rowSize            
               write(*,*) (mask(i,j), j=1, colSize)         
         ENDDO

	 print*
	 DO i=1, rowSize            
               write(*,*) (maskObj(i,j), j=1, colSize)         
         ENDDO

         

         deallocate (cfilter)
	 deallocate (convField)
	 deallocate (maskField)
	 deallocate (TESTE)



      END SUBROUTINE mode_run
    !**************************************************************************************************************************************    



    


    




END MODULE mode

     

