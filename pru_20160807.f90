PROGRAM pru

IMPLICIT NONE
   !PRIVATE


   	integer				:: i, j, ysize=15, xsize=5, treshold=1, totalObj, is_valid
	integer, allocatable 		:: precOriginalField(:,:)  ! Matrix to save precipitation observation data 
	integer, allocatable    	:: mask(:,:), maskObj(:,:)	! Mask to count objects
	! Arrays to save objects attributes 
	!integer, allocatable    	:: xcent(:), ycent(:), area(:), perimeter(:), area_hull(:)  
	!real, allocatable	    	:: angle(:), aspect_ratio(:), complexity(:) 	
	integer				:: perimeter, area, xcent, ycent
	real				:: angle, aspect_ratio
   

   	! ESTRUTURA USADA NA SUBROUTINE  singleObjIdent_Attrib   PARA IDENTIFICAR OBJETOS E ATRIBUTOS
	type mark			
      		integer			:: x, y
      		type(mark), pointer	:: next
   	end type
	   type(mark), pointer		:: points	   
      
   	type point
      		integer	:: x, y
   	end type    
      	   type(point), pointer	:: total_pts(:)  

   	type attrs
      		integer			:: id, area, xcent, ycent, perimetro
      		type(point), pointer	:: pts_per, total_pts
   	end type   
      

 

      	Allocate(precOriginalField(ysize,xsize))
	Allocate(mask(ysize,xsize))
	Allocate(maskObj(ysize,xsize))

	!Inicializando las matriz donde serán detectados los objetos
      	do i=1, ysize
           do j=1, xsize
              precOriginalField(i,j)=0
	      mask(i,j)=0
	      maskObj(i,j)=0
           enddo
      	enddo

	!precOriginalField(1,1)=1
        precOriginalField(1,1)=1
	precOriginalField(1,2)=1
	precOriginalField(2,1)=1      	
	precOriginalField(2,2)=1
      	precOriginalField(2,3)=1
	precOriginalField(2,4)=1
      	precOriginalField(3,2)=1
      	precOriginalField(3,3)=1
	!precOriginalField(3,4)=1
	precOriginalField(5,1)=1	
	precOriginalField(5,4)=1
	precOriginalField(5,5)=1
	precOriginalField(6,1)=1
	precOriginalField(6,5)=1
	precOriginalField(7,1)=1
	precOriginalField(7,2)=1
	precOriginalField(7,5)=1
	precOriginalField(8,2)=1
	precOriginalField(8,3)=1
	precOriginalField(8,5)=1
	precOriginalField(9,1)=1
	precOriginalField(9,2)=1
	precOriginalField(9,5)=1
	precOriginalField(10,5)=1
	precOriginalField(11,5)=1
	precOriginalField(12,4)=1
	precOriginalField(12,5)=1
	precOriginalField(13,4)=1
	precOriginalField(13,5)=1
	precOriginalField(14,3)=1
	precOriginalField(14,5)=1

        
	! Loop para detectar cantidad de objetos y calcular atributos de cada objeto
	print*
	 print*,' **** Object Identification ****'
         totalObj=0
	 DO j=1, xsize
	    DO i=1, ysize            
               call valid(precOriginalField, ysize, xsize, i, j, treshold, is_valid)		
	       if (is_valid .and. (mask(i,j) .EQ. 0) ) then		  
                  call singleObj_Ident_Attrib(i, j, precOriginalField, ysize, xsize, treshold, mask, totalObj, maskObj, perimeter, area, xcent, ycent, angle, aspect_ratio)
		  totalObj = totalObj + 1
	          maskObj(i,j) = totalObj	       
	       endif	       
            ENDDO
         ENDDO

	print*	

	DO i=1, ysize            
           write(*,*) (mask(i,j), j=1, xsize)         
        ENDDO

	print*
	DO i=1, ysize            
           write(*,*) (maskObj(i,j), j=1, xsize)         
        ENDDO


   Contains   
  

    !**************************************************************************************************************************************
      ! The subroutine verify valid points in the object field
      Subroutine valid(restoreField, ysize, xsize, i, j, treshold, is_valid)
	
	implicit none
	! INPUT PARAMETERS:
	integer, allocatable, intent(in)   	:: restoreField(:,:)  ! Object Field resulting from tresholding subroutine
	integer, intent(in)            		:: i, j, xsize, ysize    ! Object matrix dimensions
	integer, intent(in)	     		:: treshold

        ! OUTPUT PARAMETERS:
	integer, intent(out)		:: is_valid
        
	is_valid=0
	if ( (i .GT. 0) .And. (i .LE. ysize) .And. (j .GT. 0) .And. (j .LE. xsize) ) then 
	   if (restoreField(i,j) .GE. treshold) then 
	      is_valid = 1
	   endif
	endif

	return

      end Subroutine valid
    !**************************************************************************************************************************************



    !**************************************************************************************************************************************
	Subroutine singleObj_Ident_Attrib(row, col, restoreField, ysize, xsize, treshold, mask, totalObj, maskObj, perimeter, area, xcent, ycent, angle, aspect_ratio)

	   Implicit none
	   integer, intent(in)           	    	:: row, col		! Loop variables
	   integer, allocatable, intent(in)    		:: restoreField(:,:)	! Object Field resulting from tresholding subroutine
	   integer, intent(in)           	    	:: xsize, ysize		! Object matrix dimensions 	
	   integer, intent(in)	      	    		:: treshold, totalObj			! Rain treshold 
		
	   !OUTPUT PARAMETERS:
	   integer, allocatable, intent(inout)          :: mask(:,:), maskObj(:,:)	! Mask to count objects
	   integer, intent(out)				:: perimeter, area, xcent, ycent
	   real, intent(out)				:: angle, aspect_ratio
	   ! Arrays to save objects attributes 
	   !integer, allocatable, intent(out)   		:: xcent(:), ycent(:), area(:), perimeter(:), area_hull(:)  
	   !real, allocatable, intent(out)	    	:: angle(:), aspect_ratio(:), complexity(:)
	   

	   integer				    	:: is_valid, yini, xini, n, ip, ip2, ys, xs, ps, xds, pds, yi, xi, pi, xdi, pdi
	   integer				    	:: leftPoint, rightPoint, upperPoint, lowerPoint 
	   integer				    	:: nsx, nsy  
	   integer					:: i, nhull	   
	   integer				    	:: xleft, xright, x, temp, area_hull_int=0
           
           type(mark), pointer 		    		:: m, maux, maux1, maux2, total_points, points_aux, total_points_aux 
	   type(mark), pointer		    		:: per_linked, total_linked
	   type(point), pointer				:: pts_per(:)
	   !

	   yini=row
	   xini=col

	   nullify(m)
	   nullify(points)
	   nullify(total_points)

	   area=0
	   nsx=0
	   nsy=0
	   perimeter=0

	   print*, '************  Inicio loop exterior  **********'
	   print*, ' Objetos detectados', totalObj
      
	  Do		
 
		call valid(restoreField, ysize, xsize, yini, xini, treshold, is_valid)
		If (is_valid) then		   
		   ys=yini-1
	      	   yi=yini+1
	      	   n=-1

		   do while (n < 2)		      
		      xs=0
		      xds=0
	              xi=0
		      xdi=0

	              if (n .LT. 0) then
	                 ip=xini+1
	              else
	                 ip=xini
		      endif
		      ip2 = ip
		      
		      do
			
			 ip2=ip2+n
			 
			 call valid(restoreField, ysize, xsize, yini, ip2, treshold, is_valid)			 
			 if ( is_valid .and. (mask(yini,ip2) .EQ. 0) ) then
	                    ! Creando linked list para guardar los ptos vecinos válidos que faltan por verificar
			    call valid(restoreField, ysize, xsize, ys, ip, treshold, is_valid)			    
			    if ( is_valid .and. (mask(ys,ip) .EQ. 0) ) then			
				pds=1
				if (xds .NE. pds) then			     	   
			     	   if (.not. associated(maux)) then
			              allocate(maux)
				      maux%x = ip    
			     	      maux%y = ys 
			     	      nullify(maux%next) 		              
				   else
				         allocate(m)				   
				         m%x = ip    
			     	         m%y = ys 
			     	         m%next => maux 
				         maux => m				      
				   endif       
			  	endif
			    else
				pds=0
			    endif
			    xds=pds

			    call valid(restoreField, ysize, xsize, yi, ip, treshold, is_valid)			    
		       	    if (is_valid .and. (mask(yi,ip) .EQ. 0) ) then	          
			        pdi=1
                                if (xdi .NE. pdi) then			            
			            if (.not. associated(maux)) then
			               allocate(maux)
				       maux%x = ip
			               maux%y = yi
			               nullify(maux%next)			            			            
				    else				       
				          allocate(m)				    
				          m%x = ip
				          m%y = yi
				          m%next => maux 
				          maux => m				       				
				    endif		
			        endif
			    else
				pdi=0
		            endif
			    xdi=pdi

			    call valid(restoreField, ysize, xsize, ys, ip2, treshold, is_valid)			    
			    if ( is_valid .and. (mask(ys,ip2) .EQ. 0) ) then		        
				ps=1
				if (xs .NE. ps) then			     	   
			     	   if (.not. associated(maux)) then
			              allocate(maux)
				      maux%x = ip2    
			     	      maux%y = ys 
			     	      nullify(maux%next) 		              
				   else
				         allocate(m)				   
				         m%x = ip2    
			     	         m%y = ys 
			     	         m%next => maux 
				         maux => m
				   endif       
			  	endif
			    else
				ps=0
			    endif
			    xs=ps

			    call valid(restoreField, ysize, xsize, yi, ip2, treshold, is_valid)		    
			    if ( is_valid .and. (mask(yi,ip2) .EQ. 0) ) then			        
				pi=1
				if (xi .NE. pi) then			     	   
			     	   if (.not. associated(maux)) then
			              allocate(maux)
				      maux%x = ip2    
			     	      maux%y = yi 
			     	      nullify(maux%next) 		              
				   else
				         allocate(m)				   
				         m%x = ip2    
			     	         m%y = yi 
			     	         m%next => maux 
				         maux => m
				   endif       
			  	endif
			    else
				pi=0
			    endif
			    xi=pi

			    mask(yini,ip2) = 1
			    maskObj(yini,ip2) = totalObj + 1

			    ! Verificando si el punto es frontera para calcular el perímetro del objeto
			    call valid(restoreField, ysize, xsize, ys, ip2, treshold, is_valid)
		            leftPoint=is_valid
			    !print*, 'leftPoint', leftPoint
		       
		            call valid(restoreField, ysize, xsize, yini, ip2-1, treshold, is_valid)
		            upperPoint=is_valid
			    !print*, 'upperPoint', upperPoint
			  
		            call valid(restoreField, ysize, xsize, yi, ip2, treshold, is_valid)
		            rightPoint=is_valid
			    !print*, 'rightPoint', rightPoint
			     
		            call valid(restoreField, ysize, xsize, yini, ip2+1, treshold, is_valid)
		            lowerPoint=is_valid
			    !print*, 'lowerPoint', lowerPoint
		            if (.NOT. (leftPoint .and. upperPoint .and. rightPoint .and. lowerPoint) ) then
				!print*, '.NOT. (leftPoint .and. upperPoint .and. rightPoint .and. lowerPoint)'
			        perimeter = perimeter + 1
		 		! creando linked list para guardar los puntos del perímetro
				allocate(points_aux)
				points_aux%next => points
				points => points_aux
				points%x = ip2
				points%y = yini
		            endif

		            area=area+1
	                    nsx=nsx+ip2   !sumando las posiciones x
	                    nsy=nsy+yini  !sumando las posiciones y

			    ! creando linked list para guardar todos los puntos del objeto
			    allocate(total_points_aux)
			    total_points_aux%next => total_points
			    total_points => total_points_aux
			    total_points%x = ip2
			    total_points%y = yini
			    
			 else			    
			    exit
			 endif			 
		      enddo
		      n = n + 2
		      
		   enddo
		Endif		

		if (associated(maux)) then
		   m => maux
	           xini = maux%x
	           yini = maux%y
	           maux => maux%next		   
		   deallocate(m)		   
        	else
		   nullify(m)	
		   print*, 'perimetro', perimeter, 'area', area, 'nsx', nsx, 'nsy', nsy           
          	   exit
		Endif
	   Enddo

	   xcent = nsx/area
	   ycent = nsy/area

	   ! Object Orientation Angle 
	   call object_angle(xcent, ycent, perimeter, pts_per, angle)
	   print*, 'orientation angle', angle

	   call object_Aspect_Ratio(pts_per, xcent, ycent, perimeter, angle, aspect_ratio)
	   print*, 'Aspect Ratio', aspect_ratio

           return	   

	End Subroutine


    !**************************************************************************************************************************************
     !The subroutine calculates the orientation of the objects
      Subroutine object_angle(xcent, ycent, perimeter, pts_per, angle)        

	Implicit None
	! INPUT PARAMETERS:
	!type(mark), pointer, intent(in)		:: points
        integer, intent(in)			:: xcent, ycent, perimeter	
	!OUTPUT PARAMETERS:
	type(point), pointer, intent(out)	:: pts_per(:)
	real, intent(out)			:: angle

	!type(point), pointer	:: pts_per(:)
	integer			:: i, j
	real, parameter   	:: pi=3.141592654
	real			:: sumUp, sumDown

	sumUp=0.0
	sumDown=0.0

	allocate(pts_per(perimeter))

	Do while (associated(points))
	   do i=1, perimeter
	      pts_per(i)%x = points%x
	      pts_per(i)%y = points%y

	      sumUp = sumUp + (pts_per(i)%x - xcent) * (pts_per(i)%y - ycent)
	      sumDown = sumDown + (pts_per(i)%x - xcent)**2 - (pts_per(i)%y - ycent)**2
	
	      points => points%next
	   enddo
        Enddo

	angle = 0.5*ATAN2(2*sumUp,sumDown)

	if (angle .LT. 0.0) then
	   angle = angle + 2*pi
        else if (angle .GT. pi) then
	   angle = angle - pi
	endif

      End Subroutine 
    !**************************************************************************************************************************************




    !**************************************************************************************************************************************
     !The subroutine calculates the orientation of the objects
      Subroutine object_Aspect_Ratio(pts_pos, xcent, ycent, perimeter, angle, aspect_ratio)
      
	Implicit None
	! INPUT PARAMETERS:
	type(point), pointer, intent(in)	:: pts_pos(:)
        integer, intent(in)	:: xcent, ycent, perimeter
	real, intent(in)	:: angle
	!OUTPUT PARAMETERS:
	real, intent(out)	:: aspect_ratio

	integer			:: i
	real			:: dist, vmajor, major_up, major_down, vminor, minor_down, minor_up, major_axis, minor_axis 
	real, parameter   	:: pi=3.141592654

	Do i=1, perimeter

	   dist =  (pts_pos(i)%x - xcent)*SIN(angle) - (pts_pos(i)%y - ycent)*COS(angle)
	   vmajor = (pts_pos(i)%y - ycent) - TAN(angle)*(pts_pos(i)%x - xcent)

	   If (vmajor .GT. 0) then
	      if(dist .GT. major_up) then
		 major_up = dist
	      endif
	   Endif

	   If (vmajor .LT. 0) then
	      if(dist .GT. major_down) then
		 major_down = dist
	      endif
	   Endif

	   dist = (pts_pos(i)%x - xcent)*SIN(angle + pi/2) - (pts_pos(i)%y - ycent)*COS(angle + pi/2)
	   vminor = (pts_pos(i)%y - ycent) + (pts_pos(i)%x - xcent)/TAN(angle)

	   If (vminor .GT. 0) then
	      if(dist .GT. minor_up) then
		 minor_up = dist
	      endif
	   Endif

	   If (vminor .LT. 0) then
	      if(dist .GT. minor_down) then
		 minor_down = dist
	      endif
	   Endif

	Enddo

	major_axis = major_up + major_down
	minor_axis = minor_up + minor_down

	If (major_axis .LT. minor_axis) then
	   aspect_ratio = major_axis/minor_axis
	else 
	   aspect_ratio = minor_axis/major_axis
	Endif

      End Subroutine 
    !**************************************************************************************************************************************



    !**************************************************************************************************************************************
      ! Quickhull algorithm to find the convex hulls of the objects
      Subroutine quick_hull(array_pts, npoints, hull)
	! pts_per e hull
         Implicit None
	 type(point), pointer, intent(in)	:: array_pts(:)
	 integer, intent(in)			:: npoints

	 type(mark), pointer, intent(out)	::hull
	 type(point), pointer			:: A, B, temp, left_set(:), right_set(:)
	 type(mark), pointer			:: aux, hull_temp
	 integer				:: i, j, k, positions, nleft, nright, xmin, ymin, xmax, ymax, minpos, maxpos, side

	 nleft=0
	 nright=0
	 xmin=9999
	 xmax=-1

	 nullify(hull_temp)

	 !Loop para buscar los xmin y xmax
	 Do i=1, npoints	    
	    If (array_pts(i+1)%x .LT. xmin) then
	       xmin = array_pts(i+1)%x
	       ymin = array_pts(i+1)%y
	       minpos = i+1
	    Endif

	    If (array_pts(i+1)%x .GT. xmax) then
	       xmax = array_pts(i+1)%x
	       ymax = array_pts(i+1)%y
	       maxpos = i+1
	    Endif
	 Enddo

	 ! pts extremos de las x
	 A = array_pts(minpos)
	 B = array_pts(maxpos)

	 temp = array_pts(minpos)
	 array_pts(minpos) = array_pts(1)
	 array_pts(1) = temp
	 temp = array_pts(maxpos)
	 array_pts(maxpos) = array_pts(2)
         array_pts(2) = temp

	 ! guardo A y B como parte de la convex hull
	 allocate(aux)
	 aux%next => hull_temp
	 hull_temp => aux
	 hull_temp%x = A%x;
         hull_temp%y = A%y;

	 allocate(aux)
	 aux%next => hull_temp
	 hull_temp => aux
	 hull_temp%x = B%x;
         hull_temp%y = B%y;
	 
	 ! hallo los conjuntos a ambos lados del segmento AB
	 Do i=3, npoints	    
	    positions = (B%x - A%x)*(array_pts(i+1)%y - A%y) - (B%y - A%y)*(array_pts(i+1)%x - A%x)
	   
	    if (positions .GT. 0) then  
	       nleft = nleft + 1
	    endif

	    if (positions .LT. 0) then  
	       nright = nright + 1
	    endif
	 Enddo

	 allocate(left_set(nleft))
	 allocate(right_set(nright))

	 j=1
	 k=1

	 Do i=3, npoints	    
	    positions = (B%x - A%x)*(array_pts(i+1)%y - A%y) - (B%y - A%y)*(array_pts(i+1)%x - A%x)

	    if (positions .GT. 0) then  
	       left_set(j) = array_pts(i+1)
	       j=j+1
	    endif

	    if (positions .LT. 0) then  
	       right_set(k) = array_pts(i+1)
	       k=k+1
	    endif

	 Enddo

      End Subroutine
    !**************************************************************************************************************************************


!// llamo a la funcion que determina los pts q pertenecen a la convex hull
!   Hull_Set(A, B, right_set, nright, &hull_temp);
!   Hull_Set(A, B, left_set, nleft, &hull_temp);

!  *hull = hull_temp; 
!}
    
    !**************************************************************************************************************************************
      Subroutine Hull_Set(A, B, set, array_size, hull_temp)

	 Implicit None

	 type(point), pointer, intent(in)	:: A, B, set(:)
	 integer, intent(in)			:: array_size

	 type(mark), pointer, intent(out)	:: hull_temp

	 integer				:: dist, distmax, furthest_point, nAP, nPB, i, j, k
	 type(mark),pointer			:: aux, hull_temp2

	 distmax = -1
	 nAP = 0
	 nPB = 0

	 hull_temp2 => hull_temp

	 ! compruebo el tamaño del arreglo
	 if (array_size .EQ. 0) then
	    return
	 endif

	 if (array_size .EQ. 1) then
	    allocate(aux)
	    aux%next => hull_temp2
	    hull_temp2 => aux
	    hull_temp2%x = set(1)%x
	    hull_temp2%y = set(1)%y
	    return
	 endif

         ! determinar el pto mas lejano
	 
	 	  


      End Subroutine
    !**************************************************************************************************************************************



END PROGRAM pru
