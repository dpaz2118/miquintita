MODULE ode_path
	USE nrtype
	INTEGER(I4B) :: nok,nbad,kount
	!LOGICAL(LGT), SAVE :: save_steps=.false.
	LOGICAL(LGT), SAVE :: save_steps=.true.
	REAL(DP) :: dxsav
	REAL(DP), DIMENSION(:), POINTER :: xp
        REAL(DP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

	SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE ode_path
        !modificado
        use parametros, only : utp
        !fin modificado
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
	REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
                use parametros
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(DP), INTENT(INOUT) :: x
		REAL(DP), INTENT(IN) :: htry,eps
		REAL(DP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
                        use parametros
			IMPLICIT NONE
			REAL(DP), INTENT(IN) :: x
			REAL(DP), DIMENSION(:), INTENT(IN) :: y
			REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkqs
	END INTERFACE
	REAL(DP), PARAMETER :: TINY=1.0e-30_dp
	INTEGER(I4B), PARAMETER :: MAXSTP=10000000
	INTEGER(I4B) :: nstp
	REAL(DP) :: h,hdid,hnext,x,xsav
	REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal
	x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
	if (save_steps) then
		xsav=x-2.0_dp*dxsav
		nullify(xp,yp)
		allocate(xp(256))
		allocate(yp(size(ystart),size(xp)))
	end if
	do nstp=1,MAXSTP
		call derivs(x,y,dydx)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
		!		print*,'hola'
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

		yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
		if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
			call save_a_step
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
		!print*,'hola2'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
		call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		if ((x-x2)*(x2-x1) >= 0.0) then
			ystart(:)=y(:)
			if (save_steps) call save_a_step
                        utp=nstp
			RETURN
		end if
		if (abs(hnext) < hmin)&
			call nrerror('stepsize smaller than minimum in odeint')
		h=hnext
	end do
	call nrerror('too many steps in odeint')
	CONTAINS
!BL
	SUBROUTINE save_a_step
	kount=kount+1
	if (kount > size(xp)) then
		xp=>reallocate(xp,2*size(xp))
		yp=>reallocate(yp,size(yp,1),size(xp))
                print*,'reallocateando..',size(xp)
	end if
	xp(kount)=x
	yp(:,kount)=y(:)
	xsav=x
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	!print*,'hola4'
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	END SUBROUTINE save_a_step
	END SUBROUTINE odeint
