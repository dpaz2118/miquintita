! ifort nrtype.f90 parametros.f90 nrutil_DP.f90 nr_DP.f90 rkck_DP.f90 rkqs_DP.f90 odeint_DP.f90 oscila2.f90 -o i.x
!Hola
program oscila2
!programa principal que usa el runge-kutta del fortran
!El x es la variable independiente osea el tiempo (tau)

  use nrtype, only : DP    
  use nr,     only : odeint 
  use ode_path   ! el modulo donde guarda la solucion la subrutina odeint
  use parametros

  implicit none
  interface
     subroutine derivs(x,y,dydx)
       use nrtype
       real(DP), intent(in) :: x
       real(DP), dimension(:), intent(in) :: y
       real(DP), dimension(:), intent(out) :: dydx
     end subroutine derivs
     subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       use nrtype
       implicit none
       real(DP), dimension(:), intent(inout) :: y
       real(DP), dimension(:), intent(in) :: dydx,yscal
       real(DP), intent(inout) :: x
       real(DP), intent(in) :: htry,eps
       real(DP), intent(out) :: hdid,hnext
       interface
          subroutine derivs(x,y,dydx)
            use nrtype
            implicit none
            real(DP), intent(in) :: x
            real(DP), dimension(:), intent(in) :: y
            real(DP), dimension(:), intent(out) :: dydx
          end subroutine derivs
       end interface
     end subroutine rkqs
  end interface
  !integer               :: j,nsteps
  integer               :: j,n,i
  real(DP),dimension(2) :: ystart
  !ystart es la condicion inicial (pos, vel)
  !es de dimension 2 pq resuelvo linealizada de una
  !ODE de orden 2
  real(DP)              :: tinicial=0._dp , tfinal=6000._dp
  !los tiempos entre los que integro
  real(DP)              :: tolerancia=1.e-8 !la tolerancia
  real(DP)              :: guess_step, step_min, onda, u
  character(len=100)    :: filename_1


  filename_1='oscilador.dat'
  masa=1.
  F0=0.
  kk=1.
  amor=0.0

  !guess_step es un ancho inicial de paso que luego el programa modifica
  !step_min es el ancho minimo que le permitimos al programa achicar

  print*,'comienzo...'

  guess_step=(tfinal-tinicial)/60000.    !define el primer paso
  step_min=0._dp                         !define el minimo

  !llama al integrador dandole la CI, el intervalo de tiempo, 
  !la tolerancia, los anchos de pasos, y las subrutinas que necesita para
  !integrar
  
  ystart=(/0.65,0./)

  call odeint(ystart,tinicial,tfinal,tolerancia, &
            &    guess_step,step_min,derivs,rkqs)

  !una vez que termino de integrar, el modulo ode_path contiene ya la
  !informacion almacenada entonces podemos preguntar cuantos pasos tiene 
  !la solucion (atravez size(xp)) y podemos pedir la solucion misma (xp 
  !es un vector de los pasos que contiene los tiempos, yp es la solucion
  !para los pasos)

  !Escribimos la salida
  open(111,file=filename_1,status='unknown',form='formatted')
  !archivo de salida

  !do j=1,int(0.2_dp*utp/50._dp)
  do j=1,utp
      write(111,'(4(e16.5,1x))')xp(j),yp(1,j),yp(2,j),ystart(1)
  enddo 

  close(111)
  print*,'numero de steps=',utp


  !liberamos la memoria
  deallocate(xp)
  deallocate(yp)

  print*,'termine'

end program oscila2

!Subrutina que calcula el lado derecho de la ecuacion diff
!linealizada osea la derivad de y con respecto a x que llamamos dydx
subroutine derivs(x,y,dydx)
    use nrtype
    use parametros

    implicit none
    real(DP), intent(in)                :: x
    real(DP), dimension(:), intent(in)  :: y
    real(DP), dimension(:), intent(out) :: dydx


    dydx(1) = y(2)/masa
    dydx(2) = F0 - kk*y(1) - amor*y(2)

 
end subroutine derivs
