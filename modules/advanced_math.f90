module advanced_math
  use math, only: pi,RandomInteger
  implicit none
contains
  subroutine get_angulo_de_sincos_array3(senos,cosenos,angulos,errores)
    use heap_sort_mod, only:HEAP_SORT_ABS_real
    real(kind(1d0)), intent(in)       :: senos(:), cosenos(:)
    real(kind(1d0)), intent(out)      :: angulos(:), errores(:)
    real(kind(1d0)), allocatable      :: tsin(:),tcos(:),errors(:)
    real(kind(1d0))                   :: coseno_buscado,seno_bueno
    integer                           :: pos,pos_p,s,s_iterator
    s=size(senos);print*,"tamano del problema:",s
    if(((s.ne.size(cosenos)).or.(s.ne.size(errores))).or.(s.ne.size(angulos))) then
       print*,"algo mal en get_angulo_de_sincos_array"; stop
    end if
    angulos=0d0; errores=0d0
    allocate(tsin(s),tcos(s),errors(s))
    !do j_1=-10,10; print*,modulo(j_1,4);enddo;stop
    tcos=cosenos;tsin=senos;
    call HEAP_SORT_ABS_real(tcos);
    !print'(A10,99F9.4)',"cosenos:",tcos;print'(A10,99F9.4)',"senos:",tsin;print*
    !print*;print*,tcos**2-cshift(tcos**2,1)<0d0
    s_iterator=s;s_iterator=15;s_iterator=randominteger(s)+1
    do s_iterator=1,s
    !do s_iterator=2,2
      !print*,"iteracion ",s_iterator
      seno_bueno=tsin(s_iterator);!seno_bueno=0d0
      coseno_buscado=sqrt(1-seno_bueno**2)
      !coseno_buscado=0.1d0
      !posicion_tentativa=floor(2*s*acos(abs(seno_bueno))/pi)
      pos  =mymod(floor(s*cumulative_distribution(seno_bueno)),s); 
      if (abs(coseno_buscado)>abs(tcos(s)).or.abs(coseno_buscado)<abs(tcos(1))) then; pos=s;
      else
        !print*,"en algun pedo1:",tcos(pos)**2.,coseno_buscado**2
        do while (tcos(pos)**2.gt.coseno_buscado**2) ;pos=mymod(pos-1,s);end do
        !print*,"en algun pedo2:"
        do while (tcos(mymod(pos+1,s))**2.le.coseno_buscado**2); pos=mymod(pos+1,s); end do
        !print*,"en algun pedo3:"
      end if 
      pos_p=mymod(pos+1,s)
      ! aca tengo algo mal.. croe que puedo tener pos_p-pos diferente de uno y
      ! de s, de alguna forma tengo que ponerlos a buscar juntos
      !print*,"posiciones        ",pos,pos_p; 
      !print*,"cosenos alrededor ",tcos(pos),tcos(pos_p); 
      !print*,"errores respectivos ",tcos(pos),tcos(pos_p); 
      !print*,"coseno interesado ",coseno_buscado;print*
      if (abs(abs(tcos(pos))-abs(coseno_buscado))>abs(abs(tcos(pos_p))-abs(coseno_buscado))) then; pos=pos_p; endif 
      !print*,"posicion del bueno:",pos
      errors(s_iterator)=abs(1-tcos(pos)**2-seno_bueno**2)
      if(errors(s_iterator)>10d-10) print*,"error:",errors(s_iterator)
      !print*,"error:",j_1,errors(s_iterator)
      !print*
    enddo
    print*,count(errors>10d-4);print*,count(errors>10d-5); print*,count(errors>10d-6); print*,count(errors>10d-7); print*,count(errors>10d-8);
    print*,"unfinished, try version 2 in math package...";stop
    deallocate(tsin,tcos,errors)
    return
  contains
    function mymod(n,big_n)
      integer             :: mymod
      integer, intent(in) ::n,big_n
      mymod=modulo(n-1,big_n)+1
    end function mymod
    function cumulative_distribution(x)
      real(kind(1d0))        :: cumulative_distribution,x
      cumulative_distribution=2*acos(abs(x))/pi
    end function cumulative_distribution
  end subroutine get_angulo_de_sincos_array3
end module advanced_math

