module exact_ps
  implicit none
  interface w_goe    
    module procedure w_goe_single
    module procedure w_goe_double
  end interface
  interface wgoe
    function wgoe(s)
      real wgoe,s
    end function
  end interface
contains
  subroutine w_goe_single(points,ws)
    implicit none
    real, intent(in)      :: points(:)
    real, intent(out)     :: ws(:)
    integer               :: j_1,size_problem
     size_problem=size(points);
    if (size_problem.ne.size(ws)) then 
      print*,"error marica en w_goe_single",size_problem,size(ws)
      stop
    endif
    call initGOE;ws=0d0;ws=wgoe(0.5)
    do j_1=1,size_problem;ws(j_1)=wgoe(points(j_1));enddo
  end subroutine w_goe_single

  subroutine w_goe_double(points,ws)
    implicit none
    real(kind(1d0)), intent(in)      :: points(:)
    real(kind(1d0)), intent(out)     :: ws(:)
    integer                          :: j_1,size_problem
    size_problem=size(points);
    if (size_problem.ne.size(ws)) then 
      print*,"error marica en w_goe_single",size_problem,size(ws)
      stop
    endif
    call initGOE;ws=0d0
    do j_1=1,size_problem;ws(j_1)=wgoe(real(points(j_1),kind(1.)));enddo
  end subroutine w_goe_double

end module exact_ps

