module heap_sort_mod
  ! Dowloaded from: http://www.fortran.com/F/adsff.html
  ! An explanation of the algorithm:
  ! http://www.personal.kent.edu/~rmuhamma/Algorithms/MyAlgorithms/Sorting/heapSort.htm
  implicit  none
  !  public :: HEAP_SORT,HEAP_sort_integer,HEAP_sort_real
  !  private :: MAKE_HEAP_integer,MAKE_HEAP_real
  public :: HEAP_SORT
  private :: MAKE_HEAP_integer,MAKE_HEAP_real,HEAP_sort_integer,HEAP_sort_real
  interface MAKE_HEAP
     module procedure MAKE_HEAP_integer
     module procedure MAKE_HEAP_real
  end interface
  interface HEAP_sort
     module procedure HEAP_sort_integer
     module procedure HEAP_sort_real
  end interface
contains
  subroutine MAKE_HEAP_integer (A, Apex, No_Elements)
    integer, dimension ( : ), intent (in out) :: A
    integer, intent (in)            :: Apex, No_Elements
    integer                         :: Temp, node
    Temp = A(Apex); Node = 2*Apex
    do
       if (Node > No_Elements) then; exit; end  if
       if (Node < No_Elements) then; if (A(Node) < A(Node+1)) then; Node = Node + 1; end  if; end  if
       if (Temp >= A(Node) ) then; exit; end  if
       A(Node/2) = A(Node); Node = 2*Node
    end  do
    A(Node/2) = Temp
  end  subroutine  MAKE_HEAP_integer
  subroutine MAKE_HEAP_real (A, Apex, No_Elements)
    real(kind(1d0)), dimension ( : ), intent (in out) :: A
    integer, intent (in)                              :: Apex, No_Elements
    real(kind(1d0))                                   :: Temp
    integer                                           :: node
    Temp = A(Apex); Node = 2*Apex
    do
       if (Node > No_Elements) then; exit; end  if
       if (Node < No_Elements) then; if (A(Node) < A(Node+1)) then; Node = Node + 1; end  if; end  if
       if (Temp >= A(Node) ) then; exit; end  if
       A(Node/2) = A(Node); Node = 2*Node
    end  do
    A(Node/2) = Temp
  end  subroutine  MAKE_HEAP_real
  subroutine HEAP_SORT_integer (A)
    !! This procedure performs the overall tasks of (a) initial setup by organizing the incoming
    !! array A as a pyramid obeying the heap rule, and (b) overseeing the extraction of successive
    !! next-largest elements from the pyramid, and writing them back in order into the array A.
    !! INCOMING:   A = the array whose elements are to be sorted.
    !! OUTGOING:   A = array of sorted elements.  Elements are in ascending order.
    integer, dimension ( : ), intent (in out) :: A
    integer                                   :: Temp        !! The same type as A.
    integer                                   :: K, N, Apex
    N = ubound(A, 1);
    do Apex = N/2, 1, -1; call MAKE_HEAP (A, Apex, N); end  do
    do K=N,2,-1; Temp = A(K);  A(K) = A(1); A(1) = Temp; call MAKE_HEAP (A, 1, K-1); end  do
  end subroutine  HEAP_SORT_integer
  subroutine HEAP_SORT_real (A)
    !! INCOMING:   A = the array whose elements are to be sorted.
    !! OUTGOING:   A = array of sorted elements.  Elements are in ascending order.
    real(kind(1d0)), dimension ( : ), intent (in out) :: A
    real(kind(1d0))                                   :: Temp        !! The same type as A.
    integer                                           :: K, N, Apex
    N = ubound(A, 1);
    do Apex = N/2, 1, -1; call MAKE_HEAP (A, Apex, N); end  do
    do K=N,2,-1; Temp = A(K);  A(K) = A(1); A(1) = Temp; call MAKE_HEAP (A, 1, K-1); end  do
  end subroutine  HEAP_SORT_real

  subroutine HEAP_SORT_ABS_real (A)
    !! INCOMING:   A = the array whose elements are to be sorted.
    !! OUTGOING:   A = array of sorted elements.  Elements are in ascending order.
    real(kind(1d0)), dimension ( : ), intent (in out) :: A
    real(kind(1d0))                                   :: Temp        !! The same type as A.
    integer                                           :: K, N, Apex
    N = ubound(A, 1);
    do Apex = N/2, 1, -1; call MAKE_HEAP_ABS_real (A, Apex, N); end  do
    do K=N,2,-1; Temp = A(K);  A(K) = A(1); A(1) = Temp; call MAKE_HEAP_ABS_real (A, 1, K-1); end  do
  end subroutine  HEAP_SORT_ABS_real





  subroutine MAKE_HEAP_ABS_real (A, Apex, No_Elements)
    real(kind(1d0)), dimension ( : ), intent (in out) :: A
    integer, intent (in)                              :: Apex, No_Elements
    real(kind(1d0))                                   :: Temp
    integer                                           :: node
    Temp = A(Apex); Node = 2*Apex
    do
       if (Node > No_Elements) then; exit; end  if
       if (Node < No_Elements) then; if (abs(A(Node)) < abs(A(Node+1))) then; Node = Node + 1; end  if; end  if
       if (abs(Temp) >= abs(A(Node)) ) then; exit; end  if
       A(Node/2) = A(Node); Node = 2*Node
    end  do
    A(Node/2) = Temp
  end  subroutine  MAKE_HEAP_ABS_real

  
end  module  heap_sort_mod
