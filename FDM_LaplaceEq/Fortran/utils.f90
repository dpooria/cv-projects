
module utils

   integer, parameter:: dp=kind(0.d0) ! double precision
   real(dp), parameter :: pi=3.14159265358979_dp
   private
   public LInfNorm, L2Norm, L1Norm, write_mat, print_mat, dp, pi
contains
   function L2Norm(T1, T2, m) result(e)
      integer, intent(in) :: m
      real(dp), intent(in) :: T1(m, m), T2(m, m)
      real(dp) :: e
      e = sum((T2 - T1)**2)
   end function

   function L1Norm(T1, T2, m) result(e)
      integer, intent(in) :: m
      real(dp), intent(in) :: T1(m, m), T2(m, m)
      real(dp) :: e
      e = sum(abs(T2 - T1))
   end function
   
   function LInfNorm(T1, T2, m) result(e)
      integer, intent(in) :: m
      real(dp), intent(in) :: T1(m, m), T2(m, m)
      real(dp) :: e
      e = maxval(abs(T2 - T1))
   end function
   


   subroutine print_mat(A)
      real(dp), intent(in) :: A(:,:)
      integer :: n, i
      n = size(A, 1)
      do i = 1, n
         print *, A(i, :)
      enddo
   end subroutine


   subroutine write_mat(A, fname)
      real(dp), intent(in) :: A(:,:)
      character(len=*), intent(in)::fname
      integer::u, n, i
      open(newunit=u, file=fname, status='replace')
      n = size(A, 1)
      do i = 1, n
         write(u, *) A(i, :)
      enddo
      close(u)
   end subroutine
end module utils
