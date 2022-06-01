program cfd1proj1laplace

   use utils
   implicit none

   real(dp), parameter:: L2_1=1e-10_dp, L2_2=1e-12_dp, w1=1.1_dp, w2=1.3_dp, w3=1.5_dp
   integer, parameter :: M1=10, M2=20, M4=40

   real(dp) :: tic, toc, err_pgs_l2, err_lgs_l2, err_pj_l2, err_lj_l2, err_sor_l2 &
      , err_pgs_l1, err_lgs_l1, err_pj_l1, err_lj_l1, err_sor_l1 &
      , err_pgs_linf, err_lgs_linf, err_pj_linf, err_lj_linf, err_sor_linf &
      , t_pgs, t_lgs, t_pj, t_lj, t_sor, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   
   real(dp) :: T1(M1, M1), T1_implicit(M1, M1), T2_implicit(M2, M2),&
      T2(M2, M2), T4(M4, M4), T1_exact(M1, M1), &
      T2_exact(M2, M2), T4_exact(M4, M4), T4_implicit(M4, M4)

   integer ::  u, u6,iter_pgs, iter_lgs, iter_pj, iter_lj, iter_sor

   !Part 4: Run your program for Mesh-1 10*10 using Gauss-Seidel routine using initial condition zero everywhere and let it converge to L2=1e-10
   call BoundaryCondition(T1)
   call ExactSolution(T1_exact)
   call cpu_time(tic)
   call PointGaussSeidel(T1, '4-M1_pgs_l21e-10_converge_history.txt', iter_pgs, L2_1)
   call cpu_time(toc)
   t_pgs = toc - tic
   err_pgs_l2 = L2Norm(T1_exact, T1, M1)
   print *, '************************************************************'
   print *, 'Part 4:'
   print *, 'Point Gauss-Seidel method results for 10x10 mesh (L2 = 1e-10):'
   print *, 'l2=', err_pgs_l2, 'iter_pgs=', iter_pgs, 't_pgs=', t_pgs
   open(newunit=u, file='part4.txt', status='replace')
   write(u, *) err_pgs_l2, t_pgs, iter_pgs
   call write_mat(T1_exact, '4-T1_exactsolution.txt')
   call write_mat(T1, '4-T1_gaussseidel.txt')
   close(u)
   !Part 5: Run your program for Mesh-2: 20*20 up to Mesh-4: 40*40 using Gauss-Seidel and
   !let them converge to L2=1e-12. Report the convergence history (both in terms of CPU
   ! time and iteration number) for M-1->M-4.
   !============
   !Part 6: Verify the accuracy of your solution against the exact solution and show your code
   !is a 2nd-order code by plotting the error norms versus mesh. Report L1, L2 and Linf
   !for M1-> M4. Now, manipulate your the boundary condition subroutine in such a
   !way that for only one boundary cell of your choice the boundary condition is
   !implemented by a first-order discretization. Do the accuracy analysis again and draw
   !a conclusion about your result. Visualize the solution in both cases in different plots
   !and compare them together.

   print *, '************************************************************'
   print *, 'Part 5 and 6:'

   open(newunit=u, file='part5.txt', status='replace')
   open(newunit=u6, file='part6.txt', status='replace')

   !M1
   call BoundaryCondition(T1)
   ! call ExactSolution(T1_exact)
   call cpu_time(tic)
   call PointGaussSeidel(T1, '5-M1_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   call write_mat(T1, '5-T1_pgs_l21e-12_result.txt')
   call write_mat(T1_exact, '5-T1_exactsolution.txt')
   print *, 'Point Gauss-Seidel method results for 10x10 mesh with upwind bc(L2 = 1e-12):'
   write(*, *) M1, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u, *) M1, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M1, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs

   call FirstOrderUpwindBC(T1, 0.0_dp)
   ! call ExactSolution(T1_exact)
   call cpu_time(tic)
   call PointGaussSeidel(T1, '6-upwind_M1_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   call write_mat(T1, '6-upwind_T1_pgs_l21e-12_result_upwind.txt')
   t_pgs = toc - tic
   err_pgs_l1 = L1Norm(T1_exact, T1, M1)
   err_pgs_l2 = L2Norm(T1_exact, T1, M1)
   err_pgs_linf = LInfNorm(T1_exact, T1, M1)
   print *, 'Point Gauss-Seidel method results for 10x10 mesh with upwind bc(L2 = 1e-12):'
   write(*, *) M1, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M1, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs

   !M2
   call BoundaryCondition(T2)
   call ExactSolution(T2_exact)
   call cpu_time(tic)
   call PointGaussSeidel(T2, '5-M2_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   t_pgs = toc - tic
   err_pgs_l2 = L2Norm(T2_exact, T2, M2)
   err_pgs_l1 = L1Norm(T2_exact, T2, M2)
   err_pgs_linf = LInfNorm(T2_exact, T2, M2)
   print *, 'Point Gauss-Seidel method results for 20x20 mesh (L2 = 1e-12):'
   write(*, *) M2, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u, *) M2, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M2, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs

   call write_mat(T2, '5-T2_pgs_l21e-12_result.txt')
   call write_mat(T2_exact, '5-T2_exactsolution.txt')

   call FirstOrderUpwindBC(T2, 0.0_dp)
   call cpu_time(tic)
   call PointGaussSeidel(T2, '6-upwind_M2_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   t_pgs = toc - tic
   err_pgs_l2 = L2Norm(T2_exact, T2, M2)
   err_pgs_l1 = L1Norm(T2_exact, T2, M2)
   err_pgs_linf = LInfNorm(T2_exact, T2, M2)
   print *, 'Point Gauss-Seidel method results for 20x20 mesh with upwind bc(L2 = 1e-12):'
   write(*, *) M2, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M2, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   call write_mat(T2, '6-upwind_T2_pgs_l21e-12_result.txt')

   !M4
   call BoundaryCondition(T4)
   call ExactSolution(T4_exact)
   call cpu_time(tic)
   call PointGaussSeidel(T4, '5-M4_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   t_pgs = toc - tic
   err_pgs_l2 = L2Norm(T4_exact, T4, M4)
   err_pgs_l1 = L1Norm(T4_exact, T4, M4)
   err_pgs_linf = LInfNorm(T4_exact, T4, M4)
   print *, 'Point Gauss-Seidel method results for 40x40 mesh (L2 = 1e-12):'
   write(*, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   
   call write_mat(T4, '5-T4_pgs_l21e-12_result.txt')
   call write_mat(T4_exact, '5-T4_exactsolution.txt')
   close(u)

   call FirstOrderUpwindBC(T4, 0.0_dp)
   call cpu_time(tic)
   call PointGaussSeidel(T4, '6-upwind_M4_pgs_l21e-12_converge_history.txt', iter_pgs, L2_2)
   call cpu_time(toc)
   t_pgs = toc - tic
   err_pgs_l2 = L2Norm(T4_exact, T4, M4)
   err_pgs_l1 = L1Norm(T4_exact, T4, M4)
   err_pgs_linf = LInfNorm(T4_exact, T4, M4)
   print *, 'Point Gauss-Seidel method results for 40x40 mesh with upwind bc(L2 = 1e-12):'
   write(*, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   write(u6, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   call write_mat(T4, '6-upwind_T4_pgs_l21e-12_result.txt')
   close(u6)



   !================
   !Part 7:
   !Use SOR for Mesh-4 with w=1.1, 1.3, 1.5 as over relaxation factors and compare
   !your convergence history with Jacobi and Gauss-Seidel (both in terms of CPU time
   !and iteration number). What is your preferred choice for convergence method and
   !why it is so?
   print *, '************************************************************'
   print *, 'Part 7:'
   open(newunit=u, file="part7.txt", status="replace")

   !gauss seidel
   write(u, *) M4, err_pgs_l2, err_pgs_l1, err_pgs_linf, t_pgs, iter_pgs
   !jacobi
   call BoundaryCondition(T4)
   call cpu_time(tic)
   call PointJacobi(T4, '7-M4_jacobi_l21e-12_converge_history.txt', iter_pj, L2_2)
   call cpu_time(toc)
   t_pj = toc - tic
   err_pj_l2 = L2Norm(T4_exact, T4, M4)
   err_pj_l1 = L1Norm(T4_exact, T4, M4)
   err_pj_linf = LInfNorm(T4_exact, T4, M4)
   print *, 'Point Jacobi method results for 40x40 mesh (L2 = 1e-12):'
   write(*, *) M4, 1000, err_pj_l2, err_pj_l1, err_pgs_linf, t_pj, iter_pj
   write(u, *) M4, err_pj_l2, err_pj_l1, err_pj_linf, t_pj, iter_pj
   call write_mat(T4, '7-T4_jacobi_l21e-12_result.txt')
   !w1
   call BoundaryCondition(T4)
   call cpu_time(tic)
   call PointSOR(T4, w1, '7-w1_M4_sor_l21e-12_converge_history.txt', iter_sor, L2_2)
   call cpu_time(toc)
   t_sor= toc - tic
   err_sor_l2 = L2Norm(T4_exact, T4, M4)
   err_sor_l1 = L1Norm(T4_exact, T4, M4)
   err_sor_linf = LInfNorm(T4_exact, T4, M4)
   print *, 'Point SOR method results for 40x40 mesh and w1 = 1.1 (L2 = 1e-12):'
   write(*, *) M4, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M4, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T4, '7-T4_sor_w1_l21e-12_result.txt')
   !w2
   call BoundaryCondition(T4)
   call cpu_time(tic)
   call PointSOR(T4, w2, "7-w2_M4_sor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T4_exact, T4, M4)
   err_sor_l1 = L1Norm(T4_exact, T4, M4)
   err_sor_linf = LInfNorm(T4_exact, T4, M4)
   print *, "Point SOR method results for 40x40 mesh and w2 = 1.3 (L2 = 1e-12):"
   write(*, *) M4, w2, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M4, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T4, "7-T4_sor_w2_l21e-12_result.txt")
   !w3
   call BoundaryCondition(T4)
   call cpu_time(tic)
   call PointSOR(T4, w3, "7-w3_M4_sor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T4_exact, T4, M4)
   err_sor_l1 = L1Norm(T4_exact, T4, M4)
   err_sor_linf = LInfNorm(T4_exact, T4, M4)
   print *, "Point SOR method results for 40x40 mesh and w3 = 1.5 (L2 = 1e-12):"
   write(*, *) M4, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M4, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T4, "7-T4_sor_w3_l21e-12_result.txt")

   close(u)

   print *, '************************************************************'
   print *, "Part 8-1:"
   !Part 8-1: Solve the system for line SOR with w=1.1 for M1, M2 & M4. You can use Tri-
   !Diagonal solver for line iteration. Make sure that your solution still is correct and is
   !similar to the explicit solver result for the same mesh. Report the CPU-time and
   !compare it with the CPU-time of explicit SOR for the same mesh and over relaxation
   !factor on one graph. Draw your own conclusion at the end.
   !M1
   open(newunit=u, file="part8-1.txt", status="replace")

   call BoundaryCondition(T1_implicit)
   call cpu_time(tic)
   call LineSOR(T1_implicit, w1, "8-1-w1_M1_linesor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T1_exact, T1_implicit, M1)
   err_sor_l1 = L1Norm(T1_exact, T1_implicit, M1)
   err_sor_linf = LInfNorm(T1_exact, T1_implicit, M1)
   print *, "line SOR method results for 10x10 mesh and w1 = 1.1 (L2 = 1e-12):"
   write(*, *) M1, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M1, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor

   call BoundaryCondition(T1)
   call cpu_time(tic)
   call PointSOR(T1, w3, "8-1-w3_M1_sor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T1_exact, T1, M1)
   err_sor_l1 = L1Norm(T1_exact, T1, M1)
   err_sor_linf = LInfNorm(T1_exact, T1, M1)
   print *, "Point SOR method results for 10x10 mesh and w3 = 1.5 (over relaxation) (L2 = 1e-12):"
   write(*, *) M1, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M1, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T1, "8-1-T1_pointsor_w3_l21e-12_result.txt")


   !M2
   call BoundaryCondition(T2_implicit)
   call cpu_time(tic)
   call LineSOR(T2_implicit, w1, "8-1-w1_M2_linesor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T2_exact, T2_implicit, M2)
   err_sor_l1 = L1Norm(T2_exact, T2_implicit, M2)
   err_sor_linf = LInfNorm(T2_exact, T2_implicit, M2)
   print *, "line SOR method results for 20x20 mesh and w1 = 1.1 (L2 = 1e-12):"
   write(*, *) M2, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M2, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T2_implicit, "8-1-T2_linesor_w1_l21e-12_result.txt")


   call BoundaryCondition(T2)
   call cpu_time(tic)
   call PointSOR(T2, w3, "8-1-w3_M2_linesor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T2_exact, T2, M2)
   err_sor_l1 = L1Norm(T2_exact, T2, M2)
   err_sor_linf = LInfNorm(T2_exact, T2, M2)
   print *, "PointSOR method results for 20x20 mesh and w3 = 1.5(overrelaxation) (L2 = 1e-12):"
   write(*, *) M2, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M2, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T2_implicit, "8-1-T2_pointsor_w3_l21e-12_result.txt")

   !M4
   call BoundaryCondition(T4_implicit)
   call cpu_time(tic)
   call LineSOR(T4_implicit, w1, "8-1-w1_M4_linesor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T4_exact, T4_implicit, M4)
   err_sor_l1 = L1Norm(T4_exact, T4_implicit, M4)
   err_sor_linf = LInfNorm(T4_exact, T4_implicit, M4)
   print *, "line SOR method results for 40x40 mesh and w1 = 1.1 (L2 = 1e-12):"
   write(*, *) M4, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M4, w1, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor

   call write_mat(T4_implicit, "8-1-T4_linesor_w1_l21e-12_result.txt")

   call BoundaryCondition(T4)
   call cpu_time(tic)
   call PointSOR(T4, w3, "8-1-w3_M2_linesor_l21e-12_converge_history.txt", iter_sor, L2_2)
   call cpu_time(toc)
   t_sor = toc - tic
   err_sor_l2 = L2Norm(T4_exact, T4, M4)
   err_sor_l1 = L1Norm(T4_exact, T4, M4)
   err_sor_linf = LInfNorm(T4_exact, T4, M4)
   print *, "PointSOR method results for 40x40 mesh and w3 = 1.5(overrelaxation) (L2 = 1e-12):"
   write(*, *) M4, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   write(u, *) M4, w3, err_sor_l2, err_sor_l1, err_sor_linf, t_sor, iter_sor
   call write_mat(T4_implicit, "8-1-T4_pointsor_w3_l21e-12_result.txt")
   call BoundaryCondition(T4)
   close(u)
   !======================
   print *, '************************************************************'
   print *, "Part 8-2:"
   open(newunit=u, file='part8-2.txt', status='replace')
   call BoundaryCondition(T1_implicit)
   call cpu_time(tic)
   call GaussElimination(T1_implicit)
   call cpu_time(toc)
   t_ge = toc - tic
   err_ge_l2 = L2Norm(T1_exact, T1_implicit, M1)
   err_ge_l1 = L1Norm(T1_exact, T1_implicit, M1)
   err_ge_linf = LInfNorm(T1_exact, T1_implicit, M1)
   print *, "Gauss Elimination (10x10):"
   write(*, *) M1, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   write(u, *) M1, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   call write_mat(T1_implicit, "8-2-T1_GaussElimination_result.txt")

   !M2
   call BoundaryCondition(T2_implicit)
   call cpu_time(tic)
   call GaussElimination(T2_implicit)
   call cpu_time(toc)
   t_ge = toc - tic
   err_ge_l2 = L2Norm(T2_exact, T2_implicit, M2)
   err_ge_l1 = L1Norm(T2_exact, T2_implicit, M2)
   err_ge_linf = LInfNorm(T2_exact, T2_implicit, M2)
   print *, "Gauss Elimination (20x20):"
   write(*, *) M2, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   write(u, *) M2, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   call write_mat(T2_implicit, "8-2-T2_GaussElimination_result.txt")
   !M4
   call BoundaryCondition(T4_implicit)
   call cpu_time(tic)
   call GaussElimination(T4_implicit)
   call cpu_time(toc)
   t_ge = toc - tic
   err_ge_l2 = L2Norm(T4_exact, T4_implicit, M4)
   err_ge_l1 = L1Norm(T1_exact, T1_implicit, M4)
   err_ge_linf = LInfNorm(T4_exact, T4_implicit, M4)
   print *, "Gauss Elimination (40x40):"
   write(*, *) M4, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   write(u, *) M4, err_ge_l2, err_ge_l1, err_ge_linf, t_ge
   call write_mat(T4_implicit, "8-2-T4_GaussElimination_result.txt")

   close(u)
   !==================
   print *, '************************************************************'
   print *, "Part 8-3:"
   open(newunit=u, file='part8-3.txt', status='replace')
   call BoundaryCondition(T1_implicit)
   call cpu_time(tic)
   call LineJacobi(T1_implicit, "8-3-M1_linejacobi_l21e-12_converge_history.txt", iter_lj, L2_2)
   call cpu_time(toc)
   t_lj = toc - tic
   err_lj_l2 = L2Norm(T1_exact, T1_implicit, M1)
   err_lj_l1 = L1Norm(T1_exact, T1_implicit, M1)
   err_lj_linf = LInfNorm(T1_exact, T1_implicit, M1)
   print *, "Line Jacobi method results for 10x10 mesh(L2 = 1e-12):"
   write(*, *) M1, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj
   write(u, *) M1, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj
   call write_mat(T1_implicit, "8-3-T1_linejacobi_l21e-12_result.txt")


   call cpu_time(tic)
   call LineJacobi(T2_implicit, "8-3-M2_linejacobi_l21e-12_converge_history.txt", iter_lj, L2_2)
   call cpu_time(toc)
   t_lj = toc - tic
   err_lj_l2 = L2Norm(T2_exact, T2_implicit, M2)
   err_lj_l1 = L1Norm(T2_exact, T2_implicit, M2)
   err_lj_linf = LInfNorm(T2_exact, T2_implicit, M2)
   print *, "Line Jacobi method results for 20x20 mesh(L2 = 1e-12):"
   write(*, *) M2, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj
   write(u, *) M2, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj

   call write_mat(T2_implicit, "8-3-T2_linejacobi_l21e-12_result.txt")

   call cpu_time(tic)
   call LineJacobi(T4_implicit, "8-3-_M4_linejacobi_l21e-12_converge_history.txt", iter_lj, L2_2)
   call cpu_time(toc)
   t_lj = toc - tic
   err_lj_l2 = L2Norm(T4_exact, T4_implicit, M4)
   err_lj_l1 = L1Norm(T4_exact, T4_implicit, M4)
   err_lj_linf = LInfNorm(T4_exact, T4_implicit, M4)
   print *, "Line Jacobi method results for 40x40 mesh(L2 = 1e-12):"
   write(*, *) M4, w1, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj
   write(u, *) M4, w1, err_lj_l2, err_lj_l1, err_lj_linf, t_lj, iter_lj
   call write_mat(T4_implicit, "8-3-T4_linejacobi_l21e-12_result.txt")


   !lg
   call BoundaryCondition(T1_implicit)
   call cpu_time(tic)
   call LineGaussSeidel(T1_implicit, "8-3-M1_linegaussseidel_l21e-12_converge_history.txt", iter_lgs, L2_2)
   call cpu_time(toc)
   t_lgs = toc - tic
   err_lgs_l2 = L2Norm(T1_exact, T1_implicit, M1)
   err_lgs_l1 = L1Norm(T1_exact, T1_implicit, M1)
   err_lgs_linf = LInfNorm(T1_exact, T1_implicit, M1)
   print *, "Line GaussSeidel method results for 10x10 mesh(L2 = 1e-12):"
   write(*, *) M1, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs
   write(u, *) M1, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs

   call write_mat(T1_implicit, "8-3-T1_linegaussseidel_l21e-12_result.txt")


   call cpu_time(tic)
   call LineGaussSeidel(T2_implicit, "8-3-M2_linegaussseidel_l21e-12_converge_history.txt", iter_lgs, L2_2)
   call cpu_time(toc)
   t_lgs = toc - tic
   err_lgs_l2 = L2Norm(T2_exact, T2_implicit, M2)
   err_lgs_l1 = L1Norm(T2_exact, T2_implicit, M2)
   err_lgs_linf = LInfNorm(T2_exact, T2_implicit, M2)
   print *, "Line GaussSeiedel method results for 20x20 mesh(L2 = 1e-12):"
   write(*, *) M2, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs
   write(u, *) M2, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs
   call write_mat(T2_implicit, "8-3-T2_linegaussseidel_l21e-12_result.txt")

   call cpu_time(tic)
   call LineGaussSeidel(T4_implicit, "8-3-_M4_linegaussseidel_l21e-12_converge_history.txt", iter_lgs, L2_2)
   call cpu_time(toc)
   t_lgs = toc - tic
   err_lgs_l2 = L2Norm(T4_exact, T4_implicit, M4)
   err_lgs_l1 = L1Norm(T4_exact, T4_implicit, M4)
   err_lgs_linf = LInfNorm(T4_exact, T4_implicit, M4)
   print *, "Line GaussSeidel method results for 40x40 mesh(L2 = 1e-12):"
   write(*, *) M4, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs
   write(u, *) M4, err_lgs_l2, err_lgs_l1, err_lgs_linf, t_lgs, iter_lgs
   call write_mat(T4_implicit, "8-3-T4_linegaussseidel_l21e-12_result.txt")
   close(u)

   !++++++++++
contains

   subroutine ExactSolution(T)
      real(dp), intent(out) :: T(:, :)
      real(dp) :: x, y
      integer :: m, i, j
      m = size(T, 1)
      do j = 1, m
         y = (j - 1) * 1.0_dp / (m - 1)
         do i = 1, m
            x = (i - 1) * 1.0_dp / (m - 1)
            T(i, j) = sin(pi * x) * sinh(pi * y) / sinh(pi)
         enddo
      enddo
   end subroutine


   subroutine BoundaryCondition(T)
      real(dp), intent(out) :: T(:, :)
      real(dp) :: x
      integer :: i, m
      m = size(T, 1)

      T(:, :) = 0.0
      do i = 1, m
         x = (i - 1) * 1.0_dp / (m - 1)
         T(i, m) = sin(pi * x)
      enddo
   end subroutine

   subroutine FirstOrderUpwindBC(T, Tw)
      real(dp), intent(in) :: Tw
      real(dp), intent(out) :: T(:, :)
      real(dp) :: x
      integer :: i, j, m
      m = size(T, 1)
      x = 0.5_dp / m
      T(:, :) = 0.0
      do j = 2, m
         T(1, j) = 2*Tw - T(2, j);
      enddo
      do j = 2, m
         T(m, j) = 2*sin(pi*x) - T(m - 1, j);
         x = x + 1.0_dp / m
      enddo

      do i = 2, m
         T(1, i) = 2*Tw - T(2, i)
      enddo
      do i = 2, m
         T(m, i) = 2*Tw - T(m - 1, i)
      enddo
   end subroutine


   subroutine PointJacobi(T, history, iter, L2)
      real(dp), intent(in) :: L2
      character(len=*), intent(in) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp), allocatable :: T_old(:, :)
      real(dp) :: err
      integer :: i, j, m, k, u
      m = size(T, 1)
      allocate(T_old(m, m))
      err = 1000.0
      k = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)
         T_old = T
         do j = 2, m - 1
            do i = 2, m - 1
               T(i, j) = (1.0_dp / 4) * (T_old(i, j + 1) + T_old(i, j - 1) + T_old(i + 1, j) + T_old(i - 1, j))
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         k = k + 1
         write(u, *) k, err
      enddo
      close(u)
      iter = k
   end subroutine

   subroutine LineJacobi(T, history, iter, L2)
      real(dp), intent(in) :: L2
      character(len=*), intent(in) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp), allocatable :: T_old(:, :), a(:), b(:), c(:), d(:), tmp(:)
      real(dp) :: err
      integer :: i, j, m, k, u
      m = size(T, 1)
      allocate(T_old(m, m))
      allocate(a(m - 1))
      a = -1
      allocate(b(m))
      b = 4.0
      allocate(c(m - 1))
      c = -1

      allocate(d(m))
      allocate(tmp(m))

      err = 1000.0
      k = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)
         T_old = T
         do i = 2, m - 1
            do j = 1, m
               d(j) = (T_old(i + 1, j) + T_old(i - 1, j))
            enddo
            call solve_tridiag(a, b, c, d, tmp, m)
            do j = 2, m -1
               T(i, j) = tmp(j)
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         k = k + 1
         write(u, *) k, err
      enddo
      close(u)
      iter = k
   end subroutine


   subroutine PointGaussSeidel(T, history, iter, L2)
      real(dp), intent(in) :: L2
      character(len=*), intent(in) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp) :: err
      real(dp), allocatable :: T_old(:, :)
      integer :: i, j, m, u
      m = size(T, 1)
      allocate(T_old(m, m))
      err = 1000.0
      iter = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)
         T_old = T
         do j = 2, m - 1
            do i = 2, m - 1
               T(i, j) = (1.0_dp / 4) * (T_old(i, j + 1) + T(i, j - 1) + T_old(i + 1, j) + T(i - 1, j))
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         iter = iter + 1
         write(u, *) iter, err

      enddo
      close(u)

   end subroutine

   subroutine LineGaussSeidel(T, history, iter, L2)
      real(dp), intent(in) :: L2
      character(len=*), intent(in) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp) :: err
      real(dp), allocatable :: T_old(:, :), a(:), b(:), c(:), d(:), tmp(:)
      integer :: i, j, m, u
      m = size(T, 1)
      allocate(T_old(m, m))
      allocate(a(m))
      a = -1.0
      allocate(b(m))
      b = 4.0
      allocate(c(m))
      c = -1.0
      allocate(d(m))
      allocate(tmp(m))

      err = 1000.0
      iter = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)
         T_old = T
         do i = 2, m-1
            do j = 1, m
               d(j) = T_old(i+1, j) + T(i - 1, j)
            enddo
            call solve_tridiag(a, b, c, d, tmp, m)
            do j = 2, m - 1
               T(i, j) = tmp(j)
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         iter = iter + 1
         write(u, *) iter, err

      enddo
      close(u)

   end subroutine

   subroutine PointSOR(T, w, history, iter, L2)
      real(dp), intent(in) :: L2, w
      character(len=*) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp) :: err
      real(dp), allocatable :: T_old(:, :)
      integer :: i, j, m, u
      m = size(T, 1)
      allocate(T_old(m, m))
      err = 1000.0
      iter = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)
         T_old = T
         do j = 2, m - 1
            do i = 2, m - 1
               T(i, j) = (1.0 - w) * T_old(i, j) + (w / 4.0) *&
                  (T_old(i, j + 1) + T(i, j - 1) + T_old(i + 1, j) + T(i - 1, j))
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         iter = iter + 1
         write(u, *) iter, err

      enddo
      close(u)
   end subroutine


   subroutine LineSOR(T, w, history, iter, L2)
      real(dp), intent(in) :: L2, w
      character(len=*) :: history
      real(dp), intent(out) :: T(:, :)
      integer, intent(out) :: iter
      real(dp) :: err, h
      real(dp), allocatable :: T_old(:, :), a(:), b(:), c(:), d(:), tmp(:)
      integer :: i, j, m, u
      m = size(T, 1)
      h = 1.0_dp / (m - 1)
      allocate(T_old(m, m))
      allocate(a(m))
      allocate(b(m))
      allocate(c(m))
      allocate(d(m))
      allocate(tmp(m))
      a = -w
      b = 4.0
      c = -w
      err = 1000.0
      iter = 0
      open(newunit=u, file=history, status='replace')
      do while(err .gt. L2)

         T_old = T
         do i = 2, m - 1
            do j = 1, m
               d(j) = 4.0 * (1.0 - w) * T(i, j) + w * (T(i - 1, j) + T(i + 1, j))
            enddo
            call solve_tridiag(a, b, c, d, tmp, m)
            do j = 2, m - 1
               T(i, j) = tmp(j)
            enddo
         enddo
         err = L2Norm(T, T_old, m)
         iter = iter + 1
         write(u, *) iter, err
      enddo
      close(u)

   end subroutine

   subroutine solve_tridiag(a,b,c,d,x,n)

      !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
      !	 b - the main diagonal
      !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
      !	 d - right part
      !	 x - the answer
      !	 n - number of equations
      implicit none
      integer,intent(in) :: n
      real(dp),dimension(n),intent(in) :: a,b,c,d
      real(dp),dimension(n),intent(out) :: x
      real(dp),dimension(n) :: cp, dprime
      real(dp) :: m
      integer:: i

      ! initialize c-prime and d-prime
      cp(1) = c(1)/b(1)
      dprime(1) = d(1)/b(1)
      ! solve for vectors c-prime and d-prime
      do i = 2,n
         m = b(i)-cp(i-1)*a(i)
         cp(i) = c(i)/m
         dprime(i) = (d(i)-dprime(i-1)*a(i))/m
      end do
      ! initialize x
      x(n) = dprime(n)
      ! solve for x from the vectors c-prime and d-prime
      do i = n-1, 1, -1
         x(i) = dprime(i)-cp(i)*x(i+1)
      end do

   end subroutine solve_tridiag


   subroutine GaussElimination(T)
      real(dp), intent(out) :: T(:, :)
      real(dp), allocatable :: A(:, :), Ty(:), Tx(:), x(:)
      real(dp) :: mul
      integer :: m, i, j, k
      m = size(T, 1)
      allocate(x(m))
      allocate(Tx(m**2))
      allocate(Ty(m**2))
      allocate(A(m**2, m**2))

      do i = 1, m
         x(i) = (i - 1.0) / (m - 1)
      enddo
      A = 0.0_dp
      Ty = 0.0_dp
      do i = 1, m
         do j = 1,m
            k = j + m * (i - 1)

            if (i .eq. 1) then

               if (j .eq. 1) then
                  A(k, k) = -(1.0 + 1.0/2.0)
                  A(k, k + 1) = 1.0 / 4.0
                  A(k, k + m) = 1.0 / 4.0
               else if (j .eq. m) then
                  A(k, k) = -(1.0 + 1.0/2.0)
                  A(k, k -1) = 1.0 / 4.0
                  A(k, k + m) = 1.0/ 4.0
               else
                  A(k, k) = -(1.0 + 1.0/4.0)
                  A(k, k + 1) = 1.0/4.0
                  A(k, k - 1) =  1.0 /4.0
                  A(k, k + m) = 1.0/4.0
               endif
            else if (i .eq. m) then

               if (j .eq. 1) then
                  Ty(k) = -sin(pi*x(j))
                  A(k, k) = -(1.0 + 1.0/2.0)
                  A(k, k + 1) = 1.0 / 4.0
                  A(k, k - m) = 1.0 / 4.0
               else if (j .eq. m) then
                  Ty(k) = -sin(pi*x(j))
                  A(k, k) = -(1.0 + 1.0/2.0)
                  A(k, k -1) = 1.0 / 4.0
                  A(k, k - m) = 1.0/ 4.0
               else
                  Ty(k) = -0.5 * sin(pi*x(j))

                  A(k, k) = -(1.0 + 1.0/4.0)
                  A(k, k + 1) = 1.0/4.0
                  A(k, k - 1) =  1.0 /4.0
                  A(k, k - m) = 1.0/4.0
               endif
            else
               if (j .eq. 1) then
                  A(k, k) = -(1.0 + 1.0/4.0)
                  A(k, k + 1) = 1.0 / 4.0
                  A(k, k - m) = 1.0 / 4.0
                  A(k, k + m) = 1.0 / 4.0

               else if (j .eq. m) then
                  A(k, k) = -(1.0 + 1.0/4.0)
                  A(k, k -1) = 1.0 / 4.0
                  A(k, k - m) = 1.0/ 4.0
                  A(k, k + m) = 1.0/ 4.0
               else
                  A(k, k) = -1.0
                  A(k, k + m) = 1.0/4.0
                  A(k, k + 1) = 1.0/4.0
                  A(k, k - 1) =  1.0 /4.0
                  A(k, k - m) = 1.0/4.0
               endif
            endif
         enddo
      enddo
      do j = 1, m**2 - 1
         do i = j+1, m**2
            mul = A(i, j) / A(j, j)
            do k = j, m**2
               A(i, k) = A(i, k) - mul * A(j, k)
            enddo
            Ty(i) = Ty(i) - mul * Ty(j)
         enddo
      enddo

      Tx = Ty
      do i = 1, m**2
         if (i .eq. m**2) then
            Tx(m**2) = A(m**2, m**2)
         endif
         do j = m**2, i, -1
            Ty(i) = Ty(i) - A(i, j) * Tx(j)
         enddo
         Tx(i) = Ty(i) / A(i, i)
      enddo

      do i = 1, m
         do j = 1,m
            k = j + m * (i - 1)
            T(i, j) = Tx(k)
         enddo
      enddo

   end subroutine


end program cfd1proj1laplace
