function [ A ] = Create_Poisson_problem_A( N )

  % Create the archetypical matrix A for an N x N Poisson problem (5-point
  % stencil).  
  
  % Set the diagonal
  d = 4*ones(N^2,1);
  
  % Set the entries of the first sub and super diagonals
  sd = -1*ones(N^2,1);
  
  % Set the other off-diagonal entries
  d(N+1:N:end) = -1;
  d(N:N:end-1) = -1;
  
  % Assemble the matrix A
  A = spdiags([sd d sd], [-1 0 1], N^2, N^2);

  % Visualize the matrix A
  figure;
  spy(A);
  title('Visualization of matrix A');
  
%end