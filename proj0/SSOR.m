%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Ax=b by SSOR                                               %
% (D-\omega C_L)x_{m-1/2}=(\omega C_U+(1-\omega)D)x_{m-1}+\omega b %
% (D-\omega C_U)x_m=(\omega C_L+(1-\omega)D)x_{m-1/2}+\omega b     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iter = SSOR (A, b, omega) %omega relaxation factor
  iter = 0;
  eps = 1e-8; %||b-Ax||_2<=eps
  N = length(b);
  D = spdiags(diag(A), 0, N, N);
  C_L = tril(D-A);
  C_U = triu(D-A);
  M1 = D - omega .* C_L;
  M2 = D - omega .* C_U;
  N1 = omega .* C_U + (1 - omega) .* D;
  N2 = omega .* C_L + (1 - omega) .* D;
  x = randn(N, 1);
  r = b - A * x;
  while r' * r > eps.^2
    x = M1 \ (N1 * x + omega .* b);
    x = M2 \ (N2 * x + omega .* b);
    r = b - A * x;
    iter = iter + 1;
  end
end
