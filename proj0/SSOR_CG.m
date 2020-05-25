%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Ax=b by PCG with                 %
% M = (D-\omega C_L)D^{-1}(D-\omega C_U) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iter = SSOR_CG (A, b, omega) %omega relaxation factor
  iter = 0;
  eps = 1e-8; %||b-Ax||_2<=eps
  N = length(b);
  D = spdiags(diag(A), 0, N, N);
  C_L = tril(D-A);
  C_U = triu(D-A);
  M1 = D - omega .* C_L;
  M2 = D - omega .* C_U;
  x = randn(N, 1);
  r = b - A * x;
  z = M2 \ (D * (M1 \ r));
  p = z;
  rho_old = r' * z; %\rho_{k-1}
  while r' * r > eps.^2
    w = A * p;
    alpha = rho_old ./ (p' * w);
    x = x + alpha .* p;
    r = r - alpha .* w;
    z = M2 \ (D * (M1 \ r));
    rho_new = r' * z; %\rho_k
    beta = rho_new ./ rho_old;
    p = z + beta .* p;
    rho_old = rho_new;
    iter = iter + 1;
  end
end
