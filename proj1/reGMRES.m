%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Ax=b for sparse matrix A by restarted GMRES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iter = reGMRES (A, b, m)
  iter = 0;
  eps = 1e-8; %||b-Ax||_2<=eps
  N = length(b);
  x = randn(N, 1);
  V = zeros(N, m);
  H = zeros(m+1, m);
  e = zeros(m+1, 1);
  e(1) = 1;
  r = b - A * x;
  while r' * r > eps.^2
    beta = norm(r, 2);
    V(:, 1) = r ./ beta;
    for j = 1: m
      w = A * V(:, j);
      for i = 1: j
        H(i, j) = w' * V(:, i);
        w = w - H(i, j) .* V(:, i);
      end
      H(j+1, j) = norm(w, 2);
      if j < m
        V(:, j+1) = w ./ H(j+1, j);
      end
    end
    y = H \ (beta .* e);
    x = x + V * y;
    r = b - A * x;
    iter = iter + 1;
  end
end
