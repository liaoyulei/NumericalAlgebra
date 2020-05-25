%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Bx=b by SSOR and SSOR_CG    %
% choose the solver in line 7 and 8 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 1; %relaxation factor
%format = 'SSOR';
format = 'SSOR_CG';
str_eval = [format, '(B, b, omega)'];
if strcmp(format, 'SSOR')
  ne = 3;
  N = [128, 256, 512];
elseif strcmp(format, 'SSOR_CG')
  ne = 6;
  N = [128, 256, 512, 1024, 2048, 4096];
end
iter = zeros(1, ne);
time = zeros(1, ne);
order = zeros(1, ne);
for in = 1: ne
  M = N(in) .* 4;
  A = randn(M, N(in));
  x_0 = randn(N(in), 1);
  z = randn(N(in), 1);
  d = abs(A * x_0);
  D = spdiags(d, 0, M, M);
  B = A' * D * A;
  y = abs(A * z);
  b = A' * D * (d.^2 - y.^2) ./ 2;
  tic;
  iter(in) = eval(str_eval);
  time(in) = toc;
  if in > 1
    order(in) = log2(time(in)./time(in-1));
  end
end
subplot(1, 2, 1);
hold on;
plot(N, time);