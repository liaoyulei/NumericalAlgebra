%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Bx=b by GMRES without or with preconditioner %
% choose the solver in line 9 and 10                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = [8; 16; 32; 64; 128]; %mesh size
m = 30; %dimention of K_m
format = 'reGMRES';
%format = 'preGMRES';
str_eval = [format, '(A, b, m)'];
ne = 5;
iter = zeros(1, ne);
time = zeros(1, ne);
order = zeros(1, ne);
for in = 1: ne
  [A, b] = init(n(in));
  tic;
  iter(in) = eval(str_eval);
  time(in) = toc;
  if in > 1
    order(in) = log2(time(in)./time(in-1));
  end
end
hold on;
plot(3.*(n+1).^2, time);

n = 64;
m = 30;
u = zeros(n+1);
v = zeros(n+1);
[A, b] = init(n);
[iter, x] = preGMRES(A, b, m);
tt = 0: (1/n): 1;
[xx, yy] = meshgrid(tt, tt);
for j = 0: n
  for i = 0: n
    u(i+1, j+1) = x(idx('u', j, i, n));
    v(i+1, j+1) = x(idx('v', j, i, n));
  end
end
subplot(1, 2, 1);
surf(xx, yy, u);
xlabel('x');
ylabel('y');
zlabel('u');
title('u');
subplot(1, 2, 2);
surf(xx, yy, v);
xlabel('x');
ylabel('y');
zlabel('v');
title('v');