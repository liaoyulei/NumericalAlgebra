%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init sparse matrix A and vector b %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, b] = init (n) %n mesh size
  h = 1 ./ n;
  N = 3 .* (n + 1).^2;
  b = zeros(N, 1);
  index = idx(['u', 'u', 'u', 'u', 'v', 'v', 'v', 'v', 'w', 'w', 'w', 'w'], [0, 0, n, n, 0, 0, n, n, 0, 0, n, n], [0, n, 0, n, 0, n, 0, n, 0, n, 0, n], n);
  A = sparse(index, index, 1, N, N, N.*7); %u_{0,0}=u_{0,n}=u_{n,0}=u_{n,n}=v_{0,0}=v_{0,n}=v_{n,0}=v_{n,n}=w_{0,0}=w_{0,n}=w_{n,0}=w_{n,n}=0
  for j = 1: n-1 %x=0
    index = idx('u', 0, j, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %u_{0,j}=0
    index = idx('v', 0, j, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %v_{0,j}=0
    index = idx(['w', 'v', 'v', 'v'], [0, 0, 1, 2], [j, j, j, j], n);
    A(idx('w', 0, j, n), :) = sparse(1, index, [2, n.*3, -n.*4, n], 1, N, 4); %2w_{0,j}+3nv_{0,j}-4nv_{1,j}+nv_{2,j}=0
  end
  for j = 1: n-1 %x=1
    index = idx('u', n, j, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %u_{n,j}=0
    index = idx('v', n, j, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %v_{n,j}=0
    index = idx(['w', 'v', 'v', 'v'], [n, n-2, n-1, n], [j, j, j, j], n);
    A(idx('w', n, j, n), :) = sparse(1, index, [2, -n, n.*4, -n.*3], 1, N, 4); %2w_{n,j}-nv_{n-2,j}+4nv_{n-1,j}-3nv_{n,j}=0
  end
  for i = 1: n-1 %y=0
    index = idx('u', i, 0, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %u_{i,0}=0
    index = idx('v', i, 0, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %v_{i,0}=0
    index = idx(['w', 'u', 'u', 'u'], [i, i, i, i], [0, 0, 1, 2], n);
    A(idx('w', i, 0, n), :) = sparse(1, index, [2, -n.*3, n.*4, -n], 1, N, 4); %2w_{i,0}-3nu_{i,0}+4nu_{i,1}-nu_{i,2}=0
  end
  for i = 1: n-1 %y=1
    index = idx('u', i, n, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %u_{i,n}=1
    b(index) = 1;
    index = idx('v', i, n, n);
    A(index, :) = sparse(1, index, 1, 1, N, 1); %v_{i,n}=0
    index = idx(['w', 'u', 'u', 'u'], [i, i, i, i], [n, n-2, n-1, n], n);
    A(idx('w', i, n, n), :) = sparse(1, index, [2, n, -n.*4, n.*3], 1, N, 4); %2w_{i,n}+nu_{i,n-2}-4nu_{i,n-1}+3nu_{i,n}=0
  end
  for i = 1: n-1 %inner
    for j = 1: n-1
      index = idx(['u', 'u', 'u', 'u', 'u', 'w', 'w'], [i, i-1, i+1, i, i, i, i], [j, j, j, j-1, j+1, j-1, j+1], n);
      A(idx('u', i, j, n), :) = sparse(1, index, [4, -1, -1, -1, -1, h./2, -h./2], 1, N, 7); %4u_{i,j}-u_{i-1,j}-u_{i+1,j}-u_{i,j-1}-u_{i,j+1}+h/2w_{i,j-1}-h/2w_{i,j+1}=0
      index = idx(['v', 'v', 'v', 'v', 'v', 'w', 'w'], [i, i-1, i+1, i, i, i-1, i+1], [j, j, j, j-1, j+1, j, j], n);
      A(idx('v', i, j, n), :) = sparse(1, index, [4, -1, -1, -1, -1, -h./2, h./2], 1, N, 7); %4v_{i,j}-v_{i-1,j}-v_{i,j-1}-v_{i,j-1}-v_{i,j+1}-h/2w_{i-1,j}+h/2w_{i+1,j}=0
      index = idx(['w', 'w', 'w', 'w', 'w'], [i, i-1, i+1, i, i], [j, j, j, j-1, j+1], n);
      A(idx('w', i, j, n), :) = sparse(1, index, [4, -1, -1, -1, -1], 1, N, 5); %4w_{i,j}-w_{i-1,j}-w_{i+1,j}-w_{i,j-1}-w_{i,j+1}=0
    end
  end
end
