%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the index in x for element ch_{x,y} %
% u_{x,y} 3((n+1)y+x)+1               %
% v_{x,y} 3((n+1)y+x)+2               %
% w_{x,y} 3((n+1)y+x)+3               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function index = idx (ch, x, y, n)
  index = 3 .* ((n + 1) .* y + x) + ch - 't';
end
