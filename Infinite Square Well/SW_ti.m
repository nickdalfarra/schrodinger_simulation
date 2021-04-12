function phi = SW_ti(x, n, a)
%Time dependent solution to the infinite square well.
%   Inputs
%       x: evaluation points
%       n: energy level
%       a: well length
%   Outputs
%       phi: time independent solution


phi = sqrt(2/a)*sin(n*pi*x/a).*(x >= 0 & x <= 1);

end