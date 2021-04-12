function f = SW_td(t, n, a, m, hbar)
%Time dependent solution to infinite square well
%   Inputs
%       t: time
%       n: energy level
%       a: well length
%       m: mass
%       hbar: angular Planck's constant
%   Outputs
%       f: time dependent solution

j = sqrt(-1);
f = exp(-j*((n^2*pi^2*hbar)/(2*m*a^2))*t);

end
