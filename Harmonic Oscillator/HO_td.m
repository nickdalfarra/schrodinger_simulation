function f = HO_td(t, n, omega)
%Time independent solution to the harmonic oscillator.
%   Inputs
%       t: time
%       n: energy level
%       omega: angular frequency
%   Outputs
%       f: time dependent solution

j = sqrt(-1);
f = exp(-j*(n+1/2)*omega*t);

end


