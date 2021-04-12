function phi = HO_ti(x, n, m, hbar, omega)
%Time independent solution to the harmonic oscillator.
%   Inputs
%       x: support points where the function will be evaluated.
%       n: energy level of the exact solution (0 is ground state)
%       m: mass (?) of the particle
%       hbar: angular Planck's constant
%       omega: angular frequency
%   Outputs
%       phi: analytical solution

eta = sqrt(m*omega/hbar)*x;

phi = (m*omega/(pi*hbar))^(1/4)*(2^n * factorial(n))^(-1/2)*...
    hermiteH(n, eta).*exp(-eta.^2/2);

end

