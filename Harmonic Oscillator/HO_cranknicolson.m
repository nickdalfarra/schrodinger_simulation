% Crank-Nicolson method implemented on the Harmonic Oscillator

m=1;
energy_level = 0;
hbar = 1;
omega_vib = 1;

omega = (energy_level+0.5)*omega_vib; % found this formula 
Vfunc=@(x,t) (omega^2*m)/2*x.^2;

% simulation parameters
N = 200;
sim_r = 5; % simulation radius is over domain [-sim_r, sim_r]
h = 2*sim_r/N; 
dt = 0.2; % arbitrary for CN
t_end = 4*pi;

n_iter = t_end/dt;

x = -sim_r:h:sim_r-h; % remove last point because of periodicity
t = dt:dt:t_end;

v_old = HO_ti(x,energy_level,m,hbar,omega)*...
    HO_td(0,energy_level,omega);
v_old = v_old.';

%plot(x,v_old)

clear figure(1)
%figure(1)

errors = zeros(ceil(n_iter),1);

e = ones(N, 1);
energy_mat = spdiags([e -2*e e], -1:1, N, N); 
energy_mat(N,1) = 1; energy_mat(1,N) = 1;
energies = zeros(ceil(n_iter),1);

probs = zeros(ceil(n_iter),1);

for iter=1:n_iter
    V = Vfunc(x, t(iter)).';
    [A,B] = schro_CN(N, h, dt, m, hbar, V);
    v_new = A\(B*v_old);
    v_old = v_new;

    v_exact = HO_ti(x,energy_level,m,hbar,omega)*...
        HO_td(t(iter),energy_level,omega);

    % document energy and errors
    errors(iter) = sum(abs(v_exact.' - v_new))*h;
    energies(iter) = abs(-hbar^2/(2*m*h)*v_new'*energy_mat*v_new);
    probs(iter) = sum(abs(v_new).^2)*h;

    result_plot(x, v_new, v_exact, sim_r)
%     error_plot(x, v_new, v_exact)
%     cumulative_error_plot(iter, errors)

    figure(4)
    plot(1:iter, energies(1:iter)); title('energies')
    figure(5)
    plot(1:iter, probs(1:iter)); title('total probability')

end

cumulative_error_plot(iter, errors)


function cumulative_error_plot(iter, errors)
    figure(3)
    plot(1:iter, errors(1:iter))
    title('L1 error by iteration')
end

function error_plot(x, v_new, v_exact)
    
    err = v_exact.' - v_new;
    figure(2)
    plot(x, real(err)); hold on;
    plot(x, imag(err)); hold off;
    xlabel('x'); ylabel('error');
    title('exact minus simulated')
    legend('real', 'imag')
end

function result_plot(x, v_sim, v_exact, sim_r)
    figure(1)
    plot(x, real(v_sim), 'r-'); hold on
    plot(x, imag(v_sim), 'xr-')
    plot(x, real(v_exact), 'b-')
    plot(x, imag(v_exact), 'xb-')
    axis([-sim_r sim_r -1 1])
    xlabel('x'); ylabel('v');
    legend('sim real', 'sim imag', 'exact real', 'exact imag')
    hold off
end
    
function [A,B] = schro_CN(N, h, dt, m, hbar, V)
    % A*v_new = B*v_old
    % V is a N by 1 vector of potentials along the support
    
    e = ones(N,1); j = sqrt(-1); K = dt*hbar/(4*m*h^2);
    K_times_e = K*e;
    H = spdiags([K_times_e -2*K_times_e K_times_e], -1:1, N, N);
    H(1,N) = K; H(N,1) = K; % periodic boundary conditions
    
    V_mat = dt/(2*hbar)*spdiags(V, 0, N, N);
    I = speye(N);
    
    A = I - j*(H - V_mat);
    B = I + j*(H - V_mat);

end

