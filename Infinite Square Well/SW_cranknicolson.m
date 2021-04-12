clear;clc;
% Infinite square well test
n = 1;
m = 1;
a = 1;
hbar = 1;
dt = 0.1;
t_end = 2;
t=dt:dt:t_end;
n_iter = length(t);

x_start = -0.1; x_end = 1.1;
N = 100; h = (x_end - x_start)/N;
x = x_start:h:x_end - h;

V_func = @(x,t) realmax('single')*(x < 0 | x > 1);

v_old = SW_ti(x,n,a)*SW_td(0,n,a,m, hbar); % start w/ exact
v_old = v_old.';

V = V_func(x,0).';
[A, B] = schro_CN(N, h, dt, m, hbar, V);

probs = zeros(n_iter,1);
errors = zeros(n_iter,1);

clear figure(1)
figure(1)
for iter=1:n_iter
    
    exact = SW_ti(x,n,a)*SW_td(t(iter),n,a,m, hbar);
    v_new = A\B*v_old;
    
    v_old = v_new;
    
%     plot(x, real(v_new)); hold on
%     plot(x,imag(v_new))
%     plot(x, abs(v_new).^2); hold off
%     legend('Real(sim)', 'Imag(sim)','Prob')
%     axis([-0.1 1.1 -2 2])
%     pause(0.05)
    
    probs(iter) = sum(abs(v_new).^2)*h;
    errors(iter) = sum(abs(exact.'-v_new))*h;
end
figure(2)
plot(1:n_iter,probs); title('probability')

figure(3)
plot(1:n_iter, errors)

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
