clear;clc;
x_0=6;
sigma_squared=1e-2;
delta_x=1e-2;
x=(0:delta_x:14);
delta_t=2.5e-5;
% delta_t=2.5e-4;
N = length(x); m=1; hbar=1;

k_0=30;
times = 0:delta_t:0.1;

v_old=exp(-(x-x_0).^2/sigma_squared).'.*exp(1i*k_0*x).';

integrated = delta_x*sum(abs(v_old).^2);
v_old = v_old/integrated; % normalize

sigma_barrier = 0.01;
V = 1/(2*pi*sigma_barrier^2)*exp(-(x-7).^2/(2*sigma_barrier^2)).';

[A,B] = schro_CN(N, delta_x, delta_t, m, hbar, V);

probs = zeros(length(times),1);

times_list = [1 950 2000 3000 4000];
% times_list = [1 95 200 300 400];

counter = 1;

for iter=1:length(times)
%     string(iter)
    v_new = A\B*v_old;
    v_old = v_new;
    
    probs(iter) = delta_x*sum(abs(v_new).^2);
    
    if mod(iter,50)==0
       iter 
    end
    
    if ismember(iter, times_list)
        subplot(1,5,counter); hold on;
        prob = abs(v_new).^2;
        plot(x, prob)
        plot(x, V); 
        axis([0 14 0 9]); hold off;
%         axis([4 10 -9 9])
        counter = counter + 1;
%         figure(2)
%         plot(1:iter, probs(1:iter)); title('probability')
        pause(0.1)
    end
    
%     figure(3)
%     plo
    
end

%%
figure(2)
plot(1:iter, probs(1:iter))

%%
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

% subplot(2,2,1);
% 
% plot(x,real(psi), 'b');
% title('Real part of wavefunction ');
% xlabel('distance');
% ylabel('Re(wavefunction)');
% subplot(2,2,2);
% plot(x,imag(psi),'r');
% title('Imaginary part of wavefunction');
% xlabel('distance');
% ylabel('Im(wavefunction)');
% subplot(2,2,3);
% plot(x,(conj(psi).*psi),'k');
% title('Probability of finding particle ');
% xlabel('distance');
% ylabel('psi*conj(psi)');
