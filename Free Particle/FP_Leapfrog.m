% Free Particle using leapfrog method

clear;clc;
% Initilize parameters  
x_0=6;
sigma_squared=1e-2;
k_0=30;
% Discretisation parameters
dx=1e-3;
x = (0:dx:14-dx);
dt=5e-7;
t_end = 0.1;
n_it=t_end/dt;

m = 1;
hbar = 1;
omega = 1;
a = 1;
j = sqrt(-1);
K = hbar^2/(2*m*dx^2);



% Generate an intial wavepacket
Psi_initial=exp(-(x-x_0).^2/sigma_squared).'.*exp(1i*k_0*x).';
integrated = dx*sum(abs(Psi_initial).^2);
Psi_initial = Psi_initial/integrated; % normalize
N = length(x);

sigma_barrier = 0.01;
V = 1/(2*pi*sigma_barrier^2)*exp(-(x-6.5).^2/(2*sigma_barrier^2)).';

R_initial=real(Psi_initial);
I_initial=imag(Psi_initial);
Psi_old = R_initial+j*I_initial;
%I_current=I_initial;
%R_current=R_initial;


% first time step using forward time central space method
% first and last points treated seperately for peridoicity
Psi = zeros(length(Psi_old),1);
Psi(2:N-1) = Psi_old(2:N-1) - (j*dt/hbar)*(K*(Psi_old(3:N)-2*Psi_old(2:N-1)+Psi_old(1:N-2)) + V(2:N-1).*Psi_old(2:N-1));
Psi(1) = Psi_old(1) - (j*dt/hbar)*(K*(Psi_old(2)-2*Psi_old(1)+Psi_old(N)) + V(1).*Psi_old(1));
Psi(N) = Psi_old(N) - (j*dt/hbar)*(K*(Psi_old(1)-2*Psi_old(N)+Psi_old(N-1)) + V(N).*Psi_old(N));

R_current = real(Psi);
I_current = imag(Psi);




% Initialize arrays for errors %%%%%
%probs = zeros(ceil(n_it),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Main iteration loop
for iter = 1:n_it
    
    [R_next]=real_psi(N, R_current, I_current, dt, dx, V);
    R_current=R_next;
    
    [I_next] = imag_psi(N, I_current, R_current, dt, dx, V);

    I_current=I_next;
    prob_density = R_current.^2+I_current.^2;
    
    % Document Errors %%%%%%%%%%%%%%
    %probs(iter) = sum(prob_density)*dx;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the results
    % Plot only every 10 calculatoions for speed
    if rem(iter, 10)== 0
        plot(x, prob_density);
        hold on
        plot(x,V)
        hold off
        title('Free Particle Problem');
        axis([0 14 0 100]);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end

%%%%%%% Error Plots %%%%%%%%%%
%figure()
%plot(1:length(probs),probs)
%title('Probability Conservation');
%axis([0 iter 0 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Functions to compute Psi
%Calculates the real part of wavefunction
function [R_next]= real_psi(N, R_current, I_current, dt, dx, V)
    R_next= zeros(1,N);
    K=dt/(2*dx^2);
    for x=2:N-1
        % perform calculation at t+dt given time t
        R_next(x)=R_current(x) - K*(I_current(x+1)-2*I_current(x)+I_current(x-1))+dt.*V(x).*I_current(x);
    end
    % Boundary conditions
    R_next(1)=R_next(2);
    R_next(N)=R_next(N-1);
end

% Calculates the imaginary part of the wavefunction
function [I_next]= imag_psi(N, I_current, R_current, dt, dx, V)
    I_next= zeros(1,N);
    K=dt/(2*dx^2);
    for x=2:N-1
         %Perform calculation at time t+dt given time t
         I_next(x)=I_current(x) +K*(R_current(x+1)-2*R_current(x)+R_current(x-1))-dt.*V(x).*R_current(x);
    end
    % Boundary conditions
    I_next(1)=I_next(2);
    I_next(N)=I_next(N-1);

end

