% Harmonic Osocillator using leapfrog method

clear;clc;
% Initilize parameters  
N=50;
dx =8/N;
dt=0.01;

% Quantum Number n
n = 1;

m = 1;
hbar = 1;
omega = 1;

j = sqrt(-1);
K = hbar^2/(2*m*dx^2);

% time at which we want to end the simulation
t_end=2*pi;
% number of timesteps to be taken
n_it=t_end/dt;

% initialize spatial array
x=(-4:dx:4-dx);
% temporal grid
t=(0:dt:n_it*dt);

% Potential function
V = 0.5*m*omega.^2.*x.^2;

% the initial condition
Psi = HO_ti(x,n,m,hbar,omega).*HO_td(0,n,omega);

R_initial=real(Psi);
I_initial=imag(Psi);
Psi_old = R_initial+j*I_initial;

%I_current=I_initial;
%R_current=R_initial;


% first time step using forward time central space method
% first and last points treated seperately for peridoicity
Psi = zeros(1,length(Psi));
Psi(2:N-1) = Psi_old(2:N-1) - (j*dt/hbar)*(K*(Psi_old(3:N)-2*Psi_old(2:N-1)+Psi_old(1:N-2)) + V(2:N-1).*Psi_old(2:N-1));
Psi(1) = Psi_old(1) - (j*dt/hbar)*(K*(Psi_old(2)-2*Psi_old(1)+Psi_old(N)) + V(1).*Psi_old(1));
Psi(N) = Psi_old(N) - (j*dt/hbar)*(K*(Psi_old(1)-2*Psi_old(N)+Psi_old(N-1)) + V(N).*Psi_old(N));

R_current = real(Psi);
I_current = imag(Psi);




% Initialize arrays for errors %%%%%
probs = zeros(ceil(n_it),1);
errors = zeros(ceil(n_it),1);

e = ones(N, 1);
energy_mat = spdiags([e -2*e e], -1:1, N, N); 
energy_mat(N,1) = 1; energy_mat(1,N) = 1;
energies = zeros(ceil(n_it),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Main iteration loop
for iter = 1:n_it
    
    [R_next]=real_psi(N, R_current, I_current, dt, dx, V);
    R_current=R_next;
    
    [I_next] = imag_psi(N, I_current, R_current, dt, dx, V);

    % calculate psi*psi with R(t) and I(t+delta_t/2) and I(t-delta_t/2)

    I_current=I_next;
    prob_density = R_current.^2+I_current.^2;
    
    %Exact Solution
    Psi_exact = HO_ti(x,n,m,hbar,omega)*HO_td(t(iter),n,omega);

    % Document Errors %%%%%%%%%%%%%%
    probs(iter) = sum(prob_density)*dx;
    Psi_current = R_current+j*I_current;
    errors(iter) = sum(abs(Psi_exact- Psi_current))*dx;
    %errors(iter) = sum(abs(Psi_exact.^2- Psi_current.^2))*dx;
    energies(iter) = abs(-hbar^2/(2*m*dx)*prob_density*energy_mat*prob_density');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot the results
    % Plot only every 10 calculatoions for speed
    if rem(iter, 10)== 0

        plot(x, prob_density,'-x');
        hold on
        plot(x,R_current)
        plot(x,I_current)
        plot(x,real(Psi_exact).^2+imag(Psi_exact).^2)
        hold off
        legend("Numerical Probability","Numerical Real","Numerical Imaginary","Exact Probability")
        title('Quantum Oscillator');
        axis([-4 4 -1 1]);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end
%%%%%%% Error Plots %%%%%%%%%%
figure()
plot(1:length(probs),probs)
title('Probability Conservation');
axis([0 iter 0 2]);
xlabel("Iteration")
ylabel("Probability")

figure()
plot(1:length(errors)-1,errors(1:end-1))
title('L1 Error by Iteration');
%axis([0 iter 0 0.05])
xlabel("Iteration")
ylabel("Error Value")

figure()
plot(1:length(energies)-1,energies(1:end-1))
title("Energies")
xlabel("Iteration")
ylabel("Energy")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
% Functions to compute Psi
%Calculates the real part of wavefunction
function [R_next]= real_psi(N, R_current, I_current, dt, dx, V)
    R_next= zeros(1,N);
    K=dt/(2*dx^2);
    for x=2:N-1
        % perform calculation at t+dt given time t
        R_next(x)=R_current(x) - K*(I_current(x+1)-2*I_current(x)+I_current(x-1))+dt*V(x).*I_current(x);
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
         I_next(x)=I_current(x) +K*(R_current(x+1)-2*R_current(x)+R_current(x-1))-dt*V(x).*R_current(x);
    end
    % Boundary conditions
    I_next(1)=I_next(2);
    I_next(N)=I_next(N-1);

end

