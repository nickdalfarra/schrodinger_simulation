% 4th Order Runge-Kutta method for the Schrodinger Equation of a Particle
% in a Harmonic Oscillator Potential
%
% The solution for a potential that only depends on space, e.g. the HO, can
% be separated in a spatial and a time-dependent part, that is, can be
% solve independently.
%
% In this code, we are using such separation just for the exact and (or)
% initial condition. Therefore we are solving the time dependent version of
% the Schrodinger equation.

% Number of sample points in space
N=50;
% Min and max of point of space
x_min=-4; x_max=4;
% Space step size
dx=(x_max-x_min)/N;

% Time at which we want to end the simulation
t_end=1;
% Time step size
dt=10^-3;
% Number of timesteps to be taken
n_it=t_end/dt;

% Mass of particle
m=1;
% Reduced Planck constant
hbar=1;
% Frequency of the HO potential
omega=1;
% Enery level (why do we need this if we are solving the time dependent version?)
% ok, we need it for the inital condition
n=0;

% Define imaginary unit
%j=sqrt(-1);
% Define constant for kinetic energy term
K=hbar^2/(2*m);

% Spatial grid
x=(x_min:dx:x_max);
% Temporal grid
t=(dt:dt:t_end);
% HO potential at the points of the spatial grid
%V=0.5*m*omega^2*x.^2;

% The initial condition (t=0)
Psi_ini=HO_ti(x,n,m,hbar,omega).*HO_td(0,n,omega);
%
Psi_old=Psi_ini;
dPsidx2=zeros(size(x));
k1=zeros(size(x));
k2=k1;
k3=k1;
k4=k1;

%second derivative operator
D=second_derivative(length(Psi_old),dx);

%errors
errors=zeros(ceil(n_it)-1,1);
energies=errors;
probs=errors;

for iter = 1:n_it
    
    %We calculate the second derivative in space for tn
    % Using periodic boundary conditions
    dPsidx2=(D*Psi_old')';
    
    k1=f(dPsidx2,Psi_old,x);
    k2=f(dPsidx2,Psi_old+(dt/2)*k1,x);
    k3=f(dPsidx2,Psi_old+(dt/2)*k2,x);
    k4=f(dPsidx2,Psi_old+dt*k3,x);
    
    %psi_new
    Psi_new=Psi_old+(dt/6)*(k1+2*k2+2*k3+k4);
    
   %graphical output
    %Plotting every 1000(100) steps for speed
    %if rem(iter, 100) == 0 %if N=20
%     if rem(iter, 1000) == 0 %if N=100
%         plot(x,real(Psi_old),'xb-')
%         hold on
%         plot(x,imag(Psi_old),'or-')
%         plot(x,real(Psi_exact),'b-')
%         plot(x,imag(Psi_exact),'r-')
%         hold off
%         legend("Numerical Real","Numerical Imaginary","Exact Real","Exact Imaginary")
%         axis([-4 4 -1 1])
%         title(['t=',num2str(t(iter))])
%         pause(0.001)
%     end
%    
    %Exact Solution
    Psi_exact = HO_ti(x,n,m,hbar,omega)*HO_td(t(iter),n,omega);
    
    % Document energy and errors
    errors(iter) = dx*norm(Psi_exact-Psi_new,1);
    energies(iter) = dx*abs(-K*Psi_new*D*Psi_new');
    probs(iter) = dx*sum(abs(Psi_new).^2);
    
    Psi_old=Psi_new;
    
end

% Exact Solution
%Psi_exact = HO_ti(x,n,m,hbar,omega)*HO_td(t(iter),n,omega);

% Plots
f1=figure;
plot(x,real(Psi_new),'x-')
hold on
plot(x,imag(Psi_new),'o-')
plot(x,real(Psi_exact),'-')
plot(x,imag(Psi_exact),'-')
hold off
legend("Numerical Re[\Psi]","Numerical Im[\Psi]","Exact Re[\Psi]","Exact Im[\Psi]")
axis([-4 4 -1 1])
xlabel('\it{x}')
title('Quantum Oscillator')
%title(['t=',num2str(t(iter))])

filename=strcat('RK_HO_n_',num2str(n),'_tend_',num2str(t_end),'_N_',num2str(N),'_dx_',num2str(dx),'_dt_',num2str(dt));
print(strcat(filename,'.eps'),'-depsc')
save(strcat(filename,'.mat'),'errors','energies','probs')

% f2=figure;
% subplot(3,1,1),plot(t,errors,'b-')
% title('Error (\it{t})')
% subplot(3,1,2),plot(t,energies,'b-')
% title('Energy (\it{t})')
% subplot(3,1,3),plot(t,probs,'b-')
% title('Probability (\it{t})')
% xlabel('\it{t}')

function [f]=f(dPsidx2,Psi,x)
hbar=1;
m=1;
omega=1;
j=sqrt(-1);
K=hbar^2/(2*m);
V=0.5*m*omega^2*x.^2;
f=-j/hbar*(-K*dPsidx2+ V.*Psi);
end

%matrix to obtain the second derivative. x''=D*x
%using periodic boundary conditions
function D=second_derivative(N,dx)
%N=length(x);
e = ones(N,1);
D = (1/dx^2)*spdiags([e -2*e e], -1:1, N, N);
%periodic boundary conditions
D(1,N)=1;
D(N,1)=1;
end