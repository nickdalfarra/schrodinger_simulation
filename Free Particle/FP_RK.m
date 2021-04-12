% Free Particle using RK4 method
%clear;clc;
% Initilize parameters
x_0=6;
%sigma_squared=1e-2;
sigma_squared=1e-2;
k_0=30;

% Number of sample points in space
N=1400;
% Min and max point in space
x_min=0; x_max=14;
% Space step size
a=x_max-x_min;
dx=a/N;
x=(x_min:dx:x_max);

dt=10^-6;
t_end = 0.1;
n_it=t_end/dt;
t=(dt:dt:0.1);

m = 1;
hbar = 1;
omega = 1;
a = 1;
j = sqrt(-1);
K = hbar^2/(2*m);

% Generate an intial wavepacket
Psi_initial0=exp(-(x-x_0).^2/sigma_squared).*exp(1i*k_0*x);
%Normalizing_parameter = sqrt(dx*sum(conj(Psi_initial0).*Psi_initial0));
Psi_initial = Psi_initial0/sqrt(sqrt(pi*sigma_squared/2)); % normalized

sigma_barrier = 0.01;
V = 1/(2*pi*sigma_barrier^2)*exp(-(x-7).^2/(2*sigma_barrier^2));

Psi_old=Psi_initial;
dPsidx2=zeros(size(x));
k1=zeros(size(x));
k2=k1;
k3=k1;
k4=k1;

%second derivative operator
D=second_derivative(length(Psi_old),dx);

% % Document energy and errors
% errors=zeros(ceil(n_it)-1,1);
% energies=errors;
% probs=errors;

times_list = [1 25000 50000 75000 100000];
counter = 1;

for iter = 1:n_it
    
    %We calculate the second derivative in space for tn
    dPsidx2=(D*Psi_old')';
    
    k1=f(dPsidx2,Psi_old,x);
    k2=f(dPsidx2,Psi_old+(dt/2)*k1,x);
    k3=f(dPsidx2,Psi_old+(dt/2)*k2,x);
    k4=f(dPsidx2,Psi_old+dt*k3,x);
    
    %psi_new
    Psi_new=Psi_old+(dt/6)*(k1+2*k2+2*k3+k4);
    % fixed boundary conditions
    %Psi_new(1)=0;
    %Psi_new(N+1)=0;
    
    prob_density=conj(Psi_new).*Psi_new;
    
    if ismember(iter, times_list)
        %subplot(1,5,counter); hold on;
        plot(x, prob_density)
        hold on
        plot(x, V);
        hold off
        axis([0 14 0 9]);
        title(['Free Particle Problem. RK4 Method. ','t=',num2str(t(iter))])
        xlabel('\it{x}')
        counter = counter + 1;
        %         figure(2)
        %         plot(1:iter, probs(1:iter)); title('probability')
        %pause(0.1)
        filename=strcat('RK_FP_tend_',num2str(t_end),...
            '_N_',num2str(N),'_dx_',num2str(dx),'_dt_',num2str(dt),...
            '_iter_',num2str(iter));
        print(strcat(filename,'.eps'),'-depsc')
        save(strcat(filename,'.mat'),'prob_density')
    end    
    Psi_old=Psi_new;
end



function [f]=f(dPsidx2,Psi,x)
hbar=1;
m=1;
j=sqrt(-1);
K=hbar^2/(2*m);
sigma_barrier = 0.01;
V = 1/(2*pi*sigma_barrier^2)*exp(-(x-7).^2/(2*sigma_barrier^2));
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
%plot(x,V,x,abs(Psi_initial).^2)