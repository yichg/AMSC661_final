clc
clear
close all
fsz = 20; % fontsize
% solve schrodinger equation \psi_t = i/2*\psi_xx with spectral method
k_0 = 10;
sigma_0 = 0.1;
T_max = 0.4;
%N = 128;
%N= 4096;
L_x = 20;
N = 256;
x = linspace(-L_x,L_x,N+1);
x(N + 1) = [];

%u = zeros(1,N);
% initial data
A = 1/(2*pi*sigma_0^2)^(1/4);
psi_0 = A*exp(-x.^2/4/sigma_0^2+1j*k_0*x);
% exact solution
tmax = 0.4;
tt = [0,1,2,3,4,5]*tmax/5;
psi_exact = exact_solution(x,tmax,A,sigma_0,k_0);

% %% check the initial condition
% figure
% hold on
% plot(x,abs(psi_exact))
% plot(x,abs(psi_0),"*")
%% spectral method
k = -N/2 : (N/2 - 1);
v0 = fftshift(fft(psi_0));
fac=k.*(2*pi/2/L_x);
fac2=fac.^2;
%eL=exp((-1j/2*fac2)*t);
%psi_spectral = ifft(ifftshift(eL.*v0));

%%
int_spectral = [];
for kk=1:6
    ttt = tt(kk);
    psi_exact = exact_solution(x,ttt,A,sigma_0,k_0);
    eL=exp((-1j/2*fac2)*ttt);
    psi_spectral_tmp = ifft(ifftshift(eL.*v0));
    err = norm(psi_spectral_tmp-psi_exact)/norm(psi_exact);

    figure
    hold on
    plot(x,abs(psi_spectral_tmp))
    plot(x,abs(psi_exact),'--',LineWidth=1)
    legend("Spectral Method","exact solution")
    title(["t = "+num2str(ttt),"  relative error = "+num2str(err*100)+"%"])
    
    hold off
    int_tmp = myint(x,psi_spectral_tmp);
    int_spectral = [int_spectral;int_tmp];

end

figure
hold on
plot(int_spectral)
% plot(ones(length(int_RK),1))
% plot(0*ones(length(int_RK),1))
% plot(2*ones(length(int_RK),1))
%legend(["integral of |\psi^2|","y = 1","y=0","y=2"])
title("check the integral of |\psi^2|")

%%%

% figure
% hold on
% plot(x,abs(psi_exact))
% plot(x,abs(psi_spectral))
% legend("exact solution","Spectral method")

%% D0^2+RK4
tmax=1;
t = 0;
psi_RK = [psi_0];
dt = 0.0001;
u = psi_0;
while (t<tmax) 
    disp(t)
    t=t+dt;
    
    % RK4 step in the Fourier space
    k1=rhs(x,u);
    k2=rhs(x,u+0.5*dt*k1);
    k3=rhs(x,u+0.5*dt*k2);
    k4=rhs(x,u+dt*k3);
    
    unew=u+dt/6*(k1+2*k2+2*k3+k4);
   
    u=unew;
    psi_RK = [psi_RK;u];
    %drawnow
end
%%
tt = [0,1,2,3,4,5]*tmax/5;
for kk=1:6
    ttt = tt(kk);
    psi_exact = exact_solution(x,ttt,A,sigma_0,k_0);
    psi_RK_tmp = psi_RK(1+round(ttt/dt),:);
    err = norm(psi_RK_tmp-psi_exact)/norm(psi_exact);

    figure
    hold on
    plot(x,abs(psi_RK_tmp))
    plot(x,abs(psi_exact),'--',LineWidth=1)
    legend("D0^2+RK4","exact solution")
    title(["t = "+num2str(ttt),"  relative error = "+num2str(err*100)+"%"])
    
    hold off
end
%%
int_RK=[];
for i=1:size(psi_RK,1)
    disp(i)
    int_tmp = myint(x,psi_RK(i,:));
    int_RK = [int_RK;int_tmp];
end
figure
hold on
plot(int_RK)
% plot(ones(length(int_RK),1))
% plot(0*ones(length(int_RK),1))
% plot(2*ones(length(int_RK),1))
%legend(["integral of |\psi^2|","y = 1","y=0","y=2"])
title("check the integral of |\psi^2|")
%%
err = norm(u-psi_exact)
%% check (3) probability density
check_int_spectral = myint(x,psi_spectral)
check_int_exact = myint(x,psi_exact)
%%
check_int_RK = myint(x,u)
%%
err = norm(psi_spectral-psi_exact)
%%





%%
function RHS=rhs(x,psi)
% v should be a row vector
% RHS = i psi_xx
psi=[0,psi,0];
dx = x(2)-x(1);
RHS = 1j/2/dx^2*( psi(1:end-2)+psi(3:end)-2*psi(2:end-1) );

end

function u = exact_solution(x,t,A,sigma,k)
    C = A/sqrt(1+1j*t/2/sigma^2);
    e_up = -x.^2+4*sigma^2*k*1j*x-2*sigma^2*k^2*1j*t;
    e_down = 4*sigma^2*(1+1j*t/2/sigma^2);
    u = C*exp(e_up/e_down);
end

function output = myint(x,psi)
    dx = x(2)-x(1);
    output = sum((psi.*conj(psi))*dx);
    
end