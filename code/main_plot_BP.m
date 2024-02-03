% This Matlab script can be used to generate Fig. 2 in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28

clear
clc
%%% system settings
Prms.M = 16;  M = Prms.M; %%% number of transmit/receive antennas
Prms.N = 100; N = Prms.N; %%% number of RIS elements
Prms.K = 2;  K = Prms.K;  %%% number of users
Prms.sigmar2 = 10^(-12); sigmar2 = Prms.sigmar2; %%% radar noise   -90dBm
Prms.sigmat2 = 1; sigmat2 = 1;  %%% target RCS
Prms.sigmak2 = 10^(-12); sigmak2 = Prms.sigmak2; %%% communication noise
Prms.L = 1024; L = Prms.L; %%% number of collected samples
Prms.Nmax = 100; Nmax = Prms.Nmax;  %%% maximum iterations
Prms.res_th = 5e-4;  %%% convergence tolerance
Prms.Crb = 0.01;  %%% CRB threshold

%%%% channel settings
drt = 3; %%% distance of RIS-target
dg = 30; %%% distance of BS-RIS
location_users = [-6 -2;6 -2];  %%% locations of the users
location_RIS = [0,0];  %%% location of RIS
location_BS = [-dg/sqrt(2), -dg/sqrt(2)];  %%% location of BS
drk = zeros(K,1); %%% distance of RIS-user
dk = zeros(K,1);  %%% distance of BS-user
theta_k = zeros(K,1); %%% DoA of the user w.r.t. BS
theta_rk = zeros(K,1); %%% DoA of the user w.r.t. RIS
for k = 1:1:K
    drk(k) = sqrt((location_users(k,1))^2+(location_users(k,2))^2);
    dk(k) = sqrt((location_users(k,1)+dg/sqrt(2))^2+(location_users(k,2)+dg/sqrt(2))^2);
    theta_rk(k) = atan(-location_users(k,1)/location_users(k,2));
    theta_k(k) = asin((dg/sqrt(2)+location_users(k,2))/dk(k));
end
alpha_t = 3.5; %%% path-loss of BS-target
alpha_rt = 2.3; %%% path-loss of RIS-target
alpha_k = 3.5; %%% path-loss of BS-user
alpha_rk = 2.3; %%% path-loss of RIS-user
alpha_g = 2.2; %%% path-loss of BS-RIS

theta2 = pi/4; %%% DoA of the target w.r.t. RIS
thetar = pi/4; %%% DoA of the RIS-BS
theta1 = atan( (dg*sin(thetar)-drt*cos(theta2)) / (dg*cos(thetar)+drt*sin(theta2)) ); %%% DoA of the target w.r.t. BS
dt = (dg*sin(thetar)-drt*cos(theta2) )/sin(theta1); %%% distance of the BS-target

%%% the channel of BS-target
Channel.hdt = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)'*pi*sin(theta1))/sqrt(M); hdt = Channel.hdt;
%%% the partial derivative of hdt w.r.t. theta1
Channel.hdt_der = sqrt(10^(-3)*dt^(-alpha_t))/sqrt(M)*-1j*pi*(0:1:M-1)'*cos(theta1).*exp(-1j*(0:1:M-1)'*pi*sin(theta1));hdt_der = Channel.hdt_der;
%%% the channel of RIS-target
Channel.hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)'*pi*sin(theta2)); hrt = Channel.hrt;
%%% the partial derivative of hrt w.r.t. theta2
Channel.hrt_der = sqrt(10^(-3)*drt^(-alpha_rt))*-1j*pi*(0:1:N-1)'*cos(theta2).*exp(-1j*(0:1:N-1)'*pi*sin(theta2));hrt_der=Channel.hrt_der;
%%%%% assume LoS in Fig. 2
GLos = sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N-1)'*pi*sin(thetar))*exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);
Channel.G = GLos; G = Channel.G;

Prms.gammat = 1; gammat = Prms.gammat; %%% detection SNR threshold
P = 1;  Prms.P = P; %%% power budget

Hu = zeros(K,M); %%% channel of BS-users
Hru = zeros(K,N); %%% channel of RIS-users
for k = 1:1:K
    Hu(k,:) = sqrt(10^(-3)*dk(k)^(-alpha_k))*exp(-1j*(0:1:M-1)*pi*sin(theta_k(k)))/sqrt(M);
    Hru(k,:) = sqrt(10^(-3)*drk(k)^(-alpha_rk))*exp(-1j*(0:1:N-1)*pi*sin(theta_rk(k)));
end
Channel.Hu = Hu;  Channel.Hru = Hru;

A = diag(hrt')*conj(G)*G.'*diag(hrt);
b = 2*diag(hrt')*conj(G)*hdt;
for k = 1:1:K
    A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
    b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
end
phi_CG0 = get_initial_phi(A/norm(b),b/norm(b));  %%% initialize phi
phi = phi_CG0;
C1 = [hdt.'+hrt.'*diag(phi)*G;(Hu+Hru*diag(phi)*G)];
W_CG0 = get_initial_W(C1/norm(C1,'fro'),M,K); %%% initialize W

%%% solve the SNR-constrained joint design problem
[W,phi] = get_W_phi_SNR(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);

%%%% calculate the enhanced beampattern
x_range = (-dg/sqrt(2)-2:0.1:10);
y_range = (-dg/sqrt(2)-2:0.1:0.5);
BP = zeros(length(x_range),length(y_range));
for i = 1:1:length(x_range)
    x = x_range(i);
    for j = 1:1:length(y_range)
        y = y_range(j);
        dd = sqrt((-dg/sqrt(2)-x)^2+(-dg/sqrt(2)-y)^2);
        dr = sqrt(x^2+y^2);
        ang1 = asin(abs((dg/sqrt(2)+y)/dd));
        ang2 = atan(-x/y);
        hd = sqrt(10^(-3)*dd^(-alpha_t))*exp(-1j*(0:1:M-1)'*pi*sin(ang1))/sqrt(M);
        hr = sqrt(10^(-3)*dr^(-alpha_rt))*exp(-1j*(0:1:N-1)'*pi*sin(ang2));
        y = (hd.'+hr.'*diag(phi)*G)*W;
        BP(i,j) = norm(y,2)^2;
    end
end
BP = 10*log10(BP);  %%% in dB
figure;
pcolor(x_range,y_range,BP.');
shading interp;
colorbar;colorbar;
colormap(jet);
hold on;
plot(-dg/sqrt(2),-dg/sqrt(2),'kd','markersize',9.5,'linewidth',1.5)
plot(0,0,'ks','markersize',9.5,'linewidth',1.5)
plot(drt/sqrt(2),-drt/sqrt(2),'kp','markersize',9.5,'linewidth',1.5)
plot(location_users(1,1),location_users(1,2),'ko','markersize',9.5,'linewidth',1.5)
plot(location_users(2,1),location_users(2,2),'ko','markersize',9.5,'linewidth',1.5)
hold off
xlabel('{\it x}-axis (m)'  )
ylabel('{\it y}-axis (m)')

