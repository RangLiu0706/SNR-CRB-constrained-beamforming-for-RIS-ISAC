% This Matlab script can be used to generate Fig. 3 in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28

clear
clc
%%% system settings
Prms.M = 6;  M = Prms.M; %%% number of transmit/receive antennas
Prms.K = 4;  K = Prms.K; %%% number of users
Prms.sigmar2 = 10^(-12); sigmar2 = Prms.sigmar2; %%% radar noise   -90dBm
Prms.sigmat2 = 1; sigmat2 = 1;  %%% target RCS
Prms.sigmak2 = 10^(-12); sigmak2 = Prms.sigmak2; %%% communication noise
Prms.L = 1024; L = Prms.L; %%% number of collected samples
Prms.Nmax = 100; Nmax = Prms.Nmax;  %%% maximum iterations
Prms.res_th = 0;  %%% convergence tolerance

%%%% channel settings
drt = 3; %%% distance of RIS-target
dg = 50; %%% distance of BS-RIS
drk = 8; %%% distance of RIS-user
alpha_t = 2.4; %%% path-loss of BS-target
alpha_rt = 2.2; %%% path-loss of RIS-target
alpha_k = 3.5; %%% path-loss of BS-user
alpha_rk = 2.3; %%% path-loss of RIS-user
alpha_g = 2.2; %%% path-loss of BS-RIS

theta2 = pi/4; %%% DoA of the target w.r.t. RIS
thetar = pi/4; %%% DoA of the RIS-BS
%%% DoA of the target w.r.t. BS
theta1 = atan( (dg*sin(thetar)-drt*cos(theta2)) / (dg*cos(thetar)+drt*sin(theta2)) );
dt = (dg*sin(thetar)-drt*cos(theta2) )/sin(theta1);  %%% distance of the BS-target

%%% the channel of BS-target
Channel.hdt = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)'*pi*sin(theta1))/sqrt(M); hdt = Channel.hdt;
%%% the partial derivative of hdt w.r.t. theta1
Channel.hdt_der = sqrt(10^(-3)*dt^(-alpha_t))/sqrt(M)*-1j*pi*(0:1:M-1)'*cos(theta1).*exp(-1j*(0:1:M-1)'*pi*sin(theta1));hdt_der = Channel.hdt_der;

N_range = [30 60];
%%% the channel of RIS-target
hrt0 = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(theta2));
%%% the partial derivative of hrt w.r.t. theta2
hrt_der0 = sqrt(10^(-3)*drt^(-alpha_rt))*-1j*pi*(0:1:N_range(end)-1)'*cos(theta2).*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(theta2));
kappa = 10^(3/10);  %%%% 3dB Rician factor
%%% the LoS component of G
GLos0 = sqrt(kappa/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(thetar))*exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);

N_sim = 1000;
%%% sum-rate of the SNR-constrained design
SR_my1 = zeros(1,Nmax);
SR_my2 = zeros(1,Nmax);
SR_my3 = zeros(1,Nmax);

%%% sum-rate of the CRB-constrained design
SR_my1_CRB = zeros(1,Nmax);
SR_my2_CRB = zeros(1,Nmax);
SR_my3_CRB = zeros(1,Nmax);

for sim = 1:1:N_sim
    tic
    sim

    G0 = GLos0 + sqrt(1/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*(randn(N_range(end),M)+1i*randn(N_range(end),M))/sqrt(2*M);
    theta_ru = pi/2*rand(K,1);  %%%% DoA of users w.r.t. RIS
    Hu = zeros(K,M);  %%% channel of BS-users
    Hru0 = zeros(K,N_range(end)); %%% channel of RIS-users
    for k = 1:1:K
        theta_rk = theta_ru(k);
        theta_k = atan( (dg*sin(thetar)-drk*cos(theta_rk))/(dg*cos(thetar)+drk*sin(theta_rk)) );
        dk = (dg*sin(thetar)-drk*cos(theta_rk))/sin(theta_k);
        Hu(k,:) = sqrt(10^(-3)*dk^(-alpha_k))*( sqrt(kappa/(1+kappa))*exp(-1j*(0:1:M-1)*pi*sin(theta_k))/sqrt(M) ...
            + sqrt(1/(1+kappa))*(randn(1,M)+1i*randn(1,M))/sqrt(2*M) );
        Hru0(k,:) = sqrt(10^(-3)*drk^(-alpha_rk))*( sqrt(kappa/(1+kappa))*exp(-1j*(0:1:N_range(end)-1)*pi*sin(theta_rk)) ...
            + sqrt(1/(1+kappa))*(randn(1,N_range(end))+1i*randn(1,N_range(end)))/sqrt(2) );
    end
    Channel.Hu = Hu;

    N = 30; Prms.N = N;
    hrt = hrt0(1:N); Channel.hrt = hrt;
    hrt_der = hrt_der0(1:N); Channel.hrt_der = hrt_der;
    G = G0(1:N,:); Channel.G = G;
    Hru = Hru0(:,1:N); Channel.Hru = Hru;

    A = diag(hrt')*conj(G)*G.'*diag(hrt);
    b = 2*diag(hrt')*conj(G)*hdt;
    for k = 1:1:K
        A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
        b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
    end
    phi_CG0 = get_initial_phi(A/norm(b),b/norm(b));
    C1 = [hdt.'+hrt.'*diag(phi_CG0)*G;(Hu+Hru*diag(phi_CG0)*G)];
    W_CG0 = get_initial_W(C1/norm(C1,'fro'),M,K);

    P = 10^(0.1*32-3); Prms.P = P;

    Prms.gammat = 1; gammat = Prms.gammat;
    [~,~,vsr1] = get_W_phi_SNR(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);
    SR_my1 = SR_my1 + vsr1;

    P = 10^(0.1*28-3); Prms.P = P;
    [~,~,vsr2] = get_W_phi_SNR(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);
    SR_my2 = SR_my2 + vsr2;

    Prms.Crb = 0.01;
    [~,~,vsr21] = get_W_phi_CRB(Prms,Channel,phi_CG0,[W_CG0 zeros(M,M)]*sqrt(P),1);
    SR_my1_CRB = SR_my1_CRB + vsr21;

    N = 60; Prms.N = N;
    hrt = hrt0(1:N); Channel.hrt = hrt;
    hrt_der = hrt_der0(1:N); Channel.hrt_der = hrt_der;
    G = G0(1:N,:); Channel.G = G;
    Hru = Hru0(:,1:N); Channel.Hru = Hru;

    A = diag(hrt')*conj(G)*G.'*diag(hrt);
    b = 2*diag(hrt')*conj(G)*hdt;
    for k = 1:1:K
        A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
        b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
    end
    phi_CG0 = get_initial_phi(A/norm(b),b/norm(b));
    C1 = [hdt.'+hrt.'*diag(phi_CG0)*G;(Hu+Hru*diag(phi_CG0)*G)];
    W_CG0 = get_initial_W(C1/norm(C1,'fro'),M,K);

    P = 10^(0.1*28-3); Prms.P = P;

    Prms.gammat = 10; gammat = Prms.gammat;
    [~,~,vsr3] = get_W_phi_SNR(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);
    SR_my3 = SR_my3 + vsr3;

    [~,~,vsr22] = get_W_phi_CRB(Prms,Channel,phi_CG0,[W_CG0 zeros(M,M)]*sqrt(P),0.5);
    SR_my2_CRB = SR_my2_CRB + vsr22;

    Prms.Crb = 0.02;
    [~,~,vsr23] = get_W_phi_CRB(Prms,Channel,phi_CG0,[W_CG0 zeros(M,M)]*sqrt(P),0.5);
    SR_my3_CRB = SR_my3_CRB + vsr23;

    toc
end

SR_my1 = SR_my1/sim;
SR_my2 = SR_my2/sim;
SR_my3 = SR_my3/sim;

SR_my1_CRB = SR_my1_CRB/sim;
SR_my2_CRB = SR_my2_CRB/sim;
SR_my3_CRB = SR_my3_CRB/sim;

figure;
plot(SR_my1,'-','color',[0.8,0,0],'LineWidth',1.5)
hold on
plot(SR_my2,'-','color',[0,0.5,0],'LineWidth',1.5)
plot(SR_my3,'-','color',[0,0,0.8],'LineWidth',1.5)
hold off
xlabel('Number of iterations');
ylabel('Sum-rate {\it R} (bps/Hz)');
grid on
legend('{\it P} = 32dBm, {\it N} = 30, \Gamma_t = 0dB','{\it P} = 28dBm, {\it N} = 30, \Gamma_t = 0dB', ...%'{\it P} = 28dBm, {\it N} = 60, \Gamma_t = 0dB',
    '{\it P} = 28dBm, {\it N} = 60, \Gamma_t = 10dB');
axis([1 100 4 28])

figure
plot(SR_my1_CRB,'-','color',[0,0,0.8],'LineWidth',1.5)
hold on
plot(SR_my2_CRB,'-','color',[0,0.5,0],'LineWidth',1.5)
plot(SR_my3_CRB,'-','color',[0.8,0,0],'LineWidth',1.5)
hold off
xlabel('Number of iterations');
ylabel('Sum-rate {\it R} (bps/Hz)');
grid on
legend('{\it P} = 28dBm, {\it N} = 30, \epsilon=0.01', '{\it P} = 28dBm, {\it N} = 60, \epsilon=0.01',...
    '{\it P} = 28dBm, {\it N} = 60, \epsilon=0.02');
axis([1 100 2 21])

