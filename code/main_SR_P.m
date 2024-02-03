% This Matlab script can be used to generate Fig. 4(a) in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28

clear
clc
%%% system settings
Prms.M = 6;  M = Prms.M; %%% number of transmit/receive antennas
Prms.N = 49; N = Prms.N; %%% number of RIS elements
Prms.K = 4;  K = Prms.K; %%% number of users
Prms.sigmar2 = 10^(-12); sigmar2 = Prms.sigmar2; %%% radar noise   -90dBm
Prms.sigmat2 = 1; sigmat2 = 1;  %%% target RCS
Prms.sigmak2 = 10^(-12); sigmak2 = Prms.sigmak2; %%% communication noise
Prms.L = 1024; L = Prms.L; %%% number of collected samples
Prms.Nmax = 100; Nmax = Prms.Nmax;  %%% maximum iterations
Prms.res_th = 1e-4;  %%% convergence tolerance

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
%%% the channel of RIS-target
Channel.hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)'*pi*sin(theta2)); hrt = Channel.hrt;
%%% the partial derivative of hrt w.r.t. theta2
Channel.hrt_der = sqrt(10^(-3)*drt^(-alpha_rt))*-1j*pi*(0:1:N-1)'*cos(theta2).*exp(-1j*(0:1:N-1)'*pi*sin(theta2));hrt_der=Channel.hrt_der;
kappa = 10^(3/10);  %%%% 3dB Rician factor
%%% the LoS component of G
GLos = sqrt(kappa/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N-1)'*pi*sin(thetar))*exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);

Prms.gammat = 10^0.7; gammat = Prms.gammat;  %%% detection SNR threshold 7dB
P_dbm_range = [27:1:32]; %%% in dBm
P_W_range = 10.^(P_dbm_range/10-3); %%% in W

N_sim = 10000;
%%% achieved sum-rate
SR_my = zeros(1,length(P_W_range));  %%% proposed
SR_max_CG_RIS = zeros(1,length(P_W_range)); %%% BF only
SR_comm_only = zeros(1,length(P_W_range)); %%% Comm only
SR_sep = zeros(1,length(P_W_range)); %%% Separate

%%% achieved SNR
% Gammat_my = zeros(1,length(P_W_range));
% Gammat_max_CG = zeros(1,length(P_W_range));
% Gammat_comm_only = zeros(1,length(P_W_range));
% Gammat_sep  = zeros(1,length(P_W_range));

for sim = 1:1:N_sim
    tic
    sim

    Channel.G = GLos + sqrt(1/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*(randn(N,M)+1i*randn(N,M))/sqrt(2*M);
    G = Channel.G;

    theta_ru = pi/2*rand(K,1);  %%%% DoA of users w.r.t. RIS
    Hu = zeros(K,M);  %%% channel of BS-users
    Hru = zeros(K,N); %%% channel of RIS-users
    for k = 1:1:K
        theta_rk = theta_ru(k);
        theta_k = atan( (dg*sin(thetar)-drk*cos(theta_rk))/(dg*cos(thetar)+drk*sin(theta_rk)) );
        dk = (dg*sin(thetar)-drk*cos(theta_rk))/sin(theta_k);
        Hu(k,:) = sqrt(10^(-3)*dk^(-alpha_k))*( sqrt(kappa/(1+kappa))*exp(-1j*(0:1:M-1)*pi*sin(theta_k))/sqrt(M) ...
            + sqrt(1/(1+kappa))*(randn(1,M)+1i*randn(1,M))/sqrt(2*M) );
        Hru(k,:) = sqrt(10^(-3)*drk^(-alpha_rk))*( sqrt(kappa/(1+kappa))*exp(-1j*(0:1:N-1)*pi*sin(theta_rk)) ...
            + sqrt(1/(1+kappa))*(randn(1,N)+1i*randn(1,N))/sqrt(2) );
    end
    Channel.Hu = Hu;  Channel.Hru = Hru;

    %%% initialization
    A = diag(hrt')*conj(G)*G.'*diag(hrt);
    b = 2*diag(hrt')*conj(G)*hdt;
    for k = 1:1:K
        A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
        b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
    end
    phi_CG0 = get_initial_phi(A/norm(b),b/norm(b));
    phi = phi_CG0;
    C1 = [hdt.'+hrt.'*diag(phi)*G;(Hu+Hru*diag(phi)*G)];
    W_CG0 = get_initial_W(C1/norm(C1,'fro'),M,K);

    A = zeros(N,N);
    b = zeros(N,1);
    for k = 1:1:K
        A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
        b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
    end
    phi_c = get_initial_phi(A/norm(b),b/norm(b));
    C1 = Hu + Hru*diag(phi_c)*G;
    Wc = get_initial_W(C1/norm(C1,'fro'),M,K);

    for P_index = 1:1:length(P_W_range)
        P_index
        P = P_W_range(P_index); Prms.P = P;

        [W_my,phi_my,sr_my,gammat_my] = get_W_phi_SNR(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);
        SR_my(P_index) = SR_my(P_index) + sr_my(end);
        % Gammat_my(P_index) = Gammat_my(P_index) + gammat_my;        %

        [W_CG,sr_CG,gammat_CG] = get_W_with_phi(Prms,Channel,phi_CG0,sqrt(P)*[W_CG0 zeros(M,M)]);
        SR_max_CG_RIS(P_index) = SR_max_CG_RIS(P_index) + sr_CG(end);
        % Gammat_max_CG(P_index) = Gammat_max_CG(P_index) + gammat_CG;

        [W_comm,phi_comm,sr_comm,gammat_comm] = get_W_phi_comm(Prms,Channel,phi_CG0,sqrt(P)*W_CG0);
        SR_comm_only(P_index) = SR_comm_only(P_index) + sr_comm(end);
        % Gammat_comm_only(P_index) = Gammat_comm_only(P_index) + gammat_comm;

        [W_sep,phi_sep,sr_sep,gammat_sep] = get_W_phi_sep(Prms,Channel);
        SR_sep(P_index) = SR_sep(P_index) + sr_sep(end);
        % Gammat_sep(P_index) = Gammat_sep(P_index) + gammat_sep;

    end
    toc
end

SR_my = SR_my/sim;
SR_max_CG_RIS = SR_max_CG_RIS/sim;
SR_comm_only = SR_comm_only/sim;
SR_sep = SR_sep/sim;

% Gammat_my = Gammat_my/sim;
% Gammat_max_CG = Gammat_max_CG/sim;
% Gammat_comm_only = Gammat_comm_only/sim;
% Gammat_sep  = Gammat_sep/sim;

figure
plot(P_dbm_range,SR_my,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on
plot(P_dbm_range,SR_max_CG_RIS,'-d','color',[0,0.5,0],'LineWidth',1.5)
plot(P_dbm_range,SR_comm_only,'-x','color',[0,0,0],'LineWidth',1.5)
plot(P_dbm_range,SR_sep,'-s','color',[0,0,0.8],'LineWidth',1.5)
hold off
xlabel('Transmit power {\it P}_t (dBm)');
ylabel('Sum-rate {\it R} (bps/Hz)');
grid on
legend('Proposed','BF only', 'Comm only','Separate');
