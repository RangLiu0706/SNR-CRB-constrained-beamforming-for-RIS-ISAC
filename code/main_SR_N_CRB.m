% This Matlab script can be used to generate Fig. 5(b) in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28

clear
clc
%%% system settings
Prms.M = 6;  M = Prms.M; %%% number of transmit/receive antennas
% Prms.N = 49; N = Prms.N; %%% number of RIS elements
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

Prms.Crb = 0.02; %%%% CRB threshold
Prms.P = 10.^(32/10-3); P = Prms.P; %%% power budget

N_range = (24:8:64); %%% number of antennas
rho1_range = [1 2 4 6 8 10]; %%% initial rho
%%% the channel of RIS-target
hrt0 = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(theta2));
%%% the partial derivative of hrt w.r.t. theta2
hrt_der0 = sqrt(10^(-3)*drt^(-alpha_rt))*-1j*pi*(0:1:N_range(end)-1)'*cos(theta2).*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(theta2));
kappa = 10^(3/10);  %%%% 3dB Rician factor
%%% the LoS component of G
GLos0 = sqrt(kappa/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(thetar))*exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);

N_sim = 10000;
%%% achieved sum-rate
SR_my = zeros(1,length(N_range));
SR_max_CG_RIS = zeros(1,length(N_range));
SR_comm_only = zeros(1,length(N_range));
SR_sep = zeros(1,length(N_range));
% CRB_my = zeros(1,length(N_range));
% CRB_max_CG = zeros(1,length(N_range));
% CRB_comm_only = zeros(1,length(N_range));
% CRB_sep = zeros(1,length(N_range));
for sim = 1:1:N_sim
    tic
    sim
    G0 = GLos0 + sqrt(1/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*(randn(N_range(end),M)+1i*randn(N_range(end),M))/sqrt(2*M);
    theta_ru = pi/2*rand(K,1);  %%%% DoA of users w.r.t. RIS
    Hu = zeros(K,M); %%% channel of BS-users
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

    for N_index = 1:1:length(N_range)
        N_index
        N = N_range(N_index); Prms.N = N;

        hrt = hrt0(1:N); Channel.hrt = hrt;
        hrt_der = hrt_der0(1:N); Channel.hrt_der = hrt_der;
        G = G0(1:N,:); Channel.G = G;
        Hru = Hru0(:,1:N); Channel.Hru = Hru;

        %%% initialization
        A = diag(hrt')*conj(G)*G.'*diag(hrt);
        b = 2*diag(hrt')*conj(G)*hdt;
        for k = 1:1:K
            A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
            b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
        end
        phi_CG0 = get_initial_phi(A/norm(b),b/norm(b));
        C1 = [hdt.'+hrt.'*diag(phi_CG0)*G;(Hu+Hru*diag(phi_CG0)*G)];
        W_CG0 = get_initial_W(C1/norm(C1,'fro'),M,K);

        A = zeros(N,N);
        b = zeros(N,1);
        for k = 1:1:K
            A = A + diag(conj(Hru(k,:)))*conj(G)*G.'*diag(Hru(k,:).');
            b = b + 2*diag(conj(Hru(k,:)))*conj(G)*Hu(k,:).';
        end
        phi_c = get_initial_phi(A/norm(b),b/norm(b));
        C1 = Hu+Hru*diag(phi_c)*G;
        Wc = get_initial_W(C1/norm(C1,'fro'),M,K);

        [W_my,phi_my,sr,crb] = get_W_phi_CRB(Prms,Channel,phi_CG0,[W_CG0 zeros(M,M)]*sqrt(P),rho1_range(N_index));
        SR_my(N_index) = SR_my(N_index) + sr(end);
        % CRB_my(N_index) = CRB_my(N_index) + crb;

        [W_CG,sr_CG,crb_CG] = get_W_with_phi_CRB(Prms,Channel,phi_CG0,[W_CG0 zeros(M,M)]*sqrt(P),rho1_range(N_index));
        SR_max_CG_RIS(N_index) = SR_max_CG_RIS(N_index) + sr_CG(end);
        % CRB_max_CG(N_index) = CRB_max_CG(N_index) + crb_CG;

        [W_comm,phi_comm,sr_comm,~,crb_comm] = get_W_phi_comm(Prms,Channel,phi_CG0,W_CG0*sqrt(P));
        SR_comm_only(N_index) = SR_comm_only(N_index) + sr_comm(end);
        % CRB_comm_only(N_index) = CRB_comm_only(N_index) + crb_comm;

        [W_sep,phi_sep,sr_sep,crb_sep] = get_W_phi_sep_CRB(Prms,Channel,[W_CG0 zeros(M,M)]*sqrt(P));
        SR_sep(N_index) = SR_sep(N_index) + sr_sep(end);
        % CRB_sep(N_index) = CRB_sep(N_index) + crb_sep;
    end
    toc
end

SR_my = SR_my/sim;
SR_max_CG_RIS = SR_max_CG_RIS/sim;
SR_comm_only = SR_comm_only/sim;
SR_sep = SR_sep/sim;
% CRB_my = CRB_my/sim;
% CRB_max_CG = CRB_max_CG/sim;
% CRB_comm_only = CRB_comm_only/sim;
% CRB_sep = CRB_sep/sim;

figure;
plot(N_range,SR_my,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on;
plot(N_range,SR_max_CG_RIS,'-d','color',[0,0.5,0],'LineWidth',1.5)
plot(N_range,SR_comm_only,'-x','color',[0,0,0],'LineWidth',1.5)
plot(N_range,SR_sep,'-s','color',[0,0,0.8],'LineWidth',1.5)
hold off
xlabel('Number of reflecting elements {\it N}');
ylabel('Sum-rate {\it R} (bps/Hz)');
grid on
legend('Proposed','BF only', 'Comm only','Separate');
axis([N_range(1) N_range(end) min(SR_sep)-1 max(SR_comm_only)+1])