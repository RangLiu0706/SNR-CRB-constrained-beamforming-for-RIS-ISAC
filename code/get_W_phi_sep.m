% Optimize the SNR-constrained joint beamforming and reflection design separately.
% Generate the comparison "Separate".
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; gammat_r: the achieved radar SNR
function [W,phi,Vsr,gammat_r] = get_W_phi_sep(Prms,Channel)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
gammat = Prms.gammat; P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;

A = diag(hrt')*conj(G)*G.'*diag(hrt);
b = 2*diag(hrt')*conj(G)*hdt;
phi = get_initial_phi(A/norm(b),b/norm(b));
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);

Wr = get_initial_W(Ht/norm(Ht,'fro'),M,M);
w = vec(Wr);
u = kron(eye(M),Ht)*w/(w'*(kron(eye(M),Ht'*Ht))*w);
vSNRt = nL*sigmat2*abs(u'*kron(eye(M),Ht)*w)^2/(sigmar2*u'*u);
w = sqrt(gammat/vSNRt)*w;

Wr = reshape(w,M,M);
P = P - norm(w,2)^2;

Hk = Hu + Hru*diag(phi)*G;
Wc = sqrt(P)*get_initial_W(Hk/norm(Hk,'fro'),M,K);
W = [Wc Wr];

r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end
%%%% start iteration
iter = 1;
res = 1;
Vsr = zeros(Nmax,1);
while iter <= Nmax && res >= res_th
    A = zeros(K,M);
    B = zeros(K,M);
    for k = 1:1:K
        A(k,:) = 2*sqrt(1+r(k))*c(k)'*Hk(k,:);
        B(k,:) = abs(c(k))*Hk(k,:);
    end
    cvx_begin quiet
    variable Wc(M,K) complex
    minimize square_pos(norm(B*Wc,'fro'))-real(trace(A*Wc))
    subject to
    norm(Wc,'fro') <= sqrt(P);
    cvx_end

    W = [Wc Wr];
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    if iter > 1
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter);
    end
    iter = iter + 1;
end
Vsr(iter:end) = [];
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
u = kron(eye(K+M),Ht)*vec(W)/((vec(W))'*(kron(eye(K+M),Ht'*Ht))*vec(W));
gammat_r = real(nL*sigmat2*abs(u'*kron(eye(K+M),Ht)*vec(W))^2/(sigmar2*u'*u));
% figure;grid on
% plot(VSR);title('sum-rate');legend('Sep')
% figure;grid on
% plot(VSR);title('sum-rate');legend('Separate')
