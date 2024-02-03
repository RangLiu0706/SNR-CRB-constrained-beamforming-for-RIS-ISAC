% Optimize the CRB-constrained joint beamforming and reflection design separately.
% Generate the comparison "Separate".
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels; W: the initial W
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; crb: the achieved CRB
function [W,phi,Vsr,crb] = get_W_phi_sep_CRB(Prms,Channel,W)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
Nmax = Prms.Nmax; res_th = Prms.res_th; L = Prms.L;
P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
hdt_der = Channel.hdt_der; hrt_der = Channel.hrt_der;
CRB = Prms.Crb;

W(:,M+1:end) = [];
A = diag(hrt')*conj(G)*G.'*diag(hrt);
b = 2*diag(hrt')*conj(G)*hdt;
phi = get_initial_phi(A/norm(b),b/norm(b));

Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;

A1 = Ht1'*Ht1/sigmar2*2; A2 = Ht1'*Ht2/sigmar2*2; A3 = Ht1'*Ht/sigmar2*2;
A4 = Ht2'*Ht2/sigmar2*2; A5 = Ht2'*Ht/sigmar2*2;  A6 = Ht'*Ht/sigmar2*2;

cvx_begin quiet
variable R(M,M) complex semidefinite hermitian
variable J(2,2) semidefinite
expression F1(2,2)
expression F2(2,2)
expression F4(2,2)
minimize real(trace(R))
subject to
real(trace_inv(J)) <= CRB*L;
F1 = real([trace(A1*R)  trace(A2*R); trace(A2*R) trace(A4*R)]);
F2 = real([trace(A3*R);trace(A5*R)]*[1 1j]);
F4 = real(trace(A6*R))*eye(2);
[F1-J  F2; F2' F4] == semidefinite(4);
cvx_end

Wr = (chol(R))';
P = P-real(trace(R));
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
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;
Rw = W*W';
F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1') trace(Ht2*Rw*Ht1');trace(Ht2*Rw*Ht1') trace(Ht2*Rw*Ht2')]);
F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1')*[1 1j] );real( trace(Ht*Rw*Ht2')*[1 1j] )];
F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'))*eye(2);
J = F1-F2*inv(F4)*F2.';
crb = trace(inv(J));

% figure;grid on
% plot(Vsr);title('sum-rate');legend('Sep')
