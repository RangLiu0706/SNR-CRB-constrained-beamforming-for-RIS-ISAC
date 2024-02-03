% Optimize the CRB-constrained joint beamforming and reflection design
% separately with self interference.
% Generate the comparison "Separate".
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels; W: the initial W
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; crb: the achieved CRB
function [W,phi,Vsr,crb] = get_W_phi_sep_CRB_SI(Prms,Channel,W)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
Nmax = Prms.Nmax; res_th = Prms.res_th; L = Prms.L;
P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
hdt_der = Channel.hdt_der; hrt_der = Channel.hrt_der;
CRB = Prms.Crb; Hsi = Channel.Hsi;

W(:,M+1:end) = [];
A = diag(hrt')*conj(G)*G.'*diag(hrt);
b = 2*diag(hrt')*conj(G)*hdt;
phi = get_initial_phi(A/norm(b),b/norm(b));

Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;

R = W*W';
iter = 1; res = 1;
Vp = zeros(1,Nmax);
while iter < 5 && res > res_th

    Rn = eye(M) + Hsi*R*Hsi'/sigmar2;
    A1 = 2/sigmar2*Ht1'*inv(Rn)*Ht1; A2 = 2/sigmar2*Ht1'*inv(Rn)*Ht2;
    A3 = 2/sigmar2*Ht1'*inv(Rn)*Ht; A4 = 2/sigmar2*Ht2'*inv(Rn)*Ht2;
    A5 = 2/sigmar2*Ht2'*inv(Rn)*Ht;  A6 = 2/sigmar2*Ht'*inv(Rn)*Ht;

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

    Vp(iter) = real(trace(R));
    if iter > 1
        res = abs(1-Vp(iter)/Vp(iter-1));
    end
    iter = iter + 1;
end

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

    %%%% update dual varialbes
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
Rn = eye(M) + Hsi*Rw*Hsi'/sigmar2;
F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht1'*inv(Rn));trace(Ht2*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht2'*inv(Rn))]);
F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1'*inv(Rn))*[1 1j] );real( trace(Ht*Rw*Ht2'*inv(Rn))*[1 1j] )];
F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'*inv(Rn)))*eye(2);
J = F1-F2*inv(F4)*F2.';
crb = trace(inv(J));
