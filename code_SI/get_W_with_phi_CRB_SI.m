% Optimize the CRB-constrained beamforming design given phi with self interference.
% Generate the comparison "BF only".
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         phi: the phi; W: the initial W; rho: the initial rho
% Outputs: W: transmit beamforming;
%          Vsr: the achieved sum-rate; crb: the achieved CRB
function [W,Vsr,crb] = get_W_with_phi_CRB_SI(Prms,Channel,phi,W,rho)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
Nmax = Prms.Nmax; res_th = Prms.res_th; L = Prms.L;
P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
hdt_der = Channel.hdt_der; hrt_der = Channel.hrt_der;
Hsi = Channel.Hsi;
CRB = Prms.Crb;

Hk = Hu + Hru*diag(phi)*G;
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;
w = vec(W);
Rw = W*W';
Rn = eye(M) + Hsi*Rw*Hsi'/sigmar2;
F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht1'*inv(Rn));trace(Ht2*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht2'*inv(Rn))]);
F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1'*inv(Rn))*[1 1j] );real( trace(Ht*Rw*Ht2'*inv(Rn))*[1 1j] )];
F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'*inv(Rn)))*eye(2);
J = F1-F2*inv(F4)*F2.';
Crbo = trace(inv(J));

ft = zeros(6,1);
ft(1) = 2/sigmar2*trace(Ht1*Rw*Ht1'*inv(Rn)); ft(2) = 2/sigmar2*trace(Ht2*Rw*Ht1'*inv(Rn));
ft(3) = 2/sigmar2*trace(Ht*Rw*Ht1'*inv(Rn)); ft(4) = 2/sigmar2*trace(Ht2*Rw*Ht2'*inv(Rn));
ft(5) = 2/sigmar2*trace(Ht*Rw*Ht2'*inv(Rn)); ft(6) = 2/sigmar2*trace(Ht*Rw*Ht'*inv(Rn));
f = ft;
r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end
zeta = zeros(6,1);

a = zeros(M*(K+M),1);
B = zeros(K*(K+M),M*(K+M));
for k = 1:1:K
    a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
    for j = 1:1:K+M
        Tj = zeros(M,M*(K+M)); Tj(:,(j-1)*M+1:j*M) = eye(M);
        B((k-1)*(K+M)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
    end
end

%%%% start iteration
iter = 1;
res = 1;
rho1 = norm(f,2)^2/2/abs(norm(B*w,2)^2-real(a'*w))*rho;

Vres = zeros(Nmax,1);
Vsr = zeros(Nmax,1);
VCRB = zeros(Nmax,2);
while iter <= Nmax && res >= res_th

    cvx_begin quiet
    variable J(2,2) semidefinite
    variable f(6,1) complex
    expression F1(2,2)
    expression F2(2,2)
    expression F4(2,2)
    minimize square_pos(norm(ft-f+rho1*zeta,2))
    subject to
    trace_inv(J) <= CRB*L;
    F1 = real([f(1) f(2);f(2) f(4)]);
    F2 = real([f(3);f(5)]*[1 1j]);
    F4 = f(6)*eye(2);
    [F1-J  F2; F2.' F4] == semidefinite(4);
    cvx_end

    VCRB(iter,1) = trace(inv(J))/L;
    %%%% update w
    Ai = zeros(M,M,6);
    Ai(:,:,1) = 2/sigmar2*Ht1'*inv(Rn)*Ht1; Ai(:,:,2) = 2/sigmar2*Ht1'*inv(Rn)*Ht2;
    Ai(:,:,3) = 2/sigmar2*Ht1'*inv(Rn)*Ht; Ai(:,:,4) = 2/sigmar2*Ht2'*inv(Rn)*Ht2;
    Ai(:,:,5) = 2/sigmar2*Ht2'*inv(Rn)*Ht;  Ai(:,:,6) = 2/sigmar2*Ht'*inv(Rn)*Ht;

    a = zeros(M*(K+M),1);
    B = zeros(K*(K+M),M*(K+M));
    for k = 1:1:K
        a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
        for j = 1:1:K+M
            Tj = zeros(M,M*(K+M)); Tj(:,(j-1)*M+1:j*M) = eye(M);
            B((k-1)*(K+M)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
        end
    end

    B1 = zeros(M*(K+M),M*(K+M));  b1 = zeros(1,M*(K+M));
    for i = 1:1:6
        T1 = Ai(:,:,i)*W;  T2 = Ai(:,:,i)'*W;
        B1 = B1 + vec(T1)*vec(T1)' + vec(T2)*vec(T2)';
        di = rho1*zeta(i)-f(i);
        T3 = di'*Ai(:,:,i) + di*Ai(:,:,i)';
        [~,ev1] = eigsort(T3);
        [~,ev2] = eigsort(Ai(:,:,i));
        lam1 = max(diag(ev1));
        lam2 = max( abs(diag(ev2)).^2 );
        b1 = b1 + 2*vec(T3*W)' - w'*(2*lam1*eye(M*(M+K))+4*lam2*w*w');
    end
    B1 = B1/(2*rho1);
    b1 = b1/(2*rho1);

    cvx_begin quiet
    variable w(M*(M+K),1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w) + real(b1*w) + real(w'*B1*w)
    subject to
    norm(w,2) <= sqrt(P);
    cvx_end

    W = reshape(w,M,M+K);
    Rw = W*W';
    Rn = eye(M) + Hsi*Rw*Hsi'/sigmar2;
    F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht1'*inv(Rn));trace(Ht2*Rw*Ht1'*inv(Rn)) trace(Ht2*Rw*Ht2'*inv(Rn))]);
    F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1'*inv(Rn))*[1 1j] );real( trace(Ht*Rw*Ht2'*inv(Rn))*[1 1j] )];
    F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'*inv(Rn)))*eye(2);
    J = F1-F2*inv(F4)*F2.';
    VCRB(iter,2) = trace(inv(J));

    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    %%%% update dual varialbes
    zeta = zeta + (ft-f)/rho1;
    ft = zeros(6,1);
    ft(1) = 2/sigmar2*trace(Ht1*Rw*Ht1'*inv(Rn));  ft(2) = 2/sigmar2*trace(Ht2*Rw*Ht1'*inv(Rn));
    ft(3) = 2/sigmar2*trace(Ht*Rw*Ht1'*inv(Rn));   ft(4) = 2/sigmar2*trace(Ht2*Rw*Ht2'*inv(Rn));
    ft(5) = 2/sigmar2*trace(Ht*Rw*Ht2'*inv(Rn));   ft(6) = 2/sigmar2*trace(Ht*Rw*Ht'*inv(Rn));

    Vsr(iter) = sum(log2(1+r));
    Vres(iter) = norm(f-ft,2);
    if iter > 1
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter);
    end
    rho1 = rho1*0.9;
    if VCRB(iter,2) > 0.8*CRB
        rho1 = rho1*0.7;
    elseif VCRB(iter,2) > 0.9*CRB
        rho1 = rho1*0.1;
    end
    if VCRB(iter,2) > CRB
        res = 0;
        iter = iter - 1;
    end
    iter = iter + 1;
end
Vsr(iter:end) = [];
VCRB(iter:end,:) = [];
Vres(iter:end,:) = [];
crb = VCRB(end,2);





