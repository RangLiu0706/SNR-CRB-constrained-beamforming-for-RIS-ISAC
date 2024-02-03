% Optimize the SNR-constrained joint beamforming and reflection design problem with self-interference.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         phi: the initial phi; W: the initial W
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; gammat_r: the achieved radar SNR
function [W,phi,Vsr,gammat_r] = get_W_phi_SNR_SI(Prms,Channel,phi,W)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
gammat = Prms.gammat; P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
Hsi = Channel.Hsi;

%%%% variable follow phi
Hk = Hu + Hru*diag(phi)*G;
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end

D = zeros(N,N);
g = zeros(N,1);
for k = 1:1:K
    g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
    for j = 1:1:K+M
        temp = diag(W(:,j)'*G')*Hru(k,:)';
        D = D + abs(c(k))^2*(temp*temp');
        g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
    end
end

%%%% start iteration
iter = 1;
res = 1;
Vsr = zeros(1,Nmax);
Vres = zeros(1,Nmax);
Vgammat = zeros(1,Nmax);
mu = zeros(N,1);
varphi = phi;
rho = N/abs(phi'*D*phi-real(g'*phi));
while iter <= Nmax && res >= res_th

    Mw = kron(eye(M+K),Hsi)*vec(W)*(vec(W))'*kron(eye(M+K),Hsi') + sigmar2*eye(M*(M+K));
    s_phi_w = kron(eye(M+K),Ht)*vec(W);
    Dt = kron(eye(K+M),Hsi')*inv(Mw)*s_phi_w*s_phi_w'*inv(Mw)*kron(eye(K+M),Hsi);
    gt = inv(Mw)*s_phi_w;
    ct = real(trace(sigmar2*inv(Mw)*s_phi_w*s_phi_w'*inv(Mw)));
    e3 = real(0.5*(gammat/nL/sigmat2+(vec(W))'*Dt*vec(W)+ct));

    %%%% update reflecting coefficients phi
    D = zeros(N,N);
    g = zeros(N,1);
    for k = 1:1:K
        g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
        for j = 1:1:K+M
            temp = diag(W(:,j)'*G')*Hru(k,:)';
            D = D + abs(c(k))^2*(temp*temp');
            g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
        end
    end
    F = kron(W.'*hdt,G.'*diag(hrt)) + kron(W.'*G.'*diag(hrt),hdt);
    L = kron(W.'*G.'*diag(hrt),G.'*diag(hrt));
    Ltilde = reshape(L.'*conj(gt),N,N);
    Lbar = [-real(Ltilde) imag(Ltilde);imag(Ltilde) real(Ltilde)];
    [aa,bb] = eigsort(Lbar+Lbar.');
    lambda = bb(1,1);
    phibar = [real(phi);imag(phi)];

    uth = -gt'*F + phibar.'*(Lbar+Lbar.'-lambda*eye(2*N))*[eye(N);-1i*eye(N)];
    e4 = real(-e3 + phibar.'*Lbar.'*phibar + real(gt'*kron(eye(K+M),hdt*hdt.')*vec(W))-lambda*N);

    cvx_begin quiet
    cvx_solver SeDuMi
    variable phi(N,1) complex
    minimize real(phi'*D*phi)-real(g'*phi) + 0.5/rho*square_pos(norm(phi-varphi+mu*rho,2))
    subject to
    real(uth*phi) <= e4;
    abs(phi) <= 1;
    cvx_end

    varphi = exp(1i*angle(phi+mu*rho));
    mu = mu + (phi-varphi)/rho;

    %%%% variable follow phi
    Hk = Hu + Hru*diag(phi)*G;
    Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
    s_phi_w = kron(eye(M+K),Ht)*vec(W);
    Dt = kron(eye(K+M),Hsi')*inv(Mw)*s_phi_w*s_phi_w'*inv(Mw)*kron(eye(K+M),Hsi);
    gt = inv(Mw)*s_phi_w;
    ct = real(trace(sigmar2*inv(Mw)*s_phi_w*s_phi_w'*inv(Mw)));
    e3 = real(0.5*(gammat/nL/sigmat2+(vec(W))'*Dt*vec(W)+ct));

    %%%%% update transmit beamformer w
    a = zeros(M*(K+M),1);
    B = zeros(K*(K+M),M*(K+M));
    for k = 1:1:K
        a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
        for j = 1:1:K+M
            Tj = zeros(M,M*(K+M));
            Tj(:,(j-1)*M+1:j*M) = eye(M);
            B((k-1)*(K+M)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
        end
    end

    dt = kron(eye(K+M),Hsi')*inv(Mw)*s_phi_w;
    kHt = gt'*kron(eye(K+M),Ht);
    cvx_begin quiet
    cvx_solver SeDuMi
    variable w(M*(K+M),1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w)
    subject to
    -square_pos(norm(dt'*w)) + 2*real(kHt*w) - ct >= gammat/nL/sigmat2;
    norm(w,2) <= sqrt(P);
    cvx_end
    W = reshape(w,M,K+M);

    %%%% update auxiliary variables
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    Vres(iter) = norm(phi-varphi,2)^2;
    Vgammat(iter) = real(nL*sigmat2*s_phi_w'*inv(Mw)*s_phi_w);
    if iter > 10
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter-1);
    end
    iter = iter + 1;
    rho = rho*0.8;
end
Vsr(iter:end) = [];
Vres(iter:end) = [];
Vgammat(iter:end) = [];
gammat_r = real(nL*sigmat2*s_phi_w'*inv(Mw)*s_phi_w);





