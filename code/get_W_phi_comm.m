% Optimize the joint beamforming and reflection design problem without considering sensing.
% Generate the comparison "Comm only"
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         phi: the initial phi; W: the initial W
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; gammat_r: the achieved radar SNR;
%          crb: the achieved CRB
function [W,phi,Vsr,gammat_r,crb] = get_W_phi_comm(Prms,Channel,phi,W)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
P = Prms.P; hdt = Channel.hdt; hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
hdt_der = Channel.hdt_der; hrt_der = Channel.hrt_der;

%%%% variable follow phi
Hk = Hu + Hru*diag(phi)*G;
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
    for j = 1:1:K
        temp = diag(W(:,j)'*G')*Hru(k,:)';
        D = D + abs(c(k))^2*(temp*temp');
        g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
    end
end

%%%% start iteration
iter = 1;
res = 1;
Vsr = zeros(1,Nmax);
while iter <= Nmax && res >= res_th
    %%%% update reflecting coefficients phi
    D = zeros(N,N);
    g = zeros(N,1);
    for k = 1:1:K
        g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
        for j = 1:1:K
            temp = diag(W(:,j)'*G')*Hru(k,:)';
            D = D + abs(c(k))^2*(temp*temp');
            g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
        end
    end
    phi = get_initial_phi(-D/norm(g),g/norm(g));

    %%%% variable follow phi
    Hk = Hu + Hru*diag(phi)*G;
    %%%%% update transmit beamformer w
    a = zeros(M*(K),1);
    B = zeros(K*(K),M*(K));
    for k = 1:1:K
        a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
        for j = 1:1:K
            Tj = zeros(M,M*(K));
            Tj(:,(j-1)*M+1:j*M) = eye(M);
            B((k-1)*(K)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
        end
    end
    cvx_begin quiet
    cvx_solver SeDuMi
    variable w(M*(K),1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w)
    subject to
    norm(w,2) <= sqrt(P);
    cvx_end
    W = reshape(w,M,K);
    %%%% update auxiliary variables
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    if iter > 10
        res = (Vsr(iter)-Vsr(iter-1))/Vsr(iter-1);
    end
    iter = iter + 1;
end
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
u = kron(eye(K),Ht)*vec(W)/((vec(W))'*(kron(eye(K),Ht'*Ht))*vec(W));
Vsr(iter:end) = [];
gammat_r = real(nL*sigmat2*abs(u'*kron(eye(K),Ht)*w)^2/(sigmar2*u'*u));
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;
Rw = W*W';
F1 = 2*nL/sigmar2*real([trace(Ht1*Rw*Ht1') trace(Ht2*Rw*Ht1');trace(Ht2*Rw*Ht1') trace(Ht2*Rw*Ht2')]);
F2 = 2*nL/sigmar2*[real( trace(Ht*Rw*Ht1')*[1 1j] );real( trace(Ht*Rw*Ht2')*[1 1j] )];
F4 = 2*nL/sigmar2*real(trace(Ht*Rw*Ht'))*eye(2);
J = F1-F2*inv(F4)*F2.';
crb = trace(inv(J));
% figure;grid on
% plot(Vsr);title('sum-rate');legend('Comm only')





