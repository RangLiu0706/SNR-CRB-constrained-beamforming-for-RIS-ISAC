% Optimize the SNR-constrained beamforming design given phi.
% Generate the comparison "BF only".
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         phi: the phi; W: the initial W
% Outputs: W: transmit beamforming;
%          Vsr: the achieved sum-rate; gammat_r: the achieved radar SNR
function [W,Vsr,gammat_r] = get_W_with_phi(Prms,Channel,phi,W)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
gammat = Prms.gammat; P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;

%%%% variable follow phi
Hk = Hu + Hru*diag(phi)*G;
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end
u = kron(eye(K+M),Ht)*vec(W)/((vec(W))'*(kron(eye(K+M),Ht'*Ht))*vec(W));
e3 = sqrt(real(gammat*sigmar2*u'*u)/sigmat2/nL);

%%%% start iteration
iter = 1;
res = 1;
Vsr = zeros(1,Nmax);
Vgammat = zeros(1,Nmax);
while iter <= Nmax && res >= res_th

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
    kHt = kron(eye(K+M),Ht);
    cvx_begin quiet
    cvx_solver SeDuMi
    variable w(M*(K+M),1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w)
    subject to
    real(u'*kHt*w) >= e3;
    norm(w,2) <= sqrt(P);
    cvx_end
    W = reshape(w,M,K+M);

    %%%% update receive filter u
    u = kron(eye(K+M),Ht)*vec(W)/((vec(W))'*(kron(eye(K+M),Ht'*Ht))*vec(W));
    e3 = sqrt(real(gammat*sigmar2*u'*u)/sigmat2/nL);

    %%%% update auxiliary variables
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    Vgammat(iter) = real(nL*sigmat2*abs(u'*kron(eye(K+M),Ht)*w)^2/(sigmar2*u'*u));
    if iter > 10
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter-1);
    end
    iter = iter + 1;
end
Vsr(iter:end) = [];
gammat_r = real(nL*sigmat2*abs(u'*kron(eye(K+M),Ht)*w)^2/(sigmar2*u'*u));
% figure;grid on
% subplot(2,1,1)
% plot(Vsr);title('sum-rate');legend('BF only')
% subplot(2,1,2)
% semilogy(Vgammat);title('SNR');axis([1 iter_o-1 gammat max(Vgammat)]);legend('BF only')




