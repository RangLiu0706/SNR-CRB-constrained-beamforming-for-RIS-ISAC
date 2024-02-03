% Optimize the CRB-constrained joint beamforming and reflection design problem (19).
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         phi: the initial phi; W: the initial W; rho: the initial rho
% Outputs: W: transmit beamforming; phi: RIS reflection coefficients
%          Vsr: the achieved sum-rate; crb: the achieved CRB
function [W,phi,Vsr,crb] = get_W_phi_CRB(Prms,Channel,phi,W,rho)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
Nmax = Prms.Nmax; res_th = Prms.res_th; L = Prms.L;
P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
hdt_der = Channel.hdt_der; hrt_der = Channel.hrt_der;
CRB = Prms.Crb;

Hk = Hu + Hru*diag(phi)*G;
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
    + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;
w = vec(W);
Rw = W*W';
F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1') trace(Ht2*Rw*Ht1');trace(Ht2*Rw*Ht1') trace(Ht2*Rw*Ht2')]);
F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1')*[1 1j] );real( trace(Ht*Rw*Ht2')*[1 1j] )];
F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'))*eye(2);
J = F1-F2*inv(F4)*F2.';
Crbo = trace(inv(J));

ft = zeros(6,1);
ft(1) = 2/sigmar2*trace(Ht1*Rw*Ht1'); ft(2) = 2/sigmar2*trace(Ht2*Rw*Ht1');
ft(3) = 2/sigmar2*trace(Ht*Rw*Ht1'); ft(4) = 2/sigmar2*trace(Ht2*Rw*Ht2');
ft(5) = 2/sigmar2*trace(Ht*Rw*Ht2'); ft(6) = 2/sigmar2*trace(Ht*Rw*Ht');
f = ft;
r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end
zeta = zeros(6,1);
mu = zeros(N,1);
lambda = zeros(N^2,1);
omega = zeros(N^2,1);
varphi = phi;
v = kron(phi,phi);
nu = v;

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
rho2 = N/4/abs(norm(B*w,2)^2-real(a'*w))*rho;
rho3 = N*rho2;
rho4 = rho3;
Vres = zeros(Nmax,4);
Vsr = zeros(Nmax,1);
VCRB = zeros(Nmax,3);
while iter <= Nmax && res >= res_th

    %%%% update J and f
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
    %%%%% update phi
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
    %%%% update variables
    qk = zeros(M,K+M); q1k = zeros(M,K+M);
    Uk = zeros(M,N,K+M);U1k = zeros(M,N,K+M);U2k = zeros(M,N,K+M);
    Zk = zeros(M,N^2,K+M);Z2k = zeros(M,N^2,K+M);
    for k = 1:1:K+M
        qk(:,k) = hdt*hdt.'*W(:,k); q1k(:,k) = (hdt*hdt_der.'+hdt_der*hdt.')*W(:,k);
        Uk(:,:,k) = W(:,k).'*hdt*G.'*diag(hrt)+kron(W(:,k).'*G.'*diag(hrt),hdt);
        U1k(:,:,k) = W(:,k).'*hdt_der*G.'*diag(hrt)+kron(W(:,k).'*G.'*diag(hrt),hdt_der);
        U2k(:,:,k) = W(:,k).'*hdt*G.'*diag(hrt_der)+kron(W(:,k).'*G.'*diag(hrt_der),hdt);
        Zk(:,:,k) = kron(W(:,k).'*G.'*diag(hrt),G.'*diag(hrt));
        Z2k(:,:,k) = kron(W(:,k).'*G.'*diag(hrt_der),G.'*diag(hrt)) + kron(W(:,k).'*G.'*diag(hrt),G.'*diag(hrt_der));
    end

    q1 = trace(q1k*q1k');q3 = 0;q6 = trace(qk*qk');
    u1 = zeros(N,1);u3 = zeros(N,1);u5 = zeros(N,1);u6 = zeros(N,1);
    U1 = zeros(N,N);U2 = zeros(N,N);U3 = zeros(N,N);U4 = zeros(N,N);U5 = zeros(N,N);U6 = zeros(N,N);
    u2b = zeros(N,1);u3b = zeros(N,1);
    z2bM = zeros(N,N);z3bM = zeros(N,N);
    C2 = zeros(N,N^2);C3 = zeros(N,N^2);C5 = zeros(N,N^2);
    C4b = zeros(N^2,N);C5b = zeros(N^2,N);C6b = zeros(N^2,N);
    Z4 = zeros(N^2,N^2);Z5 = zeros(N^2,N^2);Z6 = zeros(N^2,N^2);
    z5M = zeros(N,N); z6M = zeros(N,N);
    for k = 1:1:K+M
        u1 = u1 + U1k(:,:,k)'*q1k(:,k);
        U1 = U1 + U1k(:,:,k)'*U1k(:,:,k);
        u2b = u2b + U2k(:,:,k)'*q1k(:,k);
        U2 = U2 + U1k(:,:,k)'*U2k(:,:,k);
        z2bM = z2bM + diag(hrt')*conj(G)*q1k(:,k)*W(:,k)'*G'*diag(hrt_der') + diag(hrt_der')*conj(G)*q1k(:,k)*W(:,k)'*G'*diag(hrt');
        C2 = C2 + U1k(:,:,k)'*Z2k(:,:,k);
        q3 = q3 + q1k(:,k)'*qk(:,k);
        u3 = u3 + U1k(:,:,k)'*qk(:,k);
        u3b = u3b + Uk(:,:,k)'*q1k(:,k);
        U3 = U3 + U1k(:,:,k)'*Uk(:,:,k);
        z3bM = z3bM + diag(hrt')*conj(G)*q1k(:,k)*W(:,k)'*G'*diag(hrt');
        C3 = C3 + U1k(:,:,k)'*Zk(:,:,k);
        U4 = U4 + U2k(:,:,k)'*U2k(:,:,k);
        C4b = C4b + Z2k(:,:,k)'*U2k(:,:,k);
        Z4 = Z4 + Z2k(:,:,k)'*Z2k(:,:,k);
        u5 = u5 + U2k(:,:,k)'*qk(:,k);
        U5 = U5 + U2k(:,:,k)'*Uk(:,:,k);
        z5M = z5M + diag(hrt')*conj(G)*qk(:,k)*W(:,k)'*G'*diag(hrt_der') + diag(hrt_der')*conj(G)*qk(:,k)*W(:,k)'*G'*diag(hrt');
        Z5 = Z5 + Z2k(:,:,k)'*Zk(:,:,k);
        C5b = C5b + Z2k(:,:,k)'*Uk(:,:,k);
        C5 = C5 + U2k(:,:,k)'*Zk(:,:,k);
        u6 = u6 + Uk(:,:,k)'*qk(:,k);
        U6 = U6 + Uk(:,:,k)'*Uk(:,:,k);
        z6M = z6M + diag(hrt')*conj(G)*qk(:,k)*W(:,k)'*G'*diag(hrt');
        Z6 = Z6 + Zk(:,:,k)'*Zk(:,:,k);
        C6b = C6b + Zk(:,:,k)'*Uk(:,:,k);
    end
    a_phi = zeros(6,1);  b_phi = zeros(N,6);
    a_phi(1) = 2/sigmar2*( q1 + u1'*varphi )-f(1)+rho1*zeta(1);
    b_phi(:,1) = 2/sigmar2*(u1 + U1*varphi );
    a_phi(2) = 2/sigmar2*(u2b'*varphi+(vec(z2bM))'*nu )-f(2)+rho1*zeta(2);
    b_phi(:,2) = 2/sigmar2*(U2*varphi+C2*nu );
    a_phi(3) = 2/sigmar2*(q3 + u3b'*varphi + (vec(z3bM))'*nu  )-f(3)+rho1*zeta(3);
    b_phi(:,3) = 2/sigmar2*(u3 + U3*varphi + C3*nu );
    a_phi(4) = 2/sigmar2*(v'*C4b*varphi + v'*Z4*nu  )-f(4)+rho1*zeta(4);
    b_phi(:,4) = 2/sigmar2*(U4*varphi + C4b'*nu );
    a_phi(5) = 2/sigmar2*(v'*vec(z5M) + v'*Z5*nu + v'*C5b*varphi  )-f(5)+rho1*zeta(5);
    b_phi(:,5) = 2/sigmar2*(u5 + U5*varphi + C5*nu );
    a_phi(6) = 2/sigmar2*(q6 + u6'*varphi + v'*vec(z6M) + (vec(z6M))'*nu + v'*C6b*varphi + v'*Z6*nu )-f(6)+rho1*zeta(6);
    b_phi(:,6) = 2/sigmar2*(u6 + U6*varphi + C6b'*nu );

    Dphi = D + b_phi*b_phi'/(2*rho1);
    gphi = -varphi/rho2 + mu - kron(eye(N),varphi')*(v/rho3+lambda);
    for i = 1:1:6
        gphi = gphi + a_phi(i)'*b_phi(:,i)/rho1;
    end
    gphi = gphi-g;
    %%% update phi
    phi = cal_phi_manopt(Dphi/norm(gphi),gphi/norm(gphi));
    % gtphi = gphi + 2*(Dphi-norm(Dphi,'fro')^2*eye(N))*phi;
    % phi = exp(1i*(pi+angle(gtphi)));

    %%% update variables for varphi
    a_varphi = zeros(6,1);   b_varphi = zeros(N,6);
    a_varphi(1) = 2/sigmar2*( q1 + phi'*u1  )- f(1) + rho1*zeta(1);
    b_varphi(:,1) = 2/sigmar2*( u1 + U1'*phi );
    a_varphi(2) = 2/sigmar2*( (vec(z2bM))'*nu + phi'*C2*nu ) - f(2) + rho1*zeta(2);
    b_varphi(:,2) = 2/sigmar2*( u2b + U2'*phi );
    a_varphi(3) = 2/sigmar2*( q3 + phi'*u3 + (vec(z3bM))'*nu + phi'*C3*nu  )- f(3) + rho1*zeta(3);
    b_varphi(:,3) = 2/sigmar2*( u3b + U3'*phi );
    a_varphi(4) = 2/sigmar2*( phi'*C4b'*nu + v'*Z4*nu  )- f(4) + rho1*zeta(4);
    b_varphi(:,4) = 2/sigmar2*( U4'*phi + C4b'*v );
    a_varphi(5) = 2/sigmar2*( phi'*u5 + v'*vec(z5M) + v'*Z5*nu + phi'*C5*nu ) - f(5) + rho1*zeta(5);
    b_varphi(:,5) = 2/sigmar2*( U5'*phi + C5b'*v );
    a_varphi(6) = 2/sigmar2*( q6 + phi'*u6 + v'*vec(z6M) + (vec(z6M))'*nu + phi'*C6b'*nu + v'*Z6*nu  )- f(6) + rho1*zeta(6);
    b_varphi(:,6) = 2/sigmar2*( u6 + U6'*phi + C6b'*v );
    Dvarphi = b_varphi*b_varphi'/(2*rho1);
    gvarphi = -phi/rho2-mu-kron(phi',eye(N))*(v/rho3+lambda);
    for i = 1:1:6
        gvarphi = gvarphi + a_varphi(i)*b_varphi(:,i)/rho1;
    end
    %%%% update varphi
    varphi = cal_phi_manopt(Dvarphi/norm(gvarphi),gvarphi/norm(gvarphi));
    % gtvarphi = gvarphi + 2*(Dvarphi-norm(Dvarphi,'fro')^2*eye(N))*varphi;
    % varphi = exp(1i*(pi+angle(gtvarphi)));

    a_v = zeros(6,1);   b_v = zeros(N^2,6);
    a_v(1) = 2/sigmar2*( q1 + phi'*u1 + u1'*varphi + phi'*U1*varphi  )- f(1) + rho1*zeta(1);
    a_v(2) = 2/sigmar2*( u2b'*varphi + phi'*U2*varphi + (vec(z2bM))'*nu + phi'*C2*nu  )- f(2) + rho1*zeta(2);
    a_v(3) = 2/sigmar2*( q3 + phi'*u3 + u3b'*varphi + phi'*U3*varphi + (vec(z3bM))'*nu + phi'*C3*nu  )- f(3) + rho1*zeta(3);
    a_v(4) = 2/sigmar2*( phi'*U4*varphi + phi'*C4b'*nu  )- f(4) + rho1*zeta(4);
    b_v(:,4) = 2/sigmar2*( C4b*varphi + Z4*nu );
    a_v(5) = 2/sigmar2*( phi'*u5 + phi'*U5*varphi + phi'*C5*nu  )- f(5) + rho1*zeta(5);
    b_v(:,5) = 2/sigmar2*( vec(z5M) + Z5*nu + C5b*varphi );
    a_v(6) = 2/sigmar2*( q6 + phi'*u6 + u6'*varphi + phi'*U6*varphi + (vec(z6M))'*nu + phi'*C6b'*nu )- f(6) + rho1*zeta(6);
    b_v(:,6) = 2/sigmar2*( vec(z6M) + C6b*varphi + Z6*nu );
    Dv = b_v*b_v'/(2*rho1);
    gv = lambda - kron(phi,varphi)/rho3 - omega -nu/rho4;
    for i = 1:1:6
        gv = gv + a_v(i)'*b_v(:,i)/rho1;
    end
    %%%% update v
    v = cal_phi_manopt(Dv/norm(gv),gv/norm(gv));
    % gtv = gv + 2*(Dv-norm(Dv,'fro')^2*eye(N^2))*v;
    % v = exp(1i*(pi+angle(gtv)));

    a_nu = zeros(6,1);
    b_nu = zeros(N^2,6);
    a_nu(1) = 2/sigmar2*( q1 + phi'*u1 + u1'*varphi + phi'*U1*varphi  )- f(1) + rho1*zeta(1);
    a_nu(2) = 2/sigmar2*( u2b'*varphi + phi'*U2*varphi  )- f(2) - rho1*zeta(2);
    b_nu(:,2) = 2/sigmar2*( vec(z2bM) + C2'*phi );
    a_nu(3) = 2/sigmar2*( q3 + phi'*u3 + u3b'*varphi + phi'*U3*varphi  )- f(3) + rho1*zeta(3);
    b_nu(:,3) = 2/sigmar2*( vec(z3bM) + C3'*phi );
    a_nu(4) = 2/sigmar2*( phi'*U4*varphi + v'*C4b*varphi  )- f(4) + rho1*zeta(4);
    b_nu(:,4) = 2/sigmar2*( C4b*phi + Z4'*v );
    a_nu(5) = 2/sigmar2*( phi'*u5 + phi'*U5*varphi + v'*vec(z5M) + v'*C5b*varphi  )- f(5) + rho1*zeta(5);
    b_nu(:,5) = 2/sigmar2*( Z5'*v + C5'*phi );
    a_nu(6) = 2/sigmar2*( q6 + phi'*u6 + u6'*varphi + phi'*U6*varphi + v'*vec(z6M) + v'*C6b*varphi  )- f(6) + rho1*zeta(6);
    b_nu(:,6) = 2/sigmar2*( vec(z6M) + C6b*phi + Z6'*v );
    Dnu = b_nu*b_nu'/(2*rho1);
    gnu = -v/rho4 + omega;
    for i = 1:1:6
        gnu = gnu + a_nu(i)*b_nu(:,i)/rho1;
    end
    %%% update nu
    nu = cal_phi_manopt(Dnu/norm(gnu),gnu/norm(gnu));
    % gtnu = gnu + 2*(Dnu-norm(Dnu,'fro')^2*eye(N^2))*nu;
    % nu = exp(1i*(pi+angle(gtnu)));

    %%% update dual variables
    mu = mu + (phi-varphi)/rho2;
    lambda = lambda + (v-kron(phi,varphi))/rho3;
    omega = omega + (nu-v)/rho4;

    %%%% variable follow phi
    Hk = Hu + Hru*diag(phi)*G;
    Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
    Ht1 = hdt*hdt_der.' + hdt_der*hdt.' + G.'*diag(hrt)*phi*hdt_der.' + hdt_der*phi.'*diag(hrt)*G;
    Ht2 = G.'*diag(hrt_der)*phi*hdt.' + hdt*phi.'*diag(hrt_der)*G + G.'*diag(hrt)*(phi*phi.')*diag(hrt_der)*G...
        + G.'*diag(hrt_der)*(phi*phi.')*diag(hrt)*G;
    F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1') trace(Ht2*Rw*Ht1');trace(Ht2*Rw*Ht1') trace(Ht2*Rw*Ht2')]);
    F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1')*[1 1j] );real( trace(Ht*Rw*Ht2')*[1 1j] )];
    F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'))*eye(2);
    J = F1-F2*inv(F4)*F2.';
    VCRB(iter,2) = trace(inv(J));
    %%%% update w
    Ai = zeros(M,M,6);
    Ai(:,:,1) = 2/sigmar2*Ht1'*Ht1; Ai(:,:,2) = 2/sigmar2*Ht1'*Ht2;
    Ai(:,:,3) = 2/sigmar2*Ht1'*Ht; Ai(:,:,4) = 2/sigmar2*Ht2'*Ht2;
    Ai(:,:,5) = 2/sigmar2*Ht2'*Ht;  Ai(:,:,6) = 2/sigmar2*Ht'*Ht;

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
    F1 = 2*L/sigmar2*real([trace(Ht1*Rw*Ht1') trace(Ht2*Rw*Ht1');trace(Ht2*Rw*Ht1') trace(Ht2*Rw*Ht2')]);
    F2 = 2*L/sigmar2*[real( trace(Ht*Rw*Ht1')*[1 1j] );real( trace(Ht*Rw*Ht2')*[1 1j] )];
    F4 = 2*L/sigmar2*real(trace(Ht*Rw*Ht'))*eye(2);
    J = F1-F2*inv(F4)*F2.';
    VCRB(iter,3) = trace(inv(J));

    %%% update auxillary variables
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    %%%% update dual varialbes
    zeta = zeta + (ft-f)/rho1;
    ft = zeros(6,1);
    ft(1) = 2/sigmar2*trace(Ht1*Rw*Ht1');  ft(2) = 2/sigmar2*trace(Ht2*Rw*Ht1');
    ft(3) = 2/sigmar2*trace(Ht*Rw*Ht1');   ft(4) = 2/sigmar2*trace(Ht2*Rw*Ht2');
    ft(5) = 2/sigmar2*trace(Ht*Rw*Ht2');   ft(6) = 2/sigmar2*trace(Ht*Rw*Ht');

    Vsr(iter) = sum(log2(1+r));
    Vres(iter,:) = [norm(f-ft,2)  norm(phi-varphi,2) norm(v-kron(phi,varphi),2)  norm(nu-v,2)];
    if iter > 1
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter);
    end
    rho1 = rho1*0.9;
    rho2 = rho2*0.9;
    rho3 = rho3*0.9;
    rho4 = rho4*0.9;
    if VCRB(iter,3) > 0.8*CRB
        rho1 = rho1*0.7;
        rho2 = rho2*0.7;
        rho3 = rho3*0.7;
        rho4 = rho4*0.7;
    elseif VCRB(iter,3) > 0.9*CRB
        rho1 = rho1*0.1;
        rho2 = rho2*0.1;
        rho3 = rho3*0.1;
        rho4 = rho4*0.1;
    end
    if VCRB(iter,3) > CRB
        res = 0;
        iter = iter - 1;
    end

    % [iter Vsr(iter) VCRB(iter,3) Vres(iter,:)]
    iter = iter + 1;
end
Vsr(iter:end) = [];
Vsr = Vsr.';
VCRB(iter:end,:) = [];
Vres(iter:end,:) = [];
crb = VCRB(end,end);
% figure;grid on
% subplot(3,1,1)
% plot(Vsr);title('sum-rate');legend('Proposed')
% subplot(3,1,2)
% semilogy(VCRB(:,3).');title('CRB');axis([1 iter-1 min(VCRB(:,3)) CRB]);legend('Proposed')
% subplot(3,1,3)
% semilogy(sum(Vres,2));title('res');legend('Proposed')








