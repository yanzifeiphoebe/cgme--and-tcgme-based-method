function [B, Bbar, U, Uhat, V, bbeta]=JointBid(A,L,b,k,reorth)

[m,n] = size(A); p=size(L,1);
beta=norm(b);
bbeta=beta;
u=b/beta;
U(:,1)=u;

utilde=[u; zeros(p,1)];
x = lsqr(@(z,tflag)afun(z,A,L,tflag),utilde,1e-6,n);
ss = A*x; tt = L*x;
v = [ss;tt];
alpha=norm(v);
v = v/alpha;
B(1,1) = alpha;
V(:,1) = v;

uhat = v(m+1:m+p);
alphahat = norm(uhat);
uhat = uhat/alphahat;
Bbar(1,1) = alphahat;
Uhat(:,1) = uhat;

if (reorth == 0)
    u = v(1:m) - alpha * u;
elseif (reorth == 1)
    u = v(1:m) - alpha * u;
    u = u - U * (U' * u);
elseif (reorth == 2)
    u = v(1:m) - alpha * u; 
    u = u - U * (U' * u);
    u = u - U * (U' * u);
end
beta = norm(u);
u = u/beta;
B(2,1) = beta;
U(:,2) = u;

for i = 2:k
    utilde = [U(:,i); zeros(p,1)];
    x = lsqr(@(z,tflag)afun(z,A,L,tflag), utilde,1e-6,n);
    ss = A*x; tt = L*x;
    Qu = [ss;tt];
    if (reorth == 0)
        v = Qu - B(i, i-1)*V(:,i-1);
    elseif(reorth == 1)
        v = Qu - B(i, i-1)*V(:,i-1);
        for j=1:i-1, v = v - (V(:,j)'*v)*V(:,j); end
%         v = v - V * V' * v;
    elseif (reorth == 2)
        v = Qu - B(i, i-1)*V(:,i-1);
        for j=1:i-1, v = v - (V(:,j)'*v)*V(:,j); end
        for j=1:i-1, v = v - (V(:,j)'*v)*V(:,j); end
%         v = v - V * V' * v;
%         v = v - V * V' * v;
    end
    alpha = norm(v);
    v = v/alpha;
    B(i,i) = alpha;
    V(:,i) = v;
    
    betahat=(alpha*B(i,i-1))/alphahat;
    if(mod(i,2)==0)
        Bbar(i-1,i) = -betahat;
    else
        Bbar(i-1,i) = betahat;
    end
    
    
    if(mod(i,2)==0)
        vv = -v(m+1:m+p);
    else
        vv = v(m+1:m+p);
    end
    
    
    if (reorth == 0)
        uhat = vv - betahat * Uhat(:,i-1);
    elseif (reorth == 1)
        uhat = vv - betahat * Uhat(:,i-1);
        for j=1:i-1, uhat = uhat - (Uhat(:,j)'*uhat)*Uhat(:,j); end
%         uhat = uhat - Uhat * Uhat' * uhat;
    elseif (reorth == 2)
        uhat = vv - betahat * Uhat(:,i-1);
        for j=1:i-1, uhat = uhat - (Uhat(:,j)'*uhat)*Uhat(:,j); end
        for j=1:i-1, uhat = uhat - (Uhat(:,j)'*uhat)*Uhat(:,j); end
%         uhat = uhat - Uhat * Uhat' * uhat;
%         uhat = uhat - Uhat * Uhat' * uhat;
    end
    alphahat = norm(uhat);
    if(mod(i,2)==0)
        Bbar(i,i) = -alphahat;
    else
        Bbar(i,i) = alphahat;
    end
    uhat = uhat/alphahat;
    Uhat(:,i) = uhat;
    
    if (reorth == 0)
        u = v(1:m) - alpha * u;
    elseif (reorth == 1)
        u = v(1:m) - alpha * u;
        for j=1:i-1, u = u - (U(:,j)'*u)*U(:,j); end
%         u = u - U * U' * u;
    elseif (reorth == 2)
        u = v(1:m) - alpha * u;
        for j=1:i-1, u = u - (U(:,j)'*u)*U(:,j); end
        for j=1:i-1, u = u - (U(:,j)'*u)*U(:,j); end
%         u = u - U * U' * u;
%         u = u - U * U' * u;
    end
    beta = norm(u);
    u = u/beta;
    B(i+1,i) = beta;
    U(:,i+1) = u;
end
end

function y =afun(z,A,B,transp_flag)
if strcmp(transp_flag,'transp')   % y = (A(I_n-BB^T))' * z;
    m = size(A,1);
    p = size(B,1);
    s = A'*z(1:m);
    t = B'*z(m+1:m+p);
    y = s + t;
elseif strcmp(transp_flag,'notransp') % y = (A(I_n-BB^T)) * z;
    s = A*z;
    t= B*z;
    y = [s;t];
    
end
end


