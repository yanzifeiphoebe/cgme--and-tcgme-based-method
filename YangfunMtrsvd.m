
function [x_rm,x_s,z,f,flag,relres,iter] = YangfunMtrsvd(U,S,V,L,b,k)

%%%%%% Compute the MTRSVD and the relative errors %%%%%%%%%%%%%%%%%%
for i = 1:k
    U_i = U(:,1:i); S_i = S(1:i,1:i); V_i = V(:,1:i);
    p = U_i'*b;
    x_s(:,i) = V_i*(diag(diag(S_i))\p); 
    d = L*x_s(:,i);
    [z(:,i),flag(i),relres(i),iter(i)] = lsqr(@(z,tflag)afun(z,L,V_i,tflag),d,1e-6,10000);
    x_rm(:,i) = x_s(:,i)-z(:,i);
    f(:,i) = U_i*(U_i'*b);
end
end

%%%%%%%%%%%%%%Handler function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y =afun(z,A,B,transp_flag)
if strcmp(transp_flag,'transp')   % y = M' * z;
    s = A'*z;
    t = B'*s;
    t = B*t;
    y = s-t;
elseif strcmp(transp_flag,'notransp') % y = M * z;
    s = A*z;
    t = B'*z;
    t = B*t;
    t = A*t;
    y = s-t;
    
end
end