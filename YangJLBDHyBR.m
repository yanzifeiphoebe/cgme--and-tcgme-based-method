function [x_out, output] = YangJLBDHyBR(A, L, b, k, options)
%
% [x_out,flag, relres, iter, output] = YangIterJLBD(A, b, P, options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call by means of the following way
% options== HyBRset('InSolv', 'tikhonov', 'RegPar', 'gcv', 'x_true', x_true(:),'Iter', maxit);
%
% [x_out, flag, relres, iter, output] = YangJLBDHyBR(A, L, b, [], options)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs:
%      x_out : computed solution
%     output : structure with the following fields:
%      iterations - stopping iteration (options.Iter | GCV-determined)
%         GCVstop - GCV curve used to find stopping iteration
%            Enrm - relative error norms (requires x_true)
%            Rnrm - relative residual norms
%            Xnrm - relative solution norms
%             U,V - Lanczos basis vectors
%               B - bidiagonal matrix from LBD
%            flag - a flag that describes the output/stopping condition:
%                       1 - flat GCV curve
%                       2 - min of GCV curve (within window of 4 its)
%                       3 - performed max number of iterations

%% Initialization

% A = sparse(A);
[~,n] = size(L);
defaultopt = struct('InSolv','tikhonov','RegPar','wgcv','Omega',...
    'adapt', 'Iter', [] , 'Reorth', 'off', 'x_true', 'off', 'BegReg', 2,...
    'Vx' , [], 'FlatTol',  10^-6);

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x_out = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 2
    error('HyBR: Not Enough Inputs')
elseif nargin < 3
    P = []; options = [];
elseif nargin < 4
    options = [];
end
if isempty(options)
    options = defaultopt;
end

% Get options:
% [m,n] = size(A);
% defaultopt.Iter = min([m, n, 100]);
options = HyBRset(defaultopt, options);

solver = HyBRget(options,'InSolv',[],'fast');
regpar = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');
maxiter = HyBRget(options,'Iter',[],'fast');
x_true = HyBRget(options,'x_true',[],'fast');
regstart = HyBRget(options,'BegReg',[],'fast');
degflat = HyBRget(options,'FlatTol',[],'fast');

adaptWGCV = strcmp(regpar, {'wgcv'}) && strcmp(omega, {'adapt'});
notrue = strcmp(x_true,{'off'});

%--------------------------------------------
%  The following is needed for RestoreTools:
%
% if isa(A, 'psfMatrix')
%     bSize = size(b);
%     b = b(:);
%     A.imsize = bSize;
%     if ~notrue
%         xTrueSize = size(x_true);
%         x_true = x_true(:);
%     end
% end
%
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if ~notrue
    nrmtrue = norm(L*x_true);
    nrmtrue2 = norm(x_true);
end

% Set-up output parameters:
outputparams = nargout>1;
if outputparams
    output.iterations = maxiter;
    output.GCVstop = [];
%     output.err = ones(maxiter,1);
%     output.Rnrm = ones(maxiter,1);
%     output.Xnrm = (maxiter,1);
%     output.U = [];
%     output.Uhat = [];
%     output.V = [];
%     output.B = [];
%     output.Bhat = [];
%     output.flag = 3;
end

% Test for a preconditioner:
% if isempty(P)
%     beta = norm(b); 
%     U(:,1) = b / beta;
%     handle = @YangJLBD;
% else
%     U = P\b;
%     beta = norm(U); U = U / beta;
%     handle = @PLBD;
% end

%% Main Code Begins Here
% B = []; Bbar = []; Bhat = []; Uhat = [];  V = [];   x_out = []; Omega= [];
insolve = 'none';

h = waitbar(0, 'Beginning iterations: please wait ...');


[BB, BBbar, UU, ~, VV, bbeta]=JointBid(A,L,b,k,1);

for i = 1:maxiter %Iteration (i=1) is just an initialization
    
%     [U, Uhat, B, Bbar, Bhat, V] = feval(handle, A, L, U, Uhat, B, Bbar, Bhat, V, 1);

%    [U, Uhat, B, Bbar, Bhat, V, flag, relres, iter] = YangJLBD(A, L, U, Uhat, B, Bbar, Bhat, V, 0);

    U = UU(:,1:i+1); V = VV(:,1:i);
    B = BB(1:i+1,1:i); Bbar = BBbar(1:i,1:i);
    vector = (bbeta*eye(size(U,2),1));
    
    
    
    if i >= 2 %Begin Lanczos iterations
        if i >= regstart %Begin to regularize projected problem
            insolve = solver;
        end
        switch insolve
            case 'tikhonov'
                [Ub,ss,Xb]=cgsvd(B,Bbar);
                
%                 [Ub,~,Xb,Cb,Sb] = gsvd(B,Bbar);
%                 for j=1:i
%                     ss(j,1) = Cb(j,j);
%                     ss(j,2) = Sb(j,j);
%                 end
%                 Xb = inv(Xb');
%                 f = tgsvd(Ub,ss,Xb,vector,i);
%                
                
                
                if adaptWGCV %Use the adaptive, weighted GCV method
                   o(i-1)= min(1, findomega(Ub'*vector, ss));
                   omega = mean(o);
                   lambda = yangwgcv(Ub,ss,vector,'Tikh',omega);
                else
                    lambda = gcv(Ub,ss,vector);
                end
                
                f = tikhonov(Ub,ss,Xb,vector,lambda);
                d = V*f;
                x = lsqr(@(z,tflag)afun(z,A,L,tflag),d,1e-6,n);
                x_out(:,i) = x;
       
                
            case 'none'
                
                f = B\vector;
                d = V*f;
                x = lsqr(@(z,tflag)afun(z,A,L,tflag),d,1e-6,n);
                x_out(:,i) = x;
                output.Rnrm(i-1,1) = norm(B*f-vector);
                output.Xnrm(i-1,1) = norm(Bbar*f);
                
%                 [reg_c,rho_c,eta_c] = l_corner(output.Rnrm,output.Xnrm);
%                 plot_lc(output.Rnrm,output.Xnrm,'o',2);
%                 
%                 axis([1e-3,1,1e-3,1]);
%                 ax = axis;hold on
%                 loglog([min(output.Rnrm)/100,rho_c],[eta_c,eta_c],':r',...
%                 [rho_c,rho_c],[min(output.Xnrm)/100,eta_c],':r')
%                 title(['L-curve, ',' corner at ',num2str(reg_c)]);
%                 axis(ax);hold off

                
                


            otherwise
                error('HyBR error: No inner solver!')
        end
        

%         x = lsqr(M, V*f, 1e-8, n);
%         x_out(:,i) = x;
        
        if outputparams
            if ~notrue
                output.err(i-1,1) = norm(L*(x-x_true))/nrmtrue; 
                output.err2(i-1,1) = norm(x-x_true)/nrmtrue2;
            end
%             output.Rnrm(i-1,1) = norm(b - A*x);
%             output.Xnrm(i-1,1) = norm(L*x);
%             p = size(output.Rnrm,1);
%             [reg_c,rho_c,eta_c] = l_corner(output.Rnrm,output.Xnrm,(1:p)');
%             plot_lc(output.Rnrm,output.Xnrm,'o',1,(1:p)');
%             ax = axis;hold on
%             loglog([min(output.Rnrm)/100,rho_c],[eta_c,eta_c],':r',...
%               [rho_c,rho_c],[min(output.Xnrm)/100,eta_c],':r')
%             title(['L-curve, ','tsvd',' corner at ',num2str(reg_c)]);
%             axis(ax)
%             output.k(i-1) = reg_c;
%             hold off
            
        end
        
    end
    waitbar(i/(maxiter+1), h)
end
close(h)
% [reg_c,rho_c,eta_c] = l_corner(output.Rnrm,output.Xnrm);
% plot_lc(output.Rnrm,output.Xnrm,'o',2);
% ax = axis;hold on
% loglog([min(output.Rnrm)/100,rho_c],[eta_c,eta_c],':r',...
% [rho_c,rho_c],[min(output.Xnrm)/100,eta_c],':r')
% title(['L-curve, ',' corner at ',num2str(reg_c)]);
% axis(ax);hold off
% axis([1e-3,1,1,1e3])



                



% if outputparams
%     output.U = U;
%     output.V = V;
%     output.B = B;
%     output.Bbar = Bbar;
% end



%% -----------------------SUBFUNCTION---------------------------------------
function omega = findomega(bhat, s)
%
%   omega = findomega(bhat, s)
%
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the
%  derivative equal to zero and then solve for omega.
%
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%
alpha = s(1,1)/s(1,2);

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = size(s,1);
t0 = sum(abs(bhat(n+1:m)).^2);

s2 = abs(s(:,1)) .^ 2;
alpha2 = alpha^2*(abs(s(:,2)).^2);

tt = 1 ./ (s2 + alpha2);

t1 = sum(s2 .* tt);
t2 = abs(bhat(1:n).*alpha.*s(:,1).*s(:,2)) .^2;
t3 = sum(t2 .* abs((tt.^3)));

t4 = sum((s(:,1).*tt) .^2);
t5 = sum((abs(alpha2.*bhat(1:n).*tt)).^2);

v1 = abs(bhat(1:n).*s(:,1).*s(:,2)).^2;
v2 = sum(v1.* abs((tt.^3)));

%
% Now compute omega.
%
omega = (m*alpha^2*v2)/(t1*t3 + t4*(t5 + t0));
function omega = findomega222(bhat, s)
%
%   omega = findomega(bhat, s)
%
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the 
%  derivative equal to zero and then solve for omega.
%  
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%
alpha = s(end);

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);
t0 = sum(abs(bhat(n+1:m)).^2);

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

tt = 1 ./ (s2 + alpha2);

t1 = sum(s2 .* tt);
t2 = abs(bhat(1:n).*alpha.*s) .^2;
t3 = sum(t2 .* abs((tt.^3)));

t4 = sum((s.*tt) .^2);
t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);

v1 = abs(bhat(1:n).*s).^2;
v2 = sum(v1.* abs((tt.^3)));

%
% Now compute omega.
%
omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));


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



