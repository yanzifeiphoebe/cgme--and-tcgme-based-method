clc;      
clear all;
close all;


maxit = 200;
%%%%%%%%%%%%%% 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 10000;
% [K,b,x_true] = heat(n);
% noise  = randn(size(b));
% delta=0.01;
% sigma = delta*(norm(b)/norm(noise));          
% N = sigma*noise;
% b = b + N;
% L = get_l(n,1);
%                                            
% % %%%%%%%%%%%%%%%%%% 2-D  RestoreTools %%%%%%%
N = 256;
% K = blur(N,16,2);
% load mri
% x_true = double(D(1:N,1:N,15));
load GaussianBlur440
% load GaussianBlur422
% load GaussianBlur420
% load VariantMotionBlur_small
% load satellite
% load grain
% load AtmosphericBlur50
% load AtmosphericBlur30
% load AtmosphericBlur10
% load VariantGaussianBlur1
% load VariantGaussianBlur2
% load VariantGaussianBlur3
K = psfMatrix(PSF);
x_true=im2double(f_true);
% x_true=x_true(1:N,1:N);
 x_true=x_true(:);
% b=g(:);
btrue = K*x_true;
delta=0.05;
e =  randn(size(btrue(:)))/norm(randn(size(btrue(:))));
b = btrue(:) + norm(btrue(:))*delta*e;
%  %%%%%%%%%%%%%%%% 2-D IRtools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A, btrue, x_true, ProbInfo] = PRblurspeckle;
% [b, NoiseInfo] = PRnoise(btrue, 'gauss', delta);
%%%%%%%%%%%%%%%%%%%%%%%% L matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = speye(N);
L1=get_l(N,1);
L = [kron(I,L1); kron(L1,I)];

%  %%%%%%%%%%%%%%%%%%%%%%%% JBDQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = 0;
a_1=0;
b_1=0;
for cycl_no1=1:50
    tic;
    options = HyBRset('InSolv', 'none', 'x_true', x_true(:),'Iter', maxit);
    [x_out, output] = YangJLBDHyBR(K, L, b, maxit, options);
    t1=t1+toc;
    [a_jbdqr,b_jbdqr] = min(output.err);
    a_1 = a_1+a_jbdqr;
    b_1 = b_1+b_jbdqr;
end
% result1 = {'heat' delta t1/50};
% xlswrite('C:\Users\Yang Yanfei\Desktop\heat.xlsx',result1)
%%%%%% LSMR %%%%%%%%%%%
% tic;
% [X_lsmr,V_lsmr] = yanglsmr(A,b,maxit,1);
% [x_lsmr,err_lsmr,rho_lsmr,eta_lsmr] = YangfunTik(A, b,x_true,X_lsmr,V_lsmr,L,maxit);
% toc;
% result1 = {'heat' delta t1/50};
% xlswrite('C:\Users\Yang Yanfei\Desktop\heat.xlsx',result1)
%%%% CGME %%%%%%%%%%%
t2=0;
a_2=0;
b_2=0;
for cycl_no2=1:50
    tic;
    [X,V_cgme] = yangcgme(K,b,maxit,1);
    [x_cgme,err_cgme] = YangfunTik(K, b,x_true,X,V_cgme,L,maxit);
    t2 = t2+toc;
    [a_cgme,b_cgme] = min(err_cgme);
    a_2 = a_2+a_cgme;
    b_2 = b_2+b_cgme;
end
% result2 = {'heat' delta t2/50};
% xlswrite('C:\Users\Yang Yanfei\Desktop\heat.xlsx',result2)
%%%%%%%% TCGME %%%%%%%%%%%%%%%%%
t3 = 0;
a_3=0;
b_3=0;
for cycl_no3=1:50
    tic
    [X_t,V_tcgme] = yangtcgme(K,b,maxit,1);
    [x_tcgme,err_tcgme,rho_tcgme,eta_tcgme] = YangfunTik(K,b,x_true,X_t,V_tcgme,L,maxit);
    t3 = t3+toc;
    [a_tcgme,b_tcgme] = min(err_tcgme);
    a_3 = a_3+a_tcgme;
    b_3 = b_3+b_tcgme;
end

% %%%%%%% LSQR %%%%%%%%%%
% tic;
% [X_lsqr,V_lsqr] = Yanglsqr_stop(K,b,x_true, L, maxit+1,1);
% [x_lsqr,err_lsqr,rho_lsqr,eta_lsqr] = YangfunTik(K, b,x_true,X_lsqr,V_lsqr,L,maxit);
% toc;
% 
% 
% %%%%%%% the best error and the associated number of the iteration %%%%%%%%%%
% for i=1:maxit
%     err(i) = norm(L*(X(:,i)-x_true))/norm(L*x_true);
%     err_t(i) = norm(L*(X_t(:,i)-x_true))/norm(L*x_true);
% end
% [a_jbdqr,b_jbdqr] = min(output.err);
% % [a_lsqr,b_lsqr] = min(err_lsqr);
% [a_cgme,b_cgme] = min(err_cgme);
% [a_tcgme,b_tcgme] = min(err_tcgme);
% [a_c,b_c] = min(err);
% [a_t,b_t] = min(err_t);
result = {'heat' delta t1/50 t2/50 t3/50;
   '' '' '' a_1/50 b_1/50;
   '' '' '' a_2/50 b_2/50;
  '' '' '' a_3/50 b_3/50};
xlswrite('C:\Users\Yang Yanfei\Desktop\heat.xlsx',result)

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
semilogy(output.err,'.r-');hold on
% semilogy(err_lsmr,'.g-');hold on
semilogy(err_cgme,'.b-');hold on
semilogy(err_tcgme,'.g-');hold on
legend('JBDQR','hyb-cgme','hyb-tcgme')
title('relative error with noise level 0.05')

% axes('position',[0.4 0.2 0.5 0.4]);
% semilogy(err_tcgme,'.g-');hold on
% semilogy(output.err,'.r-');hold on
% semilogy(err_cgme,'.b-');hold on
% legend('hyb-tcgme','jbdqr')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% [reg_c,rho_c,eta_c] = l_corner(output.Rnrm,output.Xnrm);
% plot_lc(output.Rnrm,output.Xnrm,'o',2);hold on
%                
% ax = axis;hold on
% loglog([min(output.Rnrm)/100,rho_c],[eta_c,eta_c],':r',...
% [rho_c,rho_c],[min(output.Xnrm)/100,eta_c],':r')
% title(['L-curve ',' corner at ', num2str(reg_c), ' the relative error is ',num2str(output.err(reg_c))]);
% axis(ax);hold off
% 
% % figure;
% % [reg_c,rho_c,eta_c] = l_corner(rho_lsmr,eta_lsmr);
% % plot_lc(rho_lsmr,eta_lsmr,'o',2);hold on
% %                
% % ax = axis;hold on
% % loglog([min(rho_lsmr)/100,rho_c],[eta_c,eta_c],':r',...
% % [rho_c,rho_c],[min(eta_lsmr)/100,eta_c],':r')
% % title(['L-curve ',' corner at ', num2str(reg_c), ' the relative error is ',num2str(err_lsmr(reg_c))]);
% % axis(ax);hold off
% 
% figure;
% [reg_c,rho_c,eta_c] = l_corner(rho_lsmr,eta_lsmr);
% plot_lc(rho_lsmr,eta_lsmr,'o',2);hold on
%                
% ax = axis;hold on
% loglog([min(rho_lsmr)/100,rho_c],[eta_c,eta_c],':r',...
% [rho_c,rho_c],[min(eta_lsmr)/100,eta_c],':r')
% title(['L-curve ',' corner at ', num2str(reg_c), ' the relative error is ',num2str(err_lsmr(reg_c))]);
% axis(ax);hold off
% 
% % figure;
% % [reg_c,rho_c,eta_c] = l_corner(rho_tcgme,eta_tcgme);
% % plot_lc(rho_tcgme,eta_tcgme,'o',2);hold on
% %                
% % ax = axis;hold on
% % loglog([min(rho_tcgme)/100,rho_c],[eta_c,eta_c],':r',...
% % [rho_c,rho_c],[min(eta_tcgme)/100,eta_c],':r')
% % title(['L-curve ',' corner at ', num2str(reg_c), ' the relative error is ',num2str(err_tcgme(reg_c))]);
% % axis(ax);hold off
% 
% % 
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(2,2,1)
% imagesc(reshape(x_true,N,N))
% colormap gray, axis image off
% title('original image')
% subplot(2,2,2)
% imagesc(reshape(x_cgme(:,b_cgme),N,N))
% colormap gray, axis image off
% title('hyb-cgme')
% subplot(2,2,3)
% imagesc(reshape(x_tcgme(:,b_tcgme),N,N))
% colormap gray, axis image off
% title('hyb-tcgme')
% subplot(2,2,4)
% imagesc(reshape(x_out(:,b_jbdqr),N,N))
% colormap gray, axis image off
% title('JBDQR')
