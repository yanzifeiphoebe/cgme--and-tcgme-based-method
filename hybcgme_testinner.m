clc;
clear all;
close all;


maxit = 200;
%%%%%%%%%%%%%% 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 1000;
% [A,b,x_true] = shaw(n);
% noise  = randn(size(b));
% delta=0.1;
% sigma = delta*(norm(b)/norm(noise));
% N = sigma*noise;
% b = b + N;
% % L = get_l(n,1);
% %                                            
% % %%%%%%%%%%%%%%%%%% 2-D  RestoreTools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 N = 256;
% A = blur(N,16,2);
% x_true=imread('mri.png'); 
% load mri
% x_true = double(D(1:N,1:N,15));
% load GaussianBlur440
% load GaussianBlur422
% load GaussianBlur420
% load VariantMotionBlur_small
% load satellite
load grain
% load AtmosphericBlur50
% load AtmosphericBlur30
% load AtmosphericBlur10
% load VariantGaussianBlur1
% load VariantGaussianBlur2
% load VariantGaussianBlur3
K = psfMatrix(PSF);
% x_true=im2double(f_true);
% x_true=x_true(1:N,1:N);
x_true=x_true(:);
% b=g(:);
btrue = K*x_true;
delta=0.001;
e =  randn(size(btrue(:)))/norm(randn(size(btrue(:))));
b = btrue(:) + norm(btrue(:))*delta*e;
%  %%%%%%%%%%%%%%%% 2-D IRtools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A, btrue, x_true, ProbInfo] = PRblurspeckle;
% [b, NoiseInfo] = PRnoise(btrue, 'gauss', delta);
%%%%%%%%%%%%%%%%%%%%%%%% L matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
I = speye(N);
L1=get_l(N,2);
L = [kron(I,L1); kron(L1,I)];

%%%% CGME %%%%%%%%%%%
t2=tic;
[X,V_cgme] = yangcgme(K,b,maxit,1);
[x_04,x_05,x_06,err_04,err_05,err_06] = YangfunTik_inner(K, x_true,X,V_cgme,L,maxit);
toc

% %%%%%%%% TCGME %%%%%%%%%%%%%%%%%

% [X_t,V_tcgme] = yangtcgme(K,b,maxit,1);
% [x_04,x_05,x_06,err_04,err_05,err_06] = YangfunTik_inner(K, x_true,X_t,V_tcgme,L,maxit);

% % %%%%%%% LSQR %%%%%%%%%%
% tic;
% [X_lsqr,V_lsqr] = Yanglsqr_stop(K,b,x_true, L, maxit+1,1);
% [x_lsqr,err_lsqr,rho_lsqr,eta_lsqr] = YangfunTik(K, b,x_true,X_lsqr,V_lsqr,L,maxit);
% toc;
% 
% 
% %%%%%%% the best error and the associated number of the iteration %%%%%%%%%%

[a_04,b_04] = min(err_04);
[a_05,b_05] = min(err_05);
[a_06,b_06] = min(err_06);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
semilogy(err_04,'.r-');hold on
semilogy(err_05,'*b-');hold on
semilogy(err_06,'og-');hold on
legend('tol=1e-4','tol=1e-5','tol=1e-6')
title('relative error')
% axes('position',[0.4 0.2 0.5 0.4]);
% semilogy(err_tcgme,'.g-');hold on
% semilogy(err_cgme,'.b-');hold on
% legend('ptcgme','pcgme')
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
% title('pcgme')
% subplot(2,2,3)
% imagesc(reshape(x_tcgme(:,b_tcgme),N,N))
% colormap gray, axis image off
% title('pcgme')
% subplot(2,2,4)
% imagesc(reshape(x_out(:,b_jbdqr),N,N))
% colormap gray, axis image off
% title('JBDQR')