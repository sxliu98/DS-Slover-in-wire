%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Torsional wave          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
XiMAX = 140;
Nk = 1000;
XiX = linspace(0,XiMAX,Nk);
funT = @(beta) beta*besselj(0,beta)-2*besselj(1,beta);
mode = 80;
funT_new = @(beta) funT(beta);
betaX = zeros(1,mode);
betaX(1) = 0;
betaX(2) = fzero(funT_new,5);
for i = 3:mode
    funT_new = @(beta) funT_new(beta)./(beta-betaX(i-1));
    betaX(i) = fzero(funT_new,betaX(i-1)+0.1);
end
WX = zeros(mode,Nk);
for i = 1:mode
    for j = 1:Nk
        WX(i,j) = (betaX(i)^2 + XiX(j)^2)^0.5;
    end
end
figure('OuterPosition',[0 0 800 600])
plot(XiX,WX,"LineWidth",1.5)
xlim([0 120])
ylim([0 120])
set(gca,'linewidth',1.5,'FontSize',14);
xlabel('$ \xi=k\cdot R $','Interpreter','latex','FontSize',20,'FontWeight','bold');
ylabel('$ \Omega = \omega\cdot R/c_T $','Interpreter','latex','FontSize',20,'FontWeight','bold');
