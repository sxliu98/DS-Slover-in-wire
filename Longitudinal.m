%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Longitudinal wave         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
cL = 8.47*10^3;
cT = 5.34*10^3;
S = cL/cT;

beta = @(Xi,alpha) sqrt((S.^2-1).*Xi.^2 + S.^2.*alpha.^2);
Fun_real =@(Xi,alpha) (Xi.^2-beta(Xi,alpha).^2).^2.*alpha.*besselj(0,alpha).*besselj(1,beta(Xi,alpha))-2.*alpha.^2.*(Xi.^2+beta(Xi,alpha).^2).*besselj(1,alpha).*besselj(1,beta(Xi,alpha))+4.*Xi.^2.*alpha.^2.*beta(Xi,alpha).*besselj(0,beta(Xi,alpha)).*besselj(1,alpha);

beta_i = @(Xi,alpha) sqrt((S.^2-1).*Xi.^2 - S.^2.*alpha.^2);
Fun_image =@(Xi,alpha) -(Xi^2-beta_i(Xi,alpha)^2)^2*alpha*besseli(0,alpha)/besseli(1,alpha)-2*alpha^2*(Xi^2+beta_i(Xi,alpha)^2)+4*Xi^2*alpha^2*beta_i(Xi,alpha)*besselj(0,beta_i(Xi,alpha))/besselj(1,beta_i(Xi,alpha));

Nk = 40000;
XiMAX = 140;
XiMIN = 0.2;
XiX = linspace(XiMIN,XiMAX,Nk);
delta = 1e-4;

Fun_real_1 = @(alpha) Fun_real(XiX(1),alpha);
Fun_real_2 = @(alpha) Fun_real(XiX(2),alpha);

alpha_mesh = (0:0.01:80);

im = 1;
for ia = 1:length(alpha_mesh)-1
    if Fun_real_1(alpha_mesh(ia))*Fun_real_1(alpha_mesh(ia+1)) < 0
       alpha01(im) = fzero(Fun_real_1,[alpha_mesh(ia) alpha_mesh(ia+1)]);
       im = im+1;
    end
end

mode = length(alpha01);

im = 1;
for ia = 1:length(alpha_mesh)-1
    if Fun_real_2(alpha_mesh(ia))*Fun_real_2(alpha_mesh(ia+1)) < 0
       alpha02(im) = fzero(Fun_real_2,[alpha_mesh(ia) alpha_mesh(ia+1)]);
       im = im+1;
    end
end

[alphaX] = FindrealrootsL(Fun_real,XiX,alpha01,alpha02);

alphaX(abs(alphaX)<1e-4) = NaN;
betaX = beta(XiX,alphaX);

XiXX = zeros(1,mode);
indK = zeros(1,mode);

for im = 1:mode
    for ik = 3:Nk
        if isnan(alphaX(im,ik))
            XiXX(im) = XiX(ik);
            indK(im) = ik;
            break
        end
    end
end

XiXX(XiXX==0) = [];

Fun_image_1 = @(alpha) Fun_image(XiX(1),alpha);
alphai0(1) = fzero(Fun_image_1,0.05);

[NUMi,alphaiX] = FindimagerootsL(Fun_image,alphai0,XiX,XiXX);


alphaiX(alphaiX==0) = NaN;

alphaXX = zeros(mode,Nk);

alphaXX(1,:) = alphaiX(1,:);

for im = 2:length(XiXX)+1
    for ik = 1:Nk-1
        if isnan(alphaX(im-1,ik))
            alphaXX(im,ik) = alphaiX(im,ik);
        else
            alphaXX(im,ik) = alphaX(im-1,ik);
        end   
    end
end

for im = length(XiXX)+2:mode
    for ik = 1:Nk-1
            alphaXX(im,ik) = alphaX(im-1,ik);   
    end
end

alphaXX(1,1)=0;
alphaXX(alphaXX<0) = alphaXX(alphaXX<0)*(-1i);
OMEGAX = sqrt(alphaXX.^2+XiX.^2);
S = cL/cT;
OMEGAX_F4 = OMEGAX*S;

XiX1 = linspace(0,140,1000);
for im = 1:mode
    OMEGAX_F(im,:) = interp1(XiX,OMEGAX_F4(im,:),XiX1,'linear','extrap');
end

figure('OuterPosition',[0 0 800 600])
plot (XiX1,OMEGAX_F,"LineWidth",1.5)
xlim([0 120])
ylim([0 120])
set(gca,'linewidth',1.5,'FontSize',14);
xlabel('$ \xi =k\cdot R $','Interpreter','latex','FontSize',20,'FontWeight','bold');
ylabel('$ \Omega = \omega\cdot R/c_T $','Interpreter','latex','FontSize',20,'FontWeight','bold');

function [alphaX] = FindrealrootsL(Fun_real,XiX,alpha01,alpha02)
modei = length(alpha01);
N = length(XiX);
alphaX = zeros(modei,N);
exreal = zeros(modei,N);
alphaX(:,1) = alpha01';
alphaX(:,2) = alpha02';
NUM = zeros(modei,N);
for im = 1:modei
    delta = 1e-5;
    for ik = 3:N
        Fun_real_N = @(alpha) Fun_real(XiX(ik),alpha);
        num = 0;
        alphaleft(im) = 2*alphaX(im,ik-1) - alphaX(im,ik-2) + delta;
        alpharight(im) = 2*alphaX(im,ik-1) - alphaX(im,ik-2) - delta;
        while Fun_real_N(alphaleft(im))*Fun_real_N(alpharight(im)) > 0
            alphaleft(im) = alphaleft(im)+delta*(2*num);
            alpharight(im) = alpharight(im)+delta*(2*num);
            if Fun_real_N(alphaleft(im))*Fun_real_N(alpharight(im)) < 0
                break
            end
            alphaleft(im) = alphaleft(im)-delta*(2*num+1);
            alpharight(im) = alpharight(im)-delta*(2*num+1);
            num = num + 1;
        end
        NUM(im,ik) = num;
        [x, ~, exitflag, ~] = fzero(Fun_real_N,[alphaleft(im) alpharight(im)]);
        alphaX(im,ik) = x;
        exreal(im,ik) = exitflag;
        if exitflag == -5
            num = 0;
            Fun_real_N = @(alpha) Fun_real_N(alpha)/(alpha-x);
            while Fun_real_N(alphaleft(im))*Fun_real_N(alpharight(im)) > 0
                alphaleft(im) = alphaleft(im)+delta*(2*num);
                alpharight(im) = alpharight(im)+delta*(2*num);
                if Fun_real_N(alphaleft(im))*Fun_real_N(alpharight(im)) < 0
                    break
                end
                alphaleft(im) = alphaleft(im)-delta*(2*num+1);
                alpharight(im) = alpharight(im)-delta*(2*num+1);
                num = num + 1;
            end
            NUM(im,ik) = num;
            [x, ~, exitflag, ~] = fzero(Fun_real_N,[alphaleft(im) alpharight(im)]);
            alphaX(im,ik) = x;
            exreal(im,ik) = exitflag;
        end
    end
end
end

function [NUMi,alphaiX] = FindimagerootsL(Fun_image,alphai0,XiX,XiXX)
XiXX = [0 XiXX];
modei = length(XiXX);
N = length(XiX);
delta =1e-4;
alphaleft = alphai0(1)*ones(1,modei)-delta;
alpharight = alphai0(1)*ones(1,modei)+delta;
alphaiX = zeros(modei,N);
for im = 1:modei
    for ik = 1:N
        if XiX(ik) >= XiXX(im)
            Fun_image_N = @(alpha) Fun_image(XiX(ik),alpha);
            num = 0;
            while Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) >= 0 
                alphaleft(im) = alphaleft(im)+delta*(num-1);
                alpharight(im) = alpharight(im)+delta*(num-1);
                if Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) < 0
                    break
                end
            end
            NUMi(im,ik) = num;
            alphaiX(im,ik) = fzero(Fun_image_N,[alphaleft(im) alpharight(im)]);
        end
    end  
end
end