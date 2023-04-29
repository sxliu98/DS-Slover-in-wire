%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Flexural wave           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
cL = 8.47*10^3;
cT = 5.34*10^3;
S = cL/cT;

Chi_a = @(alpha) alpha*besselj(0,alpha)./besselj(1,alpha);
beta = @(Xi,alpha) sqrt((cL.^2-cT.^2)./cT.^2.*Xi.^2 + cL.^2./cT.^2.*alpha.^2);
Chi_b = @(Xi,alpha) beta(Xi,alpha)*besselj(0,beta(Xi,alpha))./besselj(1,beta(Xi,alpha));
f1 = @(Xi,alpha) 2*(beta(Xi,alpha).^2-Xi.^2).^2;
f2 = @(Xi,alpha) 2*beta(Xi,alpha).^2.*(5*Xi.^2+beta(Xi,alpha).^2);
f3 = @(Xi,alpha) beta(Xi,alpha).^6-10*beta(Xi,alpha).^4-2*beta(Xi,alpha).^4*Xi.^2+2*beta(Xi,alpha).^2*Xi.^2+beta(Xi,alpha).^2*Xi.^4+4*Xi.^4;
f4 = @(Xi,alpha) 2*beta(Xi,alpha).^2.*(2*beta(Xi,alpha).^2.*Xi.^2-beta(Xi,alpha).^2-9*Xi.^2);
f5 = @(Xi,alpha) beta(Xi,alpha).^2.*(-beta(Xi,alpha).^4+8*beta(Xi,alpha).^2-2*beta(Xi,alpha).^2*Xi.^2+8*Xi.^2-Xi.^4);
Phi = @(Xi,alpha) f1(Xi,alpha).*Chi_b(Xi,alpha).^2+f2(Xi,alpha).*Chi_a(alpha).*Chi_b(Xi,alpha)+f3(Xi,alpha).*Chi_b(Xi,alpha)+f4(Xi,alpha).*Chi_a(alpha)+f5(Xi,alpha);
Fun_real =@(Xi,alpha) besselj(1,alpha).*besselj(1,beta(Xi,alpha)).^2.*Phi(Xi,alpha);

Chi_ai = @(alpha) alpha.*besseli(0,alpha)./besseli(1,alpha);
beta_i = @(Xi,alpha) sqrt((cL.^2-cT.^2)./cT.^2.*Xi.^2 - cL.^2./cT.^2.*alpha.^2);
Chi_bi = @(Xi,alpha) beta_i(Xi,alpha).*besselj(0,beta_i(Xi,alpha))./besselj(1,beta_i(Xi,alpha));
f1i = @(Xi,alpha) 2*(beta_i(Xi,alpha).^2-Xi.^2).^2;
f2i = @(Xi,alpha) 2*beta_i(Xi,alpha).^2.*(5*Xi.^2+beta_i(Xi,alpha).^2);
f3i = @(Xi,alpha) beta_i(Xi,alpha).^6-10*beta_i(Xi,alpha).^4-2*beta_i(Xi,alpha).^4*Xi.^2+2*beta_i(Xi,alpha).^2*Xi.^2+beta_i(Xi,alpha).^2*Xi.^4+4*Xi.^4;
f4i = @(Xi,alpha) 2*beta_i(Xi,alpha).^2.*(2*beta_i(Xi,alpha).^2.*Xi.^2-beta_i(Xi,alpha).^2-9*Xi.^2);
f5i = @(Xi,alpha) beta_i(Xi,alpha).^2.*(-beta_i(Xi,alpha).^4+8*beta_i(Xi,alpha).^2-2*beta_i(Xi,alpha).^2*Xi.^2+8*Xi.^2-Xi.^4);
Phii = @(Xi,alpha) f1i(Xi,alpha).*Chi_bi(Xi,alpha).^2+f2i(Xi,alpha).*Chi_ai(alpha).*Chi_bi(Xi,alpha)+f3i(Xi,alpha).*Chi_bi(Xi,alpha)+f4i(Xi,alpha).*Chi_ai(alpha)+f5i(Xi,alpha);
Fun_image =@(Xi,alpha) besseli(1,alpha).*besselj(1,beta_i(Xi,alpha)).^2.*Phii(Xi,alpha);

Nk = 40000;
XiMAX = 140;
XiMIN = 0.2;
XiX = linspace(XiMIN,XiMAX,Nk);
delta = 1e-4;

Fun_image_1 = @(alpha) Fun_image(XiX(100),alpha);

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

[alphaX] = Findrealroots(Fun_real,XiX,alpha01,alpha02);

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

Fun_image_a = @(Xi) Fun_image(Xi,0.7);
Ximesh = linspace(0,140,1e5);

ip = 1;
for i = 2:1e5
    if ip > mode+1
        break
    end
    if Fun_image_a(Ximesh(i))*Fun_image_a(Ximesh(i-1))<0
        Xipoint(ip) = fzero(Fun_image_a,[Ximesh(i) Ximesh(i-1)]);
        ip = ip+1;
    end
end

[alphaiX] = Findimageroots(Fun_image,mode,XiX,Xipoint);

alphaXX = zeros(mode,Nk);

alphaXX(1,:) = -alphaiX(1,:);
for im = 2:mode
    for ik = 1:Nk-1
        if isnan(alphaX(im-1,ik))
            alphaXX(im,ik) = -alphaiX(im,ik);
        else
            alphaXX(im,ik) = alphaX(im-1,ik);
        end   
    end
end

alphaXX(alphaXX<0) = alphaXX(alphaXX<0)*(-1i);
OMEGAX = sqrt(alphaXX.^2+XiX.^2);
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


function [alphaX] = Findrealroots(Fun_real,XiX,alpha01,alpha02)
mode = length(alpha01);
N = length(XiX);
alphaX = zeros(mode,N);
exreal = zeros(mode,N);
alphaX(:,1) = alpha01';
alphaX(:,2) = alpha02';
NUM = zeros(mode,N);
for im = 1:mode
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

function [alphaiX] = Findimageroots(Fun_image,mode,XiX,Xipoint)
modei = mode + 1;
delta = 1e-4;
alphai0 = [0.4*ones(1,3) 0.7*ones(1,17) 1*ones(1,5) 2*ones(1,16) 3*ones(1,modei-41)];
alphaleft = alphai0-delta;
alpharight = alphai0+delta;
Nk = length(XiX);

alphaiX = zeros(modei,Nk);
for im = 1:modei
    ip = 1;
    for ik = 1:Nk
        if XiX(ik) > Xipoint(im) && ip<=10
            Fun_image_N = @(alpha) Fun_image(XiX(ik),alpha);
            num = 0;
            while Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) >= 0 
                alphaleft(im) = alphaleft(im)+delta*(2*num-1);
                alpharight(im) = alpharight(im)+delta*(2*num-1);
                if Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) < 0
                    break
                end
                alphaleft(im) = alphaleft(im)-delta*(2*num);
                alpharight(im) = alpharight(im)-delta*(2*num);
                num = num + 1;
            end
            NUMi(im,ik) = num;
            [x, ~, exitflag, ~] = fzero(Fun_image_N,[alphaleft(im) alpharight(im)]);
            alphaiX(im,ik) = x;
            eX(im,ik) = exitflag;
            ip = ip + 1;
        end
        if ip > 10
            Fun_image_N = @(alpha) Fun_image(XiX(ik),alpha);
            num = 0;
            alphaleft(im) = 2*alphaiX(im,ik-1) - alphaiX(im,ik-2) + delta;
            alpharight(im) = 2*alphaiX(im,ik-1) - alphaiX(im,ik-2) - delta;
            while Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) >= 0
                alphaleft(im) = alphaleft(im)+delta*(2*num-1);
                alpharight(im) = alpharight(im)+delta*(2*num-1);
                if Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) < 0
                    break
                end
                alphaleft(im) = alphaleft(im)-delta*(2*num);
                alpharight(im) = alpharight(im)-delta*(2*num);
                num = num + 1;
            end
            NUMi(im,ik) = num;
            [x, ~, exitflag, ~] = fzero(Fun_image_N,[alphaleft(im) alpharight(im)]);
            alphaiX(im,ik) = x;
            eX(im,ik) = exitflag;
        end
    end  
end

for im = 1:modei
    for ik = Nk-2:-1:1
        Fun_image_N = @(alpha) Fun_image(XiX(ik),alpha);
        num = 0;
        alphaleft(im) = 2*alphaiX(im,ik+1) - alphaiX(im,ik+2) + delta;
        alpharight(im) = 2*alphaiX(im,ik+1) - alphaiX(im,ik+2) - delta;
        while Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) >= 0
            alphaleft(im) = alphaleft(im)+delta*(2*num-1);
            alpharight(im) = alpharight(im)+delta*(2*num-1);
            if Fun_image_N(alphaleft(im))*Fun_image_N(alpharight(im)) < 0
                break
            end
            alphaleft(im) = alphaleft(im)-delta*(2*num);
            alpharight(im) = alpharight(im)-delta*(2*num);
            num = num + 1;
        end
        NUMi(im,ik) = num;
        [x, ~, exitflag, ~] = fzero(Fun_image_N,[alphaleft(im) alpharight(im)]);
        alphaiX(im,ik) = x;
        eX(im,ik) = exitflag;
    end  
end

alphaiX(abs(alphaiX)<1e-4) = NaN;
alphaiX = alphaiX(2:modei,:);
end