% initialize
clear;clc;
% flow condition
u_s=1.41; 
delta=3.2; % turbulent boundary thickness
h=5; % water depth
kappa=0.41; % karman constant
A=5.3; % parameter in mean velocity profile formula
%dp=0.25/10;U_s=1.23;
%dp=0.37/10;U_s=1.34;
dp=0.5/10;U_s=1.53; % particle diameter cm, shear velocity
% dp=0.5/10;U_s=1.54;
%dp=1/10;U_s=1.64;

%nu=dp*U_s/2.86;max_h=0.9494*h;
%nu=dp*U_s/4.59;max_h=0.9310*h;
nu=dp*U_s/6.79;max_h=0.884*h; % kinematic viscosity
% nu=dp*U_s/7.71;max_h=0.884*h;
%nu=dp*U_s/14.92;max_h=0.6064*h;

%up&intensity from NN2009
xlsfile='normalized reynolds stress.xlsx';xlsfile0='up.xlsx';xlsfile1='normalized u intensity.xlsx';xlsfile2='length scale.xlsx';xlsfile3='normalized v intensity.xlsx';
%[y_n]=xlsread(xlsfile,'A1:A32');[uv_n]=xlsread(xlsfile,'B1:B32');[yp]=xlsread(xlsfile0,'A1:A21');[up_uf]=xlsread(xlsfile0,'B1:B21');[y_d]=xlsread(xlsfile1,'A1:A22');[int_us]=xlsread(xlsfile1,'B1:B22');[y_l]=xlsread(xlsfile2,'A1:A23');[L_x]=xlsread(xlsfile2,'B1:B23');[y_dd]=xlsread(xlsfile3,'A1:A23');[int_vs]=xlsread(xlsfile3,'B1:B23');
%[y_n]=xlsread(xlsfile,'C1:C41');[uv_n]=xlsread(xlsfile,'D1:D41');[yp]=xlsread(xlsfile0,'C1:C32');[up_uf]=xlsread(xlsfile0,'D1:D32');[y_d]=xlsread(xlsfile1,'C1:C23');[int_us]=xlsread(xlsfile1,'D1:D23');[y_l]=xlsread(xlsfile2,'C1:C23');[L_x]=xlsread(xlsfile2,'D1:D23');[y_dd]=xlsread(xlsfile3,'C1:C23');[int_vs]=xlsread(xlsfile3,'D1:D23');
[y_n]=xlsread(xlsfile,'E1:E33');[uv_n]=xlsread(xlsfile,'F1:F33');[yp]=xlsread(xlsfile0,'E1:E31');[up_uf]=xlsread(xlsfile0,'F1:F31');[y_d]=xlsread(xlsfile1,'E1:E23');[int_us]=xlsread(xlsfile1,'F1:F23');[y_l]=xlsread(xlsfile2,'E1:E23');[L_x]=xlsread(xlsfile2,'F1:F23');[y_dd]=xlsread(xlsfile3,'E1:E23');[int_vs]=xlsread(xlsfile3,'F1:F23');
%[y_n]=xlsread(xlsfile,'G1:G42');[uv_n]=xlsread(xlsfile,'H1:H42');[yp]=xlsread(xlsfile0,'G1:G39');[up_uf]=xlsread(xlsfile0,'H1:H39');[y_d]=xlsread(xlsfile1,'G1:G23');[int_us]=xlsread(xlsfile1,'H1:H23');[y_l]=xlsread(xlsfile2,'G1:G23');[L_x]=xlsread(xlsfile2,'H1:H23');[y_dd]=xlsread(xlsfile3,'G1:G23');[int_vs]=xlsread(xlsfile3,'H1:H23');

%ws=0.165;wp=0.3123;h_r=dp;
%ws=0.37;wp=0.2973;h_r=dp;
ws=0.65;wp=0.2743;h_r=dp;
% ws=0.65;wp=0.2743;h_r=dp;
% ws=0.73;
%ws=1.4;wp=0.3627;h_r=dp;

rho_p=1.2;rho_f=1;Res=U_s*h/nu;y100=100*nu/U_s; % particle density, fluid density, shear Reynolds number, y^+ = 100

rtb=10;

deltat=0.05;n=500;m=1;simu=1;
%% model
timer1=tic;
for k=1:simu
    clear X Y meet u0 tt tn rt thetao drx dry ddx ddy dr lx Lx sigm_u sig_v theta Ts ue uv uv_A uv_B vb Y_plus Yi
    X=zeros(n,m);Y=zeros(n,m);Yi=zeros(n,m);Y_plus=zeros(n,m);Y(1,:)=delta-dp;X(1,:)=0;i=1;
    for j=1:m
        i=1;
        for I=1:n-1
            eta=Y(i,j)/delta;
            ypl=Y(i,j)*U_s/nu;

            Lx(i,j)=interp1(y_l,L_x,ypl,'linear');
            ll=linspace(6.9115*(100*nu/U_s)^0.74,3.5*3.2,100);
            r=rtb*(nu^3/(U_s^2/sqrt(0.09)*exp(-(Y(i,j)*U_s/nu-(0.3*Res-100))/(0.55*Res-13)))^(3/2).*ll).^0.25;
            p=r.^(4/3);
        
            uv(i,j)=interp1(y_n,uv_n,eta,'linear')*U_s^2;
            uv_A(i,j)=(1-eta+eta*log(eta))*U_s^2;
            if uv(i,j)>0&&uv_A(i,j)>uv(i,j)
                    uv_A(i,j)=uv(i,j);
            end
            uv_B(i,j)=uv(i,j)-uv_A(i,j);
            if uv(i,j)<0
                uv_B(i,j)=uv(i,j);uv_A(i,j)=0;
            end
            TA(i,j)=uv_A(i,j)/uv(i,j);TB(i,j)=uv_B(i,j)/uv(i,j);

            if Y(i,j)==h_r
               vb(i,j)=normrnd(log(5.52*U_s),sqrt(0.123));%wu and lin
               B=sqrt(3*dp*(rho_p-rho_f)*981/(3*rho_f*0.178));
               if vb(i,j)>log(B)
                   con=1;
               else
                   con=0;
               end
            else
                con=1;
            end
            if con==1
                a=rand;
                if Y(i,j)>=0.7*h
                    meet(i,j)=2;
%                     rt(i,j)=50*nu/U_s;
                    lx(i,j)=3.5*delta;
                    rt(i,j)=rtb*(nu^3/(U_s^2/sqrt(0.09)*exp(-(Y(i,j)*U_s/nu-(0.3*Res-100))/(0.55*Res-13)))^(3/2)*lx(i,j))^0.25;
                    ue(i,j)=U_s*(1/0.41*log(Y(i,j)*U_s/nu)+5.1)+U_s*wp/0.41*(2*eta^2*(3-2*eta)-eta^2/wp*(1-eta)*(1-2*eta));
                    Ts(i,j)=lx(i,j)/ue(i,j);
                    tt(i,j)=rand*Ts(i,j);tn(i,j)=ceil(tt(i,j)/deltat);
                elseif Y(i,j)>=delta&&Y(i,j)<=0.7*h
                    meet(i,j)=3;
%                     rt(i,j)=(10+rand*40)*nu/U_s;
                    rt(i,j)=randpdf(p,r,[1,1]);
                    lx(i,j)=(rt(i,j)/rtb)^4*(U_s^2/sqrt(0.09)*exp(-(Y(i,j)*U_s/nu-(0.3*Res-100))/(0.55*Res-13)))^(3/2)/nu^3;
                    int=2.3*exp(-eta)*U_s;
                    ubar=U_s*(1/kappa*log(ypl)+A);
                    Ts(i,j)=Y(i,j)/int^2/U_s*int*U_s/ubar*0.41/0.36;
                    tt(i,j)=rand*Ts(i,j);tn(i,j)=ceil(tt(i,j)/deltat);
                elseif a<=TA(i,j)
                    meet(i,j)=1;
%                     rt(i,j)=(10+rand*40)*nu/U_s;
%                     lx(i,j)=6.9115*((rt(i,j)-10*nu/U_s)/(40*nu/U_s)*(0.6*delta-y100)+y100)^0.74;
                    rt(i,j)=randpdf(p,r,[1,1]);
                    lx(i,j)=(rt(i,j)/rtb)^4*(U_s^2/sqrt(0.09)*exp(-(Y(i,j)*U_s/nu-(0.3*Res-100))/(0.55*Res-13)))^(3/2)/nu^3;
                    ue(i,j)=U_s*(1/0.41*log(Y(i,j)*U_s/nu)+5.1);
                    Ts(i,j)=lx(i,j)/ue(i,j);
                    tt(i,j)=rand*Ts(i,j);tn(i,j)=ceil(tt(i,j)/deltat);
                else
                    meet(i,j)=2;
%                     rt(i,j)=50*nu/U_s;
                    lx(i,j)=3.5*delta;
                    rt(i,j)=rtb*(nu^3/(U_s^2/sqrt(0.09)*exp(-(Y(i,j)*U_s/nu-(0.3*Res-100))/(0.55*Res-13)))^(3/2)*lx(i,j))^0.25;
                    ue(i,j)=U_s*(1/0.41*log(Y(i,j)*U_s/nu)+5.1)+U_s*wp/0.41*(2*eta^2*(3-2*eta)-eta^2/wp*(1-eta)*(1-2*eta));
                    Ts(i,j)=lx(i,j)/ue(i,j);
                    tt(i,j)=rand*Ts(i,j);tn(i,j)=ceil(tt(i,j)/deltat);
                end
                meet(i+1:i+tn(i,j)-1,j)=meet(i,j);
                rt(i+1:i+tn(i,j)-1,j)=rt(i,j);
                tt(i+1:i+tn(i,j)-1,j)=tt(i,j);
                c_theta=0;Ubarc=0;Vbarc=0;dirc=rand;tdx=-1;tdy=1;dircx=-1;dircy=1;
                if dirc>0.25&&dirc<=0.5
                        tdx=-1;tdy=-1;dircx=-1;dircy=-1;
                    elseif dirc>0.5&&dirc<=0.75
                        tdx=1;tdy=1;dircx=1;dircy=1;
                    elseif dirc>0.75&&dirc<=1
                        tdx=1;tdy=-1;dircx=1;dircy=-1;
                end
%                 if Y(i,j)==h_r
%                     tdy=1;
%                 end

                for ii=0:tn(i,j)-1
                    z=i+ii;
                    eta=Y(z,j)/delta;
                    Y_plus(z,j)=Y(z,j)/nu*U_s;
                    if Y_plus(z,j)>16
                        U_bar(z,j)=U_s*(1/kappa*log(Y_plus(z,j))+A);
                    else
                        U_bar(z,j)=Y_plus(z,j)*U_s;
                    end
                    sigmau0(z,j)=U_bar(i,j)/sqrt(pi/2);
                    u0(z,j)=raylrnd(sigmau0(z,j));
                    sigmauf=2.3*exp(-eta)*U_s;
                    uf=normrnd(0,sigmauf);
                    u0(z,j)=U_bar(z,j)+uf;
                    u0(z,j)=U_bar(z,j);
                    
                    
                    thetao(z,j)=u0(z,j)*deltat/rt(z,j);
                    theta(z,j)=rem(thetao(z,j),2*pi);
                    c_theta=rem(c_theta+theta(z,j),2*pi);
                    if c_theta>=pi && c_theta<=2*pi && Y(i,j)~=h_r
                        tdx=-tdx;
%                     elseif Y(i,j)==h_r
%                         c_theta=rem(c_theta+theta(z,j),pi);
%                         if c_theta>=pi/2
%                             tdx=-tdx;
%                         end
                    end


                    drr(z,j)=2*rt(z,j)*abs(sin(c_theta));
                    drxx(z,j)=tdx*drr(z,j)*abs(sin(c_theta/2));
                    dryy(z,j)=tdy*drr(z,j)*abs(cos(c_theta/2));
%                     ddx(z,j)=4*rt(z,j)*abs(cos(c_theta))*sin(c_theta)^2*sin(c_theta/2)^2/deltat/(ii+1)-2*rt(z,j)^2/deltat/(ii+1)*(rt(z,j)*cos(c_theta)*dircx*u0(z,j)*deltat/(ii+1)*(rt(z,j)^2)^(-1.5))*abs(cos(c_theta/2)*sin(c_theta/2)*abs(sin(c_theta))^2+2*sin(c_theta/2)^2*abs(sin(c_theta))*sign(sin(c_theta))*cos(c_theta));
%                     ddy(z,j)=4*rt(z,j)*abs(sin(c_theta))*sin(c_theta)^2*cos(c_theta/2)^2/deltat/(ii+1)-2*rt(z,j)^2/deltat/(ii+1)*(rt(z,j)*sin(c_theta)*dircy*u0(z,j)*deltat/(ii+1)*(rt(z,j)^2)^(-1.5))*abs(2*cos(c_theta/2)^2*abs(sin(c_theta))*sign(sin(c_theta))*cos(c_theta)-cos(c_theta/2)*sin(c_theta/2)*abs(sin(c_theta))^2);
                    
                    if ii==0
                        drx(z,j)=drxx(z,j);
                        dry(z,j)=dryy(z,j);
                    else
                        drx(z,j)=drxx(z,j)-drxx(z-1,j);
                        dry(z,j)=dryy(z,j)-dryy(z-1,j);
                    end
                    
                    ddx(z,j)=0;
%                     ddy(z,j)=U_s/kappa/Y(z,j)*deltat/rt(i,j)*(2*rt(i,j)^2/deltat*(2*abs(sin(theta(z,j)))*sign(sin(theta(z,j)))*cos(theta(z,j))*(cos(theta(z,j)/2))^2-(sin(theta(z,j)))^2*sin(theta(z,j)/2)*abs(cos(theta(z,j)/2))*sign(cos(theta(z,j)/2))));
%                     ddy(z,j)=(U_s^2*Y(z,j)*exp(-((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50))*(1/((U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)*(rho_f/rho_p-1))-1)-(U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)*(rho_f/rho_p-1)^2)))/(nu*((23*Res)/50-299/50))-U_s*Y(z,j)*exp(-((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50))*(((U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)^2*(rho_f/rho_p-1))+(U_s^2*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(z,j)*(rho_f/rho_p-1)*((11*Res)/20-13))-(U_s^2*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(z,j)*(rho_f/rho_p-1)*((23*Res)/50-299/50)))/((U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)*(rho_f/rho_p-1))-1)^2+(U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)^2*(rho_f/rho_p-1)^2)+(U_s^2*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(z,j)*(rho_f/rho_p-1)^2*((11*Res)/20-13))-(U_s^2*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(z,j)*(rho_f/rho_p-1)^2*((23*Res)/50-299/50)))-U_s*exp(-((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50))*(1/((U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)*(rho_f/rho_p-1))-1)-(U_s*ws*exp(-((U_s*Y(z,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(z,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(z,j)*(rho_f/rho_p-1)^2));
                    ddy(z,j)=(U_s^2*Y(i,j)*exp(-((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50))*(1/((U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)*(rho_f/rho_p-1))-1)-(U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)*(rho_f/rho_p-1)^2)))/(nu*((23*Res)/50-299/50))-U_s*Y(i,j)*exp(-((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50))*(((U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)^2*(rho_f/rho_p-1))+(U_s^2*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(i,j)*(rho_f/rho_p-1)*((11*Res)/20-13))-(U_s^2*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(i,j)*(rho_f/rho_p-1)*((23*Res)/50-299/50)))/((U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)*(rho_f/rho_p-1))-1)^2+(U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)^2*(rho_f/rho_p-1)^2)+(U_s^2*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(i,j)*(rho_f/rho_p-1)^2*((11*Res)/20-13))-(U_s^2*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*nu*Y(i,j)*(rho_f/rho_p-1)^2*((23*Res)/50-299/50)))-U_s*exp(-((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50))*(1/((U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)*(rho_f/rho_p-1))-1)-(U_s*ws*exp(-((U_s*Y(i,j))/nu-(3*Res)/10+100)/((11*Res)/20-13))*exp(((17*Res)/50+(U_s*Y(i,j))/nu-23/2)/((23*Res)/50-299/50)))/(3270*Y(i,j)*(rho_f/rho_p-1)^2));
                    
%                     X(z+1,j)=X(i,j)+(Ubarc+U_bar(z,j))*deltat+drx(z,j);
%                     Y(z+1,j)=Y(i,j)+(Vbarc+ddy(z,j)-ws)*deltat+dry(z,j);
                    
                    X(z+1,j)=X(z,j)+(U_bar(i,j))*deltat+drx(z,j);
                    Y(z+1,j)=Y(z,j)+(ddy(z,j)-ws)*deltat+dry(z,j);
                    
                    Ubarc=Ubarc+U_bar(z,j);    
                    Vbarc=Vbarc+ddy(z,j)-ws;

                    
                    if Y(z+1,j)>h
                        Y(z+1,j)=h;
                    elseif Y(z+1,j)<h_r
                        Y(z+1,j)=h_r;i=z+1;break
                    end
                end
            elseif con==0
                X(i+1,j)=X(i,j);Y(i+1,j)=h_r;
                z=i;
            end
            i=z+1;
            if i>n+1
                break
            end
        end

    end
    
    XX{k}=X;YY{k}=Y;meett{k}=meet;ttt{k}=tt;rtt{k}=rt;tnn{k}=tn;thetaoo{k}=thetao;ddrx{k}=drx;ddry{k}=dry;u00{k}=u0;
    disp(k)
end
toc(timer1)
ex=mean(X(1:n,:),2);ey=mean(Y(1:n,:),2);plot(ex,ey);
% save('D:\eddytrap50_rt4_10000p_3_bar.mat','XX','YY','meett','ttt','rtt','thetaoo','tnn','ddrx','ddry','u00','-v7.3')
% save('D:\eddytrap50_rt9_5_10000p_3_bar.mat','XX','YY','-v7.3');