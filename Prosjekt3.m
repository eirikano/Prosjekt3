%Prosjekt 3 - Eirik Andre Nordbø, Tobias Mohn Werner
%

clear all
clc
close all
%---------------------------Constants----------------------------------
c.Tf = 660.4+273;   %[K]        Melting point of pure Al
c.Te = 577+273;     %[K]        Eutectic point Al-Si system
c.Ce = 12.6;        %[wt%]      Eutectic comp. Al-Si system
c.Cs = 1.8;         %[wt%]      Max solubility of Si in Al(s) Al-Si system
c.k = 0.13;         %           Partition coefficient (=Cl/Cs)
c.pc = 2.58*10^6;   %[J/K*m^3]  Volume heat capacity
c.dHf =9.5*10^8;    %[J/m^3]    Latent heat of fusion
c.ks = 210;         %[W/K*m]    Thermal conductivity (solid)
c.kl = 95;          %[W/K*m]    Thermal conductivity (liquid)
c.m=(c.Tf-c.Te)/c.Ce; %[K/wt%]
%----------------------------------------------------------------------


%---------------------Assumptions for comparison-----------------------
C0_1=1; %[wt%Si]
C0_2=8; %[wt%Si]

TL(1)=c.Tf -c.m*C0_1;
TL(2)=c.Tf -c.m*C0_2;
%----------------------------------------------------------------------

%--------------Lever rule equilibrium fraction solid ------------------
fs_eq=@(T,TL) (1/(1-c.k))*((TL-T)/(c.Tf-T));
dfs_eq=@(T,TL) (1/(c.k-1))*((TL-c.Tf)/(c.Tf-T)^2);
%----------------------------------------------------------------------

%c)
dT=0.5; %Temperature step
%initial values
T(1,1)=TL(1);
T(1,2)=TL(2);
fs(1,1)=0; 
dfs(1,1)=0; %df/dT
fs(1,2)=0;
dfs(1,2)=0;
j=1;
while (fs(j,1)<1) || (fs(j,2)<1)
    T(j+1,1)=T(j,1)-dT;
    T(j+1,2)=T(j,2)-dT;
    fs(j+1,1)=fs_eq(T(j+1,1),TL(1));
    dfs(j+1,1)=dfs_eq(T(j+1,1),TL(1));
    fs(j+1,2)=fs_eq(T(j+1,2),TL(2));
    dfs(j+1,2)=dfs_eq(T(j+1,2),TL(2));
    j=j+1;
end

%-----------------------Eutectic fraction------------------------------
fe_eq=@(TL) 1-(1/(1-c.k))*((TL-c.Te)/(c.Tf-c.Te));
%----------------------------------------------------------------------

figure
subplot(2,1,1)

for i=1:length(TL)
    plot(T(:,i)-273,fs(:,i))
    hold on
end

grid
xlabel('T(�C)')
ylabel('Solid fraction, f');

subplot(2,1,2)
for i=1:length(TL)
    plot(T(:,i)-273,dfs(:,i))
    hold on
end
axis([0 c.Tf -inf inf]);
set(gca,'xdir','reverse');
title('Solid fraction: Lever rule vs Scheil')
grid
xlabel('T(°C)')
ylabel('Change in solid fraction, df/dT');


% oppgave 1d

fs_star_eq=@(T,TL) 1-((c.Tf-T)/(c.Tf-TL))^(1/(c.k-1));
clear T 
T(1,1)=TL(1);
T(1,2)=TL(2);

fs_star(1,2)=0;
fs_star(1,1)=0;
j=1;
while T(j,1)>273 || T(j,2)>273
    T(j+1,1)=T(j,1)-dT;
    T(j+1,2)=T(j,2)-dT;
    fs_star(j+1,1)=fs_star_eq(T(j+1,1),TL(1));
    fs_star(j+1,2)=fs_star_eq(T(j+1,2),TL(2));
    j=j+1;
end

figure(1)
subplot(2,1,1)
plot(T(:,1)-273,fs_star(:,1),'--b',T(:,2)-273,fs_star(:,2),'--r')

%------------equilibrium fraction scheil----
fe_eq_scheil=@(TL) ((c.Tf-c.Te)/(c.Tf-TL))^(1/(c.k-1));
%------------------------------------
%Creating the perfect plot------------------------------------------
axis([0 c.Tf 0 1]);
set(gca,'xdir','reverse');
title('Solid fraction: Lever rule vs Scheil')
fe1_lever=fe_eq(TL(1));
fe2_lever=fe_eq(TL(2));
fe1_scheil=fe_eq_scheil(TL(1));
fe2_scheil=fe_eq_scheil(TL(2));
if fe1_lever<0
    fe1_lever=0;
elseif fe2_lever<0
    fe2_lever=0;
end
strstart1=['Lever rule, ',num2str(C0_1)];
strend1=['wt%Si, f_e = ', num2str(fe1_lever,'%.3f')];
strstart2=['Lever rule, ',num2str(C0_2)];
strend2=['wt%Si, f_e = ', num2str(fe2_lever,'%.3f')];
str1=strcat(strstart1,strend1);
str2=strcat(strstart2,strend2);
strstart3=['Scheil, ',num2str(C0_1)];
strend3=['wt%Si, f_e = ', num2str(fe1_scheil,'%.3f')];
strstart4=['Scheil, ',num2str(C0_2)];
strend4=['wt%Si, f_e = ', num2str(fe2_scheil,'%.3f')];
str1=strcat(strstart1,strend1);
str2=strcat(strstart2,strend2);
str3=strcat(strstart3,strend3);
str4=strcat(strstart4,strend4);
legend(str1,str2,str3,str4,'Location','southeast')
%---------------------------------------------------------------------


clear T 
T(1,1)=TL(1);
T(1,2)=TL(2);
dfs_star_eq = @(T,TL) (1)/(c.k-1)*((1)/(c.Tf-TL))^((1)/(c.k-1))*(c.Tf-T)^((2-c.k)/(c.k-1));

j=1;
while T(j,1)>273 || T(j,2)>273
    T(j+1,1)=T(j,1)-dT;
    T(j+1,2)=T(j,2)-dT;
    dfs_star(j+1,1)=dfs_star_eq(T(j+1,1),TL(1));
    dfs_star(j+1,2)=dfs_star_eq(T(j+1,2),TL(2));
    j=j+1;
end
figure(1)
subplot(2,1,2)
plot(T(:,1)-273,-dfs_star(:,1),'--b',T(:,2)-273,-dfs_star(:,2),'--r')

%e) 

t_star=1; 
%----------------------Johnson-Mehl-Avrami equation------------------------
X_JMA=@(Xc,t,n) 1-(1-Xc)^((t/t_star)^n);
%--------------------------------------------------------------------------


% n=1;
% %Asking user for n value
% prompt = {'Enter dimention size, n=[1,2,3]:'};
% dlg_title = 'Input';
% num_lines = 1;
% defaultans = {num2str(n)};
% answer=inputdlg(prompt,dlg_title,num_lines,defaultans);
% if str2num(answer{1})>3 || str2num(answer{1})<1
%     display('Not a valid dimention!')
%     return
% end
% n = str2num(answer{1});
% %-----------------------
Xc=[0.05, 0.15];
n=[1,2,3];
dt=0.6;
X=ones(length(n),length(Xc),1000); %Creating a 3D matrix with 1000 time values.
for k=1:length(n)
    for i=1:length(Xc)
        j=1;
        X(k,i,j)=0;
        t(j)=0;
        while X(k,i,j)<0.97
            t(j+1)=t(j)+dt;
            X(k,i,j+1) = X_JMA(Xc(i),t(j+1),n(k));
            j=j+1;
        end
        t_plot=t; %saves the longest time vector for plotting
    end
end

figure
str1=['Xc = ', num2str(Xc(1))];
str2=['Xc = ', num2str(Xc(2))];  

for k=1:length(n)
    subplot(3,1,k)
    str_n=[', n = ', num2str(n(k))];
    for i=1:length(Xc)
        X_n(i,:)=X(k,i,:);
        plot(t_plot./t_star,X_n(i,1:length(t_plot)));
        hold on
    end
    if k==1
       title('JMA eq. X vs t/t*') 
    elseif k==2
        ylabel('X')
    end
    legend(strcat(str1,str_n),strcat(str2,str_n));  
end

xlabel('t/t*')

%h)
%dX/dt vs X
%--------------------------------------------------------------------------
dXdt_eq = @(t_star,Xc,n,X) (n.*(1-X).*log(1-X))./(t_star.*(log(1-X)./log(1-Xc)).^(1/n));
%--------------------------------------------------------------------------

n=[1,2,3];
Xc=[0.05, 0.15];
t_star=1;
Xgrid=1000;
dXdt=ones(length(n),length(Xc),Xgrid);

for k=1:length(n)
    for i=1:length(Xc)
        X(k,i,:)=linspace(0,1,Xgrid); %Creating a 3D matrix with 1000 time values.
        dXdt(k,i,:) = dXdt_eq(t_star,Xc(i),n(k),X(k,i,:));
    end
end

figure
str1=['Xc = ', num2str(Xc(1))];
str2=['Xc = ', num2str(Xc(2))];  

for k=1:length(n)
    subplot(3,1,k)
    str_n=[', n = ', num2str(n(k))];
    for i=1:length(Xc)
        Xplot(:)=X(k,i,:);
        dXdtplot(:)=dXdt(k,i,:);
        plot(Xplot,dXdtplot);
        hold on
    end
    if k==1
       title('X vs dX/dt') 
    elseif k==2
        ylabel('dX/dt')
    end
    legend(strcat(str1,str_n),strcat(str2,str_n));  
end
xlabel('X')

%%
%Heat flow model

clear n t X TL Xc fs
%-----------------------------Input parameters-----------------------------
t_r=6;               %s
C0_r=4;              %wt%Si
TL=c.Tf -c.m*C0_r;%K
T_n=TL-2;        %K
Xc=0.05;      
a=1;                 %K/s
Neq=1;               %Neq=Nr/N
%--------------------------------------------------------------------------



%=========================Equilibrium Lever rule===========================
n=1; %can be changed to n=[1,2,3]
%Finding value for fm_r
dT_r=TL-T_n;
fm_r=fs_eq(T_n,TL);


%----------------------------Equation for t*-------------------------------
t_star_eq=@(dT,C0,fm,n) t_r*((dT_r/dT)^2)*(C0/C0_r)*((Neq)^(1/n))*(((fm/fm_r))^(1/n));
%--------------------------------------------------------------------------

%Starting value-----
t_star_i=t_star_eq(dT_r,C0_r,fm_r,n);
t=1;
dt=1;
dX=-dt*((1-Xc)^((t/t_star_i)^n))*log(1-Xc)*((t/t_star_i)^n)*(n/t);
X(1)=0+dX;
%-------------------
j=1;
k=1;
T_lev(1)=TL+5;
t(1)=1;
fs(1)=0;

while X(k)<=1
    if T_lev(j)<TL
        dT(k)=TL-T_lev(j);
        fm(k)=fs_eq(T_lev(j),TL);
        t_star(k)=t_star_eq(dT(k),C0_r,fm(k),n);
        X(k+1)=X(k)-dt*(n*(1-X(k))*log(abs(1-X(k))))/(t_star(k)*(log(abs(1-X(k)))/log(1-Xc))^(1/n));
        fs(k+1)=fs(k)+fm(k)*(X(k+1)-X(k));
        T_lev(j+1)=T_lev(j)-a+(c.dHf/c.pc)*(fs(k+1)-fs(k))/dt;
        k=k+1;
    else
        T_lev(j+1)=T_lev(j)-dt*a;
    end
    t(j+1)=t(j)+dt;
    j=j+1;
end
j=j-1;
k=k-1;
limit=T_lev(j)-10;
a_star=a*(c.ks/c.kl);
fm(k+1)=fs_eq(T(j+1),TL);
while T_lev(j)>limit
    dT(k)=TL-T_lev(j);
    dfm=abs(fm(k+1)-fm(k));
    T_lev(j+1)=T_lev(j)+a_star*((c.dHf/(c.pc*dT(k)))*dfm-1)^-1;
    t(j+1)=t(j)+dt;
    k=k+1;
    fm(k+1)=fs_eq(T(j+1),TL);
    j=j+1;
end

figure
plot(t,T_lev);




