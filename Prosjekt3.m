%Prosjekt 3 - Eirik Andre Nordb�, Tobias Mohn Werner
%



clear all
clc
close all
loop=1;
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

while loop ==1
clear n t X T Xc fs TLplot TEutPlot fm fs T 
%clearing all saved variables
%---------------------Assumptions for comparison-----------------------
C0_1=1; %[wt%Si]
C0_2=8; %[wt%Si]

TL(1)=c.Tf -c.m*C0_1;
TL(2)=c.Tf -c.m*C0_2;
%----------------------------------------------------------------------

%--------------Lever rule equilibrium fraction solid ------------------
fs_eq_lever=@(T,TL) (1/(1-c.k))*((TL-T)/(c.Tf-T));
dfs_eq=@(T,TL) (1/(c.k-1))*((TL-c.Tf)/(c.Tf-T)^2);
%----------------------------------------------------------------------

%-----------------------Lever rule Eutectic fraction-------------------
fe_eq=@(TL) 1-(1/(1-c.k))*((TL-c.Te)/(c.Tf-c.Te));
%----------------------------------------------------------------------

%--------------Scheil non-equilibrium fraction solid ------------------
fs_star_eq=@(T,TL) 1-((c.Tf-T)/(c.Tf-TL))^(1/(c.k-1));

%-----------------------Equilibrium fraction scheil-----------------
fe_eq_scheil=@(TL) ((c.Tf-c.Te)/(c.Tf-TL))^(1/(c.k-1));
%-------------------------------------------------------------------

o=menu('Velg oppgave:', 'Scheil vs Equilibrium','X vs t/t*','X vs dX/dt','Heat flow model', 'Avslutt program');
%c)
switch o
    case 1
        
clear t
        
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
    fs(j+1,1)=fs_eq_lever(T(j+1,1),TL(1));
    dfs(j+1,1)=dfs_eq(T(j+1,1),TL(1));
    fs(j+1,2)=fs_eq_lever(T(j+1,2),TL(2));
    dfs(j+1,2)=dfs_eq(T(j+1,2),TL(2));
    j=j+1;
end



figure(1)
subplot(2,1,1)

for i=1:length(TL)
    plot(T(:,i)-273,fs(:,i))
    hold on
end


grid;
xlabel('T(�C)')

ylabel('Solid fraction, f_s');

subplot(2,1,2)
for i=1:length(TL)
    plot(T(:,i)-273,dfs(:,i))
    hold on
end
axis([0 c.Tf -inf inf]);
set(gca,'xdir','reverse');
title('Solid fraction: Equilibrium vs Scheil')
grid;
xlabel('T(�C)')
ylabel('Change in solid fraction, df/dT');


% oppgave 1d


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


%Creating the perfect plot------------------------------------------
axis([0 c.Tf 0 1]);
set(gca,'xdir','reverse');
title('Solid fraction: Equilibrium vs Scheil')
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

    case 2
t_star=1; 
%----------------------Johnson-Mehl-Avrami equation------------------------
X_JMA=@(Xc,t,n) 1-(1-Xc)^((t/t_star)^n);

clear t, 

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

figure(2)
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
    case 3
%dX/dt vs X
%--------------------------------------------------------------------------
dXdt_eq = @(t_star,Xc,n,X) -(n.*(1-X).*log(1-X))./(t_star.*(log(1-X)./log(1-Xc)).^(1/n));
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

figure (3)
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
    case 4
    clear n t X TL Xc fs T TLplot TEutPlot fm fs dfs_dt
    %Asking user for parameters------------
    n=1;
    eq_or_non_eq=1;
    prompt = {'Enter dimention size, n=[1,2,3]:','Scheil=1, Equilibrium=2'};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {num2str(n),num2str(eq_or_non_eq)};
    answer=inputdlg(prompt,dlg_title,num_lines,defaultans);
    if str2num(answer{1})>3 || str2num(answer{1})<1
        display('Not a valid number of dimensions!')
        return
    end
    n = str2num(answer{1});
    eq_or_non_eq=str2num(answer{2});

    if eq_or_non_eq==2
        fs_eq=fs_eq_lever;
        f_eut=fe_eq;
    else
        fs_eq=fs_star_eq;
        f_eut=fe_eq_scheil;
    end
    %---------------------------------------
loop2=1;
while loop2==1

    %-----------------------------Input parameters-----------------------------
    t_r=6;               %s
    C0_r=4;              %wt%Si
    TL=c.Tf -c.m*C0_r;%K
    T_n=TL-2;        %K
    Xc=0.05;      
    a=1;                 %K/s
    Neq=1;               %Neq=Nr/N
    %--------------------------------------------------------------------------
    
    o2=menu('Velg deloppgave:', 'Full heat model','Comparison with dfs/dt','Change parameters','Tilbake');

    %Finding value for fm_r
    dT_r=TL-T_n;
    fm_r=fs_eq(T_n,TL);

    %----------------------------Equation for t*-------------------------------
    t_star_eq=@(dT,C0,fm,n) t_r*((dT_r/dT)^2)*(C0/C0_r)*((Neq)^(1/n))*(((fm/fm_r))^(1/n));
    %--------------------------------------------------------------------------

    %Starting value-----
    t_star_i=t_star_eq(dT_r,C0_r,fm_r,n);
    dt=0.01;
    t(1)=dt;
    dX=-dt*((1-Xc)^((t(1)/t_star_i)^n))*log(1-Xc)*((t(1)/t_star_i)^n)*(n/t(1));
    X(1)=dX;
    %-------------------
    j=1;
    k=1;
    T(1)=TL+5;
    fs(1)=0;
    fm(1)=fs_eq(T_n,TL);
    TLplot(1)=TL;
    TEutPlot(1)=c.Te;
    index=0;
    while X(k)<=0.99 %Transient part and solid growth (nucleation)
        if T(j)<=T_n
            dT(k)=TL-T(j);
            t_star(k)=t_star_eq(dT(k),C0_r,fm(k),n);
            dX_dt=-(n*(1-X(k))*log(abs(1-X(k))))/(t_star(k)*(log(abs(1-X(k)))/log(1-Xc))^(1/n));
            fs(k+1)=fs(k)+fm(k)*(dX_dt)*dt;
            dfs_dt(k)=(fs(k+1)-fs(k))/dt;
            T(j+1)=T(j)+dt*(-a+(c.dHf/c.pc)*(fs(k+1)-fs(k))/dt);
            fm(k+1)=fs_eq(T(j+1),TL);
            X(k+1)=fs(k+1)/fm(k+1);
            k=k+1;
        else
            T(j+1)=T(j)-dt*a;
            index=index+1;
        end
        TLplot(j+1)=TL;
        TEutPlot(j+1)=c.Te;
        t(j+1)=t(j)+dt;
        j=j+1;
    end
    a_star=a*(c.ks/c.kl);
    grid=3000;
    smooth_a=linspace(a,a_star,grid);
    l=1;
    dt=0.001;
    while T(j)>c.Te %Steady state growth
        if l<grid+1 %for smooth transition between a and a*
           a_star_adj=smooth_a(l);
           l=l+1;
        end
        dT(k)=TL-T(j); 
        dfm_dT=(fm(k)-fm(k-1))/((TL-dT(k))-(TL-dT(k-1)));
        T(j+1)=T(j)+a_star_adj*dt*((c.dHf/(c.pc)*dfm_dT)-1)^(-1);
        t(j+1)=t(j)+dt;
        fm(k+1)=fs_eq(T(j+1),TL);
        dfs_dt(k)=(fm(k+1)-fm(k))/dt;

        TLplot(j+1)=TL;
        TEutPlot(j+1)=c.Te;
        k=k+1;
        j=j+1;
    end

    %Isothermal eutectic growth
    t1=t(j-1);
    t2=t1+(c.dHf/(a_star*c.pc))*f_eut(TL);
    while t(j)<t2
        T(j+1)=c.Te;
        t(j+1)=t(j)+dt;
        TLplot(j+1)=TL;
        TEutPlot(j+1)=c.Te;
        j=j+1;
    end

        switch o2
            case 1 %Full heat model plot
                figure
                plot(t,T,t,TLplot,'--y',t,TEutPlot, '--y');
                
            case 2 %Confirming the model (dfs/dt)
                figure
                subplot(2,1,1)
                plot(t(index:index+k-2),T(index:index+k-2),t(index:index+k-2),TLplot(index:index+k-2),'--y')
                grid;
                ylabel('T (K)')
                subplot(2,1,2)
                plot(t(index:index+k-2),dfs_dt);
                grid;
                ylabel('dfs/dt')
                xlabel('t (s)')
                
            case 3 %Change parameters
                clear n t X TL Xc fs T TLplot TEutPlot fm fs dfs_dt
                %Asking user for parameters------------
                n=1;
                eq_or_non_eq=1;
                prompt = {'Enter dimention size, n=[1,2,3]:','Scheil=1, Equilibrium=2'};
                dlg_title = 'Input';
                num_lines = 1;
                defaultans = {num2str(n),num2str(eq_or_non_eq)};
                answer=inputdlg(prompt,dlg_title,num_lines,defaultans);
                if str2num(answer{1})>3 || str2num(answer{1})<1
                    display('Not a valid number of dimensions!')
                    return
                end
                n = str2num(answer{1});
                eq_or_non_eq=str2num(answer{2});

                if eq_or_non_eq==2
                    fs_eq=fs_eq_lever;
                    f_eut=fe_eq;
                else
                    fs_eq=fs_star_eq;
                    f_eut=fe_eq_scheil;
                end
                %---------------------------------------
            case 4
                loop2=0;
        end

end %while loop for heat flow model menu


    case 5

        loop=0;
        close all
end     %menyswitch
end     %restartloop




