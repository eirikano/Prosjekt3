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
c.ks = 95;           %[W/K*m]    Thermal conductivity (solid)
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
%Creating the perfect plot------------------------------------------
axis([0 c.Tf 0 1]);
set(gca,'xdir','reverse');
title('Solid fraction: Lever rule vs Scheil')
fe1=fe_eq(TL(1));
fe2=fe_eq(TL(2));
if fe1<0
    fe1=0;
elseif fe2<0
    fe2=0;
end
strstart1=['Lever rule, ',num2str(C0_1)];
strend1=['wt%Si f_e = ', num2str(fe1,'%.3f')];
strstart2=['Lever rule, ',num2str(C0_2)];
strend2=['wt%Si f_e = ', num2str(fe2,'%.3f')];
str1=strcat(strstart1,strend1);
str2=strcat(strstart2,strend2);
legend(str1,str2,'Location','southeast')
%---------------------------------------------------------------------
grid
xlabel('T(°C)')
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
fs_star_eq=@(T,TL) 1-((c.Tf-T)/(c.Tf-TL))^(1/(c.k-1))
for i=1:length(TL)
    fs_star(1,2)=0;
    fs_star(1,1)=0;
    j=1;
    while (fs_star(j,1)<1) || (fs_star(j,2)<1)
        T(j+1)=T(j)-dT;
        fs_star(j+1,1)=fs_star_eq(T(j),TL(1))
        fs_star(j+1,2)=fs_star_eq(T(j),TL(2));
        j=j+1;
    end
end
