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

%-------------------Lever rule derivated fraction----------------------
fs_eq=@(T,TL) (1/(1-c.k))*((TL-T)/(c.Tf-T));
%----------------------------------------------------------------------

%c)
dT=0.5; 
T(1)=TL(1);
for i=1:length(TL)
    fs(1,1)=0;
    fs(1,2)=0;
    j=1;
    while (fs(j,1)<1) || (fs(j,2)<1)
        T(j+1)=T(j)-dT;
        fs(j+1,1)=fs_eq(T(j),TL(1));
        fs(j+1,2)=fs_eq(T(j),TL(2));
        j=j+1;
    end
end

figure
for i=1:length(TL)
    plot(T-273,fs(:,i))
    hold on
end
axis([0 c.Tf 0 1])
title('Eq. lever rule solid fraction')
legend('1wt%Si','8wt%Si')
grid
xlabel('T(°C)')
ylabel('Solid fraction, f');

%fs_star=@(T,TL) 1-((1/(c.Tf-TL))^())*()^()


%testtest