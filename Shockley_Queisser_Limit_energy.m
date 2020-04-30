function [Jsc_SQ_m,Jo_SQ_m,Voc_SQ,FF_SQ,eff_SQ]=Shockley_Queisser_Limit_energy(E_bandgap,E)
%E_bandgap=1.61; %bandgap energy (eV)
%% Konstanten
q=1.602176462e-19;                      % [As], elementary charge
VT=25.8e-3;                             % [V], 25.8mV thermal voltage at 300K
h=6.62606876e-34;                       % [Js], Planck's constant
%h=4.135667662e-15;                     % [eVs], Planck's constant
c=299792458;                            % [m/s], speed of light c_0
vor=((h*c)/q)/(1e-9);                   % prefactor for converting between energy and wavelength in (eVnm)
P_sun=0.1;                              % (W/cm^2)Power density of the sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E=linspace(0.1,6,1000)';                % [eV], general energy axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% absorptance (stepfunction)
a_SQ=E;
for ii=1:length(E);
if a_SQ(ii,:)<E_bandgap
a_SQ(ii,:)=0;
else
a_SQ(ii,:)=1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BLACK BODY SPECTRUM BB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
BB=(1e-4*2*pi*q^3*(E).^2)./(h^3*c^2*(exp(E./VT)-1)); %black body spectrum BB (s^-1cm^-2(eV)^-1)(in photon flux!)
%% SOLAR SPECTRUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('AM15G.dat');                        % solar spectrum AM1.5 global
AM15g=interp1(AM15G(:,1), AM15G(:,2), E); % interpolation to the general energy axis (eV) 
AM15g(isnan(AM15g))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jo_SQ=trapz(E,q*a_SQ.*BB);            %saturation current density J0,SQ (A/cm^2)
Jo_SQ_m=1000*Jo_SQ;                   %saturation current density J0,SQ (mA/cm^2)
Jsc_SQ=trapz(E,q*a_SQ.*AM15g);        %short-circuit current denstiy Jsc,SQ (A/cm^2)
Jsc_SQ_m=Jsc_SQ*1000;                 %short-circuit current denstiy Jsc,SQ (mA/cm^2)
Voc_SQ=(VT)*log((Jsc_SQ./Jo_SQ)+1);   %open-circuit voltage Voc,SQ (V)
voc=Voc_SQ/VT;  
FF_SQ=(voc-log(voc+0.72))./(voc+1);                          %fill factor FF,SQ
eff_SQ=100*(Jsc_SQ.*Voc_SQ.*FF_SQ)./(P_sun);                 %efficieny eta,SQ in (%)
