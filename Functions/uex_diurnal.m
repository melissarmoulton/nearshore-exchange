function uex = uex_diurnal(ustartotal, Qmag, bslope, rhow)

% UEX_DIURNAL computes parameterized diurnal heating/cooling exchange
%
%  Syntax: uex = uex_diurnal()
%
%  Inputs:
%     ustartotal - total shear velocity (m/s)
%     Qmag - daytime-nighttime difference heating-cooling (W/m^2)
%     bslope - bottom slope
%     rhow - seawater density
%
%  Outputs:
%     uex - diurnal heating/cooling exchange velocity
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to as M2023)

% Constants
g = 9.81; % m/s^2
kappavc = 0.4;

% Thermal expansion coefficient
alphaT = 1.7*10^-4; %/degC
% Could place a seawater toolbox function here with temperature as input
% to get a more exact answer

%Specific heat of seawater
Cp = 4020; % (J/(kg K))
% Could use a temperature dependent function

Tday = 24*3600; % seconds in a day

uex = (2/3)*g*alphaT*Tday/(rhow*pi*Cp*kappavc*ustartotal)*bslope*Qmag;

end