function uex = uex_Stokes(H, omega, k, theta)

% UEX_DIURNAL computes parameterized Stokes-drift and undertow exchange
%
%  Syntax: uex = uex_Stokes(H, omega, k, theta)
%
%  Inputs:
%     H - wave height profile (m)
%     omega - angular wave frequency (rad/s)
%     k - wavenumber (1/m)
%     theta - wave angle profile (deg)
%
%  Outputs:
%     uex - Stokes drift and undertow exchange velocity
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to as M2023)

% Constants
g = 9.81; % m/s^2

uex = 0.61*(1/16)*(H).^2*omega.*k.*cosd(theta);

end