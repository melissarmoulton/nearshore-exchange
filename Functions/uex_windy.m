function uex = uex_windy(tauy, deltas, rhow, f, coeff)

% UEX_WINDX computes parameterized cross-shore wind-driven exchange
%
%  Syntax: uex = uex_windy()
%
%  Inputs:
%     tauy - magnitude of alongshore component of wind stress (kg m /s^2)
%     deltas - surface boundary layer thickness (m)
%     rhow - water density (kg/m^3)
%     f - Coriolis frequency (1/s)
%     coeff - coefficient corresponding to regime from Fig 5, M2023
%           (alternate version of this code could be written to use 
%            list of depth values h and compute the regimes)
%
%  Outputs:
%     uex - alongshore-wind-driven exchange velocity
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to as M2023)

uex = coeff*tauy/(rhow*f*deltas);

end