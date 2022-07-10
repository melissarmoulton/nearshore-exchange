function [uex, uex_szedge, uex_profile, x] = uex_trc_profile(sprd, Hbr, gammabr, Stp, Stpi, Iri, Lsz)

% UEX_TRC computes parameterized transient rip-current exchange velocity
%  maximum, surfzone-edge value, and profile versus cross-shore coord.
%
%  Syntax: uex = uex_trc(sprd, Hbr, gammabr, Stp, Lsz)
%
%  Inputs:
%     sprd - directional spread at breaking (degrees)
%     Hbr - Wave height (m) at breaking
%     gammabr - depth-limited breaking gamma
%       (Note: Hbr/gammabr = hbr, water-depth at breaking)
%     Stp - Wave steepness at breaking
%     Stpi - Offshore wave steepness
%     Iri - Offshore Irribarren number 
%     Lsz - Width of surf zone (m)
%
%  Outputs:
%     uex - Maximum transient rip-current exchange velocity (m/s)
%     uex_szedge - Transient rip-current exchange velocity at edge of surf zone (m/s)
%           Note: uex_szedge is more representative of exchange by
%                 transient rips at surfzone edge. The max uex may be
%                 more representative of surfzone eddy mixing.
%     uex_profile - Cross-shore profile of exchange (m/s)
%     x - cross shore coordinate (m)
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to below as M2023)

% Max Uex:
uex = uex_trc(sprd, Hbr, gammabr, Stp);

% Compute parameters for profile
Xo = -(2e-3*(Stp^-1) + 0.39)*Lsz;
Ld = 1.59*Iri*(sprd^-.25)*(22.16*(Stpi^0.5) + 1)*Lsz;

% Compute profile
x = -5000:.5:0; % cross-shore coordinate (m)
uex_profile = uex*exp(-((x-Xo).^2)./(Ld^2)); % uex_profile (m/s)

% Compute value at surfzone edge
br_ind = find(min(abs(abs(x)-Lsz))==(abs(abs(x)-Lsz)));

uex_szedge = uex_profile(br_ind);

end