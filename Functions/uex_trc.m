function uex = uex_trc(sprd, Hbr, gammabr, Stp)

% UEX_TRC computes parameterized maximum transient rip-current exchange velocity
%
%  Syntax: uex = uex_trc(sigmabr, Hbr, gammabr, Stp)
%
%  Inputs:
%     sprd - directional spread at breaking (degrees)
%     Hbr - Wave height (m) at breaking
%     gammabr - depth-limited breaking gamma
%     Stp - Wave steepness at breaking
%     (note: Hbr/gammabr = hbr, water-depth at breaking)
%
%  Outputs:
%     uex - Transient rip-current exchange velocity maximum
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to below as M2023)

% Max Uex:
uex = (5e-4)*sprd*((70*Stp)+1)*sqrt(9.81*(Hbr/gammabr));

end