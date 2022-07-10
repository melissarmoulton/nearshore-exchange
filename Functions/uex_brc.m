function uex = uex_brc(Dr, Hbr, thetabr, hbar, hchan, gammabr, vAlong)

% UEX_BRC computes parameterized bathymetric rip-current exchange velocity
%
%  Syntax: uex = uex_brc(Dr,Urip)
%
%  Inputs:
%     Dr - Rip-current density
%     Hbr - Wave height (m) at breaking
%     thetabr - Wave direction (deg. from shore-normal) at breaking
%     hbar - Water depth on the bar-crest (m)
%     hchan - Water depth in the channel (m)
%     gammabr - Wave breaking gamma value, typically 0.3<gammaBr<0.9
%     vAlong - Breaking-wave-driven alongshore current speed (if 1 value) OR
%              use vAlong = [dragC bSlope] (2 values) to compute with
%              quadratic drag law:
%                 dragC - Quadratic drag coefficient, ~3E-3
%                 bSlope - Approximate beach slope
%
%  Outputs:
%     uex - Bathymetric rip-current exchange velocity
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%     (referred to below as M2023)

%% Set constants

g = 9.8; % m/s^2
vAdjust = 0.5; % Option to adjust sensitivity to alongshore flow, 0 to 1
%               ('beta' in Moulton et al., 2016)

%% Determine wave-breaking regime

R = ones(size(Hbr)); % Initialize with bar-break value, 1
R(Hbr<gammabr*hbar) = 0; % Set to 0 if shore-break
R(Hbr>gammabr*hchan) = 2; % Set to 2 if saturated

%% Compute sea-level tilt

% Initialize with bar-break value,
% If waves break on the bar but not in the channel:
etay = 1/16*gammabr^2.*(cosd(abs(thetabr)).^2+1/2).*...
    (Hbr/gammabr-hbar);

% Set to zero if shore-break,
% No sea-level difference if waves break onshore of bathy variability:
etay(R==0) = 0;

% Set to saturated value if saturated,
% If waves break in channel and on bar (offshore of bathy variability),
% the sea level difference are limited by the water depths
etay(R==2) = 1/16*gammabr^2.*(cosd(abs(thetabr(R==2))).^2+1/2).*...
          (hchan(R==2)-hbar(R==2)); 

%% Compute alongshore flow

% Assign speed or drag coefficient and beach slope from vAlong

if numel(vAlong)==1
    V = vAlong;
else
   
dragC = vAlong(:,1);
bSlope = vAlong(:,2);
    
% Initialize with maximum alongshore flow for waves breaking on the bar,
% Alongshore flow speed near the bar crest is limited by depth on bar:
V = sqrt(5/32)*dragC.^(-1/2)*bSlope.^(1/2)*(gammabr)^(1/2)*...
    (g*gammabr*hbar).^(1/2).*...
    sind(abs(thetabr)).^(1/2).*sign(sind(thetabr));

% For shore-break cases, fill in value based on wave height
V(R==0) = sqrt(5/32)*dragC.^(-1/2)*bSlope.^(1/2)*(gammabr)^(1/2)*...
    (g*Hbr(R==0)).^(1/2).*...
    sind(abs(thetabr(R==0))).^(1/2).*sign(sind(thetabr(R==0)));

% Option to use value based on wave height for all regimes (uncomment if desired)
% V = sqrt(5/32)*dragC.^(-1/2)*bSlope.^(1/2)*(gammabr)^(1/2)*...
%    (g*Hbr).^(1/2).*...
%    sind(abs(thetabr(R==0))).^(1/2).*sign(sind(thetabr(R==0)));
     
end

%% Compute alongshore flow suppression factor F

Fv = ((1+vAdjust^2.*abs(V).^2./(2*g*etay)).^(1/2)-vAdjust*abs(V)./sqrt(2*g*etay));
Fv(etay==0) = NaN; % set to NaN if etay = 0

%% Compute U

U = sqrt(2*g*etay).*Fv;
U(isnan(Fv)) = sqrt(2*g*etay(isnan(Fv))); % Cases with etay = 0

%% Compute uex

uex = U.*Dr;

end