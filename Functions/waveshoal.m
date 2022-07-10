function wave = waveshoal(T, h0, H0, theta0, gamma)

% Inputs:
% T - wave period (s)
% h0 - offshore water depth (m)
% H0 - significant wave height (m) at offshore water depth
% theta0 - wave direction at offshore location (degrees)
% gamma - wave breaking gamma threshold

% Outputs:
% wave.h = h; - water depths (m)
% wave.L = L; - wave lengths (m)
% wave.Ldeep = Ldeep; - deep-water wavelength (m)
% wave.H = H; - wave heights (m)
% wave.c = c; - wave phase speeds (m/s)
% wave.cg = cg; - wave group velocities (m/s)
% wave.theta = theta; - wave directions (degrees)
% wave.breaking_depth = breaking_depth; - depth at breaking (m)
% wave.breaking_height = breaking_height; - wave height at breaking (m)
% wave.breaking_angle = breaking_angle; - wave angle at breaking (theta)

% Constants
g=9.8; % m/s^2

% Find a water depth near but offshore of breaking,
% shoal waves on finer depth grid onshore of this depth
hBr_approx = round(10*H0/gamma*2)/10;
h = [h0:-.2:hBr_approx (hBr_approx-0.02):-0.02:0]';

% Calculate wavelengths in deep water at depths h

% Deep water wavelength:
Ldeep = g*T.^2/(2*pi); 

% Wavelength, Ole Madsen approx (nice because it doesn't require iteration):
L = Ldeep.*(1-exp(-(2.*pi.*h./Ldeep).^1.25)).^0.4;

% Calculate group and phase speeds at depths h
c = L./T; % Phase speed
k = 2.*pi./L; % Wavenumber
cg = (L./(2.*T)).*(1+2.*(k).*h./sinh(2.*(k).*h)); % Group velocity

% Calculate group and phase speeds at depth h0
c0 = c(1); % Phase speed at depth h0
cg0 = cg(1); % Phase speed at depth h0

% Compute wave height and angle at depths h
theta = asind(c./c0.*sind(theta0));
H = H0*sqrt(cg0./cg).*sqrt(cosd(abs(theta0))./cosd(abs(theta)));

% Calculate breaking variables
breaking_index = abs(H./h-gamma)==min(abs(H./h-gamma));
breaking_depth = h(breaking_index);
breaking_height = H(breaking_index);
breaking_angle = theta(breaking_index);

% Store variables
wave.h = h;
wave.L = L;
wave.Ldeep = Ldeep;
wave.H = H;
wave.c = c;
wave.cg = cg;
wave.theta = theta;
wave.breaking_depth = breaking_depth;
wave.breaking_height = breaking_height;
wave.breaking_angle = breaking_angle;

end