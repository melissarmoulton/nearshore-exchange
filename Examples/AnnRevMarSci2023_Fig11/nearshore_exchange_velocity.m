% Code to generate Figure 11 of Moulton et al. (2023).
%
% The exchange velocity for each process is computed for two cases.
%
% Dependencies:
%   - uex_windx.m
%   - uex_windy.m
%   - uex_Stokes.m
%   - uex_brc.m
%   - uex_trc_profile.m, uex_trc.m
%   - KumarDiagram_IW.m, IWtracks.m, ddz.m
%   - uex_diurnal.m
%   - waveshoal.m
% 
% References:
%
% Moulton M, Suanda S, Garwood J, Kumar N, Fewings M, Pringle J. (2023)
%     Exchange of plankton, pollutants, and particles across the nearshore
%     region. Annual Review of Marine Science.
%   (referred to below as M2023)
% 
%% Add paths to needed functions

addpath('IW_scripts/')
addpath('../../Functions/')

%% Set depth profile and other parameters for both Cases
% These parameters are the same for both cases

x = 0:1:5000; % cross-shore coordinate (m)
bslope = 0.025; % bottom slope
h = bslope*x; % water depths (m)
gammabr = 0.73; % depth limited wave breaking parameter
theta0 = 0; % mean wave angle (degrees)

rhow = 1023; % water density (kg m3 /s)
kappavc = 0.4; % von Karmen's constant
f = 2*7.2921*10^-5*sind(35); % Coriolis freq at location of interest

h0 = 120; % approximate depth of 'offshore' wave height (m) 

Dr = 10/200; % Rip-density, 10-m wide rip every 200 m
hbar = 0.95; % depth on bar crest
hchan = 0.95+1; % depth in channel
thetabr = 0; % both cases shore-normal waves 
Cd = 3.3*10^(-3); % drag coefficient

%% Set additional parameters for two cases
% The rest of the code looks at two cases.
% The parameters can be modified for each
% Case 1: stratified
% Case 2: unstratified or weakly stratified

% Case 1

N_Case1 = 0.01; % 1/s, buoyance frequency
% N=0.01 corresponds to about T = 19 to 14 C over 15m depth

Hs_Case1 = 0.605; % significant wave height (offshore value)
T_Case1 = 11; % period (s)
sprd_Case1 = 8; % directional spread (degrees)

% Set wind stress
tau_Case1 = 0.03; % wind stress N/m2
% Here assume equal cross-shore and alongshore components
taux_Case1 = tau_Case1/sqrt(2); % wind stress N/m2
tauy_Case1 = tau_Case1/sqrt(2); % wind stress N/m2

% Compute ustar and deltas
ustar_Case1 = sqrt(tau_Case1/rhow);
deltas_Case1 = 1.5*ustar_Case1/sqrt(N_Case1*f); % Stratified, M2023 Eq 3

% Index of edge of surf zone based on depth-limited wave-breaking gamma
sz_ind_Case1 = find(h==(Hs_Case1/gammabr));

% Case 2

Hs_Case2 = 1.31; % significant wave height (offshore value)
T_Case2 = 7; % period (s)
sprd_Case2 = 20; % directional spread (degrees)

% Set wind stress
tau_Case2 = 0.15;
% Here assume equal cross-shore and alongshore components
taux_Case2 = tau_Case2/sqrt(2); % wind stress N/m2
tauy_Case2 = tau_Case2/sqrt(2); % wind stress N/m2

% Compute ustar and deltas
ustar_Case2 = sqrt(tau_Case2/rhow);
deltas_Case2 = kappavc*ustar_Case2/f; % Unstratified, M2023 Eq 7

% Index of edge of surf zone based on depth-limited wave-breaking gamma
sz_ind_Case2 = find(h==(Hs_Case2/gammabr));

%% Compute surface wave transformation:

wave_Case1 = waveshoal(T_Case1, h0, Hs_Case1, theta0, gammabr);
wave_Case2 = waveshoal(T_Case2, h0, Hs_Case2, theta0, gammabr);
% see waveshoal function for description of inputs/outputs

% Compute profile of wave height H with breaking
Hprofile_Case1=zeros(size(wave_Case1.h));
for ii = 1:(length(wave_Case1.h))
    h_ii = wave_Case1.h(ii);    
    Hs_ii = wave_Case1.H(ii);
    if h_ii<wave_Case1.breaking_depth
        Hs_ii = h_ii*gammabr;
    end
    Hprofile_Case1(ii) = Hs_ii;
end

% Compute profile of wave height H with breaking
Hprofile_Case2=zeros(size(wave_Case2.h));
for ii = 1:(length(wave_Case2.h))
    h_ii = wave_Case2.h(ii);    
    Hs_ii = wave_Case2.H(ii);
    if h_ii<wave_Case2.breaking_depth
        Hs_ii = h_ii*gammabr;
    end
    Hprofile_Case2(ii) = Hs_ii;
end

%% Wind-driven exchange

% Case 1 on the stratified shelf with non-overlapping boundary layers
% Intermediate stratification, using B = 0.9 in Eq 9a
uex_windx_Case1a = uex_windx(taux_Case1, deltas_Case1, rhow, f, 0.9);

% Case 1 in shallow water h<deltas
% Eq 9b
uex_windx_Case1b = uex_windx(taux_Case1, deltas_Case1, rhow, f, 1.2);

% Case 1 alongshore-wind-driven exchange for h>0.5*deltas
uex_windy_Case1 = uex_windy(tauy_Case1, deltas_Case1, rhow, f, 1);

% Case 2 - unstratified/ weakly stratified, h<deltas
uex_windx_Case2_B = uex_windx(taux_Case2, deltas_Case2, rhow, f, 1.2); 
% Case 2 - unstratified/ weakly stratified, h>deltas, weak stratification
uex_windx_Case2_C = uex_windx(taux_Case2, deltas_Case2, rhow, f, 1.6); 

% Case 2 alongshore-wind-driven exchange for h>0.5*deltas
uex_windy_Case2_B = uex_windy(tauy_Case2, deltas_Case2, rhow, f, 1);
% Case 2 alongshore-wind-driven exchange for h>2*deltas
uex_windy_Case2_C = uex_windy(tauy_Case2, deltas_Case2, rhow, f, 1); 

%% Stokes drift and undertow exchange - set params

omega_Case1 = 2*pi/T_Case1;
omega_Case2 = 2*pi/T_Case2;

k_Case1 = 2*pi./(wave_Case1.L); k_Case1(end)=k_Case1(end-1); % remove infinite value
k_Case2 = 2*pi./(wave_Case2.L); k_Case2(end)=k_Case2(end-1);

h_Case1 = wave_Case1.h;
h_Case2 = wave_Case2.h;

%% Stokes drift and undertow - Case 1

% Compute uex Stokes with Eq 14
uex_Stokes_Case1 = uex_Stokes(Hprofile_Case1, omega_Case1, k_Case1, 0); % 0 wave angle
% Locations where uex_Stokes may = 0 are estimated in plot below 

%% Stokes drift and undertow - Case 2

% Compute ex_Stokes with Eq 14
uex_Stokes_Case2 = uex_Stokes(Hprofile_Case2, omega_Case2, k_Case2, 0); % 0 wave angle
% Locations where uex_Stokes may = 0 are estimated in plot below 

%% uex Bathy rip currents

Hbr_Case1 = wave_Case1.breaking_height;
Hbr_Case2 = wave_Case2.breaking_height;

uex_brip_Case1 = uex_brc(Dr, Hbr_Case1, thetabr, hbar, hchan, gammabr, [Cd bslope]);
uex_brip_Case2 = uex_brc(Dr, Hbr_Case2, thetabr, hbar, hchan, gammabr, [Cd bslope]);

%% Transient rip current exchange

% Compute quantities for Case 1
Hsi_Case1 = Hs_Case1; % Significant wave height offshore
Li_Case1 = wave_Case1.Ldeep; % Offshore wavelength 
br_ind_Case1 = find(wave_Case1.h==wave_Case1.breaking_depth);
Lb_Case1 = wave_Case1.L(br_ind_Case1); % wavelength at breaking
dirspb_Case1 = sprd_Case1; % Directional spread at breaking
Lsz_Case1 = wave_Case1.breaking_depth/bslope; % Length of the surfzone
hb_Case1 = bslope*Lsz_Case1; % Water depth at breaking
Stp_Case1 = Hsi_Case1/Lb_Case1; % Wave steepness at breaking (using offshore wave height following Suanda)
%Stp_Case1 = wave_Case1.breaking_height/Lb_Case2; % Wave steepness at breaking (alternate using shoaled wave height, leads to about 10% difference in uex)
Stpi_Case1 = Hsi_Case1/Li_Case1; % Wave steepness offshore
Iri_Case1 = bslope/(Stpi_Case1^.5); % Irribarren number offshore

% Compute quantities for Case 2
Hsi_Case2 = Hs_Case2; % Significant wave height offshore
Li_Case2 = wave_Case2.Ldeep; % Offshore wavelength 
br_ind_Case2 = find(wave_Case2.h==wave_Case2.breaking_depth);
Lb_Case2 = wave_Case2.L(br_ind_Case2); % wavelength at breaking
dirspb_Case2 = sprd_Case2; % Directional spread at breaking
Lsz_Case2 = wave_Case2.breaking_depth/bslope; % Length of the surfzone
hb_Case2 = bslope*Lsz_Case2; % Water depth at breaking
Stp_Case2 = Hsi_Case2/Lb_Case2; % Wave steepness at breaking (using offshore wave height following Suanda)
%Stp_Case2 = wave_Case2.breaking_height/Lb_Case2; % Wave steepness at breaking (alternate using shoaled wave height, leads to about 10% difference in uex)
Stpi_Case2 = Hsi_Case2/Li_Case2; % Wave steepness offshore
Iri_Case2 = bslope/(Stpi_Case2^.5); % Irribarren number offshore

% Generate cross-shore profile of transient-rip current exchange 

[uex_trc_Case1, uex_trc_szedge_Case1, uex_trc_profile_Case1, x_Case1] = ...
    uex_trc_profile(sprd_Case1, Hbr_Case1, gammabr, ...
        Stp_Case1, Stpi_Case1, Iri_Case1, Lsz_Case1);

[uex_trc_Case2, uex_trc_szedge_Case2, uex_trc_profile_Case2, x_Case2] = ...
    uex_trc_profile(sprd_Case2, Hbr_Case2, gammabr, ...
        Stp_Case2, Stpi_Case2, Iri_Case2, Lsz_Case2);

%% Internal wave exchange:

% Script to compute internal wave exchange velocities in IW_scripts/
% from Garwood 2022. Uncomment the next 3 lines to run it.
% KumarDiagram_IW
% h_IW = KumarD_IW.H;
% uex_IW = KumarD_IW.u_ex;

% Alternatively, here are the values generated by that code:
h_IW = [20 30 40 50]; % Depths of internal wave estimates (m)
uex_IW = [.00062 .00282 .00327 .00579]; % IW exchange velocities (m/s)

%% Diurnal heating and cooling

Qmag_Case2 = 650; % W/m^2, peak, or half difference between max heating/cooling

% Write script for diurnal code
uex_diurnal_Case2 = uex_diurnal(ustar_Case2, Qmag_Case2, bslope, rhow);

%% Plot results for Case 1 and Case 2

fs=16; % Font size
lw = 1.5; % Line widths

xls1 = [0 50]; % x limits for Case 1
xls2 = [0 80]; % x limits for Case 2
yls1 = [0 0.07]; % y limits for Case 1
yls2 = [0 0.07]; % y limits for Case 2

fig1 = figure('Color','w');
set(fig1, 'Position', [1 1 1400 350])

s1 = subplot(1,3,1); hold on
set(s1, 'FontSize', fs, 'LineWidth', lw, 'Box', 'on'); hold on

% Stokes drift and undertow
tzone1 = 2/4; tzone2 = 5/4; % choose (arbitrary) fractional width of transition zones
hSt0 = 4*wave_Case1.breaking_depth;
iSt1 = find(min(abs(wave_Case1.h-hSt0*tzone1))==abs(wave_Case1.h-hSt0*tzone1));
iSt2 = find(min(abs(wave_Case1.h-hSt0*tzone2))==abs(wave_Case1.h-hSt0*tzone2));

xpatch = [flipud(wave_Case1.h(1:iSt1));(wave_Case1.h(1:iSt2));([wave_Case1.h(iSt1);wave_Case1.h(iSt2)])];
ypatch = [flipud(uex_Stokes_Case1(1:iSt1));(0*uex_Stokes_Case1(1:iSt2));([uex_Stokes_Case1(iSt1);0])];
ph0 = patch(xpatch,ypatch,[1 1 1]*.8);
set(ph0,'EdgeColor','none'); %,'FaceAlpha',0.6); % optional transparency

plot(wave_Case1.h(iSt1:end),uex_Stokes_Case1(iSt1:end),'-','LineWidth',3,'Color',[1 1 1]*.5);
plot([wave_Case1.h(iSt1) wave_Case1.h(iSt2)],[uex_Stokes_Case1(iSt1) 0],':','LineWidth',3,'Color',[1 1 1]/2)
plot(wave_Case1.h(1:iSt2),0*uex_Stokes_Case1(1:iSt2),'-','LineWidth',3,'Color',[1 1 1]*.5);
plot(wave_Case1.h(1:iSt1),uex_Stokes_Case1(1:iSt1),'-','LineWidth',3,'Color',[1 1 1]*.5);

% Rip currents
plot(-x_Case1*bslope, uex_trc_profile_Case1,'r-','LineWidth',3)
plot(wave_Case1.breaking_depth, uex_trc_szedge_Case1,'rs','MarkerSize',12,'LineWidth',3)
plot(wave_Case1.breaking_depth, uex_brip_Case1,'ro','MarkerSize',12,'LineWidth',3,'Color',[.7 0 0])

% Cross-shore wind
tzone1 = 3/5; tzone2 = 7/5; % choose (arbitrary) fractional width of transition zones
plot([0 wave_Case1.breaking_depth*2],uex_windx_Case1b*[1 1],'b:','LineWidth',3);
plot([wave_Case1.breaking_depth*2 deltas_Case1*tzone1],uex_windx_Case1b*[1 1],'b-','LineWidth',3);
plot([deltas_Case1*tzone1 deltas_Case1*tzone2],[uex_windx_Case1b uex_windx_Case1a],'b:','LineWidth',3);
plot([deltas_Case1*tzone2 120],uex_windx_Case1a*[1 1],'b-','LineWidth',3);

% Alongshore wind
tzone2 = 7/5; % choose (arbitrary) fractional width of transition zones
plot([0 deltas_Case1/2],0*[1 1],'c-','LineWidth',3,'Color',[0 1 1]*.9);
plot([deltas_Case1/2 deltas_Case1*tzone2],[0 uex_windy_Case1],'c:','LineWidth',3,'Color',[0 1 1]*.9);
plot([deltas_Case1*tzone2 120],uex_windy_Case1*[1 1],'c-','LineWidth',3,'Color',[0 1 1]*.9);

% Internal waves
plot(h_IW,uex_IW,'d-','LineWidth',3,'Color',[1 .5 0],'MarkerSize',10);

xlim(xls1)
ylim(yls1)
xlabel('water depth: $h$ $\mathrm{(m)}$','interpreter','latex')
ylabel({'cross-shore exchange velocity','$u_\mathrm{ex}$ $\mathrm{(m/s)}$'},'interpreter','latex');

s2 = subplot(1,3,[2 3]);
set(s2, 'FontSize', fs, 'LineWidth', lw, 'Box', 'on'); hold on

% Stokes drift and undertow
tzone1 = 2/4; tzone2 = 5/4; % choose (arbitrary) fractional width of transition zones
hSt0 = 4*wave_Case2.breaking_depth; % Stokes drift exchange may transition to 0 several surfzone widths from shore
iSt1 = find(min(abs(wave_Case2.h-hSt0*tzone1))==abs(wave_Case2.h-hSt0*tzone1));
iSt2 = find(min(abs(wave_Case2.h-hSt0*tzone2))==abs(wave_Case2.h-hSt0*tzone2));
iSt2=iSt2(1);

xpatch = [flipud(wave_Case2.h(1:iSt1));(wave_Case2.h(1:iSt2));([wave_Case2.h(iSt1);wave_Case2.h(iSt2)])];
ypatch = [flipud(uex_Stokes_Case2(1:iSt1));(0*uex_Stokes_Case2(1:iSt2));([uex_Stokes_Case2(iSt1);0])];
ph = patch(xpatch,ypatch,[1 1 1]*.8);
set(ph,'EdgeColor','none')% ,'FaceAlpha',0.6); % optional transparency

plot(wave_Case2.h(iSt1:end),uex_Stokes_Case2(iSt1:end),'-','LineWidth',3,'Color',[1 1 1]*.5);
plot([wave_Case2.h(iSt1) wave_Case2.h(iSt2)],[uex_Stokes_Case2(iSt1) 0],':','LineWidth',3,'Color',[1 1 1]/2)
plot(wave_Case2.h(1:iSt2),0*uex_Stokes_Case2(1:iSt2),'-','LineWidth',2,'Color',[1 1 1]*.5);
plot(wave_Case2.h(1:iSt1),uex_Stokes_Case2(1:iSt1),'-','LineWidth',3,'Color',[1 1 1]*.5);

% Diurnal heating and cooling
plot([0 7],uex_diurnal_Case2*[1 1],'LineWidth',5,'Color',[1 .85 0])
plot([7 10],0.025*[1 0],':','LineWidth',3,'Color',[1 .85 0])
plot([10 80],[0 0]+.0005,'LineWidth',4.5,'Color',[1 .85 0])

% Rip currents
plot(-x_Case2*bslope, uex_trc_profile_Case2,'r-','LineWidth',3)
plot(wave_Case2.breaking_depth, uex_trc_szedge_Case2,'rs','MarkerSize',12,'LineWidth',3)
plot(wave_Case2.breaking_depth, uex_brip_Case2,'ro','MarkerSize',12,'LineWidth',3,'Color',[.7 0 0])

% Cross-shore wind
tzone1 = 3/5; tzone2 = 6/5; %7/5; % choose (arbitrary) fractional width of transition zones
plot([0 wave_Case2.breaking_depth*2],uex_windx_Case2_B*[1 1],'b:','LineWidth',3);
plot([wave_Case2.breaking_depth*2 deltas_Case2*tzone1],uex_windx_Case2_B*[1 1],'b-','LineWidth',3); %deltas_Case2*tzone1
plot([deltas_Case2*tzone1 deltas_Case2*tzone2],[uex_windx_Case2_B uex_windx_Case2_C],'b:','LineWidth',3);
plot([deltas_Case2*tzone2 120],uex_windx_Case2_C*[1 1],'b-','LineWidth',3);

% Alongshore wind
tzone1 = 3/5; tzone2 = 7/5; % choose (arbitrary) fractional width of transition zones
plot([0 1/2*deltas_Case2*tzone1],0*[1 1],'c-','LineWidth',3,'Color',[0 1 1]*.9);
plot([1/2*deltas_Case2*tzone1 deltas_Case2*tzone1],[0 uex_windy_Case2_C],'c:','LineWidth',3,'Color',[0 1 1]*.9);
plot([deltas_Case2*tzone1 150],uex_windy_Case2_C*[1 1],'c-','LineWidth',3,'Color',[0 1 1]*.9);

set(gca,'XTick',0:10:120)
xlim(xls2)
ylim(yls2)
xlabel('water depth: $h$ $\mathrm{(m)}$','interpreter','latex')
ylabel('$u_\mathrm{ex}$ $\mathrm{(m/s)}$','interpreter','latex'); %/U_{Ek}')

set(s2,'Position',[.4108 .1323 .4942-.153 .7927])

%% Print figure

% addpath('/yourpath/export_fig'))
% export_fig('-pdf','AnnRevMarSci2023_Fig11')
% print(gcf,'-dpsc2','Figure11_Kumardiagrams')
% exportgraphics(gcf,'Figure11_Kumardiagrams.pdf','ContentType','vector','Resolution',1200)
% saveas(gcf, 'Figure11_Kumardiagrams.svg')
% saveas(gcf, 'Figure11_Kumardiagrams.pdf')
