% Code to calculate the exchange velocity due to internal waves for the 
% Kumar Diagram (Fig. 11) of Moulton et al. 2023.
% 
% Requires the function IWtracks.m described in Garwood et al. (2020) and 
% archived on Zenodo (doi.org/10.5281/zenodo.6784257), and Smyth et al. 
% (2011)'s ddz.m function. Both are included here for reference. To better
% visualize the internal wave velocity fields, it is recommended to replace
% the 'cool' colormap by one's favorite divergent colormap.
%
% References
% Garwood et al. 2020. Life in Internal Waves. Oceanography.
%
% Moulton et al. 2023. Exchange of plankton, pollutants, and particles 
%   across the nearshore region. Annual Review of Marine Science.
%
% Smyth, W.D., J.N. Moum and J.D. Nash. 2011. Narrowband, high-frequency 
% 	oscillations at the equator. Part II: Properties of shear instabilities.

%% Create a table with all of Colosi's information
site    = categorical({'S50'; 'LR50'; 'O50'; 'N50'; 'O40'; ...
        'S30'; 'LR30'; 'O30'; 'N30'; 'O20'; 'a50'; 'a30'});

nWave   = [23 22 21 24 17 18 14 14 15 12 22 15]';
eta     = [7.3 6.9 7.3 7.0 5.7 5.1 4.6 5.1 4.7 2.1 7.1 4.7]';  % m
twidth  = [7.6 8.2 8.2 8.3 7.6 6.5 6.6 6.8 7.1 5.6 8.1 6.7]' * 60;  % s
c       = [19 19 19 19 15 14 13 13 14 8.1 19 14]' ./ 100;  % m/s 
H       = [50 50 50 50 40 30 30 30 30 20 50 30]';
lambda  = twidth .* c;

colosi  = table(site, nWave, eta, twidth, c, H, lambda);

%% Parameters to generate particle tracks
z       = (-1:0.01:0)';  % vertical vector (normalized to H)
phi     = -sin(pi*z);  % wave mode

deltat  = 1;  % time steps for tracks
nsteps  = 8000;  % number of iterations
wshape  = 'sech2';  % wave shape

% Linear stratification
rho_bot = 1025;  % bottom density
rho_top = 1024;
mu      = (rho_bot/rho_top - 1)/(1);

rho0    = rho_bot * (1 - mu * z);

%% Calculate particle horizontal displacement for each wave
% Initialize matrix to store particle horizontal displacement
deltax_dk = NaN*twidth;

% Initialize figures for reference
figure(1)
clf
set(gcf, 'Paperunits', 'inches' )
set(gcf, 'PaperSize', [7.5 3.4])
set(gcf,'PaperPosition',[0 0 7.5 3.4])
set(gcf, 'Units', 'inches')

figure(2)
clf
t2 = tiledlayout(3,4);
set(gcf, 'Paperunits', 'inches' )
set(gcf, 'PaperSize', [12 10])
set(gcf,'PaperPosition',[0 0 12 10])
set(gcf, 'Units', 'inches')

% Loop through each wave
for iwave = 1:length(site)

    % Generate tracks for depth-keeping particles
    display(['Generating tracks for wave ' num2str(iwave) ' out of ' ...
        num2str(length(site)) '...'])
    [dk, ~, wave]    = IWtracks(z*H(iwave), H(iwave), phi, rho0, ...
        eta(iwave)/2, twidth(iwave), c(iwave), -H(iwave):0.1:0, deltat, ...
        nsteps, wshape);
    
    % Interpolate to higher resolution for better answer
    z_interp         = -H(iwave):0.001:0;
    dk_dx_interp     = interp1(-H(iwave):0.1:0, dk.deltax, z_interp);

    % Calculate mean delta x for dk particles with positive transport
%    deltax_dk(iwave) = nanmean(dk_dx_interp(dk_dx_interp >= 0));
    deltax_dk(iwave) = mean(dk_dx_interp(dk_dx_interp >= 0),'omitnan');

    % Plot the surface track for reference
    figure(1)
    plot(dk.xp(end,:) - dk.xp(end,1), 'linewidth', 1.5); hold on  

    % Plot the wave horizontal velocity field
    figure(2)
    nexttile(iwave)    
    pcolor(wave.xgrid, wave.zgrid, wave.u); hold on
        shading flat
    contour(wave.xgrid, wave.zgrid, wave.rho, 'linewidth', 1, 'color', 'k')
        colormap('cool')
        colorbar()

    title(colosi.site(iwave))

        caxis(0.4 * [-1 1])
end

% Fix look of figures
figure(1)
set(gca, 'linewidth', 1)
ylabel('Horizontal displacement (m)')
xlabel('Time (s)')

figure(2)
set(gca, 'linewidth', 1)
ylabel(t2, 'z (m)')
xlabel(t2, 'x (m)')
title(t2, 'Wave horizontal velocity (m/s)')

%% Plot the mean horizontal displacement for each sampling line
figure(3)
clf
set(gcf, 'Paperunits', 'inches' )
set(gcf, 'PaperSize', [4 3])
set(gcf,'PaperPosition',[0 0 4 3])
set(gcf, 'Units', 'inches')

% North line
plot(H(site == 'N30' | site == 'N50'), ...
    deltax_dk(site == 'N30' | site == 'N50') ...
    .* nWave(site == 'N30' | site == 'N50') / 1000, ...
    '.', 'markersize', 20); hold on

% O line
plot(H(site == 'O20' | site == 'O30' | site == 'O40' | site == 'O50'), ...
    deltax_dk(site == 'O20' | site == 'O30' | site == 'O40' | site == 'O50') ...
    .* nWave(site == 'O20' | site == 'O30' | site == 'O40' | site == 'O50') / 1000, ...
    '.', 'markersize', 20)

% LR line
plot(H(site == 'LR30' | site == 'LR50'), ...
    deltax_dk(site == 'LR30' | site == 'LR50') ...
    .* nWave(site == 'LR30' | site == 'LR50') / 1000, ...
    '.', 'markersize', 20)

% S line
plot(H(site == 'S30' | site == 'S50'), ...
    deltax_dk(site == 'S30' | site == 'S50') ...
    .* nWave(site == 'S30' | site == 'S50') / 1000, ...
    '.', 'markersize', 20)


    legend('N', 'O', 'LR', 'S', 'location', 'northwest')

    xlim([0 55])

    xlabel('Isobath (m)')
    ylabel('Daily internal wave transport (km)')

    set(gca, 'linewidth', 1)

%% Save mean value at each isobath for Kumar diagram
KumarD_IW.H    = unique(H);
KumarD_IW.u_ex = NaN*KumarD_IW.H;

for iH = 1:length(KumarD_IW.H)
    iwaves = H == KumarD_IW.H(iH);
    KumarD_IW.u_ex(iH) = mean(deltax_dk(iwaves) .* nWave(iwaves)) ...
        / (24*60*60);
end

save('KumarDiagram_IW.mat', 'KumarD_IW');