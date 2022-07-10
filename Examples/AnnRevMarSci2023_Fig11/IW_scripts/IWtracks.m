function [dk, p, wave] = IWtracks(z, H, phi, rho0, amp, period, cp, z_org, deltat, nsteps, wshape)
% [dk, p, wave] = IWtracks(z, H, phi, rho0, amp, period, cp, z_org, deltat, nsteps, wshape)
%
% Uses Matlab's ODE solver to generate tracks for depth-keeping (dk) and 
% passive (p) particles in a linear or weakly nonlinear internal wave. 
% Also calculates the total horizontal displacement of these particles and
% their residence time in the wave.
%
% In an oscillatory wave, the residence time of a particle is defined as
% the time required for the particle to encounter two adjacent, identical 
% wave phases (e.g., from crest to crest). In a solitary wave, the
% residence time of a particle is defined as the time required for the
% particle to travel from one wavelength ahead of the wave trough to one 
% wavelength behind the wave trough.
%
% Note that deltat sets the time interval between particle position
% estimates but is distinct from the ODE internal step size, which is set
% to be at most 1/100th of the wave period, as the default step size was 
% too coarse.
%
% Input:
% z        = vertical vector (positive up, 0 at surface)
% H        = water column height [m]
% phi      = wave vertical function
% rho0     = background density profile [kg m-3]
% amp      = wave amplitude, i.e., wave height/2 or eta_max/2 [m]
% period   = wave period [s]  
% cp       = wave propagation speed [m/s]
% z_org    = vertical position [m] of organisms
% deltat   = time interval [s]
% nsteps   = number of time steps
% wshape   = shape of wave ('cos', 'cos2', 'sech2')
%
% Output - dk, and p:
% .xp      = time series of horizontal position [m]
% .zp      = ... of vertical position [m]
% .deltat  = time step for time series [s]
% .time    = time associated with time series [s]
% .restime = residence time in the wave [s]
% .deltax  = net horizontal displacement [m]
%
% Output - wave:
% .x       = x vector [m]
% .z       = z vector [m]
% .cp      = wave propagation speed [m/s]
% .H       = water height [m]
% .Amp     = wave amplitude, i.e., wave height/2 [m]
% .period  = wave period [s]
%
% Dependencies:
% - ddz.m from Smyth et al. (2011) (included)
%   Smyth, W.D., J.N. Moum and J.D. Nash. 2011. Narrowband, high-frequency 
%        oscillations at the equator. Part II: Properties of shear 
%        instabilities. J. Phys. Oceanogr. 41(3): 412-428.
%
% Author:
% Jessica C. Garwood - for application see:
%   Garwood, J.C., R.C. Musgrave and A.J. Lucas. 2020. Life in internal 
%        waves. Oceanography. 33(3): 38-49.
%
% Archived on Zenodo at doi.org/10.5281/zenodo.6784257

%% Function parameters
omega0 = 2 * pi / period;  % wave frequency [s-1]
lambda = period * cp;  % horizontal wavelength [m]
k      = 2*pi / lambda;  % Horizontal wavenumber [m-1]
time   = 0:deltat:(nsteps-1)*deltat;  % time vector [s]

phi    = reshape(phi, numel(phi), 1);  % make sure phi is a column vector

%% Grids - Wave
x      = -3 * lambda : lambda/1000 : 3 * lambda;
z      = reshape(z, numel(z), 1);  % make sure z is a column vector

%% Starting x-position - Organisms
if strcmp(wshape, 'sech2') 
    x0  = lambda;
else
    x0  = lambda/2;
end

%% Define wave functions
phi_interp  = @(zp)(interp1(z, phi, zp));
dPhidz      = @(zp)(interp1(z, ddz(z) * phi, zp));

% Isopycnal displacement
if strcmp(wshape, 'cos')  % linear wave
    eta = @(xp,zp,tp) (amp * phi_interp(zp) .* ...
            cos(k * xp - omega0 * tp));  
    % d(eta)/d(z)    
    detadz = @(xp,zp,tp) (amp * dPhidz(zp) .* ...
            cos(k * xp - omega0 * tp));
    % d(eta)/d(x)
    detadx = @(xp,zp,tp) (-k * amp * phi_interp(zp) .* ...
            sin(k * xp - omega0 * tp));
                     
elseif strcmp(wshape, 'cos2')  % weakly nonlinear, oscillatory    
    eta = @(xp,zp,tp) (-amp*2 * phi_interp(zp) .* ...
            cos(k/2 * xp - omega0/2 * tp).^2);        
    % d(eta)/d(z)
    detadz = @(xp,zp,tp) (-amp*2 * dPhidz(zp) .* ...
            cos(k/2 * xp - omega0/2 * tp).^2);
    % d(eta)/d(x)
    detadx = @(xp,zp,tp) (k*amp * phi_interp(zp) .* ...
        sin(k * xp - omega0 * tp));
        
elseif strcmp(wshape, 'sech2')  % weakly nonlinear, solitary    
    % Isopycnal displacement
    eta     = @(xp,zp,tp) (-amp*2 * phi_interp(zp) .* ...
                sech(k/2 * xp - omega0/2 * tp).^2);        
    % d(eta)/d(z)
    detadz  = @(xp,zp,tp) (-amp*2 * dPhidz(zp) .* ...
                sech(k/2 * xp - omega0/2 * tp).^2);
    % d(eta)/d(x)
    detadx  = @(xp,zp,tp) (k*amp*2 * phi_interp(zp) .* ...
                sech(k/2 * xp - omega0/2 * tp).^2 .* ...
                tanh(k/2 * xp - omega0/2 * tp));
end

% Wave velocities
wave_u  = @(xp,zp,tp) (cp * detadz(xp,zp,tp));
wave_w  = @(xp,zp,tp) (-cp * detadx(xp,zp,tp));

% Deformation of rho
rho     = @(xp,zp,tp) (interp1(z, rho0, zp - eta(xp, zp, tp)));

%% Initialize loop conditions
% Number of particles
npart = length(z_org);

% Initialize position, residence time, and net displacement matrices
p_xpos     = NaN(npart, nsteps);
p_zpos     = NaN(npart, nsteps);
p_restime  = NaN(npart, 1);
p_deltax   = NaN(npart, 1);

dk_xpos    = NaN(npart, nsteps);
dk_zpos    = NaN(npart, nsteps);
dk_restime = NaN(npart, 1);
dk_deltax  = NaN(npart, 1);

%% Advect particles using ODE45 solver
% Time points at which positions are evaluated
tspan    = [0:deltat:nsteps];

% Pass event function to stop integration and set ODE internal step size
xoverFcn = @(tp, pp) crestEventsFcn(tp,pp,x0,cp);
opts     = odeset('MaxStep', period/10, 'Events', xoverFcn);

% Solve for passive and depth-keeping particles at every depth
for iorg = 1:npart
    % Passive 
    [T,P,restime,~,~] = ode45(@(tp,pp) odefcn_p(tp, pp, wave_u, wave_w), ...
        tspan, [x0 z_org(iorg)], opts);
    
    p_xpos(iorg, 1:length(T)) = P(:,1);
    p_zpos(iorg, 1:length(T)) = P(:,2);
    p_restime(iorg)           = restime;
    p_deltax(iorg)            = P(end,1) - P(1,1);
    
    % Depth-keeping
    [T,P,restime,~,~] = ode45(@(tp,pp) odefcn_dk(tp, pp, wave_u), ...
        tspan, [x0 z_org(iorg)], opts);
    
    dk_xpos(iorg, 1:length(T)) = P(:,1);
    dk_zpos(iorg, 1:length(T)) = P(:,2);
    dk_restime(iorg)           = restime;
    dk_deltax(iorg)            = P(end,1) - P(1,1);
    
end

% Store particle information for output
p.xp        = p_xpos;   
p.zp        = p_zpos;  
p.deltat    = deltat;  
p.time      = time;    
p.restime   = p_restime;
p.deltax    = p_deltax;

dk.xp       = dk_xpos;
dk.zp       = dk_zpos;
dk.deltat   = deltat;
dk.time     = time;  
dk.restime  = dk_restime;
dk.deltax   = dk_deltax;

% Wave information
wave.Amp    = amp;
wave.H      = H;
wave.period = period;
wave.cp     = cp;
wave.lambda = lambda;
wave.wshape = wshape;

% Wave vectors
wave.xgrid  = x;
wave.zgrid  = z;
wave.rho0   = rho0;

% Wave fields
wave.u      = wave_u(wave.xgrid, wave.zgrid, 0);
wave.w      = wave_w(wave.xgrid, wave.zgrid, 0);
wave.eta    = eta(wave.xgrid, wave.zgrid, 0);
wave.rho    = rho(wave.xgrid, wave.zgrid, 0);

end

%% Define ODE functions to integrate wave velocities
% Passive
function dppdt = odefcn_p(tp, pp, wave_u, wave_w)
    dppdt = zeros(2,1);  % pp(1) is x, pp(2) is z
    dppdt(1) = wave_u(pp(1),pp(2),tp);  
	dppdt(2) = wave_w(pp(1),pp(2),tp); 
end

% Depth-keeping
function dppdt = odefcn_dk(tp, pp, wave_u)
    dppdt = zeros(2,1);  % pp(1) is x, pp(2) is z
    dppdt(1) = wave_u(pp(1),pp(2),tp);
	dppdt(2) = 0;  % particle resists vertical velocities
end

%% Define event function to stop integration after one wave
function [value,isterminal,direction] = crestEventsFcn(tp,pp,x0,cp)
    wcrest = (-x0) + cp * tp - pp(1);  % crest catches up to particle
    value = [wcrest];  % value that we want to be zero
    isterminal = 1;  % halt integration 
    direction = 0;   % zero can be approached from either direction
end
