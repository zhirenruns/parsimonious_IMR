%                       Parsimonious IMR (pIMR) 
%           Toward rapid first-estimate of viscoelastic properties
%      
%       Main Driver Script
%
%       Zhiren Zhu (zhiren@umich.edu)
%
%       Updated: Jan. 2025
%
% =========================================================================
% Usage:
%
%   This main driver scripts reads in a matrix containing summary info from
%   multiple micro-cavitation experiments and batch fits a suitable set of
%   viscoelastic properties.
%
% =========================================================================
%
% Required input format:
%   size(data_in) = N x 3
%       where N is the # of experiments performed
%       Data in each column:
%           (1) Rmax: maximum bubble radius (at beginning of collapse)
%           (2) Lmax = Rmax/Req: amplification factor
%           (3) t1: estimated collapse time
%
% =========================================================================

clc; close all; 
clearvars;

addpath('../graphics'); % To access color map etc.

%% User input

% Input file
infile = 'data/pIMR_PA_10_Refined.mat';   % Name of file to read
load(infile);

data_in = pIMR_array;                     % Name of variables to read

% Output options
make_plot = 1;       % Set to 1 to plot cost function for KV fit, as demo.

% Physical parameters used
p_inf = 101325;         % (Pa) Atmospheric Pressure
rho = 998.2;            % (kg/m^3) Density of characterized material
gam = 0.070;            % (N/m) Surface tension
cwave = 1484;           % (m/s) Wave speed in characterized material
pvsat = 3116.7757;      % (Pa) Saturated vapor pressure at far-field temperature

% pIMR non-viscoelastic parameters
C_kap = 1.4942;

%% Read data and pre-process

nX = size(data_in,1);   % # of experiments
RX = data_in(:,1);      % All Rmax 
LX = data_in(:,2);      % All amplification Lmax
T1X = data_in(:,3);     % All collapse time t1

%% Set up model parameters

pbar = p_inf;

We = pbar*RX/(2*gam);
uc = sqrt(pbar/rho);
Ma = cwave/uc;

ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) ); % ~= 1/0.9147
f_Ma = (2/sqrt( 1 + (1 + 4*(Ma/ARC)^2 ))) * ones(nX,1); % This is constant given fixed Ma, but create vector
f_We = - pi*ARC./(sqrt(6) * We); % This varies with Rmax


Req = (RX./LX);
PG0 = (p_inf).*(LX.^(-3));  
Alpha = 1 + pvsat./PG0;

f_gas = C_kap*(LX.^(-3)).*Alpha;

Ca_scale = pbar*ones(nX,1); % Ca*G
Re_scale = rho*uc*RX; % = Re*mu
De_scale = RX/uc; % This is the characteristic time scale

% Construct new matrix to pass to solver:
data_fit = [RX,LX,f_Ma,f_We,f_gas,Ca_scale,Re_scale,De_scale]; % These are the info needed to estimate t1 for given viscoelastic parameters

% Convert T1X to dimensionless, relative to t_{RC} at each scale.
% This way we don't need to (a) play with small numbers, (b) pass the 
% characteristic scale to solver functions repetitively.

tRC = (RX/uc)/ARC;
T1_ND = T1X./tRC;

%% NHKV fit

err_NHKV = @(X) (T1_ND./fit_NHKV(X(1),X(2),data_fit)).^2 - 1;
err_fn_NHKV = @(X) log10((err_NHKV(X))'*(err_NHKV(X))/nX);

G_min = 1E0;
G_max = 1E6;
mu_min = 0.0;
mu_max = 1.0;
G_start = 1E4;
mu_start = 0.1;

opt_fit = fminsearchbnd(err_fn_NHKV,[G_start,mu_start],[G_min,mu_min],[G_max,mu_max]);

G_opt_KV = opt_fit(1);
mu_opt_KV = opt_fit(2);

% Also find single-parameter fits:
opt_fit_justNeoH = fminsearchbnd(err_fn_NHKV,[G_start,0],[G_min,0],[G_max,0]);
opt_fit_justNewt = fminsearchbnd(err_fn_NHKV,[0,mu_start],[0,mu_min],[0,mu_max]);

disp("NH Best Fit: G = " + opt_fit_justNeoH(1) + " Pa.")
disp("Newtonian Best Fit: mu = " + opt_fit_justNewt(2) + " Pa*s.")
disp("KV Best Fit: G = " + G_opt_KV + " Pa, mu = " + mu_opt_KV + " Pa*s.")



%% SLS fit

err_SLS = @(X) (T1_ND./fit_SLS(X(1),X(2),X(3),data_fit)).^2 - 1;
err_fn_SLS = @(X) log10((err_SLS(X))'*(err_SLS(X))/nX);

G_min = 1;
G_max = 1E6; 
G_start = G_opt_KV;

mu_min = 0.0;
mu_max = 1.0;
mu_start = 0.01;

tau1_min = 1E-9;
tau1_max = 1E-1;
tau1_start = 1E-7;

options = optimset('TolFun',1E-8,'MaxIter', 8000, 'MaxFunEvals', 2000);

opt_fit = fminsearchbnd(err_fn_SLS,[G_start,mu_start,tau1_start],[G_min,mu_min,tau1_min],[G_max,mu_max,tau1_max], options);

G_opt_SLS = opt_fit(1);
mu_opt_SLS = opt_fit(2);
tau1_opt_SLS = opt_fit(3);

disp("SLS Best Fit: G = " + G_opt_SLS + " Pa, mu = " + mu_opt_SLS + " Pa*s, tau1 = " + tau1_opt_SLS + "s.")

%% Plot of Kelvin-Voigt model cost function space

if make_plot == 1

    % (1) Define parameter space to sweep
    dx = 2E-3;
    G_try = 10.^(3:dx:5);
    mu_try = 10.^(-3:dx:0);
    
    nG = length(G_try);
    nmu = length(mu_try);
    
    Err_NHKV = zeros(nG,nmu);

    % (2) Calculate cost function 
    disp('Preparing contour plot ...')

    for ii = 1:nG
    
        Gi = G_try(ii);
    
        for jj = 1:nmu
    
            muj = mu_try(jj);
    
            Err_NHKV(ii,jj) = err_fn_NHKV([Gi,muj]);
        end
    end

    % (3) Normalize
    psi_min = err_fn_NHKV([G_opt_KV,mu_opt_KV]);
    Err2 = Err_NHKV - psi_min;

    % (4) Plot
    figure(100);
    hold on; box on;

    kilo = 1E3;
    ft_sz = 12;
    bnd = [0.1,0.5:0.5:10]; % Contour levels

    [C1,h1] = contour(G_try/kilo, mu_try, Err2', bnd, 'LineWidth', 1.0 );
    clabel(C1, h1, bnd, 'FontSize',ft_sz*(2/3),'Interpreter','Latex'); 

    plot(G_opt_KV/kilo, mu_opt_KV, 'Marker', 'pentagram', 'MarkerSize', 10, 'Color', 'k', 'MarkerFaceColor', 'k');

    colormap(inferno)

    cmin = 0;
    cmax = 4;
    clim([cmin cmax])

    cbar = colorbar;
    set(cbar,'TickLabelInterpreter','Latex','FontSize',ft_sz)
    cbar_ttl = get(cbar,'Title');
    set(cbar_ttl ,'String',"$\hat{\psi}$",'Interpreter',"Latex",'FontSize',ft_sz);

    set(gca,'TickLabelInterpreter','Latex','FontSize',ft_sz)
    xl = xlabel('$G$ [kPa]');
    yl = ylabel('$\mu~[\rm Pa \cdot s]$');
    set(xl,'Interpreter','Latex','FontSize',ft_sz)
    set(yl,'Interpreter','Latex','FontSize',ft_sz)

    pbaspect([1,1,1])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end