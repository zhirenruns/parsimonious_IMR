%                       Parsimonious IMR (pIMR) 
%           Toward rapid first-estimate of viscoelastic properties
%      
%       SLS Solver
%
%       Zhiren Zhu (zhiren@umich.edu)
%
%       Updated: March 2025, analytically estimate growth effect
%
% =========================================================================
% Usage:
%
%   This function estimates collapse time for a given set of input data
%   describing a finite deformation SLS material ongoing inertial
%   cavitation.
%
% =========================================================================

function t1_approx = fit_SLS(G_guess,mu_guess,tau1_guess,data_fit)

% Inputs:
%   G_guess - elastic shear modulus (Pa)
%   mu_guess - viscous shear modulus (Pa*s)
%   tau1_guess - relaxation time scale (s)

% Unwrap input
RX = data_fit(:,1);     % Not used, but let's keep for debugging
LX = data_fit(:,2);
f_Ma = data_fit(:,3);
f_We = data_fit(:,4);
f_gas = data_fit(:,5);
Ca_scale = data_fit(:,6);
Re_scale = data_fit(:,7);
De_scale = data_fit(:,8);
hypg = data_fit(:,9); % Hypergeometric function!

% Other constants to use:
ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) ); % ~= 1/0.9147

% Get fitting parameters
Ca_guess = Ca_scale./G_guess;    % Cauchy #
Re_guess = Re_scale./mu_guess;   % Reynolds #
De_guess = tau1_guess./De_scale;       % Deborah #

B_elast = 5/2 - sqrt(2/3)*pi*ARC./LX;
Beta = 1./(1 + B_elast./Ca_guess);

Y = 2*ARC * (0.4637./(Re_guess) + 0.56598./(Re_guess.^2) + 5.7331./(Re_guess.^3));
Z = 1 - ( Y + sqrt(Y.^2 + 1) ).^(-2); % This is fv
blob = ARC*De_guess;

fNH = (1 - 1./Beta);

% Get baseline effect, without growth:
fM0 = Z + (blob).*((exp(-1./blob)-1).*Z);
fgrowth = fM0 + fNH + (f_gas + f_We + f_Ma); % Reverse effect during growth

% Get growth time:
R0 = 1./LX;

% Get growth time. An array involving hypergeometric function, so run loop:
% nx = length(LX);
% TG = zeros(nx,1);
% for ii = 1:nx
%     TG(ii) = (5*sqrt(pi)*gamma(5/6)-6*R0(ii)^(5/2)*gamma(4/3)*hypergeom([1/2,5/6],11/6, R0(ii)^3))./(5*sqrt(6-6.*fgrowth(ii))*gamma(4/3));
% end

% Speed up: pass hypergeom result in, since it only depends on R0:
TG = (5*sqrt(pi)*gamma(5/6)-6*(R0.^(5/2))*gamma(4/3).*hypg)./(5*sqrt(6-6.*fgrowth)*gamma(4/3));

Z0 = - Z.*(1 - exp(-TG./De_guess)); % Negative valued

fM = Z + (blob).*((exp(-1./blob) - 1).*(Z - Z0));
fmodel = fNH + fM;

fsum = 1 - fmodel - f_We - f_Ma - f_gas; 

% Add an extra layer of check: if growth time is absurd, reject

if min(fsum) < 0 
    % warning('Negative fsum encountered for SLS fit, G = %f, mu = %f, tau1 = %f', G_guess, mu_guess, tau1_guess);
    t1_approx = nan;
else
    t1_approx = (fsum).^(-1/2);
end

end