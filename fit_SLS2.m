%                       Parsimonious IMR (pIMR) 
%           Toward rapid first-estimate of viscoelastic properties
%      
%       SLS Solver
%
%       Zhiren Zhu (zhiren@umich.edu)
%
%       Updated: Sept 2024
%
% =========================================================================
% Usage:
%
%   This function estimates collapse time for a given set of input data
%   describing a finite deformation SLS material ongoing inertial
%   cavitation.
%
% =========================================================================

function t1_approx = fit_SLS2(G_guess,mu_guess,tau1_guess,data_fit)

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

nx = length(LX);

% Other constants to use:
ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) ); % ~= 1/0.9147

% Get fitting parameters
Ca_guess = Ca_scale./G_guess;    % Cauchy #
Re_guess = Re_scale./mu_guess;   % Reynolds #
De_guess = tau1_guess./De_scale;       % Deborah #

B_elast = 5/2 - sqrt(2/3)*pi*ARC./LX;
Beta = 1./(1 + B_elast./Ca_guess);

Y = 2*ARC * (0.4637./(Re_guess) + 0.56598./(Re_guess.^2) + 5.7331./(Re_guess.^3));
Z = 1 - ( Y + sqrt(Y.^2 + 1) ).^(-2);

blob = ARC*De_guess;

% Find initial stress:
FM0 = zeros(nx,1);

for ii = 1:nx
    FM0(ii) = get_S0_SLS(LX(ii),Re_guess(ii),De_guess(ii),Ca_guess(ii));
end

DF = Z - FM0;
FMT = Z + (blob).*( DF .* exp(-blob.^(-1)) - DF );
fmodel = (1 - 1./Beta) + FMT;

fsum = 1 - fmodel - f_We - f_Ma - f_gas; 

if min(fsum) < 0 
    % warning('Negative fsum encountered for SLS fit, G = %f, mu = %f, tau1 = %f', G_guess, mu_guess, tau1_guess);
    t1_approx = nan;
else
    t1_approx = (fsum).^(-1/2);
end

end