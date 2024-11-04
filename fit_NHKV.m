%                       Parsimonious IMR (pIMR) 
%           Toward rapid first-estimate of viscoelastic properties
%      
%       NHKV Solver
%
%       Zhiren Zhu (zhiren@umich.edu)
%
%       Updated: July 2024
%
% =========================================================================
% Usage:
%
%   This function estimates collapse time for a given set of input data
%   describing a Neo-Hookean Kelvin-Voigt material ongoing inertial
%   cavitation.
%
% =========================================================================

function t1_approx = fit_NHKV(G_guess,mu_guess,data_fit)

% Inputs:
%   G_guess - elastic shear modulus (Pa)
%   mu_guess - viscous shear modulus (Pa*s)

% Unwrap input
RX = data_fit(:,1);     % Not used, but let's keep for debugging
LX = data_fit(:,2);
f_Ma = data_fit(:,3);
f_We = data_fit(:,4);
f_gas = data_fit(:,5);
Ca_scale = data_fit(:,6);
Re_scale = data_fit(:,7);

% Other constants to use:
ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) ); % ~= 1/0.9147

% Get fitting parameters
Ca_guess = Ca_scale/G_guess;    % Cauchy #
Re_guess = Re_scale/mu_guess;   % Reynolds #

B_elast = 5/2 - sqrt(2/3)*pi*ARC./LX;
Beta = 1./(1 + B_elast./Ca_guess);

Y = 2*ARC * (0.4637./(Re_guess) + 0.56598./(Re_guess.^2) + 5.7331./(Re_guess.^3));
Z = Beta.*Y;

fmodel = ( Z + sqrt(Z.^2 + Beta) ).^(-2);

fsum = fmodel - f_We - f_Ma - f_gas; 

if min(fsum) < 0 
    % Re too high for at least one scale (i.e., ~ over-damped)
    % warning('Negative fsum encountered for NHKV fit, G = %f, mu = %f', G_guess, mu_guess);
    t1_approx = nan;
else
    t1_approx = (fsum).^(-1/2);
end

% We have solution for t1_approx

% ===== End of Function =====
end 