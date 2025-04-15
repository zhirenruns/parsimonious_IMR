%                       Parsimonious IMR (pIMR) 
%           Toward rapid first-estimate of viscoelastic properties
%      
%       NHKV Solver
%
%       Zhiren Zhu (zhiren@umich.edu)
%
%       Updated: April 2025
%
% =========================================================================
% Usage:
%
%   This function estimates collapse time for a given set of input data
%   describing a quadratic Kelvin-Voigt material ongoing inertial
%   cavitation.
%
% =========================================================================

function t1_approx = fit_QKV(G_guess,mu_guess,alp_guess,data_fit)

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

trc = 1/ARC;
Ro = 1./LX;
C = 2.1844;

% Get fitting parameters
Ca_guess = Ca_scale/G_guess;    % Cauchy #
Re_guess = Re_scale/mu_guess;   % Reynolds #

fbarv = -4.*C.^2./(2.*C.^2+sqrt(4.*C.^4+C.^2.*Re_guess.^2.*trc^2));

part1 = 40*sqrt(pi).*Ro*(1-3*alp_guess);
part2 = 120*(-1+2*Ro.^3)*alp_guess*gamma(7/6) ./ (Ro*gamma(2/3));
part3 = ((-50+177*alp_guess))*gamma(5/6)/gamma(4/3);
fbare = 1 - 1/(60*Ca_guess*gamma(5/6)) * gamma(1/3) * (part1 + part2 + part3);

% Viscoelastic contribution
fmodel = fbare + fbarv;

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