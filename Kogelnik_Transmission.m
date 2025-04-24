%==========================================================================
% Kogelnik calculation of diffraction efficiency for a transmission
% hologram.
%
% Input parameters
%   d          Thickness of grating (in meters)
%   n          Average index of material 
%   n1         Peak to mean index change of hologram
%   Lambda     Period of grating (in meters)
%   lambda0    Read wavelength in free space (meters)
%   phi        Slant angle (90 deg is unslanted) [rad]
%   theta      Input values of angles over which to measure the
%               diffraction efficiency, internal angles [rad]
%
% % History
% J. Kowalski & A Sullivan  Oct 27, 2020    Originate
% R. McLeod                 Oct 28, 2020    Small changes to units (e.g. deg->rad)
%==========================================================================
function DE = Kogelnik_Transmission(d, n, n1, Lambda, lambda0, phi, theta)

thetaB  = asin(lambda0/(2*n*Lambda))+phi - pi/2; 	% Bragg angle (rad)
dtheta  = theta - thetaB;                           % Angle detuning (rad)

% Kogelnik variables
K       = 2*pi/Lambda;
Beta    = 2*pi*n./lambda0;

cR      = cos(theta);
cS      = cos(theta)-(K./Beta).*cos(phi);

nu      = pi*n1*d./(lambda0.*sqrt(cR.*cS));
xi      = dtheta.*K.*d.*sin(phi-thetaB)./(2.*cS);

% Calculate diffraction efficiency as a function of incident angle
DE      = (sin(sqrt(xi.^2+nu.^2))).^2./(1+xi.^2./nu.^2);

end