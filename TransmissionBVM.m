function DE = TransmissionBVM(d, n1, theta_B, dtheta, lambda0, n)

% July 14, 2022:  This appears to still only work for unslanted gratings!

% Transmission Holography
% Calculate diffraction efficiency of a volume transmission 
% hologram using the beta value method. 
%
% Input parameters
% d         Thickness of hologram, in meters
% n1        Peak-to-mean index variation of hologram
% theta_B   Internal incident angle at Bragg condition [rad]
% dtheta    Internal incident angular sweep variable [rad]
% lambda0   Readout wavelength, in meters, in free space
% n         Average background index of the material
%
% Output
% Diffraction efficiency (fractional) as a function of dtheta (deg)
%
% History
% Amy Sullivan                      Originate
% Robert McLeod     Oct 23, 2020    Handle slanted gratings
%                                   Input angles in radians

% Calculate Bragg angle and define angles, theta, over which the hologram
% is evaluated.
% theta_B = asind(lambda0/(2*n*period));   % Bragg angle, in material -
% UNSLANTED ONLY.  REPLACE WITH PASSING BRAGG ANGLE

theta = theta_B + dtheta;       	% Array of incident angles in material

% Define variables as defined in Fally, M., Klepp, J., & Tomita, Y. (2012).
% "An experimental study on the validity of diffraction theories for 
% off-Bragg replay of volume holographic gratings." 
% Applied Physics B, 108(1), 89-96.
Beta    = 2*pi*n./lambda0;
cR      = cos(theta);
cS      = sqrt(1-(2.*sin(theta_B)-sin(theta)).^2);
nu      = pi*n1*d./(lambda0.*sqrt(cR.*cS));
xi      = (cos(theta)-cS).*Beta.*d./2;

% Calculate diffraction efficiency as a function of incident angle
% FKT Equation 1
DE      = (cR./cS).*nu.^2.*(sinc(sqrt(nu.^2+xi.^2)./pi)).^2;

end

