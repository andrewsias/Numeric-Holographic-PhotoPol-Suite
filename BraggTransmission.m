%==========================================================================
% BraggTransmission - Calculate diffraction efficiency of a transmission
% hologram in weak or strong limit with uniform or variable amplitude and
% phase.  Valid for any grating slant.
%
% The derivation and notation follow Uchida *except* here diffraction is
% off of +K, not -K as in Uchida.  That is, here k1 = k0 + K, not k0 - K.
% Also, all optical quantities are found without small angle approximations
% via k-space geometry.  These agree (to within a sign) with Uchida for
% thick gratings (small Bragg detuning angle) but diverge noticibly for
% thin (e.g. 15 micron) gratings in the sidelobes.
%
% Input parameters: See dictionary in parser for name value pairs of optional inputs
%   dRedTheta        Scalar or vector of read angles relative to Bragg matched [deg]
%
% Per eq. 13 of Uchida, DE is c1 |A1|^2 / (c0 |A0|^2) where the c1/c0 term is the projection
% of the incident and diffracted Poynting vectors into the direction normal
% to the surface.  This terms only applies if power in the lab is measured
% appropriately and, in most experimental cases, this is not the case.  
% The input parameter PowerNormal controls this.
%
% DE can be calculated in various limits including weak or strong
% diffraction and distorted or undistorted fringes.  The simplest
% expression available will be used depending on input parameters.
%
% Rather than use grating period and tilt, the routine requests the
% experimental writing conditions (beam angles internal or external plus
% wavelength) and calculates the required quanties, returning them in the
% "out" structure.
%
% Sources 
% 1. Oculus Holographic Photopolymers/New program June 2022/Aims/Aim 3 Optical/Fitting Bragg selectivity with phase
% 2. "All derivitaions Bragg with distortion.nb" RRM June 2022
% 3. ﻿Calculation of diffraction efficiency in hologram gratings attenuated along the direction perpendicular to the grating vector, Uchida 1973
%
% Sample call comparing undistorted, polynomial and Fourier series distorted
% dRedTheta = [-6:.01:6]; % degrees
% [Eta,in,out]    = BraggTransmission(dRedTheta,'n1',.01); 
% [EtaP,inP,outP] = BraggTransmission(dRedTheta,'dPhi2',pi/2,'n1',.01); 
% [EtaF,inF,outF] = BraggTransmission(dRedTheta, 'Phi1',pi/2,'n1',.01);
% plot(dRedTheta,Eta,'k',dRedTheta,EtaP,'r',dRedTheta,EtaF,'b');
% set(gcf,'color','w');set(gca,'FontSize',18);set(gca,'FontName','Arial');
% xlabel('\delta{\theta}_{ext} [^o]');ylabel('\eta');
%
% History
% RRM         Sep 12 2022   Originate
% RRM         Oct 16 2022   V1.0: Polynomial and Fourier distortion, weak &
%                           strong, uniform and exponential amplitude.
%                           To be done: table export for weak with amplitude
%                           Control over ODE function and options
% RRM         Nov 16, 2022  Fourier series for amplitude
%==========================================================================
function [Eta, in, out] = BraggTransmission(dRedTheta,varargin)

%--------------------------------------------------------------------------
% Parse and calculate holography values
%--------------------------------------------------------------------------
in  = ParseTransBVMInputs(dRedTheta,varargin);
out = Holography(in);

% Set directory for lookup tables
in.TableDir = '/Users/mcleod/Google Drive/My documents/Professional/CU/Public/Research/Funded projects/Active/Oculus Holographic Photopolymers/New program June 2022/Aims/Aim 3 Optical/Fitting Bragg selectivity with phase';

%--------------------------------------------------------------------------
% Calculate diffraction efficiency with appropriate method
%--------------------------------------------------------------------------
% If no phase or amplitude distortion terms
if                  in.dPhi2 == 0 && in.dPhi3 == 0 && in.dPhi4 == 0 && ...
   in.Phi1  == 0 && in.Phi2  == 0 && in.Phi3  == 0 && in.Phi4  == 0 && in.Phi5  == 0 && in.Phi6  == 0 && ...
   in.A1    == 0 && in.A2    == 0 && in.A3    == 0 && in.A4    == 0 && in.A5    == 0 && in.A6    == 0

    %----------------------------------------------------------------------
    % Undistorted 
    %----------------------------------------------------------------------
    if in.Born
        % Weak limit: Undistorted Born
        out = DETransmissionBornAbsorbing(in,out);
    else
        % Strong limit: Uchida BVM
        out = DETransmissionUndistorted(in,out);
    end
else
    %----------------------------------------------------------------------
    % Distorted and/or attenuated Born (weak) diffraction
    %----------------------------------------------------------------------
    if in.Born
        if in.dPhi2 ~= 0 || in.dPhi3 ~= 0 || in.dPhi4 ~= 0      % Have polynomial phase terms been set
            if in.dPhi3 == 0 && in.dPhi4 == 0                   % Only quadratic term, dPhi2 set          
                % Weak limit and quadratic distortion
                out = DETransmissionBornDistortedAbsorbing(in,out);
            else                                                % Full polynomial solution requested
                % Weak limit and polynomial distortion
                out = DETransmissionBornPolynomialDistortedNumeric(in,out); 
            end
        else                                                    % Not polynomial, so Fourier series
            out = DETransmissionBornFourierDistortedAndAttenuatedNumeric(in,out);
        end
    else
    %----------------------------------------------------------------------
    % Distorted and/or attenuated strong diffraction
    %----------------------------------------------------------------------
        % If no Fourier series terms set
        if in.Phi1  == 0 && in.Phi2  == 0 && in.Phi3  == 0 && in.Phi4  == 0 && in.Phi5  == 0 && in.Phi6  == 0 && ...
           in.A1    == 0 && in.A2    == 0 && in.A3    == 0 && in.A4    == 0 && in.A5    == 0 && in.A6    == 0
            % Polynomial distortion
            if in.dPhi2 ~= 0 && in.dPhi3 == 0 && in.dPhi4 == 0
                % Strong limit and quadratic distortion
                out = DETransmissionDistortedTheory(in,out);    
            else
                % Strong limit and polynomical disortion, numerical 
                out = DETransmissionPolynomialDistortedNumeric(in,out);   
            end
        else
            % Strong limit and Fourier series distortion and attenuation
            out = DETransmissionFourierDistortedAndAttenuatedNumeric(in,out);   
        end
    end
end

%--------------------------------------------------------------------------
% Post processing
%--------------------------------------------------------------------------
% Correct for cosine projection
if in.PowerNormal
    out.Eta = (out.c1/out.c0) .* out.Eta;
end

% Return variable.  Insure row vector for matlab fit function
Eta = reshape(out.Eta,length(out.Eta),1);    

end % TransmissionBVMWithDistortion

%==========================================================================
% Parse inputs
%==========================================================================
function in = ParseTransBVMInputs(dRedTheta,VariablesPassed)

nano  = 10^-9;                      % units.  All calculations MKS.
micro = 10^-6;
%milli = 10^-3;

p = inputParser;                    % Parser structure

%---Parser validation functions. 
isPosScalar  = @(x) isnumeric(x) && isscalar(x) && (x >= 0);     % Scalar > 0

%------------------------------Parameter dictionary------------------------
%                       PARAM             DEFAULT  VALIDITY      DEFINITION
%--------------------------------------------------------------------------
%---Program options and controls
addParameter(p,      'External',            false, @isscalar   ); % Beam angles are external to medium
addParameter(p,          'Born',            false, @isscalar   ); % Assume Born (weak) approximation
addParameter(p,  'RadianAngles',            false, @isscalar   ); % false = angles in degrees, true = radians
addParameter(p,   'PowerNormal',            false, @isscalar   ); % TRUE: power is measured normal to surface with an overfilled 
                                                                  % detector (so Sz).  FALSE: power is measured normal to k or the 
                                                                  % detector area is larger than the beam
%---Holography experiment description
addParameter(p,   'RefWrtTheta',           -22.5,  @isscalar   ); % Reference write beam angle from surface normal [deg]
addParameter(p,   'ObjWrtTheta',            22.5,  @isscalar   ); % Object write beam angle from surface normal [deg]
addParameter(p,    'lambdaWrt0',        405*nano,   isPosScalar); % Vacuum write optical wavelength
addParameter(p,    'lambdaRed0',        630*nano,   isPosScalar); % Vacuum  read optical wavelength
addParameter(p,          'nWrt',             1.5,   isPosScalar); % Refractive index at write wavelength
addParameter(p,          'nRed',             1.5,   isPosScalar); % Refractive index at read wavelength

%---Primary fit coefficients
addParameter(p,            'n1',            0.01,   isPosScalar); % Peak-to-mean index variation at read wavelength.  
addParameter(p,             'L',        15*micro,   isPosScalar); % Material thickness [m]
addParameter(p,         'alpha',               0,  @isscalar   ); % Exponential grating decay (>0) of grating amplitude in z [1/m]

%---Polynomial phase coefficients
addParameter(p,         'dPhi2',               0,  @isscalar   ); % Peak quadratic phase distortion [rad]
addParameter(p,         'dPhi3',               0,  @isscalar   ); % Peak cubic phase distortion [rad]
addParameter(p,         'dPhi4',               0,  @isscalar   ); % Coefficient of quartic phase distortion [rad]
addParameter(p,      'z04OverL',            0.25,  @isscalar   ); % Location of symmetric nulls in quartric distortion as fraction of L.

%---Fourier series phase coefficients
addParameter(p,          'Phi1',               0,  @isscalar   ); % First  phase harmonic (Period = 2 L / 1)
addParameter(p,          'Phi2',               0,  @isscalar   ); % Second phase harmonic (Period = 2 L / 2)
addParameter(p,          'Phi3',               0,  @isscalar   ); % Third  phase harmonic (Period = 2 L / 3)
addParameter(p,          'Phi4',               0,  @isscalar   ); % Fourth phase harmonic (Period = 2 L / 4)
addParameter(p,          'Phi5',               0,  @isscalar   ); % Fifth  phase harmonic (Period = 2 L / 5)
addParameter(p,          'Phi6',               0,  @isscalar   ); % Sixth  phase harmonic (Period = 2 L / 6)

%---Fourier series amplitude coefficients
addParameter(p,            'A1',               0,  @isscalar   ); % First  amplitude harmonic (Period = 2 L / 1)
addParameter(p,            'A2',               0,  @isscalar   ); % Second amplitude harmonic (Period = 2 L / 2)
addParameter(p,            'A3',               0,  @isscalar   ); % Third  amplitude harmonic (Period = 2 L / 3)
addParameter(p,            'A4',               0,  @isscalar   ); % Fourth amplitude harmonic (Period = 2 L / 4)
addParameter(p,            'A5',               0,  @isscalar   ); % Fifth  amplitude harmonic (Period = 2 L / 5)
addParameter(p,            'A6',               0,  @isscalar   ); % Sixth  amplitude harmonic (Period = 2 L / 6)

%--------------------------------------------------------------------------
% Parse
%--------------------------------------------------------------------------
parse(p,VariablesPassed{:});                        % Parse inputs into struct p
in = p.Results;                                     % Short hand notation 

%--------------------------------------------------------------------------
% Massage inputs as needed
%--------------------------------------------------------------------------
% Add angle sweep as input field for ease of passing around
in.dRedTheta = reshape(dRedTheta,length(dRedTheta),1);  % Insure column vector

if ~in.RadianAngles  % Calculations in radians.  Convert from deg if needed
    in.dRedTheta    = deg2rad(in.dRedTheta);
    in.RefWrtTheta  = deg2rad(in.RefWrtTheta);
    in.ObjWrtTheta  = deg2rad(in.ObjWrtTheta);
end

%--------------------------------------------------------------------------
% Error checking
%--------------------------------------------------------------------------
if in.dPhi4 ~= 0 && (in.z04OverL <0 || in.z04OverL > 1), error("Quartic null location / L must be between 0 and 1"); end

end % ParseTransBVMInputs

%==========================================================================
% Do Bragg calculations.
%==========================================================================
function out = Holography(in)

NAngle = length(in.dRedTheta);                  % How many incident angles

%--------------------------------------------------------------------------
% Refract writing beams into material if required
%--------------------------------------------------------------------------
if in.External
    out.RefWrtThetaInternal = asin(in.nWrt * in.RefWrtTheta);   % Snell's law
    out.ObjWrtThetaInternal = asin(in.nWrt * in.ObjWrtTheta);
else
    out.RefWrtThetaInternal = in.RefWrtTheta;                   % Just copy
    out.ObjWrtThetaInternal = in.ObjWrtTheta;
end

%--------------------------------------------------------------------------
% Write hologram.  2D coordinate system is [z,x], see Uchida Figure 1.
%--------------------------------------------------------------------------
out.RefWrtk     = 2*pi/in.lambdaWrt0 * [cos(out.RefWrtThetaInternal),sin(out.RefWrtThetaInternal)];
out.ObjWrtk     = 2*pi/in.lambdaWrt0 * [cos(out.ObjWrtThetaInternal),sin(out.ObjWrtThetaInternal)];
out.K           = out.ObjWrtk - out.RefWrtk;
out.Lambda      = 2*pi/norm(out.K);
% Grating angle, K relative to z, Uchida phi as shown in Figure 1.
out.GratingPhi  = atan2(out.K(2),out.K(1));         

%--------------------------------------------------------------------------
% Calculate Bragg matched reference and object k vectors.  This technique
% uses vector math, not trig functions, to avoid sign ambiguities.
%--------------------------------------------------------------------------
out.RedNormk = 2*pi/in.lambdaRed0;         % |k_read| 

% Construct unit vectors parallel and normal to K
KHat = out.K/norm(out.K);
NHat = (out.RefWrtk + out.ObjWrtk) / norm(out.RefWrtk + out.ObjWrtk);   % Unit vector normal to K that bisects internal writing k vectors
    
% Project reference and object writing beams onto K to find Bragg matched reading k's
out.RefRedkB = dot(out.RefWrtk,KHat) * KHat;  % Components of Bragg matched reading k in direction of K
out.ObjRedkB = dot(out.ObjWrtk,KHat) * KHat;  % This calculation conserves momentum in KHat direction

% Keep this conserved momentum and add orthogonal component such that the final norm is |k| 
out.RefRedkB = out.RefRedkB + sqrt(out.RedNormk.^2 - norm(out.RefRedkB)^2) * NHat;
out.ObjRedkB = out.ObjRedkB + sqrt(out.RedNormk.^2 - norm(out.ObjRedkB)^2) * NHat;

%--------------------------------------------------------------------------
% Create vector of (non Bragg matched) reference k vectors based on set of angular inputs.
%--------------------------------------------------------------------------
out.RefRedThetaBInternal = atan2(out.RefRedkB(2),out.RefRedkB(1));                      % Internal Bragg angles
out.ObjRedThetaBInternal = atan2(out.ObjRedkB(2),out.ObjRedkB(1));                      

out.RefRedThetaBExternal = asin(in.nRed*sin(atan2(out.RefRedkB(2),out.RefRedkB(1))));   % External Bragg angles
out.ObjRedThetaBExternal = asin(in.nRed*sin(atan2(out.ObjRedkB(2),out.ObjRedkB(1))));  

% Add vector of reference delta angles to either external or internal reference ("0") Bragg angle
% to calculate reference ("0") read angle vector internal to the material
if in.External
    out.RefRedThetaInternal = asin(in.Red*sin(out.RefRedThetaBExternal + in.dRedTheta)); % Add delta theta to external angle and refract in
else
    out.RefRedThetaInternal = out.RefRedThetaBInternal + in.dRedTheta;     % Add delta theta to internal angle
end

%--------------------------------------------------------------------------
% Create vector of (non Bragg matched) object k vectors from references and K
% Using k space geometry.  Uchida expressions are off by sign (see figure 1
% - he is diffracting off of -K) but also slightly different.  Are his
% equations in a small angle limit?
%--------------------------------------------------------------------------
out.RefRedk = out.RedNormk * [cos(out.RefRedThetaInternal),sin(out.RefRedThetaInternal)];  % All reading reference k vectors
out.ObjRedk = out.RefRedk + repmat(out.K,NAngle,1);                 % NOT Bragg matched diffracted object
kTran       = out.ObjRedk(:,2);                                     % Mag of internal obj k transverse to surface (conserved)   
SignNorm    = sign(out.ObjRedk(:,1));                               % Dir of obj k relative to surface normal
out.ObjRedk = [SignNorm.*sqrt(out.RedNormk^2 - kTran.^2),kTran ];   % Bragg matching: Keep  kx, set z component to get |k| 
out.ObjRedThetaInternal = atan2(out.ObjRedk(:,2),out.ObjRedk(:,1));  

%--------------------------------------------------------------------------
% Cosine factors used (see Uchida 1c).  Ref = incident = 0, Obj = diffracted = 1
%--------------------------------------------------------------------------
out.c0          = cos(out.RefRedThetaInternal);
out.c1          = cos(out.ObjRedThetaInternal);
out.c1Uchida    = -cos(2*out.GratingPhi-2*out.RefRedThetaBInternal+out.RefRedThetaInternal);  % Uchida 1c

%--------------------------------------------------------------------------
% Bragg mismatch and coupling factors 
%--------------------------------------------------------------------------
out.DeltaK          = out.RefRedk - out.ObjRedk + repmat(out.K,NAngle,1);    % Uchida  2 but using +K so that k1 = k0 + K
out.DeltaKzL        = out.DeltaK(:,1) * in.L;               % z component, times L
out.DeltaKzLUchida  = norm(out.K)*in.L * (cos(out.GratingPhi-out.RefRedThetaBInternal+out.RefRedThetaInternal) - cos(out.GratingPhi)); % Uchida 1c

out.kappa0L         = in.n1*out.RedNormk *in.L/2;           % Uchida 1c
out.AlphaL          = in.alpha * in.L;
out.EtaMax          = out.kappa0L^2 ./(out.c0 .* out.c1);   % Max DE in weak limit for power normal to surface

end % Holography

%==========================================================================
% Born limit (weak) undistorted transmission hologram with exponential
% amplitude variation in z
%==========================================================================
function out = DETransmissionBornAbsorbing(~,out)

out.Eta = (out.c0./out.c1) .* out.EtaMax .* exp(-out.AlphaL) .* ...
          sincNoPi(out.DeltaKzL/2+1i*out.AlphaL/2).^2;

end % DETransmissionBornAbsorbing

%==========================================================================
% Born limit (weak) quadratically distorted transmission hologram with exponential
% amplitude variation in z
%==========================================================================
function out = DETransmissionBornDistortedAbsorbing(in,out)

out.Eta = (out.c0./out.c1) .* out.EtaMax .* pi./(32*in.dPhi2) .* ...
   exp(-(2*out.AlphaL*(in.dPhi2-out.DeltaKzL/4))/(8*in.dPhi2)) .* ...
   abs(         sqrt(2)*erfi((-1)^(1/4)*(in.dPhi2-out.DeltaKzL/4-1i*out.AlphaL/4)./sqrt(in.dPhi2)) + ...
    (1-1i).*(-1)^(1/4).*erfi((-1)^(1/4)*(in.dPhi2+out.DeltaKzL/4+1i*out.AlphaL/4)./sqrt(in.dPhi2))).^2;

end % DETransmissionBornDistortedAbsorbing

%==========================================================================
% Strong undistorted transmission hologram.  
%==========================================================================
function out = DETransmissionUndistorted(~,out)

out.Eta = (out.c0./out.c1) .* out.EtaMax .* sincNoPi(sqrt(out.DeltaKzL.^2/4+out.EtaMax)).^2;

end % DETransmissionUndistorted

%==========================================================================
% Strong quadraticaly distorted transmission hologram.  
% 
% See "Export Tables for Distorted Bragg.nb" for source calculation
%==========================================================================
function out = DETransmissionDistortedTheory(in,out)

% Declare lookup table and index vectors as persistent so only read once
persistent EtaMax  DeltakzL  DeltaPhi  A1atLSqrdc1Overc0_DP_DkzL_EM

%--------------------------------------------------------------------------
% If first call, read pre-calculated table for interpolation
%--------------------------------------------------------------------------
if isempty(EtaMax)          % Need to read?

    % Read table from Mathematica export containing solution of CPM
    % equations in presence of quadratic phase distortion
    % See Mathematica file for details 
    % Since .mat format can't handle multidimensional arrays (WTF?) this is a
    % flattened data set
    % Directory (set in main) must be given to locate this file
    A1atLSqrdc1Overc0 = importdata([in.TableDir,'/A1atLSqrdc1Overc0.mat']);

    % Vectors that index data set.
    % These must be defined consistently with the .nb file
    EtaMax   = 0.001 : 0.05 : 2;
    DeltakzL = 0     : 0.1  :(12*pi);
    DeltaPhi = 0.001 : 0.05 :(4*pi);

    % Restore dimensions.  Note that order of subscripts is reverse of the
    % indices in the Mathematica Table command
    A1atLSqrdc1Overc0_DP_DkzL_EM = reshape(A1atLSqrdc1Overc0,length(DeltaPhi), length(DeltakzL), length(EtaMax));

    % Add final entry of dPhi2 = 0 to interpolation table.  The expression
    % used in Mathematica does not converge at this limit, so instead use
    % the "strong limit with no distortion and no absorption (Uchida)"
    % expression, stripping the leading c0/c1 term as done for the export
    % of the rest of the table.

    % Create 2D coordinate arrays for DeltakzL
    [DeltakzL_DkzL_EM, EtaMax_DkzL_EM] = ndgrid(DeltakzL,EtaMax);

    % Create new entry using these values.  See DETransmissionUndistorted function
    A1atLSqrdc1Overc0_DP0_DkzL_EM = zeros(1,length(DeltakzL),length(EtaMax));  % first dim must = 1
    A1atLSqrdc1Overc0_DP0_DkzL_EM(1,:,:) = EtaMax_DkzL_EM .* ...
        abs(sincNoPi(sqrt(DeltakzL_DkzL_EM.^2/4 + EtaMax_DkzL_EM))).^2;

    % Add new entry to table and to dPhi2 vector
    A1atLSqrdc1Overc0_DP_DkzL_EM = cat(1,A1atLSqrdc1Overc0_DP0_DkzL_EM,A1atLSqrdc1Overc0_DP_DkzL_EM);
    DeltaPhi                     = [0,DeltaPhi];
end

%--------------------------------------------------------------------------
% Create normalized variables used to index table
% Currently assuming symmetry in dkzL
% Note INSANE matlab convention that must exchange first two variables
%--------------------------------------------------------------------------
EtaMaxQ   = out.kappa0L.^2 ./ (out.c0.*out.c1);   % Query vectors, all same length
DeltaPhiQ = repmat(in.dPhi2,size(EtaMaxQ));
DeltakzLQ = abs(out.DeltaKzL); 

% Interpolate
out.Eta = interp3(DeltakzL,  DeltaPhi,   EtaMax,  A1atLSqrdc1Overc0_DP_DkzL_EM,...
                  DeltakzLQ, DeltaPhiQ,  EtaMaxQ, "spline",0);

out.Eta = out.Eta ./ (out.c1./out.c0);  % Return |A1|^2

end % DETransmissionDistorted

%==========================================================================
% Strong power-law distorted transmission hologram.  
%==========================================================================
function out = DETransmissionPolynomialDistortedNumeric(in,out)

%--------------------------------------------------------------------------
% Create struct with required parameters
%--------------------------------------------------------------------------
P.kappa0L   = out.kappa0L;
P.dPhi2     = in.dPhi2;
P.dPhi3     = in.dPhi3;
P.dPhi4     = in.dPhi4;
P.z04OverL  = in.z04OverL;

%--------------------------------------------------------------------------
% For each Bragg mismatch
%--------------------------------------------------------------------------
for idKzL = 1:length(out.DeltaKzL)
    P.c0        = out.c0(idKzL);
    P.c1        = out.c1(idKzL);
    P.DeltaKzL  = out.DeltaKzL(idKzL);

    % Solve coupled mode equations for R and I of Inc and Diff fields
    % This combo of options ('AbsTol',1e-14) and solver (ODE78) yields
    % reasonbable confidence bounds on the polynomial phase coefficients.
    % All other combinations yield the same fit values but orders of
    % magnitude higher confidence bounds.
    options = odeset('AbsTol',1e-8);
    [~, A] = ode78(@(zOverL,A) RealCMEPolynomialDistorted(zOverL,A,P),[0 1],[1 0 0 0],options);

    % Efficiency is |A1|^2 at end of zOverL grid
    out.Eta(idKzL) = A(end,3)^2 + A(end,4)^2;
end

end % DETransmissionPolynomialDistortedNumeric

%==========================================================================
% RHS of coupled mode equations with polynomial phase distortion written in
% complex form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   A           Complex amplitudes of incident A(1) and diffracted A(2) 
%   P           Struct of parameters
%==========================================================================
function dA = ComplexCMEPolynomialDistorted(zOverL,A,P)

% Polynomial distortion of phase
Phi = PhiPoly(zOverL,P.dPhi2,P.dPhi3,P.dPhi4,P.z04OverL);

% RHS of coupled mode equations using complex representation
dA    = zeros(2,1);  
dA(1) = 1i * P.kappa0L * exp(-1i * P.DeltaKzL * zOverL + 1i * Phi) * A(2) / P.c0;
dA(2) = 1i * P.kappa0L * exp(+1i * P.DeltaKzL * zOverL - 1i * Phi) * A(1) / P.c1;

end % ComplexCMEPolynomialDistorted

%==========================================================================
% RHS of coupled mode equations with polynomial phase distortion written in
% real form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitudes of incident A(1)+jA(2) and diffracted A(3)+jA(4)
%   P           Struct of parameters
%==========================================================================
function dAReal = RealCMEPolynomialDistorted(zOverL,AReal,P)

% Accept real and imaginary parts of two fields, assemble into two complex variables
AComplex  = [AReal(1)+1i*AReal(2); AReal(3)+1i*AReal(4)];

% Call CME in complex form to calculate two complex results
dAComplex = ComplexCMEPolynomialDistorted(zOverL,AComplex,P);

% Separate two complex results into four real numbers
dAReal    = [real(dAComplex(1)); imag(dAComplex(1)); real(dAComplex(2)); imag(dAComplex(2))];

end % RealCMEPolynomialDistorted

%==========================================================================
% Weak (Born) power-law distorted transmission hologram.  
%==========================================================================
function out = DETransmissionBornPolynomialDistortedNumeric(in,out)

%--------------------------------------------------------------------------
% Create struct with required parameters
%--------------------------------------------------------------------------
P.kappa0L   = out.kappa0L;
P.dPhi2     = in.dPhi2;
P.dPhi3     = in.dPhi3;
P.dPhi4     = in.dPhi4;
P.z04OverL  = in.z04OverL;

%--------------------------------------------------------------------------
% For each Bragg mismatch
%--------------------------------------------------------------------------
for idKzL = 1:length(out.DeltaKzL)
    P.c0        = out.c0(idKzL);                % Parameter struct terms that change with deltaKzL
    P.c1        = out.c1(idKzL);
    P.DeltaKzL  = out.DeltaKzL(idKzL);

    % Solve coupled mode equations for R and I of Diff field.  BC is A1 = 0 + j 0
    [~, A] = ode45(@(zOverL,A) RealCMEBornPolynomialDistorted(zOverL,A,P),[0 1],[0 0]);

    % Efficiency is |A1|^2 at end of zOverL grid
    out.Eta(idKzL) = A(end,1)^2 + A(end,2)^2;
end

end % DETransmissionPolynomialDistortedNumeric

%==========================================================================
% RHS of coupled mode equations with polynomial phase distortion written in
% complex form in weak (Born) limit.
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitudes of diffrcted field A(1)+jA(2)
%   P           Struct of parameters
%==========================================================================
function dA = ComplexCMEBornPolynomialDistorted(zOverL,~,P)

% Polynomial distortion of phase
Phi = PhiPoly(zOverL,P.dPhi2,P.dPhi3,P.dPhi4,P.z04OverL);

% RHS of coupled mode equations using complex representation
dA = 1i * P.kappa0L * exp(+1i * P.DeltaKzL * zOverL - 1i * Phi) * (1) / P.c1;

end % ComplexCMEBornPolynomialDistorted

%==========================================================================
% RHS of coupled mode equations in Born limit with polynomial phase distortion written in
% real form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitude of diffracted A(1)+jA(2)
%   P           Struct of parameters
%==========================================================================
function dAReal = RealCMEBornPolynomialDistorted(zOverL,AReal,P)

% Accept real and imaginary parts of A1 field, assemble into complex variable
AComplex  = AReal(1)+1i*AReal(2);

% Call CME in complex form to calculate two complex results
dAComplex = ComplexCMEBornPolynomialDistorted(zOverL,AComplex,P);

% Separate complex results into two real numbers
dAReal    = [real(dAComplex); imag(dAComplex)];

end % RealCMEBornPolynomialDistorted

%==========================================================================
% Weak (Born) Fourier distorted and attenuated transmission hologram,
% solved numerically.
%==========================================================================
function out = DETransmissionBornFourierDistortedAndAttenuatedNumeric(in,out)

%--------------------------------------------------------------------------
% Create struct with required parameters
%--------------------------------------------------------------------------
P.kappa0L   = out.kappa0L;
P.Phi1      = in.Phi1;          P.A1    = in.A1;
P.Phi2      = in.Phi2;          P.A2    = in.A2;
P.Phi3      = in.Phi3;          P.A3    = in.A3;
P.Phi4      = in.Phi4;          P.A4    = in.A4;
P.Phi5      = in.Phi4;          P.A5    = in.A5;
P.Phi6      = in.Phi6;          P.A6    = in.A6;

%--------------------------------------------------------------------------
% For each Bragg mismatch
%--------------------------------------------------------------------------
for idKzL = 1:length(out.DeltaKzL)
    P.c0        = out.c0(idKzL);                % Parameter struct terms that change with deltaKzL
    P.c1        = out.c1(idKzL);
    P.DeltaKzL  = out.DeltaKzL(idKzL);

    % Solve coupled mode equations for R and I of Diff field.  BC is A1 = 0 + j 0
    [~, A] = ode45(@(zOverL,A) RealCMEBornFourierDistortedAndAttenuated(zOverL,A,P),[0 1],[0 0]);

    % Efficiency is |A1|^2 at end of zOverL grid
    out.Eta(idKzL) = A(end,1)^2 + A(end,2)^2;
end

end % DETransmissionBornFourierDistortedAndAttenuatedNumeric

%==========================================================================
% RHS of coupled mode equations with Fourier phase and amplitude distortion written in
% complex form in weak (Born) limit.
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitudes of diffrcted field A(1)+jA(2)
%   P           Struct of parameters
%==========================================================================
function dA = ComplexCMEBornFourierDistortedAndAttenuated(zOverL,~,P)

% Fourier series distortion of phase and amplitude
Phi = FourierSeries(zOverL,P.Phi1,P.Phi2,P.Phi3,P.Phi4,P.Phi5,P.Phi6);
Atn = FourierSeries(zOverL,  P.A1,  P.A2,  P.A3,  P.A4,  P.A5,  P.A6);

% RHS of coupled mode equations using complex representation
% Note A scales amplitude (kappa which is proportional to n1) with a fixed
% 0th order term "1" so that if A = 0, we have the usual equation
dA = 1i * P.kappa0L * (1 + Atn) * exp(+1i * P.DeltaKzL * zOverL - 1i * Phi) * (1) / P.c1;

end % ComplexCMEBornFourierDistortedAndAttenuated

%==========================================================================
% RHS of coupled mode equations in Born limit with Fourier phase and amplitude distortion written in
% real form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitude of diffracted A(1)+jA(2)
%   P           Struct of parameters
%==========================================================================
function dAReal = RealCMEBornFourierDistortedAndAttenuated(zOverL,AReal,P)

% Accept real and imaginary parts of A1 field, assemble into complex variable
AComplex  = AReal(1)+1i*AReal(2);

% Call CME in complex form to calculate two complex results
dAComplex = ComplexCMEBornFourierDistortedAndAttenuated(zOverL,AComplex,P);

% Separate complex results into two real numbers
dAReal    = [real(dAComplex); imag(dAComplex)];

end % RealCMEBornFourierDistortedAndAttenuated

%==========================================================================
% Strong Fourier series distorted and attenuated transmission hologram.  
%==========================================================================
function out = DETransmissionFourierDistortedAndAttenuatedNumeric(in,out)

%--------------------------------------------------------------------------
% Create struct with required parameters
%--------------------------------------------------------------------------
P.kappa0L   = out.kappa0L;
P.Phi1      = in.Phi1;          P.A1    = in.A1;
P.Phi2      = in.Phi2;          P.A2    = in.A2;
P.Phi3      = in.Phi3;          P.A3    = in.A3;
P.Phi4      = in.Phi4;          P.A4    = in.A4;
P.Phi5      = in.Phi4;          P.A5    = in.A5;
P.Phi6      = in.Phi6;          P.A6    = in.A6;

%--------------------------------------------------------------------------
% For each Bragg mismatch
%--------------------------------------------------------------------------
for idKzL = 1:length(out.DeltaKzL)
    P.c0        = out.c0(idKzL);
    P.c1        = out.c1(idKzL);
    P.DeltaKzL  = out.DeltaKzL(idKzL);

    % Solve coupled mode equations for R and I of Inc and Diff fields
    % This combo of options ('AbsTol',1e-14) and solver (ODE78) yields
    % reasonbable confidence bounds on the polynomial phase coefficients.
    % All other combinations yield the same fit values but orders of
    % magnitude higher confidence bounds.
    options = odeset('AbsTol',1e-14);
    [~, A] = ode78(@(zOverL,A) RealCMEFourierDistortedAndAttenuated(zOverL,A,P),[0 1],[1 0 0 0],options);

    % Efficiency is |A1|^2 at end of zOverL grid
    out.Eta(idKzL) = A(end,3)^2 + A(end,4)^2;
end

end % DETransmissionFourierDistortedAndAttenuatedNumeric

%==========================================================================
% RHS of coupled mode equations with Fourier series phase distortion written in
% complex form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   Amp         Complex amplitudes of incident A(1) and diffracted A(2) 
%   P           Struct of parameters
%==========================================================================
function dA = ComplexCMEFourierDistortedAndAttenuated(zOverL,A,P)

% Fourier series distortion of phase and amplitude
Phi = FourierSeries(zOverL,P.Phi1,P.Phi2,P.Phi3,P.Phi4,P.Phi5,P.Phi6);
Atn = FourierSeries(zOverL,  P.A1,  P.A2,  P.A3,  P.A4,  P.A5,  P.A6);

% RHS of coupled mode equations using complex representation
dA    = zeros(2,1);  
dA(1) = 1i * P.kappa0L*(1 + Atn) * exp(-1i * P.DeltaKzL * zOverL + 1i * Phi) * A(2) / P.c0;
dA(2) = 1i * P.kappa0L*(1 + Atn) * exp(+1i * P.DeltaKzL * zOverL - 1i * Phi) * A(1) / P.c1;

end % ComplexCMEFourierDistortedAndAttenuated

%==========================================================================
% RHS of coupled mode equations with polynomial phase distortion written in
% real form.  
% 
% Inputs
%   zOverL      Distance normalized by thickness, [0,1]
%   AReal       Real amplitudes of incident A(1)+jA(2) and diffracted A(3)+jA(4)
%   P           Struct of parameters
%==========================================================================
function dAReal = RealCMEFourierDistortedAndAttenuated(zOverL,AReal,P)

% Accept real and imaginary parts of two fields, assemble into two complex variables
AComplex  = [AReal(1)+1i*AReal(2); AReal(3)+1i*AReal(4)];

% Call CME in complex form to calculate two complex results
dAComplex = ComplexCMEFourierDistortedAndAttenuated(zOverL,AComplex,P);

% Separate two complex results into four real numbers
dAReal    = [real(dAComplex(1)); imag(dAComplex(1)); real(dAComplex(2)); imag(dAComplex(2))];

end % RealCMEFourierDistortedAndAttenuated

%==========================================================================
% Sinc defined as sinc(x)/x, not the matlab definition of sinc(pi x)/(pi x)
%==========================================================================
function result = sincNoPi(x)
    result = sinc(x/pi);
end

%==========================================================================
% Polynomial distortion of phase
%==========================================================================
function Phi = PhiPoly(zOverL,dPhi2,dPhi3,dPhi4,z04OverL)

% Normalization of quartic depends on z04OverL since two possible extrema
% for a 4th order polynomial
Norm4 = max([(1-2*z04OverL)^2/16,(z04OverL-1)^2 * z04OverL^2 / 4 ]);

% No linear term since that is equivalent to Bragg match shift
Phi = dPhi2 * 4 *          zOverL.*(1-zOverL) + ...
      dPhi3 * 12*sqrt(3) * zOverL.*(zOverL-.5).*(zOverL-1) + ...
      dPhi4 *              zOverL.*(zOverL-z04OverL).*(zOverL+z04OverL-1).*(zOverL-1)/Norm4;

end % PhiPoly

%==========================================================================
% Fourier series distortion of phase
%==========================================================================
function Phi = FourierSeries(zOverL,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6)

Phi = Phi1 * sin(2*pi*zOverL * 1 /2) + ...
      Phi2 * sin(2*pi*zOverL * 2 /2) + ...
      Phi3 * sin(2*pi*zOverL * 3 /2) + ...
      Phi4 * sin(2*pi*zOverL * 4 /2) + ...
      Phi5 * sin(2*pi*zOverL * 5 /2) + ...
      Phi6 * sin(2*pi*zOverL * 6 /2);

end % FourierSeries
