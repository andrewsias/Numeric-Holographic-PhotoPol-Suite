%==========================================================================
% Parse inputs for simulation of holographic photopolymers.
% Return various structures to organize inputs, outputs and coordinates
% systems.
%
% NEEDS: divisors from Symbolic Math toolbox
%
% History
% RRM   Nov 20 2019     Version 1.0
% RRM   Feb 3 2020      Version 1.1 Changed date format to avoid : character
% RRM   Feb 19, 2020 	Support for shrinakge/swelling displacement
% RRM   Sep 22, 2020    Version 2.1 Spectral scatter, start reflection, radical dynamics
% RRM   Dec 18, 2020    Version 2.2 Multiplexing, slanted gratings, read
%                                   color, new output options
% RRM   Sep  1, 2022    Version 2.3 deltan for perfect hologram to test distorted reading
%==========================================================================
function [in,out,crd] = ParseNHInputs(VariablesPassed)

crd.nano  = 10^-9;	crd.micro = 10^-6;	crd.milli = 10^-3;            

%--------------------------------------------------------------------------
% Set up to parse
%--------------------------------------------------------------------------
%---Parser validation functions
validPosScalar  = @(x) isnumeric(x) && isscalar(x) && (x >= 0);     % Scalar > 0
validScalar     = @(x) isnumeric(x) && isscalar(x);                 % Scalar
validEven       = @(x) isnumeric(x) && isscalar(x) && (x > 0) && 2*round(x/2)==x; % Even integers
validOdd        = @(x) isnumeric(x) && isscalar(x) && (x > 0) && 2*round(x/2)~=x; % Odd integers

ExpectedScatType= {'RayleighV','RayleighS','Isotropic','None'};                  % Scatter angular distributions
validScatType 	= @(x) any(validatestring(x,ExpectedScatType)); 	% Check types

%---Default save file names
InFName        = strcat("HoloPolymer In ", ...
    string(datetime('now','format','uuuu-MM-dd'' ''HH-mm-ss' )));   % Default input save name, compat with Windows
OutFName       = strcat("HoloPolymer Out ", ...
    string(datetime('now','format','uuuu-MM-dd'' ''HH-mm-ss' )));   % Default input save name

%--- Parser structure
p = inputParser;                                              

%------------------------------Parameter dictionary------------------------
% Default material is bisphenol A diacrylate in polyurethane matrix
% Material parameters from Izzy "dn infinity results 92019.ppt"
%
%                       PARAM        DEFAULT  VALIDITY            DEFINITION
%--------------------------------------------------------------------------
% Save files
addParameter(p,  'InSaveFile',       InFName, @isstring);       % Input save file name.  Set to "" to supress save.
addParameter(p, 'OutSaveFile',      OutFName, @isstring);       % Output save file name.  Set to "" to supress save.
addParameter(p,  'InLoadFile',            "", @isstring);       % Load previous save. Default = "" = don't load.

% Material parameters
addParameter(p,       'nmono',         1.545, validPosScalar);  % Monomer index
addParameter(p,          'nP',        1.5789, validPosScalar);  % Polymer index
addParameter(p,          'nM',         1.476, validPosScalar);  % Matrix index
addParameter(p,       'sigma',          0.14, validPosScalar);  % Full conv shrinkage
addParameter(p,       'phim0',        0.1884, validPosScalar);  % Formulation monomer volume fraction at 20 wt% rho = 1/(.2/1.146 + .8/1.064), phim0 = .2*rho/1.146
addParameter(p,          'mu',             0, validPosScalar);  % Attenuation coefficient, [1/m]
addParameter(p,           'Z',  10*crd.micro, validPosScalar);  % Sample thickness
addParameter(p,          'Rm',           0.4, validPosScalar);  % Rm = D_monomer K^2 / R_P. See Zhao and Mouroulis definition of R
addParameter(p,          'Rr',             0, validPosScalar);  % Rr = D_radical K^2 / R_T.  Analogous quantity for radicals. 
addParameter(p,      'deltan',           NaN, validPosScalar);  % If set, overides all writing physics.  Perfect grating is created with this deltan.  Used to test readout physics.

addParameter(p,     'Distort',             0,     @islogical); 	% Apply distortion model?

% Writing parameters
% Any parameter marked MUX can be a vector which will then cause multiplex
% writing.  If 2 or more such parameters are vectors, must be = length.
addParameter(p,         'tau',           2.5, validPosScalar);  % MUX Exposure time multiplied by mean polymerization rate / exposure
addParameter(p, 'ThetaIncWrt',        0.1346,     @isnumeric);  % MUX Writing angle of incident (reference) beam positive from normal, internal [rad]
addParameter(p, 'ThetaDifWrt',       -0.1346,     @isnumeric);  % MUX Writing angle of diffracted (object) beam positive from normal, internal [rad]
addParameter(p, 'MuxSchedule',             0, @islogical);      % Adjust vector tau to achieve equal DE for all writing conditions?
addParameter(p, 'TrackAllMux',             0, @islogical);      % Read DE of each multiplexed hologram at each time step?


addParameter(p,  'lambda0Wrt',  405*crd.nano, validPosScalar);  % Optical wavelength for writing
addParameter(p,   'BeamRatio',             1, validPosScalar);  % Diffracted (object) beam power relative to incident (reference)
addParameter(p,       'w0wrt', 100*crd.micro, validPosScalar);  % Waist radius of writing beams


% Reading parameters
addParameter(p,  'lambda0Red',  405*crd.nano, validPosScalar);  % Optical wavelength for reading - NOT DONE
addParameter(p,       'w0red',  20*crd.micro, validPosScalar);  % Waist radius of reading beams
addParameter(p,   'CalcBragg',             0, @islogical);      % Compute Bragg selectivity?
addParameter(p, 'BraggTSteps',             0, @isnumeric);      % Which time steps to calculate Bragg selectivity
                                                                % Post exposure = NMux*in.Nt+1, post dark = NMux*in.Nt+2, post flood = NMux*in.Nt+3
addParameter(p,  'BraggNulls',             2, validPosScalar);  % Angular extent of theta for Bragg readout relative to first null
addParameter(p,'NBraggAngles',            41, validOdd);        % Number of theta samples for Bragg readout.  Should be odd.

% Grid parameters
addParameter(p,          'Nx',          1024, validEven);       % Number of x cells.  Most efficient if 2^^N 
addParameter(p,          'Ny',          1024, validEven);       % Number of y cells
addParameter(p,          'Nz',            36, validEven);       % Number of z cells
addParameter(p,          'dx',  203*crd.nano, validPosScalar);  % x cell size (Nyquist in air)
addParameter(p,          'dy',  203*crd.nano, validPosScalar);  % y cell size (Nyquist in air)
addParameter(p,          'Nt',            20, validPosScalar);  % Number of simulation time steps for each multiplexed write

% Scatter parameters - New system in V2
addParameter(p, 'InitialHaze',          0.01, validPosScalar);  % Haze due to static scatter centers at time = 0
addParameter(p,  'NScatInVol',         36000, validPosScalar);  % N centers in simulation volume. Increase = more structure to scatter. 
addParameter(p,  'NScatInSub',             0, validPosScalar);  % N centers in first substrate 
addParameter(p,    'ScatType',        'None',  validScatType);  % Angular distribution of single scatterer
addParameter(p,    'NHaxeVsZ',             1,  validPosScalar); % Number of thickness to measure haxze.  1 = exit only.  Will be rounded up to divisor of Nz

addParameter(p, 'PSThreshold',            -1, validScalar);     % Threshold polymer value at onset of phase sep.  Norm to [m0].  Set <0 to turn off
addParameter(p,         'PS1',           0.1, validScalar);     % First term in power series expansion of derivative of normalized chemical potential of monomer wrt normalized polymer conc
addParameter(p,         'PS2',           0.0, validScalar);     % Second term in power series expansion of derivative of normalized chemical potential of monomer wrt normalized polymer conc
addParameter(p,         'PS3',           0.0, validScalar);     % Third term in power series expansion of derivative of normalized chemical potential of monomer wrt normalized polymer conc

% ---PARSE---
parse(p,VariablesPassed{:});                                    % Parse inputs into struct p

%--------------------------------------------------------------------------
% Load old or parse new inputs.  If name+value pairs are used along with an
% input load file, the new name+value pairs will overwrite loaded values or
% create new fields if the parsed parameter did not exist at the time the
% file was saved.  Messages to that effect are displayed.
%--------------------------------------------------------------------------
if strlength(p.Results.InLoadFile) > 0                      % Is there a load file name specified?
    load(p.Results.InLoadFile);                             % Load previous "in" struct
    for ivar = 1:2:length(VariablesPassed)                  % Overwrite or add fields
        if ~strcmpi(VariablesPassed{ivar},'InLoadFile')     % Don't do anything to parameter input load file
            if isfield(in,VariablesPassed{ivar})            % Does this field exsist in the loaded struct?
                disp(strcat("Loaded parameter ",VariablesPassed{ivar}," = "));
                disp(in.(VariablesPassed{ivar}));
                disp(strcat("Replaced with input parameter "));
            else
                disp(strcat("Input parameter ",VariablesPassed{ivar}," added to loaded parameter set with value"));    
            end % Did field exist in load file
            disp(VariablesPassed{ivar+1});                  % Write field to struct
            in.(VariablesPassed{ivar}) = VariablesPassed{ivar+1};
        end % is parameter input load file ?
    end % for all passed name+variable pairs
else                                                        % No load file.  Use parsed values   
    in = p.Results;                                         % Short hand notation "in"
end % is there a load file?
  
%--------------------------------------------------------------------------
% Additional parameters not under user control
%--------------------------------------------------------------------------
out.Version     = 2.2;                      % Version number
out.VersionDate = '16-Oct-2020';
out.HazeStopDeg = 2.5;                      % Half angle [deg] excluded for haze measurement

%--------------------------------------------------------------------------
% Calculated parameters
%--------------------------------------------------------------------------
out.n0    	= LorentzLorenz(in.nmono,in.phim0,in.nM,1-in.phim0);  % Initial index (real part)
out.nprime	= in.mu/(2*pi/in.lambda0Wrt); % Imaginary refractive index = attenuation coefficient / vacuum wave number

% Index contrast under ideal recording conditions (Rm->inf, Rr->0) and
% complete monomer consumption.  Twice original monomner concentration 
% converted to polynmer at bright fringe, no polymer at dark fringe.  No
% remaining monomner anywhere (all converted to polymer).
out.dnPred	= LorentzLorenz(in.nP,in.phim0*(1-in.sigma)*2, in.nM,1-in.phim0*(1-in.sigma)*2) - ...
              LorentzLorenz(in.nP,                      0, in.nM,                        1);

out.Vr      = 1/(1+in.Rr);                  % Steady state solution for harmonic series
    
%--------------------------------------------------------------------------
% Replicate values related to multiplexing so all are arrays of same
% length.  This simplifies access throughout rest of program.
%--------------------------------------------------------------------------
out.NBraggT = length(in.BraggTSteps);              % Number of Bragg selectivity times
out.NMux    = max([length(in.tau),length(in.ThetaIncWrt),length(in.ThetaDifWrt)]); % Number of writing mux conditions

if length(in.tau) == 1
    % If scheduling, calculate tau vector to give equal dn
    if in.MuxSchedule
        dnOverDn = 1 - exp(-in.tau);        % Fractional dynamic range, first exposure
        if dnOverDn * out.NMux > 1
            error(['ParseNHInputs: Requested scheduling cannot be satisifed with ',...
                num2str(dnOverDn),' fractional dynamic range per exposure']);
        else
            % Calculate vector of exposure times based on first-order kinetics
            in.tau = -log(1-dnOverDn*(2:(out.NMux+1))) + log(1-dnOverDn*(1:out.NMux));
        end
    else
        in.tau = in.tau .* ones(out.NMux,1);
    end
elseif length(in.tau) ~= out.NMux
    error('ParseNHInputs: Number of elements in tau must = 1 or other multiplex controls');
end

if length(in.ThetaIncWrt) == 1
    in.ThetaIncWrt = in.ThetaIncWrt .* ones(out.NMux,1);
elseif length(in.ThetaIncWrt) ~= out.NMux
    error('ParseNHInputs: Number of elements in ThetaIncWrt must = 1 or other multiplex controls');
end

if length(in.ThetaDifWrt) == 1
    in.ThetaDifWrt = in.ThetaDifWrt .* ones(out.NMux,1);
elseif length(in.ThetaDifWrt) ~= out.NMux
    error('ParseNHInputs: Number of elements in ThetaDifWrt must = 1 or other multiplex controls');
end

% If Bragg spectrum requested with no time
if in.CalcBragg && length(in.BraggTSteps) == 1 && in.BraggTSteps == 0      
    in.BraggTSteps = out.NMux * in.Nt + 1;  % ...choose end of illumination period
end
          
%--------------------------------------------------------------------------
% BPM controls
%--------------------------------------------------------------------------
crd.ABCMinT	= 0;                        % Atten at grid edge
crd.ABCWid  = 10;                       % Width of absorber in cells               
out.X      	= in.Nx*in.dx;              % Size of sim
out.Y      	= in.Ny*in.dy;

% What z indices to sample haze?
NzDivisors      = divisors(in.Nz);                          % Possible sub-sampling in z = divisors of Nz
HazeSubSamp     = NzDivisors(find(NzDivisors<=in.NHaxeVsZ,1,'last'));  % Largest divisor <= requested number
out.HazeZIndxs  = flip(in.Nz:-in.Nz/HazeSubSamp:1);         % Vector of locations.  Count down from Nz so Nz always included.

end % ParseInputs