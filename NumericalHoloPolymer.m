%==========================================================================
% HoloPolymer.  Model diffraction, scatter and polymerization dynamics of
% holographic photopolymer.
% 
% NEEDS:    Name                        Origin
%           -------------------         -----------------
%           LorentzLorenz               CU 
%           BPM                         CU
%           ParseNHInputs               CU
%           Signal Processing Toolbox   MatLab
%           Curve Fitting Toolbox       Matlab
% 
% INPUTS:  See parameter dictionary in ParseNHInputs. Not case sensitive.
%
% OUTPUTS:
% in        Struct with all input values including defaults
% out       Struct with all output values
% crd       Struct with useful coordinate information.
%
% Meaning of time index "it" is:
%   Start           End
%   1                           Initial condition
%   2               Nt+1    	Illumination for hologram 1
%   Nt+2            2*Nt+1      Illumination for hologram 2
%            ...
%   (NMux-1)*Nt+2   NMux*Nt+1   Illumination for hologram NMux
%   NMux*Nt+2                   Dark development,
%   NMux*Nt+3                   Flood cure
%
% VERSION HISTORY:
%-----V1.0----
% RRM   June 3  2019  	First version HoloPolymerHaze
% RRM   July 30 2019   	Function for more general analysis + inputParser
% RRM   Sept 2, 2019   	VERSION 1.0 release
%-----V1.1----
% RRM   Feb 19, 2020  	Post-exposure coordinate displacement 
%                       Modified file names for Windows compatibility
% RRM   Apr 07, 2020    Phase separation 
% RRM   Apr 15, 2020    Normalized variables tau, Rm, Rs   
% RRM   June 3, 2020    Fixed bug involving out.HazeToStop 
%                       Added Rayleigh scatter type
% RRM   Sept 9, 2020    Added tomographic readout option
%-----V2.0----
% RRM   Jul 14, 2020    Switched to BPM_V2 in prep for reflection
%-----V2.1----
% RRM  Sept 22, 2020    Combined various improvements made to V1.1 into 2.0
%                       Scatter now described in angular domain
%                       More control over optics inc tilted holograms - unfinished
%                       Radical kinetics (not just steady state)
%-----V2.2----
% RRM   Oct 16, 2020   Distortion during exposure 
%                      Multiplexing
%-----V2.3----
% RRM   Sep  1, 2022   Option to supress writing and just load hologram
%==========================================================================
function [in,out,crd] = NumericalHoloPolymer(varargin)

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
[in,out,crd]            = ParseNHInputs(varargin);          % Read user input
[out,crd]               = DefineCoordinates(in,out,crd);    % Fourier coords  
[in,out,crd,mat,inc] 	= Initialize(in,out,crd);       	% Allocate and initialize
out                     = CreateScatterers(in,out,crd);     % Set up scattering 

%--------------------------------------------------------------------------
% Time loop for normal operation including writing of hologram
%--------------------------------------------------------------------------
if out.Writing
    for it = 1:(out.NMux * in.Nt+3) 
    
        %------Illumination step----- 
        if it >1 && it <= out.NMux * in.Nt+1       
            [Ewrt_xyz,crd] 	= Diffract(it,in,out,inc,mat,crd);	% Diffract writing light
            mat             = React(Ewrt_xyz,in,mat,crd);       % Photochemical reactions
            mat             = Diffuse(in,out,mat,crd);       	% Transport
        
        %------Dark development step--------
        elseif it == out.NMux * in.Nt+2        
            mat             = DarkDevelop(mat);                 % Complete diffusion
    
        %------Flood cure step--------------        
        elseif it == out.NMux * in.Nt+3              
            mat             = FloodCure(mat);                 	% Complete polymerization
        end 
       
        %------Update and read hologram every step--------------
	    out       = UpdateIndex(it,in,out,mat,crd);             % New refractive index
        [out,crd] = ReadOut(it,in,out,inc,crd);                 % Measure results
    
    end % for it

%--------------------------------------------------------------------------
% If just reading an ideal hologram (possibly with scatter and distortion)
%--------------------------------------------------------------------------   
else
    out       = UpdateIndex(1,in,out,mat,crd);                  % Calculate refractive index
    [out,crd] = ReadOut(1,in,out,inc,crd);                      % Measure results
end % writing

%--------------------------------------------------------------------------
% Save files
%--------------------------------------------------------------------------
if strlength(in.InSaveFile)  > 0, save(in.InSaveFile,  'in'); end
if strlength(in.OutSaveFile) > 0, save(in.OutSaveFile,'out'); end

end % NumericalHoloPolymer

%--------------------------------------------------------------------------
%-------------------------------LOCAL METHODS------------------------------
%--------------------------------------------------------------------------

%==========================================================================
% Set up real and Fourier space coordinates
%==========================================================================
function [out,crd] = DefineCoordinates(in,out,crd)

out.dz      = in.Z/in.Nz;

% Real space 1, 2 and 3D grids
out.x_x     = ((-in.Nx/2) : (in.Nx/2-1))*in.dx; 
out.y_y     = ((-in.Ny/2) : (in.Ny/2-1))*in.dy; 
out.z_z     = (0          : (in.Nz-1))* out.dz; 

[crd.x_xy,crd.y_xy] = ndgrid(out.x_x,out.y_y);
[crd.x_xyz,crd.y_xyz,crd.z_xyz] = ndgrid(out.x_x,out.y_y,out.z_z);

% Fourier space 1 and 3D grids
out.fx_fx  	= ((-in.Nx/2) : (in.Nx/2-1))*1/(in.Nx*in.dx);
out.fy_fy 	= ((-in.Ny/2) : (in.Ny/2-1))*1/(in.Ny*in.dy);
fz_fz       = ((-in.Nz/2) : (in.Nz/2-1))*1/(in.Nz*out.dz);
[crd.fx_fxfyfz, crd.fy_fxfyfz, crd.fz_fxfyfz] = ndgrid(out.fx_fx,out.fy_fy,fz_fz);

end % DefineCoordinates

%==========================================================================
% Initialize variables
%==========================================================================
function [in,out,crd,mat,inc] = Initialize(in,out,crd)

%--------------------------------------------------------------------------
% Control flags
%--------------------------------------------------------------------------
out.Writing     = isnan(in.deltan);                 % Simulate writing process or just load perfect grating with specified dn?
out.Scattering  = ~strcmp(in.ScatType,'None');      % Are there scattering centers?

%--------------------------------------------------------------------------
% Initialize counters
%--------------------------------------------------------------------------
crd.iBraggT	= 0;            % Bragg time stamp counter, incremented on entry to AngleSpectrum
crd.iMux    = 1;            % Hologram number, updated in Diffract
mat.u_xyz   = 0;            % So that mat is defined (with no content)

%--------------------------------------------------------------------------
% Allocate arrays as needed.  
%--------------------------------------------------------------------------
if out.Scattering
    NHazeZ      = length(out.HazeZIndxs);           % Number of haze samples in depth
    out.Haze    = zeros(out.NMux*in.Nt+3,NHazeZ);	% Haze between out.HazeStopDeg and 90 deg
    out.HazeX 	= zeros(out.NMux*in.Nt+3,NHazeZ); 	% Haze in x slice 
    out.HazeY 	= zeros(out.NMux*in.Nt+3,NHazeZ); 	% Haze in y slice 
    out.HazeRef = zeros(out.NMux*in.Nt+3,1);     	% Haze measured with reference beam
    out.HazeXRef= zeros(out.NMux*in.Nt+3,1);     	% Haze in x slice measured with reference
end

if out.Writing                                      % Concentrations only used in writing
    out.Eta    	= zeros(out.NMux*in.Nt+3,out.NMux); % Diffraction efficiency vs. time for each hologram
    out.DnNum   = zeros(out.NMux*in.Nt+3,1);        % Largest peak to mean polymer index change

    out.m0Num  	= zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated monomer harmonic 0
    out.m1Num   = zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated monomer harmonic 1
    out.P0Num   = zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated polymer harmonic 0
    out.P1Num   = zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated polymer harmonic 1
    out.P2Num   = zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated polymer harmonic 2
    out.P3Num   = zeros(out.NMux*in.Nt+3,1);        % Numericaly calculated polymer harmonic 3
    
    mat.m_xyz 	=  ones(in.Nx,in.Ny,in.Nz);         % Initial normalized monomer concentration
    mat.P_xyz 	= zeros(in.Nx,in.Ny,in.Nz);         % Initial normalized polymer concentration
    if in.Rr ~= 0                                   % Tracking radical dynamics?
        mat.r_xyz   = zeros(in.Nx,in.Ny,in.Nz);     % Initial normalized radical concentration
    end
end

%--------------------------------------------------------------------------
% Spatial frequencies needed for haze calculations
%--------------------------------------------------------------------------
crd.dfx         = 1/(in.Nx*in.dx);          % Frequency spacings
crd.dfy         = 1/(in.Ny*in.dy);	

% Integration area for diffraction efficiency measurements
crd.iFxDet      = round(-log(.001)/(pi*in.w0red)/crd.dfx);     	% x index of half angle to capture 99.9% of read beam
crd.jFyDet      = round(-log(.001)/(pi*in.w0red)/crd.dfy);     	% y index of half angle to capture 99.9% of read beam

% Indices of spatial frequencies used in signal and haze calculations
crd.iFxHazeHA   = round(sind(out.HazeStopDeg)/in.lambda0Red/crd.dfx); 	% Frequency index half width of haze cone in fx 
crd.jFyHazeHA   = round(sind(out.HazeStopDeg)/in.lambda0Red/crd.dfy);  	% Frequency index half width of haze cone in fy 

crd.iFxHazeM90  = find(out.fx_fx*in.lambda0Red <= -1,1,'last');     % X freq index of -90 degrees in air   
if isempty(crd.iFxHazeM90), crd.iFxHazeM90 = 1; end                 % Set = 1 if-90 outside of simulation
crd.iFxHazeP90  = find(out.fx_fx*in.lambda0Red >=1,1,'first');      % X freq index of +90 degrees in air
if isempty(crd.iFxHazeP90), crd.iFxHazeP90 = in.Nx; end             % Set = Nx if +90 outside of simulation

crd.jFyHazeM90  = find(out.fy_fy*in.lambda0Red <= -1,1,'last');     % Y freq index of -90 degrees in air
if isempty(crd.jFyHazeM90), crd.jFyHazeM90 = 1; end                 %  Set = 1 if-90 outside of simulation
crd.jFyHazeP90  = find(out.fy_fy*in.lambda0Red >=1,1,'first');      % Y freq index of +90 degrees in air
if isempty(crd.jFyHazeP90), crd.jFyHazeP90 = in.Ny; end             % Set = Ny if +90 outside of simulation          

%--------------------------------------------------------------------------
% Normal incident field and power in Fourier space for haze measurement.  
%--------------------------------------------------------------------------
inc             = struct;                           % Incident field structure
inc.Enrm_xy0    = exp(-(crd.x_xy.^2+crd.y_xy.^2)/in.w0red^2);   % Gaussian E of read beam diameter
out.Pnrm        = sum(sum(abs( fftshift(fft2(fftshift(inc.Enrm_xy0))) ).^2));   

%--------------------------------------------------------------------------
% If distorting due to shrinkage, calculate normalized shrinkage distortion
%--------------------------------------------------------------------------
if in.Distort
    mat = IncompressibleU(crd,in,mat);      % Normalized distortion field mat.u_xyz
end

%--------------------------------------------------------------------------
% Configure variables that describe writing and reading of possibly
% multiplexed holograms.
%--------------------------------------------------------------------------
inc.Ewrt_xy0    = zeros(in.Nx,in.Ny,out.NMux);      % Writing and reading fields
inc.Ered_xy0    = zeros(in.Nx,in.Ny,out.NMux);
out.Pref        = zeros(1,out.NMux);                % Power in read beam (in Fourier space)

for iMux = 1:out.NMux                   % For all holograms...
           
    %----------------------------------------------------------------------
    % Bragg condition
    %----------------------------------------------------------------------
    out.Kvec(:,iMux)= 2*pi*(out.n0/in.lambda0Wrt) * ...   	% Grating wave vector
                      ([sin(in.ThetaDifWrt(iMux)),cos(in.ThetaDifWrt(iMux)),0] - ...
                       [sin(in.ThetaIncWrt(iMux)),cos(in.ThetaIncWrt(iMux)),0]);
    out.K(iMux)     = norm(out.Kvec(:,iMux));               % Magnitude of grating vector
    out.Lambda(iMux)= 2*pi/out.K(iMux);                     % Hologram period
    out.phi(iMux)   = acos(dot([0,-1,0],out.Kvec(:,iMux)/out.K(iMux)));   % Slant angle, 90 deg = unslanted

    % Bragg matched angle for reference beams.  
    % Find bisector angle of writing inc and diff k vecrtors which is orthogonal to K
    ThetaBisect = (in.ThetaIncWrt(iMux) + in.ThetaDifWrt(iMux))/2;
    % Usual Bragg calculation relative to bisector
    out.ThetaIncRed(iMux) = asin((in.lambda0Red/in.lambda0Wrt) * ...
               sin(in.ThetaIncWrt(iMux)-ThetaBisect))+ThetaBisect; 
	out.ThetaDifRed(iMux) = asin((in.lambda0Red/in.lambda0Wrt) * ...
               sin(in.ThetaDifWrt(iMux)-ThetaBisect))+ThetaBisect; 

    %----------------------------------------------------------------------
    % Calculate various spatial frequencies needed to write and read this
    % mux condition
    %----------------------------------------------------------------------
    % Spatial frequencies of incident fields and hologram
    FxIncWrt   	= sin(in.ThetaIncWrt(iMux))/(in.lambda0Wrt / out.n0);  	% X spatial frequency of incident   writing wave [1/m]
    FxDifWrt  	= sin(in.ThetaDifWrt(iMux))/(in.lambda0Wrt / out.n0);  	% X spatial frequency of diffracted writing wave [1/m]
    FxIncRed   	= sin(out.ThetaIncRed(iMux))/(in.lambda0Red / out.n0);	% X spatial frequency of incident   reading wave [1/m]
    FxDifRed  	= sin(out.ThetaDifRed(iMux))/(in.lambda0Red / out.n0);	% X spatial frequency of diffracted reading wave [1/m]

    crd.FxHolo(iMux) = FxDifWrt - FxIncWrt ;                           	% Spatial frequency of grating [1/m]

    % Convert spatial frequencie of diffracted read beam to array indices for convenience
    crd.iFxDifRed(iMux)	= in.Nx/2+1 + round(FxDifRed/crd.dfx);          % Index of center of reading dif (object) wave 

    %----------------------------------------------------------------------
    % Incident fields for writing and reading.
    % Center transverse profile at Z/2, not front surface.
    %----------------------------------------------------------------------
    inc.Ewrt_xy0(:,:,iMux)	= (exp(-((crd.x_xy+in.Z/2*tan(in.ThetaIncWrt(iMux))).^2+crd.y_xy.^2)/in.w0wrt^2) ...
                        .* exp(+1i*2*pi*FxIncWrt*crd.x_xy) + ...        % Write reference beam
                           exp(-((crd.x_xy+in.Z/2*tan(in.ThetaDifWrt(iMux))).^2+crd.y_xy.^2)/in.w0wrt^2) ...
                        .* exp(1i*2*pi*FxDifWrt*crd.x_xy)) / sqrt(2);   % Write object beam   

    inc.Ered_xy0(:,:,iMux)	= exp(-((crd.x_xy-in.Z/2*tan(out.ThetaIncRed(iMux))).^2+crd.y_xy.^2)/in.w0red^2) ...
                       .* exp(+1i*2*pi*FxIncRed*crd.x_xy);              % Read reference beam

    Ered_fxfy           = fftshift(fft2(fftshift(inc.Ered_xy0(:,:,iMux))));   	% Read field in Fourier space
    out.Pref(iMux)    	= sum(sum(abs(Ered_fxfy).^2));                  % Ref beam power in Fourier space

end % For all multiplex conditions

%--------------------------------------------------------------------------
% Fick's law transfer function.  See "Diffusion model of polymer
% development.ppt"  Normalized spatial transfer function is H(f,t) =
% exp(-(f/f0)^2) where f0 = 1/(2 pi sqrt(D t)).  Here t = one time step dt
% and f magnitude of 3D spatial frequency vector.
% Using definitions from "Harmonic analysis.ppt", 
% Note same fFick for monomer and radicals because we assume same D.
% Reaction rates can be different (see React)
% Note that D is same for mux'ed holograms so fFick and thus H is same
% Must = 1 at f = 0 to conserve mass
% Must calculate each mux condition since dtau = in.tau(iMux)/in.Nt changes
%--------------------------------------------------------------------------
if out.Writing && (in.Rm ~= 0 || in.Rr ~= 0)            % Diffusion happening?
    fFick           = 1/(out.Lambda(1)*sqrt(in.Rm*in.tau(1)/in.Nt)); 
    mat.H_fxfyfz    = exp(-(crd.fx_fxfyfz.^2+crd.fy_fxfyfz.^2+crd.fz_fxfyfz.^2)/fFick^2);  % Transfer func
end   

end % Initialize

%==========================================================================
% Create scatter distribution. 
% Updated for V1.1 to allow scatter centers to be added via phase sep.
%==========================================================================
function out = CreateScatterers(in,out,crd)

if out.Scattering               % Are there any scatterers?
    NsX = 8;                    % Transverse grid size for scatter center
    NsY = 8;                    % Should this be measured from Escat_xy?

    %--------------------------------------------------------------------------
    % Define scatter distribution in far field. 
    %--------------------------------------------------------------------------
    switch in.ScatType
        case 'RayleighV'                % Vector              
        % E proportional to cos(alpha) where alpha is angle out of plane of
        % polarization.  In other words, inner product of inc pol and radiated
        % pol.  Hologram K in xz plane, so E polarization is y. alpha thus is
        % arcsin(ky/k).  cos(arcsin(ky/k)) = sqrt(1-(ky/k)^2)
            ESingleScat_fxfy = sqrt(1- ((in.lambda0Wrt/out.n0) * crd.fy_fxfyfz(:,:,1)).^2);
        case 'RayleighS'                % Scalar
            ESingleScat_fxfy = sqrt(1- ((in.lambda0Wrt/out.n0) * (crd.fx_fxfyfz(:,:,1)).^2+crd.fy_fxfyfz(:,:,1)).^2);
        case 'Isotropic'
            ESingleScat_fxfy = ones(in.Nx,in.Ny);
    end

    % Evanescent fields = 0
    ESingleScat_fxfy((crd.fx_fxfyfz(:,:,1).^2 + crd.fy_fxfyfz(:,:,1).^2) > ...
        (out.n0/mean(in.lambda0Wrt))^2) = 0;

    %--------------------------------------------------------------------------
    % Normalize to meet initial scatter requirement
    %--------------------------------------------------------------------------
    % How must power of single scatterer be scaled to yield desired initial haze
    IntensScaling = in.InitialHaze / CalcHaze(ESingleScat_fxfy,in,out,crd) / ...
        (in.NScatInVol + in.NScatInSub);

    % This scatterer is on axis so intersects probe beam.  Probe beam is small
    % so must increase strength of scatterers so finite area probed by haze
    % beam has desired scatter.  That is, NScat will not be seen by probe.
    IntensScaling = IntensScaling * (out.X * out.Y) / (in.w0red^2*pi/2);  % Increase by ratio of areas

    ESingleScat_fxfy = ESingleScat_fxfy * sqrt(IntensScaling);  % Scale scatter strength

    % ======NOTE: Initial haze is off by factor < 2.  Don't see why===========

    %IFT to find spatial distribution of E corresponding to this angular dist
    ESingleScat_xy  = ifftshift(ifftn(ifftshift(ESingleScat_fxfy)));

    %--------------------------------------------------------------------------
    % Calculate dn for a scatter center from this desired scattered field
    % dn relates to scattered field through series expansion of real-space
    % propagator exp(i k0 dn dz) = 1 + i k0 dn dz + O(dn^2) = Einc + Escat.
    % Therefore dn = Escat / (k0 dz)
    %--------------------------------------------------------------------------
    dnS_xy = ESingleScat_xy(in.Nx/2+1+((-NsX/2):(NsX/2-1)),in.Ny/2+1+((-NsY/2):(NsY/2-1))) / ...
            (2*pi/mean(in.lambda0Wrt) * out.dz);

    %--------------------------------------------------------------------------
    % Create 3D index array
    %--------------------------------------------------------------------------
    NScatPerLayer = round(in.NScatInVol/(in.Nz-1),0);           % Scatter centers / BPM z step

    for iz = 1:in.Nz
        if iz == 1                                              % First layer scattering = substrate
            XInds   = round(1+rand(in.NScatInSub,1)*(in.Nx-1)); % Random coordinates
            YInds   = round(1+rand(in.NScatInSub,1)*(in.Ny-1));
        else                                                    % Volume layers
            XInds   = round(1+rand(NScatPerLayer,1)*(in.Nx-1)); % Random coordinates
            YInds   = round(1+rand(NScatPerLayer,1)*(in.Ny-1));
        end

        dn_xy   = zeros(in.Nx,in.Ny);
        dn_xy(sub2ind(size(dn_xy),XInds,YInds)) = 1;     	    % Binary, 1 at scat loc
        dn_xy   = conv2(dn_xy,dnS_xy,'same');                   % Add shape of scatterers     
        out.dnS_xyz(:,:,iz)  = dn_xy;                           % Accumulate in final array
    end

else                                                            % No scattering.  
    out.dnS_xyz	= 0;
end % if scattering

end % CreateScatterers

%==========================================================================
% Diffract writing beam.  Use incident fields and material descriptions
% passed through structures to calcuclate E fields in material at this t.
% Note that index has three parts:
%           1) homogeneous, 2) dynamic dn, 3) static (scattering) dn
%==========================================================================
function [Ewrt_xyz,crd] = Diffract(it,in,out,inc,mat,crd)

% Choose writing condition from set of NMux possible
crd.iMux = max(1,min(out.NMux,floor((it-2)/in.Nt)+1));  	
    
%--------------------------------------------------------------------------
% Diffract light with BPM
%--------------------------------------------------------------------------
[~, BPMout] = BPM_v2(inc.Ewrt_xy0(:,:,crd.iMux),out.z_z,'dx',in.dx,'dy',in.dy,...
            'lambda0',in.lambda0Wrt,'nTrans',(out.n0+1i*out.nprime)+out.dn_xyz+out.dnS_xyz,...
            'ABCMinTrans',crd.ABCMinT,'ABCCellWidth',crd.ABCWid,...
            'RetXInd1',1:in.Nx,'RetYInd1',1:in.Ny,'RetZInd1',1:in.Nz);
        
Ewrt_xyz = BPMout.Eret1;

%--------------------------------------------------------------------------
% If modelling distortion, warp diffracted fields in the actual distorted
% coordinate system into the undistorted Cartesian coordinate system used
% for reaction and diffusion.  Note that sign of dVOverV is opposite here
% to that used in UpdateIndex.  This is the opposite transformation.
% Note swap of x and y coordinate in call to interp3, as requried by
% assumed mesghrid indexing (y,x,z). Finally note that only shrinkage is
% modelled (no swelling) and flood cure at it = NMux * Nt + 3 is assumed to
% release radial stress, thus returning material to original condition.
%--------------------------------------------------------------------------
if in.Distort && it >= 2 && it <= out.NMux * in.Nt+2    % Distortion during exposure & dark cure
	dVOverV = +in.phim0 * in.sigma * out.P0Num(it);     % Negative of volume shrinkage due to mean polymerization           
   
    Ewrt_xyz = interp3(out.y_y,out.x_x,out.z_z,Ewrt_xyz,...
        crd.y_xyz-dVOverV*squeeze(mat.u_xyz(2,:,:,:)),...
        crd.x_xyz-dVOverV*squeeze(mat.u_xyz(1,:,:,:)),...
        crd.z_xyz-dVOverV*squeeze(mat.u_xyz(3,:,:,:)),'cubic',0);
end

end % Diffract

%==========================================================================
% React species.
% For improved accuracy of time stepping, use solutions to first order
% equations, not simple forward finite differences.  
% For no radical diffusion:
%   Calc conversion from dm/dtau = -m I where I is normalized intensity
%   Solution is m(tau+dtau) = c exp(-I (tau+dtau)) where
%   c = m(tau) exp(+I tau) so m(tau+dtau) = m(tau) exp(-I dtau).
%   Conversion = m(tau) - m(tau+dtau) = m(tau)(1-exp(-I dtau))
% else
%   Normalized radical reaction equation w/o diffusion is
%   dr/dtau = (Rm/Rr) ( -r + I )  (See Harmonic analysis.ppt)
%   Solution is r(tau+dtau) = c exp(-(Rm/Rr)(tau+dtau)) + I
%   where c = (r(tau)-I)exp(+(Rm/Rr)tau) so
%   r(tau+dtau) = (r(tau)-I) exp(-(Rm/Rr)dtau) + I 
%               = r(tau) exp(-(Rm/Rr)dtau) + I (1-exp(-(Rm/Rr)dtau))
%==========================================================================
function [mat] = React(Ewrt_xyz,in,mat,crd)

dtau        = in.tau(crd.iMux)/in.Nt;             % Time step size for convenience

% Calculate monomer conversion this time step
if in.Rr == 0           % No radical diffusion
    Conv_xyz = mat.m_xyz .* (1 - exp(-dtau * abs(Ewrt_xyz).^2));  % conversion this time step
else                    % Radical diffusion
    % Radical dynamics
    mat.r_xyz   = mat.r_xyz * exp(- in.Rm/in.Rr * dtau) + ...
                  abs(Ewrt_xyz).^2 * (1-exp(- in.Rm/in.Rr * dtau));

    % Then calculate conversion as above with r in place of I
    Conv_xyz = mat.m_xyz .* (1 - exp(-dtau * mat.r_xyz));  % conversion this time step
end

% Update monomer and polymer
mat.m_xyz 	= mat.m_xyz - Conv_xyz;         % Conversion decreases monomer
mat.P_xyz  	= mat.P_xyz + Conv_xyz;         % and increases polymer
    
end % React

%==========================================================================
% Diffusion of mobile species  
% Concentration driven transport via convolution with Fick's law impulse response
%==========================================================================
function [mat] = Diffuse(in,out,mat,crd)
  
dtau        = in.tau(crd.iMux)/in.Nt;       % Time step for convenience

%--------------------------------------------------------------------------
% Concentration driven diffusion 
%--------------------------------------------------------------------------
if in.Rm ~= 0
    m_fxfyfz    = fftshift(fftn(fftshift(mat.m_xyz)));      % monomer
    m_fxfyfz    = m_fxfyfz .* mat.H_fxfyfz;
    mat.m_xyz   = abs(ifftshift(ifftn(ifftshift(m_fxfyfz)))); % force real
end

if in.Rr ~= 0
    r_fxfyfz    = fftshift(fftn(fftshift(mat.r_xyz)));      % radicals
    r_fxfyfz    = r_fxfyfz .* mat.H_fxfyfz;
    mat.r_xyz   = abs(ifftshift(ifftn(ifftshift(r_fxfyfz)))); % force real
end

%--------------------------------------------------------------------------
% Miscibility driven diffusion of monomer via direct FD solution
%--------------------------------------------------------------------------
if in.PSThreshold >= 0                         % Flag to invoke

    % Note that the gradient command assumes that the target array is
    % produced by meshgrid and thus indexed (y,x,z).  Since this code
    % indexes arrays (x,y,z), we must reverse the order of the output
    % variables.
    [dPdy_xyz,dPdx_xyz,dPdz_xyz] = gradient(mat.P_xyz,in.dy,in.dx,out.dz); 

    mat.m_xyz = mat.m_xyz + dtau * in.Rm/(2*pi)^2 * out.Lambda^2 * DuDp(mat,in) .* ... 
        (FiniteDiffCentral(mat.m_xyz .* dPdx_xyz,1,in.dx) + ...
         FiniteDiffCentral(mat.m_xyz .* dPdy_xyz,2,in.dy) + ...
         FiniteDiffCentral(mat.m_xyz .* dPdz_xyz,3,out.dz));
end
    
end % Diffuse

%==========================================================================
% Dark development of hologram.
% Currently only considers continued diffusion, not continuted reaction due
% to radical kinetics.
%==========================================================================
function [mat] = DarkDevelop(mat)

mat.m_xyz(:,:,:) = mean(mat.m_xyz,'all');   % Monomer becomes uniform

end

%==========================================================================
% Flood cure.
%==========================================================================
function [mat] = FloodCure(mat)

mat.P_xyz           = mat.P_xyz + mat.m_xyz; 	% All monomer -> polymer
mat.m_xyz(:,:,:)	= 0 ;                       % No monomer remains

end
%==========================================================================
% Update refractive index from constituent concentrations
%==========================================================================
function out = UpdateIndex(it,in,out,mat,crd)

%--------------------------------------------------------------------------
% If user provided a deltan value, simply fill space with a uniform grating
% of this strength.
%--------------------------------------------------------------------------
if ~out.Writing                                 % Simulating writing physics?
    
    % Fill entire grid with perfect hologram(s)
    out.dn_xyz = zeros(size(crd.x_xyz));        % Initialize
    for iMux = 1:out.NMux                       % For all holograms...
        out.dn_xyz = out.dn_xyz + in.deltan * cos(out.Kvec(1,iMux).*crd.x_xyz + ...
            out.Kvec(2,iMux).*crd.y_xyz + out.Kvec(3,iMux).*crd.z_xyz);
    end

else
    %----------------------------------------------------------------------
    % Use Lorentz Lorenz equation to calcualte index from constituents
    %----------------------------------------------------------------------
    out.dn_xyz = LorentzLorenz(in.nmono,in.phim0*mat.m_xyz,...
                 in.nP,in.phim0*(1-in.sigma)*mat.P_xyz, ...
                 in.nM,1-in.phim0*mat.m_xyz-in.phim0*(1-in.sigma)*mat.P_xyz)... 
                 - out.n0;   % Subtract background to get delta n     
    
    %----------------------------------------------------------------------
    % Fit chemical constituents and index to a harmonic series for comparison 
    % to theory and for use in shrinkage calculation
    %----------------------------------------------------------------------
    [out.m0Num(it), out.m1Num(it)] = ...
            fitHarmonics(in,out,out.Lambda(crd.iMux),mat.m_xyz(:,in.Ny/2+1,1));
    [out.P0Num(it), out.P1Num(it), out.P2Num(it), out.P3Num(it)] = ...
            fitHarmonics(in,out,out.Lambda(crd.iMux),mat.P_xyz(:,in.Ny/2+1,1));
    [ ~, out.DnNum(it)] = ...
	    fitHarmonics(in,out,out.Lambda(crd.iMux),out.dn_xyz(:,in.Ny/2+1,1));  
end

%--------------------------------------------------------------------------
% If modelling distortion, warp index in the undistorted Cartesian 
% coordinate system used for reaction and diffusion into the actual distorted
% coordinate system used for diffraction.  Note that sign of dVOverV is opposite here
% to that used in Diffract.  This is the opposite transformation.
% Note swap of x and y coordinate in call to interp3, as requried by
% assumed mesghrid indexing (y,x,z). Finally note that only shrinkage is
% modelled (no swelling) and flood cure at it = NMux * Nt + 3 is assumed to
% release radial stress, thus returning material to original condition.
%
% RRM 9/22: Wait, wouldn't it be better to warp the constituent
% concentrations so they dephase with intensity rather than just dn?
%--------------------------------------------------------------------------
if in.Distort && ((it >= 2 && it <= out.NMux * in.Nt+2) || ~out.Writing)    % If using distortion and either within writing time or using ideal grating

    if out.Writing
	    dVOverV = -in.phim0 * in.sigma * out.P0Num(it);     % Volume shrinkage due to current mean polymerization    
    else
        dVOverV = -in.phim0 * in.sigma;                     % Volume shrinkage due to complete polymerization    
    end

    out.dn_xyz = interp3(out.y_y,out.x_x,out.z_z,out.dn_xyz,...
        crd.y_xyz-dVOverV*squeeze(mat.u_xyz(2,:,:,:)),...
        crd.x_xyz-dVOverV*squeeze(mat.u_xyz(1,:,:,:)),...
        crd.z_xyz-dVOverV*squeeze(mat.u_xyz(3,:,:,:)),'cubic',0);
end

end % UpdateIndex

%==========================================================================
% DuDp: derivative of normalized chemical potential wrt normalized polymer
% concentration from power series terms provided by user
%==========================================================================
function DuDp_xyz = DuDp(mat,in)

DuDp_xyz = zeros(in.Nx,in.Ny,in.Nz);

if in.PS1 ~= 0
    DuDp_xyz = in.PS1;
end
if in.PS2 ~= 0
    DuDp_xyz = (mat.P_xyz - in.PSThreshold) * in.PS2;
end
if in.PS3 ~= 0
    DuDp_xyz = (mat.P_xyz - in.PSThreshold).^2 * in.PS3;
end

DuDp_xyz = DuDp_xyz .* (mat.P_xyz > in.PSThreshold);        % Threshold

end % DuDp

%==========================================================================
% FiniteDiffCentral - Use central finite difference to
% approximate first derivative.  Use circular shift command to keep array same
% size.  Appropriate for periodic BCs.
%==========================================================================
function result = FiniteDiffCentral(mat,dim,h)
    result = (circshift(mat,-1,dim) -  circshift(mat,+1,dim))/(2*h);
end % FiniteDiffCentral

%==========================================================================
% Read out = record various properties of system
%==========================================================================
function [out,crd] = ReadOut(it,in,out,inc,crd)
    
%--------------------------------------------------------------------------
% Bragg selectivity vs. angle measurment at end of exposure or at any time
% if just reading out ideal grating
%--------------------------------------------------------------------------
if in.CalcBragg && (ismember(it,in.BraggTSteps) || ~out.Writing)
   [out,crd] = AngleSpectrum(in,out,crd);           
end
        
%--------------------------------------------------------------------------
% Bragg matched holographic efficiency at each time step
% Return y=0 slice for later rendering
% Note that this is the Bragg matched efficiency so is reasonably accurate
% even for small read beams.  Thus we do not need the skewing calculation
% used in angular spectrum reading.
%--------------------------------------------------------------------------
if in.TrackAllMux               
    StartiMux = 1;          % Measure efficiency of all holograms
else
    StartiMux = crd.iMux;   % ...or measure efficiency only at latest hologram
end

for iMux = StartiMux:crd.iMux         % For all requested writing conditions
    % Diffract specified reading beam
    [~, BPMout]  = BPM_v2(inc.Ered_xy0(:,:,iMux),out.z_z,'dx',in.dx,'dy',in.dy,...
                'lambda0',in.lambda0Red,'nTrans',(out.n0+1i*out.nprime)+out.dn_xyz+out.dnS_xyz,...
                'ABCMinTrans',crd.ABCMinT,'ABCCellWidth',crd.ABCWid,...
                'RetXInd1',1:in.Nx,'RetYInd1',1:in.Ny,  'RetZInd1',in.Nz,...
                'RetXInd2',1:in.Nx,'RetYInd2',in.Ny/2+1,'RetZInd2',1:in.Nz);

    out.Eref_x0z = BPMout.Eret2;
    out.Eref_xyZ = BPMout.Eret1;
    Ered_fxfyZ   = fftshift(fft2(fftshift(out.Eref_xyZ)));  

    % Diffracted signal measured at spatial frequency centered on output for
    % this writing (mux) condition through crd.iFxDifRed 
    out.Eta(it,iMux) = sum(sum(abs(...
        Ered_fxfyZ(crd.iFxDifRed(iMux)+(-crd.iFxDet:crd.iFxDet), ...
        in.Ny/2+1+(-crd.jFyDet:crd.jFyDet))).^2)) / out.Pref(iMux);  
end % mux condition

%--------------------------------------------------------------------------
% Haze measurement if simulating scattering
%--------------------------------------------------------------------------  
if out.Scattering
    % "Haze" measurement with reference beam readout.  Due to strong
    % diffraction and off-axis read, the only meaningful value is Haze - HazeX
    % which blocks light in a strip centered on fy = 0 
    % Note that this uses the current hologram (mux) reading beam
    [out.HazeRef(it), out.HazeXRef(it),~] = CalcHaze(Ered_fxfyZ,in,out,crd);

    % True haze measurment with normally incident beam
    [~, BPMout]  = BPM_v2(inc.Enrm_xy0,out.z_z,'dx',in.dx,'dy',in.dy,...
                    'lambda0',in.lambda0Red,'nTrans',(out.n0+1i*out.nprime)+out.dn_xyz+out.dnS_xyz,...
                    'ABCMinTrans',crd.ABCMinT,'ABCCellWidth',crd.ABCWid,...
                    'RetXInd1',1:in.Nx,'RetYInd1',1:in.Ny,  'RetZInd1',out.HazeZIndxs);
    
    out.Enrm_xyHazeZ = BPMout.Eret1;            % Extract requested fields
    
    for iz = 1:length(out.HazeZIndxs)           % For all requested thicknesses
        [out.Haze(it,iz), out.HazeX(it,iz), out.HazeY(it,iz)] = ...
            CalcHaze(fftshift(fft2(fftshift(out.Enrm_xyHazeZ(:,:,iz)))),in,out,crd);
    end
end % if scattering
  
end % Readout

%==========================================================================
% Fit harmonics to cos modulated Gaussian to accurately estimate properties of a
% hologram.  
%
% Fit to raised cosine of finite contrast 
%
% Returns:      An = Harmonics proportional to Gaussian intensity.
% NEEDS:        Curve Fitting Toolbox
% To validate:  figure;plot(f1,out.x_x',slice_x)
%==========================================================================
function [ZeroHarmonic, OneHarmonic, A2, A3] = fitHarmonics(in,out,Lambda,slice_x)

FracWidth = 0.5;                                % Fraction of slice to fit
ixLow     = round(in.Nx*(1-FracWidth)/2);       % Low index
ixHi      = round(in.Nx*(1+FracWidth)/2);       % High index

LambdaStr = num2str(Lambda);                    % Convert to string

EqnForm = ['C0 + C1 * cos(2*pi*x/',LambdaStr,') + exp(-(x/w0IntStr)^2) * (A0 + A1 * cos(2*pi*x/',LambdaStr,') ', ...
   ' + A2 * cos(2*pi*x*2/',LambdaStr,') + A3 * cos(2*pi*x*3/',LambdaStr,'))'];


% For speed and to guarantee convergence, guess initial values
estA0 = min(min(min(slice_x)));                 % Background shift 
estA1 = (max(max(max(slice_x))) - estA0)/2;     % First harmonic
estw0 = in.w0wrt/sqrt(2);                       % sqrt(2) = Intensity

% Fit
f1 = fit(out.x_x(ixLow:ixHi)',squeeze(slice_x(ixLow:ixHi)),EqnForm,'Start', ...
    [0 0 estA0 estA1 0 0 estw0]);

ZeroHarmonic = f1.C0 + f1.A0;
OneHarmonic  = f1.C1 + f1.A1;
A2 = f1.A2;
A3 = f1.A3;

end % fitHarmonics

%==========================================================================
% Calculate haze (sum of scattered p in all SF except near origin) / inc P
% Also calc haze in X and Y slices such that haze - hazeX is haze with
% linear stop in x with width = haze stop.
%==========================================================================
function [Haze, HazeX, HazeY] = CalcHaze(E_fxfy,in,out,crd)

% Efficiency of scattering into on axis stop.  Approximated as square.
HazeToStop  = sum(sum(abs( ...
    E_fxfy(in.Nx/2+1+(-crd.iFxHazeHA:crd.iFxHazeHA), ...
	in.Ny/2+1+(-crd.jFyHazeHA:crd.jFyHazeHA))).^2))/out.Pnrm;      
          
Haze = sum(sum(abs( ...
    E_fxfy(crd.iFxHazeM90:crd.iFxHazeP90,crd.jFyHazeM90:crd.jFyHazeP90)).^2))/out.Pnrm ...
    - HazeToStop;     

HazeX  = sum(sum(abs( ...
    E_fxfy(crd.iFxHazeM90:crd.iFxHazeP90,in.Ny/2+1+(-crd.jFyHazeHA:crd.jFyHazeHA))).^2))/out.Pnrm ...
    - HazeToStop; 
HazeY  = sum(sum(abs( ...
    E_fxfy(in.Nx/2+1+(-crd.iFxHazeHA:crd.iFxHazeHA),crd.jFyHazeM90:crd.jFyHazeP90)).^2))/out.Pnrm ...
    - HazeToStop; 

end % CalcHaze 

%==========================================================================
% AngleSpectrum - Bragg selectivity calculation
%==========================================================================
function [out,crd] = AngleSpectrum(in,out,crd)

crd.iBraggT = crd.iBraggT + 1;          % Handle multiple measurement times

%--------------------------------------------------------------------------
% Bragg selectivity setup
%-------------------------------------------------------------------------- 
if crd.iBraggT == 1

    % Probe angles that span all read angles plut requested Bragg selectivity  
    HoloDTheta              = max(out.Lambda) / in.Z;   % Angular selectivity peak to first null
    out.HoloThetaRed_theta  = linspace(min(out.ThetaIncRed)-HoloDTheta*in.BraggNulls,...
                                       max(out.ThetaIncRed)+HoloDTheta*in.BraggNulls,...
                                       in.NBraggAngles); % Incident angles for reading

    % Spatial frequency range of detector in fx grid coordinates.
    % For multiplexed holograms, extend detector area to cover possible outputs
    crd.iFxDifDet = ceil((max(crd.iFxDifRed) - min(crd.iFxDifRed))/2) + crd.iFxDet;

    % Allocate space to record diffracted fields at exit (Z) during Bragg readout
    % Detected frequencies span diffracted frequencies of all mux cases plus
    % detector bandwidth calculated via size of beam 
    crd.SubSamp = 4;                            % Subsample interpolation for improved accuracy
    out.Ediff_fxfyZtheta = zeros(2*crd.SubSamp*crd.iFxDifDet+1,...
                                 2*crd.SubSamp*crd.jFyDet+1,...
                                 in.NBraggAngles,out.NBraggT); % Diffracted fields
    % Grating spatial frequencies used in small read beam correction & tomography                         
    out.Fx_fxfyZtheta    = zeros(2*crd.SubSamp*crd.iFxDifDet+1,...
                                 2*crd.SubSamp*crd.jFyDet+1,...
                                 in.NBraggAngles);
    out.Fy_fxfyZtheta    = out.Fx_fxfyZtheta;
    out.Fz_fxfyZtheta    = out.Fx_fxfyZtheta;
    out.Eta_theta        = zeros(in.NBraggAngles,out.NBraggT);
    out.Eta_Fz           = zeros(in.NBraggAngles,out.NBraggT);
end

%--------------------------------------------------------------------------
% Angle sweep.  Since the read beam is small and the hologram is large, the
% uncertainly P(kx,ky,kz) that is sampled by the diffracted k surface is
% sheared by the finite kz(kx) of the incident beam.  This (e.g.) fills in 
% the nulls upon Bragg readout of an otherwise ideal hologram. To simulate an
% infinite hologram, this shear must be removed by calculating the
% diffraction efficiency in a corrected (unsheared) k space before
% integrating.  That is, one cannot just integrate the total power at each
% incident angle but must instead calculate the power diffracted vs. k,
% then sum along transverse k.
%-------------------------------------------------------------------------- 
for itheta = 1:in.NBraggAngles
    % Spatial frequency of center of incident read beam.  
    fxInc = sin(out.HoloThetaRed_theta(itheta))/(in.lambda0Red / out.n0);
    
    % Incident probe beam tilted by this spatial frequency and centered on z = Z/2
    Ered_xy0    = exp(-((crd.x_xy-in.Z/2*tan(out.HoloThetaRed_theta(itheta))).^2+...
                             crd.y_xy.^2)/in.w0red^2) .* exp(+1i*2*pi*fxInc*crd.x_xy);  % Tilted incident beam for reading

    % Diffract off of hologram          
    [~, BPMout]  = BPM_v2(Ered_xy0,out.z_z,'dx',in.dx,'dy',in.dy,...
        'lambda0',in.lambda0Red,'nTrans',(out.n0+1i*out.nprime)+out.dn_xyz+out.dnS_xyz, ...
        'ABCMinTrans',crd.ABCMinT,'ABCCellWidth',crd.ABCWid,...
        'RetXInd1',1:in.Nx,'RetYInd1',1:in.Ny,'RetZInd1',in.Nz);          
    
    % Calculate the transverse angular spectrum at the output face z = Z
    Ered_fxfyZ      = fftshift(fft2(fftshift(BPMout.Eret1)));   
    
    % Calcualte central spatial frequency center of diffracted beam off Bragg
    % When multiplexing, use mean hologram spatial frequency
	fxDiff          = fxInc + mean(crd.FxHolo);        

    % Extract peak of diffracted fields.  I tried to just sample the grid
    % and keep track of shifts but couln'd make it work, so am
    % interpolating here to avoid this. There is some danger in
    % interpolation of complex numbers.
    % First construct noninteger grid coordinates centered at zero (not
    % grid ceter of N/2+1) centered on Bragg peak in frequency space
    ifxD_fx                 = (-crd.iFxDifDet:1/crd.SubSamp:crd.iFxDifDet) + fxDiff/crd.dfx;
    jfyD_fy                 = (-crd.jFyDet:1/crd.SubSamp:crd.jFyDet);
    [ifxD_fxfy,jfyD_fxfy]   = meshgrid(ifxD_fx,jfyD_fy);  
    % Interpolate to these points using center of grid N/2+1
    % Transposes deal with stupid matlab meshgrid index reversal
    Ediff_fxfyZ             = interp2(Ered_fxfyZ',(in.Nx/2+1) + ifxD_fxfy,...
                                                  (in.Ny/2+1) + jfyD_fxfy)';
    
    % And save for later calculation (or tomographic reconstruction)
    out.Ediff_fxfyZtheta(:,:,itheta,crd.iBraggT) = Ediff_fxfyZ;
    
    % Construct spatial frequencies of diffracted fields over the reduced
    % range, centered again on zero.  Transpose to index x,y not y,x
    fxD_fxfy                = crd.dfx * ifxD_fxfy';  
    fyD_fxfy                = crd.dfy * jfyD_fxfy';
    
    % Save grating spatial frequencies associated with each transverse 
    % spatial frequency and incident angle       
    out.Fx_fxfyZtheta(:,:,itheta) = fxD_fxfy - fxInc;   
	out.Fy_fxfyZtheta(:,:,itheta) = fyD_fxfy;	% minus fyInc = 0
    
    % Calculate fz of diffracted fields, then subtract incident to get Fz
	out.Fz_fxfyZtheta(:,:,itheta) = sqrt( (out.n0 / in.lambda0Red)^2 - fxD_fxfy.^2 - fyD_fxfy.^2 );
    
    % fzInc is not single number but varies across incident beam, shearing
    % apparent uncertainty.  To correct for this, subtract fzInc from fzDif
    % where fzInc is function of fx.  Follows calculation of diffracted
    % frequencies above.
    ifxI_fx                 = (-crd.iFxDifDet:1/crd.SubSamp:crd.iFxDifDet) + fxInc/crd.dfx;
    [ifxI_fxfy,jfyI_fxfy]   = ndgrid(ifxI_fx,jfyD_fy);  
    fxI_fxfy                = crd.dfx*ifxI_fxfy;  
    fyI_fxfy                = crd.dfy*jfyI_fxfy;
    
    % Find fz for entire incident beam 
    fzI_fxfy                = sqrt((out.n0 / in.lambda0Red)^2 - fxI_fxfy.^2 - fyI_fxfy.^2);  
        
    % and subtract from diffracted fz to get grating Fz
    out.Fz_fxfyZtheta(:,:,itheta) = out.Fz_fxfyZtheta(:,:,itheta) - fzI_fxfy;                                                                   

    % Integrate all diffracted spatial frequencies at this theta inc
    % This does NOT correct for tilt of uncertainty due to finite incident
    % beam.  That's done next. This result, therefore, will have reduced
    % null depth etc.
    out.Eta_theta(itheta,crd.iBraggT)  = sum(sum(abs(Ediff_fxfyZ).^2)) / out.Pref(crd.iMux) / crd.SubSamp^2;
                                                  
end % for itheta

%--------------------------------------------------------------------------
% Now that all diffracted fields are calculated, we can resample the skewed
% uncertainty to a regular grid and "unskew" the results.  This can be
% summed along the transverse spatial frequencies to simulate plane-wave
% readout.  This code was taken from HolographicTomographyV2.m
%--------------------------------------------------------------------------
% Fz samples at center of unskewed uncertainty grid
FzMin               = out.Fz_fxfyZtheta(crd.SubSamp*crd.iFxDifDet+1,crd.SubSamp*crd.jFyDet+1+1,1);
FzMax               = out.Fz_fxfyZtheta(crd.SubSamp*crd.iFxDifDet+1,crd.SubSamp*crd.jFyDet+1+1,end);

% Extents in Fx and Fy
FxMin               = min(out.Fx_fxfyZtheta,[],[1,2,3]);
FxMax               = max(out.Fx_fxfyZtheta,[],[1,2,3]);
FyMin               = min(out.Fy_fxfyZtheta,[],[1,2,3]);
FyMax               = max(out.Fy_fxfyZtheta,[],[1,2,3]);

% 2D grids
Fx_Fx               = linspace(FxMin,FxMax,2*crd.SubSamp*crd.iFxDifDet+1);
Fy_Fy               = linspace(FyMin,FyMax,2*crd.SubSamp*crd.jFyDet+1);
Fz_Fz               = linspace(FzMin,FzMax,in.NBraggAngles);
[Fx_FxFyFz,Fy_FxFyFz,Fz_FxFyFz]   = ndgrid(Fx_Fx,Fy_Fy,Fz_Fz);

% Create interpolating object from irregular sampled grid
% Strange indexing because scatInterp requires collumn vectors
Nelements = (2*crd.SubSamp*crd.iFxDifDet+1)*(2*crd.SubSamp*crd.jFyDet+1)*...
	in.NBraggAngles;
InterpToRegular = scatteredInterpolant(out.Fx_fxfyZtheta(:),...
    out.Fy_fxfyZtheta(:),out.Fz_fxfyZtheta(:),...
    reshape(out.Ediff_fxfyZtheta(:,:,:,crd.iBraggT),[Nelements,1]),...
    'linear','none');

% ...and interpolate to the regular grid
EDif_FxFyFz = InterpToRegular(Fx_FxFyFz,Fy_FxFyFz,Fz_FxFyFz);

% Set all extrapolated values to zero.  With ExtrapMethod = 'none', these
% will be NaNs
EDif_FxFyFz(isnan(EDif_FxFyFz)) = 0;

% Integrate over Fx and Fy, divide by incident power
out.Eta_Fz(:,crd.iBraggT) = squeeze(sum(abs(EDif_FxFyFz).^2,[1,2])) / out.Pref(crd.iMux) / crd.SubSamp^2;

end % AngleSpectrum

%==========================================================================
% IncompressibleU - Simple displacment model based on incompressibility.
% See "Bragg selectivity of distorted holograms.ppt"
%==========================================================================
function mat = IncompressibleU(crd,in,mat)

mat.u_xyz   = zeros(3,in.Nx,in.Ny,in.Nz);

r_xyz       = sqrt(crd.x_xyz.^2 + crd.y_xyz.^2);

% Displacement magnitude scaled to unit shrinkage, represented by (1)
umag_xyz    = (1) * 3 * in.w0wrt^2./r_xyz .* (1-exp(-(r_xyz/in.w0wrt).^2)) ...
                 .* (crd.z_xyz/in.Z) .* (1-crd.z_xyz/in.Z);
         
mat.u_xyz(1,:,:,:) = (crd.x_xyz./r_xyz) .* umag_xyz;    % mag of u times x unit vector
mat.u_xyz(2,:,:,:) = (crd.y_xyz./r_xyz) .* umag_xyz;    % mag of u times y unit vector
% displacement in this model is assumed to be radial so z component of u=0

% Deal with limit at origin
mat.u_xyz(1,crd.x_xyz == 0 & crd.y_xyz == 0) = 0;  % limit of umag as r->0 is 0
mat.u_xyz(2,crd.x_xyz == 0 & crd.y_xyz == 0) = 0;  

end % IncompressibleU