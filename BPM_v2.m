%==========================================================================
% BPM: Use beam propagation to calc E(x,y,z) given n(r), E(x,y,0)         
%      
% Required inputs:
%   E:      Incident e field vs transverse coord at z(1)
%           If first dimension = 2, taken to be x and y components of E
%           Dimensionality of E sets problem (x or x,y)
%   z       Z coordinate vector.  z(1) is plane of Einc
%           If z is vector, regular BPM to xy planes at locations z(iz)
%           If ndims(z)>1, XY dims must match E. BPM to z(x,y,iz)
% 
% Optional inputs via name value pairs (see ParseBPMInputs.m for details)
%   dx, dy: Cell size transversely
%   lambda0:Vacuum wavelength [same units as dx, z]
%   n:      Possibly complex homo+inhomo index. Can be z invarient or not.
%           If n=scalar, invarient in x,y,z
%           If dims(n)=dims(Einc), invarient in z. 
%           If dims(n)=1+dims(Einc), varies in z. 
%           At each z, avg (DC,homo) index removed and used for dif prop.
%           Dims must agree with E.  For 1D, E=[Nx,1],n=[Nx,1,Nz] or [Nx,1]
%   ABCMinTrans, ABCCellWidth:  Absorbing boundary data
%   RetXInd(1,2,3),RetYInd(1,2,3),RetZInd(1,2,3)  Indices of E returned
%   LocalnTransAvg              How to find uniform value of n for refraction.
%   
% Outputs:
%   E:      E field at requested indices.                                 
%           Vector problems return 3 vector components as first dimen.    
%                                                                         
% Original version:   3/06/09                                             
%   4/13/15:  Added vector isotropic via AOP paper                        
%   6/25/15:  Calculated average index via weighting with intensity.      
%             This handles refraction better than just an average.      
%   1/22/20   Input parsing.
%   1/29/20   Cleaned up structure, local funcitons, in/out/local structs
%                                                                         
%==========================================================================
function [in, out] = BPM_v2(Einc0,z,varargin)

%--------------------------------------------------------------------------
% Read inputs
%--------------------------------------------------------------------------
[Einc0,in]   	= ParseBPMInputs(Einc0,z,varargin);  
[in,out,local] 	= SetUpBPM(in);
                    
%--------------------------------------------------------------------------
% BPM z loop
%--------------------------------------------------------------------------
[out,local] = RecordResults(Einc0,1,in,out,local);

for iz = 2:in.Nz
	[out, local] = CalcPropagators(Einc0,iz,in,out,local);
    Einc0        = Propagate(Einc0,in,local);
    [out,local]  = RecordResults(Einc0,iz,in,out,local);
end

end % main

%==========================================================================
%-------------------------------LOCAL METHODS------------------------------
%==========================================================================

%==========================================================================
% Parse inputs for BPM_v2.
%==========================================================================
function [Einc0,in] = ParseBPMInputs(Einc0,z,VariablesPassed)
                 
nano  = 10^-9;                          % units.  All calculations MKS.
micro = 10^-6;
milli = 10^-3;

%--------------------------------------------------------------------------
% Determine problem size
%--------------------------------------------------------------------------
[NE1, NE2, NE3]  = size(Einc0);         % Get dimensions of incident field
Vector = (NE1 == 2 || NE1 == 3);     	% Vector propagation problem?

if Vector 
    Nx = NE2;                           % First dimension is components of E
    Ny = NE3;
else
    Nx = NE1;                           % else scalar, size = dims of E
    Ny = NE2;
end

% Planar or variable thickness (in xy) BPM?
if isvector(z)                      % Extract relevant dimensions from z
    Nz          = length(z);
else                                % Variable thickness BPM
	[~, Nz2, Nz3]  = size(z);       % That is, dz is function of xy
    if ismatrix(z)                    
        Nz = Nz2;
    else
        Nz = Nz3;
    end
     % Force same dims as E. This will throw error if z does not match E.
    z = reshape(z,Nx,Ny,Nz);          
end

%--------------------------------------------------------------------------
% Parser validation functions
%--------------------------------------------------------------------------
validArray      = @(x) isnumeric(x) && length(x) >1;                % Numeric but not a scalar
validScalar     = @(x) isnumeric(x) && isscalar(x) && (x >= 0);    	% Check most scalar inputs
validInteger    = @(x) x ==round(x) && isscalar(x) && (x >= 0);     % Scalar integer > 0
validXcell      = @(x) isnumeric(x) && (min(x(:)) >= 1) && (max(x(:)) <= Nx); % X cell coord
validYcell      = @(y) isnumeric(y) && (min(y(:)) >= 1) && (max(y(:)) <= Ny); % Y cell coord
validZcell      = @(z) isnumeric(z) && (min(z(:)) >= 1) && (max(z(:)) <= Nz); % Z cell coord
validDims       = @(E) numel(E) == Vector*Nx*Ny;                    % Array has same # elements as Einc0
validIndex      = @(n) isa(n,'function_handle') || length(n)==1 || numel(n) == Nx*Ny*Nz; % n is function or scalar or same size as E

p = inputParser;                                                    % Parser structure

%------------------------------Parameter dictionary------------------------
%                        PARAM     DEFAULT   VALIDITY         DEFINITION
%--------------------------------------------------------------------------
addParameter(p,       'lambda0',   1*micro, validScalar);     % Vacuum wavelength [m]
addParameter(p,            'dx', 0.5*micro, validScalar);     % X cell size [m]
addParameter(p,            'dy', 0.5*micro, validScalar);     % Y cell size [m]
addParameter(p,   'ABCMinTrans',         0, validScalar);     % Minimum absorbing BC transmission. 0...1
addParameter(p,  'ABCCellWidth',        10, validInteger); 	  % Absorbing BC width in cells
addParameter(p,      'RetXInd1',      1:Nx, validXcell);      % X coord of returned region 1
addParameter(p,      'RetYInd1',      1:Ny, validYcell);      % Y coord of returned region 1
addParameter(p,      'RetZInd1',      1:Nz, validZcell);      % Z coord of returned region 1
addParameter(p,      'RetXInd2',         0, validXcell);      % X coord of returned region 2
addParameter(p,      'RetYInd2',         0, validYcell);      % Y coord of returned region 2
addParameter(p,      'RetZInd2',         0, validZcell);      % Z coord of returned region 2
addParameter(p,      'RetXInd3',         0, validXcell);      % X coord of returned region 3
addParameter(p,      'RetYInd3',         0, validYcell);      % Y coord of returned region 3
addParameter(p,      'RetZInd3',         0, validZcell);      % Z coord of returned region 3
addParameter(p,'LocalnTransAvg',         1, validScalar);     % Intensity weighted or global average method for nTrans
addParameter(p,        'CalcEz',         0, validScalar);     % Calculate z-directed E fields?
addParameter(p,         'EincZ',         0, validDims);       % Electric field incident from z = Z
addParameter(p,        'nTrans',         0, validIndex);      % Index at fz = 0
addParameter(p,         'nRefl',         0, validIndex);      % Index at fz = 2 k
addParameter(p,         'ERefl',         0, validArray);      % Reflected fields each z from previous iteration
addParameter(p,      'ThetaInc',         0, validScalar);     % Incident angle, used for improved accy in dn projection

%--------------------------------------------------------------------------
% Parse.  Note that all assignments to "in" must happen after this point.
%--------------------------------------------------------------------------
parse(p,VariablesPassed{:});          	% Parse inputs into struct p
in = p.Results;                       	% Short hand notation "in"

in.Nx = Nx;  in.Ny = Ny;  in.Nz = Nz; 	% add problem dimensions
in.Vector       = Vector;
in.z            = z;
in.ZPlanar      = isvector(z);         	% Same z for each step?

if Vector && NE1 == 2 && in.CalcEz    	% Add Ez component if needed
    temp            = zeros(3,Nx,Ny);
    temp(1:2,:,:)   = Einc0;
    Einc0           = temp;
end

%--------------------------------------------------------------------------
% Detect nature of n.
%   A scalar:           Index is uniform in all space.  No refraction step.
%   An array
%     Same size as E    n does not vary with z, use same n in all z steps
%     One more dim      n varies in each z step
%   A function:         Call with coordinate z to get n at each plane
%--------------------------------------------------------------------------
if ~isa(in.nTrans,'function_handle')                    % Transmission
    in.nTransFunction   = 0;                            % Not a function
    in.nTransScalar     = isscalar(in.nTrans);      	% Is n scalar?
   	in.nTransVaries     = (numel(in.nTrans) == Nx*Ny*Nz);	% Does n vary with z?
else
    in.nTransFunction   = 1;                            % Is a function
    in.nTransVaries     = 1;                            % Assume function = variable n(z) 
end

if ~isa(in.nRefl,'function_handle')                     % Reflection
    in.nReflFunction    = 0;                            % Not a function
    in.nReflScalar      = isscalar((in.nRefl));         % Is n scalar?
    in.nReflVaries      = (numel(in.nRefl) == Nx*Ny*Nz);     % Does n vary with z?
else
    in.nReflFunction    = 1;                            % Is a function
    in.nReflVaries      = 1;                            % Assume function = variable n(z) 
end

%--------------------------------------------------------------------------
% Set flags
%--------------------------------------------------------------------------
in.Return1	= ~any(in.RetXInd1 == 0);       % Return region 1 requested? 
in.Return2	= ~any(in.RetXInd2 == 0);     	% Return region 2 requested? 
in.Return3	= ~any(in.RetXInd3 == 0);     	% Return region 3 requested? 

%--------------------------------------------------------------------------
% Check for needed inputs if bidirectional
%--------------------------------------------------------------------------
if length(in.EincZ) > 1 && length(in.nRefl) > 1 && length(in.ERefl) > 1
     in.Bidirectional    = true; 
elseif length(in.EincZ) == 1 && length(in.nRefl) == 1 && length(in.ERefl) == 1
     in.Bidirectional    = false;   
else
    error('Bidirectional simulation requires EincZ, nRefl and ERefl');
end

end % ParseBPMInputs

%==========================================================================
% Set up BPM
%==========================================================================
function [in,out,local] = SetUpBPM(in)

%--------------------------------------------------------------------------
% Set up Fourier coordinate arrays
%--------------------------------------------------------------------------
out.kx_kx 	= 0;
out.ky_ky	= 0;
    
if in.Nx > 1 
    dkx     	= 2*pi/(in.Nx*in.dx);
    out.kx_kx 	= (-in.Nx/2)*dkx: dkx : (in.Nx/2-1)*dkx;
end

if in.Ny > 1 
    dky         = 2*pi/(in.Ny*in.dy);
    out.ky_ky  	= (-in.Ny/2)*dky: dky : (in.Ny/2-1)*dky;
end

%--------------------------------------------------------------------------
% Calculate ABC transmission masks
%--------------------------------------------------------------------------
if in.Nx > 1
    XABC = [(1+ in.ABCMinTrans)/2 + (1-in.ABCMinTrans)/2 * ...
        cos(pi * (in.ABCCellWidth+1-(1:in.ABCCellWidth))/in.ABCCellWidth),...
        ones([1,in.Nx-2*in.ABCCellWidth]),...
        (1+in.ABCMinTrans)/2 + (1-in.ABCMinTrans)/2 * ...
        cos(pi * (in.ABCCellWidth-in.Nx+((in.Nx-in.ABCCellWidth+1):in.Nx))/in.ABCCellWidth)];
else
    XABC = 1;
end

if in.Ny > 1
    YABC = [(1+in.ABCMinTrans)/2 + (1-in.ABCMinTrans)/2 * ...
        cos(pi * (in.ABCCellWidth+1-(1:in.ABCCellWidth))/in.ABCCellWidth),...
        ones([1,in.Ny-2*in.ABCCellWidth]),...
        (1+in.ABCMinTrans)/2 + (1-in.ABCMinTrans)/2 * ...
        cos(pi * (in.ABCCellWidth-in.Ny+((in.Ny-in.ABCCellWidth+1):in.Ny))/in.ABCCellWidth)];
else
    YABC = 1;
end

local.XYABC = XABC' * YABC;              	% Outer product to make 2D condition

%--------------------------------------------------------------------------
% Initialize scalars and counters
%--------------------------------------------------------------------------
local.nbarlast	= -101;                     % To see if nbar has changed
local.dzlast  	= -1;                    	% To see if dz has changed.
local.Retiz1    = 1;                        % Counters for return arrays
local.Retiz2    = 1;
local.Retiz3    = 1;
local.k0        = 2 * pi / in.lambda0;   	% Free-space wavenumber for conven.

%--------------------------------------------------------------------------
% Allocate field reporting array(s)
%--------------------------------------------------------------------------
if in.Return1
    if in.Vector                   
        out.Eret1 = (1+1i*1)*zeros([2+in.CalcEz,length(in.RetXInd1),...
                            length(in.RetYInd1),length(in.RetZInd1)]);
    else
        out.Eret1 = (1+1i*1)*zeros([  length(in.RetXInd1),...
                  length(in.RetYInd1),length(in.RetZInd1)]);
    end
end

if in.Return2              
    if in.Vector
        out.Eret2 = (1+1i*1)*zeros([2+in.CalcEz,length(in.RetXInd2),...
                            length(in.RetYInd2),length(in.RetZInd2)]);
    else
        out.Eret2 = (1+1i*1)*zeros([  length(in.RetXInd2),...
                  length(in.RetYInd2),length(in.RetZInd2)]);
    end
end

if in.Return3             
    if in.Vector
        out.Eret3 = (1+1i*1)*zeros([2+in.CalcEz,length(in.RetXInd3),...
                            length(in.RetYInd3),length(in.RetZInd3)]);
    else
        out.Eret3 = (1+1i*1)*zeros([  length(in.RetXInd3),...
                  length(in.RetYInd3),length(in.RetZInd3)]);
    end
end

end % SetUpBPM

%==========================================================================
% CalcPropagators
%   Calculate diff and ref propagators
%
% Grid layout:
%      z=0                                     z=Z
%      |       |       |       ...     |       |
%      | n(1)  | n(2)  | n(3)  ...     |n(Nz-1)| n(Nz) is not used
%      |       |       |               |       |
%   iz 1       2       3               Nz-1    Nz
%
%==========================================================================
function [out, local] = CalcPropagators(Einc0,iz,in,out,local)

%--------------------------------------------------------------------------
% Z step.  In version 1, Einc0 specified to be at z=0 so first step is to
%  z(1).  Changed in version 2, Einc0 is at z(1)
%--------------------------------------------------------------------------
if in.ZPlanar
 	dz      = in.z(iz)-in.z(iz-1);          % dz is same for all x,y
    dzBar   = dz;
else
    dz      = in.z(:,:,iz)-in.z(:,:,iz-1);  % dz varies with x,y
    dzBar   = mean(dz,'all');
end

%--------------------------------------------------------------------------
% Calculate dn(x,y) and nBar = average of re(nTrans) over transverse (xy) slice
% Note array index iz-1.  iz is location we will propagate TO
%--------------------------------------------------------------------------
if in.nTransVaries
    % Calculate or extract n array at this z
    if in.nTransFunction
        nT_xy = in.nTrans(z(iz-1));         % Call function, get n at this z
    else
        if in.nTransScalar
            nT_xy = in.nTrans * ones(in.Nx,in.Ny);  % Replicate scalar n
        else
            if in.Nx ==1 || in.Ny ==1
                nT_xy = in.nTrans(:,iz-1);          % Extract 1D slice
            else
                nT_xy = in.nTrans(:,:,iz-1);        % Extract 2D slice
            end
        end
    end
    
    % Calculate nbar by one of two methods
    if in.LocalnTransAvg
        % Calculate average nTrans weighted with intensity 
        % This works well for (e.g.) a narrow beam on a tilted
        % interface.
        nBar = sum(real(nT_xy).*abs(Einc0).^2,'all')/sum(abs(Einc0).^2,'all');  
    else
        % Calcualte average over all transverse space
        nBar	= sum(real(nT_xy),'all')/(in.Nx*in.Ny);
    end
    
    % Combine two in order to get delta n
    dnTrans_xy  = nT_xy - nBar;

else
    % Homogenous space.  nbar = n.
    nBar = in.nTrans;
end

%--------------------------------------------------------------------------
% Set flags to control what gets calculated.  The goal is to only
% recalculate quantities when needed.
%--------------------------------------------------------------------------
dzChanged     	= ~in.ZPlanar    || (dz ~= local.dzlast);    % Second not checked if first is true
nBarChanged     = (iz == 2)      || (in.nTransVaries && nBar ~= local.nbarlast);
dnTransChanged	= (iz == 2)      || in.nTransVaries;
CalcKz          = nBarChanged;
CalcDiffProp   	= CalcKz         || dzChanged;
CalcRefProp  	= dnTransChanged || dzChanged;

%--------------------------------------------------------------------------
% Calculate kz if needed
%--------------------------------------------------------------------------
if CalcKz
    local.kz_kxky = conj(sqrt((local.k0*nBar)^2-...
         (repmat(out.kx_kx'.^2,[1,in.Ny])+repmat(out.ky_ky.^2,[in.Nx,1]))));
     
	%----------------------------------------------------------------------
    % If vector propagation, calculate vector modification to scalar
    % transfer function (AOP paper, equation 14)
    %----------------------------------------------------------------------
    if in.Vector
        local.kxOverkz_kxky = repmat(out.kx_kx',[1,in.Ny]) ./ local.kz_kxky;
        local.kxOverkz_kxky(isinf(local.kxOverkz_kxky)) = 10^6;   	% remove div by 0
        local.kxOverkz_kxky = ifftshift(local.kxOverkz_kxky);       % Put into fft coords

        local.kyOverkz_kxky = repmat(out.ky_ky ,[in.Nx,1]) ./ local.kz_kxky; 
        local.kyOverkz_kxky(isinf(local.kyOverkz_kxky)) = 10^6;   	% remove div by 0
        local.kyOverkz_kxky = ifftshift(local.kyOverkz_kxky);      	% Put into fft coords
    end
end

%--------------------------------------------------------------------------
% Calculate diffraction and refraction propagators if needed
% Note 1/cos(ThetaInc) factor.  This can improve the accy of refraction
% since the projection should be along the direction of propagation.  The
% correction only applies for fields of low and constant angular bandwidth
% such as holography.
%--------------------------------------------------------------------------
if CalcDiffProp
    % Note use of average dz in Fourier space.
    local.diffprop    = exp(-1i * ifftshift(local.kz_kxky) * dzBar);  %  in fft coords  
end   

if CalcRefProp
    if in.nTransVaries || ~in.ZPlanar        % XY variable phase?
        local.refprop = exp(-1i*local.k0*dnTrans_xy/cos(in.ThetaInc).*dz) .* local.XYABC;  
    else
        local.refprop = local.XYABC;       	% Phase flat.  Uniform phase in kz, diffprop
    end
end

% Record this nbar and dz for next time
local.nbarlast    = nBar;    
if in.ZPlanar, local.dzlast = dz; end   	% Only save for planar BPM

%--------------------------------------------------------------------------
% Reflection and transmission coefficients if bidi propagation.
% Even for a uniform reflection hologram, the phase of nRefl will vary with
% z, so don't bother with mechanisms used above for nTrans to detect if it
% is a constant.  
% Note we need nR of iz, the slab we are going IN to to calculate r
% t depends on transmision index of next slab (iz) and last (iz-1)
%--------------------------------------------------------------------------
if in.Bidirectional
    % Calculate or extract n array at this z
    if in.nReflFunction
        nR_xy  	= in.nRefl(z(iz));          % Call function, get n at this z
    else
        if in.Nx ==1 || in.Ny ==1           % Extract 1 or 2D nR
            nR_xy = in.nRefl(:,iz);
        else
            nR_xy = in.nRefl(:,:,iz);
        end
    end
    
    % See "Bidirectional BPM.ppt" for derivation.  Note use of ThetaInc to
    % approximiate kz.  See previous comment on when this is an acceptable
    % parameter to use.
    local.r_xy = -1i * local.k0 / (2 * cos(in.ThetaInc)) * nR_xy;
    local.t_xy = sqrt(1-(local.nTlast./nT_xy) .* abs(local.r_xy.^2));
    
	local.nTlast = nT_xy;                   % Save for t calculation
end
    
end % CalcPropagator

%==========================================================================
% RecordResults
%   Record E at requested locations
%==========================================================================
function [out,local] = RecordResults(E_xy,iz,in,out,local)

if in.Return1                         % Requested?
    if iz==in.RetZInd1(local.Retiz1) 
        if in.Vector
            out.Eret1(:,:,:,local.Retiz1) = E_xy(:,in.RetXInd1,in.RetYInd1); 
        else
            out.Eret1(:,:,local.Retiz1)   = E_xy(in.RetXInd1,in.RetYInd1); 
        end
        if local.Retiz1 < length(in.RetZInd1)       % Don't go beyond end of array
            local.Retiz1 = local.Retiz1 + 1;
        end
    end
    if iz == in.RetZInd1(end), out.Eret1 = squeeze(out.Eret1); end % Remove singlton dimensions
end

if in.Return2                         % Requested?
    if iz==in.RetZInd2(local.Retiz2) 
        if in.Vector
            out.Eret2(:,:,:,local.Retiz2) = E_xy(:,in.RetXInd2,in.RetYInd2); 
        else
            out.Eret2(:,:,local.Retiz2)   = E_xy(in.RetXInd2,in.RetYInd2); 
        end
        if local.Retiz2 < length(in.RetZInd2)  	% Don't go beyond end of array
            local.Retiz2 = local.Retiz2 + 1;
        end
    end
    if iz == in.RetZInd2(end), out.Eret2 = squeeze(out.Eret2); end % Remove singlton dimensions
end

if in.Return3                         % Requested?
    if iz==in.RetZInd3(local.Retiz3) 
        if in.Vector
            out.Eret3(:,:,:,local.Retiz3) = E_xy(:,in.RetXInd3,in.RetYInd3); 
        else
            out.Eret2(:,:,local.Retiz3)   = E_xy(in.RetXInd2,in.RetYInd3); 
        end
        if local.Retiz3 < length(in.RetZInd3)  	% Don't go beyond end of array
            local.Retiz3 = local.Retiz3 + 1;
        end
    end
    if iz == in.RetZInd3(end), out.Eret3 = squeeze(out.Eret3); end % Remove singlton dimensions
end

end % RecordResults

%==========================================================================
% Propagate
%   Single z step using diffraction and refraction BPM propagators
%     Note reversal of fft and ifft.  This is to switch sign convention on
%     fft kernal so that forward transform (ifft) is exp(+j k z).
%
% % Notation: A and A' are forward going fields before and after interface
%           B and B' are reverse going fields before and after interface
% Forward and reverse here are direction of BPM, not +z and -z
% Note that B' is from previous reverse going BPM and not now available.  
% rB' was calculated and stored, however.
%
% Inputs
%   A       Forward going BPM field that will transmit across the interface
%   rBprime Previously calculated reflection.  Forward going.
%
% Outputs
%   Aprime  Forward field transmitted through interface
%   rA      Reverse field, reflection of A from interface
%==========================================================================
function E_xy = Propagate(E_xy,in,local)

if in.Vector
    Ex_kxky     = ifft2(squeeze(E_xy(1,:,:))) .* local.diffprop;
    E_xy(1,:,:) = fft2(Ex_kxky) .* local.refprop;

    Ey_kxky     = ifft2(squeeze(E_xy(2,:,:))) .* local.diffprop;
    E_xy(2,:,:) = fft2(Ey_kxky) .* local.refprop;
    
    if in.CalcEz
        E_xy(3,:,:) = - fft2(Ex_kxky.*local.kxOverkz_kxky ...
                            +Ey_kxky.*local.kyOverkz_kxky) ...
                       .* local.refprop;
    end
    
else
    Ex_kxky     = ifft2(E_xy);
    E_xy        = fft2(Ex_kxky.*local.diffprop) .* local.refprop;
end

end % Propagate