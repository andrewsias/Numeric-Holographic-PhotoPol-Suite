%==========================================================================
% TheoryHoloPolymer.  Use harmonic expansion of RD equations plus strong
% limit theory for index formation to predict diffraction efficiency vs.
% conversion a priori.
% 
% NEEDS:    Name                Origin
%           ------------------- -----------------
%           LorentzLorenz       CU 
%           P0...P3             CU
%           m0...m3             CU
%           PolymerRoots        CU
%
% INPUTS:   See parameter dictionary in ParseNHInputs. Not case sensitive.
% in        Struct with all input values including defaults
%
% OUTPUTS:
% out       Struct with all output values
%
% Time index: 1 = initial condition, 2:Nt+1 = illumination, Nt+2 = dark
% development via full diffusion, Nt+3 = flood cure.
%
% VERSION HISTORY:
% RRM   July 30 2019      Create
% RRM   Sept 2, 2019      VERSION 1.0
%==========================================================================
function out = TheoryHoloPolymer(in,out)

% Do not execute if angle between writing beams is zero
if in.ThetaIncWrt == in.ThetaDifWrt
    disp('Warning from TheoryHoloPolymer: Beam angle = 0, returning.');
    return; 
end

%--------------------------------------------------------------------------
% Calculate harmonics of monomer and polymer 
%--------------------------------------------------------------------------
% Build time array for exposure period with multiplexing
tau_tau = in.tau(1)*(0:in.Nt)/in.Nt;    % First entry is initial position
for iMux = 2:out.NMux
    tau_tau = [tau_tau, tau_tau(end) + in.tau(iMux)*(1:in.Nt)/in.Nt];
end

% Polynomial solution was defined with Gm and ExpT.  Should change over to
% Rm and tau.  For the moment, just pass Gm = Rm / (2*pi)^2 and ExpT = 2 tau
Proots = PolymerRoots(in.Rm/(2*pi)^2/2,out.Vr);        % Solve for roots of characteristic polynomial

out.m0Thy = m0(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);  	% Monomer harmonics
out.m1Thy = m1(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);
out.m2Thy = m2(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);
out.m3Thy = m3(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);

out.P0Thy = P0(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);  	% Polymer harmonics
out.P1Thy = P1(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);
out.P2Thy = P2(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);
out.P3Thy = P3(in.Rm/(2*pi)^2/2,out.Vr,2*tau_tau,Proots);

%--------------------------------------------------------------------------
%---Calculate Dn 
% Note that don't need molar volume.  
% phi_m = gamma_m [m] = phim0 (gamma_m [m] / gamma_m [m0]) = phim0 m
% Note that m1 will generally be negative and P1 positive so concentrations
% in bright/dark fringe are (m0 +/- m1) and (P0 +/- P1) respectively.
%--------------------------------------------------------------------------

nPeak_t = LorentzLorenz(in.nmono,in.phim0*(out.m0Thy + out.m1Thy),...
             in.nP,in.phim0*(1-in.sigma)*(out.P0Thy + out.P1Thy), ...
             in.nM,1-in.phim0*(out.m0Thy + out.m1Thy)-in.phim0*(1-in.sigma)*(out.P0Thy + out.P1Thy));
         
nValy_t = LorentzLorenz(in.nmono,in.phim0*(out.m0Thy - out.m1Thy),...
             in.nP,in.phim0*(1-in.sigma)*(out.P0Thy - out.P1Thy), ...
             in.nM,1-in.phim0*(out.m0Thy - out.m1Thy)-in.phim0*(1-in.sigma)*(out.P0Thy - out.P1Thy));        
         
out.DnThy = (nPeak_t - nValy_t)/2;

%--------------------------------------------------------------------------
%--- Calculate maximum possible Dn
% Polymer in bright fringe = 2 and 0 in dark
%--------------------------------------------------------------------------
m0Thy =0;  m1Thy = 0; P0Thy = 1;  P1Thy = 1; 

nPeak_t = LorentzLorenz(in.nmono,in.phim0*(m0Thy + m1Thy),...
             in.nP,in.phim0*(1-in.sigma)*(P0Thy + P1Thy), ...
             in.nM,1-in.phim0*(m0Thy + m1Thy)-in.phim0*(1-in.sigma)*(P0Thy + P1Thy));
         
nValy_t = LorentzLorenz(in.nmono,in.phim0*(m0Thy - m1Thy),...
             in.nP,in.phim0*(1-in.sigma)*(P0Thy - P1Thy), ...
             in.nM,1-in.phim0*(m0Thy - m1Thy)-in.phim0*(1-in.sigma)*(P0Thy - P1Thy));  

out.DnThyMax = (nPeak_t - nValy_t)/2;

end
         
         
         