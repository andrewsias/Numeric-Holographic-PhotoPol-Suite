%==========================================================================
% RenderHoloPolymer.  Render results of HoloPolymer function
%
% NEEDS:    Name                    Origin
%           -------------------     -----------------
%           LorentzLorenz           CU 
%           drawcircle              CU modified from Matlab file exchange
%           TransmisionBVM          CU
%           Kogelnik_Transmission   CU
%           arrow                   Matlab file exchange 
%           HolographicRenderer2D   CU
%           wavelength2color        Matlab file exchange
%
% INPUTS:   See parameter dictionary in ParseNHInputs. Not case sensitive.
% in        Struct with all input values including defaults
%
% OUTPUTS:
% out       Struct with all output values
%
% VERSION HISTORY:
% RRM   July 30 2019      Create
% RRM   Sept 2, 2019      VERSION 1.0
% RRM   June 3, 2020      Added Rayleigh scattering calculation
% RRM   Nov 1, 2020       Added k space render
%==========================================================================
function RenderHoloPolymer(in,out)

close all;

nano  = 10^-9;                          % units.  All calculations MKS.
micro = 10^-6;
milli = 10^-3;

% Build time array for exposure period with multiplexing
tau_tau = in.tau(1)*(0:in.Nt)/in.Nt;    % First entry is initial position
for iMux = 2:out.NMux
    tau_tau = [tau_tau, tau_tau(end) + in.tau(iMux)*(1:in.Nt)/in.Nt];
end

%--------------------------------------------------------------------------
% Draw k space
%--------------------------------------------------------------------------
figure(200); hold on;
for iMux = 1:out.NMux
    % Writing
    HolographicRenderer2D('lambda0',in.lambda0Wrt,...
    'RefTheta0',in.ThetaIncWrt(iMux),'ObjTheta0',in.ThetaDifWrt(iMux),... 
   'L',in.Z,'n0',out.n0,'DrawSnell',false,'DrawAngles',false,'Reading',false,...
   'RefKLabel',['$\bar{k}_{i,',num2str(iMux),'}$'], 'ObjKLabel',['$\bar{k}_{d,',num2str(iMux),'}$']);

    % Reading
    HolographicRenderer2D('lambda0',in.lambda0Red,...
    'RefTheta0',out.ThetaIncRed(iMux),'ObjTheta0',out.ThetaDifRed(iMux),...
   'L',in.Z,'n0',out.n0,'DrawSnell',false,'DrawAngles',false,'Reading',false,...
   'RefKLabel','', 'ObjKLabel','','HorizLabel',"",'VertLabel',"",...
   'DrawHorizAxis',false,'DrawVertAxis',false);  
end
xlim(1.5*xlim);         % Show whole two color diagram
ylim(1.5*ylim);

%--------------------------------------------------------------------------
% Index distribution
%--------------------------------------------------------------------------
figure(201);

% Scatter distribution
if length(out.dnS_xyz) == 1
    dnSSlice = zeros(in.Nx,in.Nz);                % No scattering
else
    dnSSlice = squeeze(out.dnS_xyz(:,in.Ny/2+1,:));
end

subplot(1,2,1);
h = pcolor(out.z_z/micro,out.x_x/micro,squeeze(abs(dnSSlice)));
LabelPlot('z [{\mu}m]','x [{\mu}m]','{\delta}n(x,y,z) scattering');
%daspect([1 1 1]);
colorbar;

% Recorded final (after flood cure) distribution.  No scatter centers
subplot(1,2,2);
holodn_x0z = squeeze(abs(squeeze(out.dn_xyz(:,in.Ny/2+1,:))-dnSSlice));
holodn_x0z = holodn_x0z - min(holodn_x0z,[],'all');  % Remove background
h = pcolor(out.z_z/micro,out.x_x/micro,holodn_x0z);

% Reading and writing beam w(z) contours
hold on;
zcen_z = out.z_z-in.Z/2;                            % Coordinates centered on material
z0wrt = pi*in.w0wrt^2/(in.lambda0Wrt/out.n0);       % Rayleigh range of writing beam
z0red = pi*in.w0red^2/(in.lambda0Wrt/out.n0);       % Rayleigh range of reading beam

% Pick colors to render write and read beams
WriteRGB = wavelength2color(in.lambda0Wrt/nano, 'maxIntensity', 1, 'colorSpace', 'rgb');
ReadRGB = wavelength2color(in.lambda0Red/nano, 'maxIntensity', 1, 'colorSpace', 'rgb');

for iMux = 1:out.NMux
    plot(out.z_z/micro,(+in.w0wrt*sqrt(1+(zcen_z/z0wrt).^2)+zcen_z*tan(in.ThetaIncWrt(iMux)))/micro,'Color',WriteRGB,'LineWidth',2);
    plot(out.z_z/micro,(-in.w0wrt*sqrt(1+(zcen_z/z0wrt).^2)+zcen_z*tan(in.ThetaIncWrt(iMux)))/micro,'Color',WriteRGB,'LineWidth',2);
    plot(out.z_z/micro,(+in.w0wrt*sqrt(1+(zcen_z/z0wrt).^2)+zcen_z*tan(in.ThetaDifWrt(iMux)))/micro,'Color',WriteRGB,'LineWidth',2);
    plot(out.z_z/micro,(-in.w0wrt*sqrt(1+(zcen_z/z0wrt).^2)+zcen_z*tan(in.ThetaDifWrt(iMux)))/micro,'Color',WriteRGB,'LineWidth',2);

    plot(out.z_z/micro,(+in.w0red*sqrt(1+(zcen_z/z0red).^2)+zcen_z*tan(out.ThetaIncRed(iMux)))/micro,'Color',ReadRGB,'LineWidth',2);
    plot(out.z_z/micro,(-in.w0red*sqrt(1+(zcen_z/z0red).^2)+zcen_z*tan(out.ThetaIncRed(iMux)))/micro,'Color',ReadRGB,'LineWidth',2);

end

LabelPlot('z [{\mu}m]','x [{\mu}m]','{\delta}n(x,y,z) recorded, final');
%daspect([1 1 1]);
colorbar;

%--------------------------------------------------------------------------
% Haze and holographic readout angular spectra
%--------------------------------------------------------------------------
dfy                 = 1/(in.Ny*in.dy);
FyHazeHA            = sind(out.HazeStopDeg)/in.lambda0Wrt;  % Haze half angle in spatial frequency
jFyHazeHA           = round(FyHazeHA/dfy);                  % Haze stop index
[fx_fxfy, fy_fxfy]  = ndgrid(out.fx_fx,out.fy_fy);

%---Holographic readout spectrum 
if out.Writing
    sf = figure(202);
    set(sf,'Position',[560   528   2*560   420]);
    subplot(1,2,1);
    Eref_fxfyZ  = fftshift(fft2(fftshift(out.Eref_xyZ)));
    Eref_fxfyZ = Eref_fxfyZ .* ((1/in.lambda0Wrt)^2 > (fx_fxfy.^2 + fy_fxfy.^2)); % Zero evanescent
    %CAxisMax = max(max(abs(Eref_fxfyZ).^2 .* (abs(fy_fxfy)>FyHazeHA))); % Max for color axis outside Bragg plane
    CAxisMax = max(max(abs(Eref_fxfyZ).^2));
    
    
    h = pcolor(out.fx_fx*micro,out.fy_fy*micro,log10(abs(Eref_fxfyZ).^2)');
    hold on;
    h = drawcircle([0,0],micro*sind(90)/in.lambda0Red,1000,'w');  h.LineWidth = 2; % Limits of diffraction in air
    %h = drawcircle([sin(in.ThetaIncWrt)*out.n0/in.lambda0Wrt*micro,0],micro*sind(out.HazeStopDeg)/in.lambda0Wrt,1000,'w'); h.LineWidth = 1.5;% Inc beam
    % h = drawcircle([sin(in.ThetaDifWrt)*out.n0/in.lambda0Wrt*micro,0],micro*sind(out.HazeStopDeg)/in.lambda0Wrt,1000,'w'); h.LineWidth = 1.5;% Inc beam
    % line(1/in.lambda0Wrt*micro*[-1,1],+FyHazeHA*micro*[1,1],'Color','w','LineWidth',1.2);
    % line(1/in.lambda0Wrt*micro*[-1,1],-FyHazeHA*micro*[1,1],'Color','w','LineWidth',1.2);
    % text(micro*sind(90)/in.lambda0Red,0,'90^o','color','w','HorizontalAlignment','right','FontSize',12); 
    % text((sin(mean(out.ThetaIncRed))*out.n0)/in.lambda0Red*micro,-2*sind(out.HazeStopDeg),'Inc','color','w','HorizontalAlignment','center','FontSize',14);
    % text((sin(mean(out.ThetaDifRed))*out.n0)/in.lambda0Red*micro,-2*sind(out.HazeStopDeg),'Dif','color','w','HorizontalAlignment','center','FontSize',14);


    LabelPlot('f_x [1/{\mu}m]','f_y [1/{\mu}m]',...
        ['Holographic replay: log_{10}(I), Haze_{ref} = ',...
        num2str(round(100*(out.HazeRef(out.NMux*in.Nt+3)-out.HazeXRef(out.NMux*in.Nt+3)),2)),'%']);
    shading('interp'); colormap('bone');
    caxis(log10(CAxisMax)+[-6, 0]);
    daspect([1 1 1]);
    colorbar;
    xlim(1/in.lambda0Red*micro*[-1.1,1.1]);
    ylim(1/in.lambda0Red*micro*[-1.1,1.1]);
    
    %---Haze spectrum 
    subplot(1,2,2);
    Enrm_fxfyZ  = fftshift(fft2(fftshift(out.Enrm_xyHazeZ(:,:,end))));  % Last z index if multiple
    Enrm_fxfyZ = Enrm_fxfyZ .* ((1/in.lambda0Wrt)^2 > (fx_fxfy.^2 + fy_fxfy.^2)); % Zero evanescent
    
    %CAxisMax = max(max(abs(Enrm_fxfyZ).^2 .* (abs(fy_fxfy)>sind(out.HazeStopDeg)/in.lambda0Wrt)));
    CAxisMax = max(max(abs(Enrm_fxfyZ).^2));


    h = pcolor(out.fx_fx*micro,out.fy_fy*micro,log10(abs(Enrm_fxfyZ).^2)');
    hold on;
    h = drawcircle([0,0],micro*sind(out.HazeStopDeg)/in.lambda0Red,1000,'w'); h.LineWidth = 1.5;% Haze stop
    h = drawcircle([0,0],micro*sind(90)/in.lambda0Red,1000,'w');  h.LineWidth = 2; % Limits of diffraction in air
    % line(1/in.lambda0Wrt*micro*[-1,1],+FyHazeHA*micro*[1,1],'Color','w','LineWidth',1.2);
    % line(1/in.lambda0Wrt*micro*[-1,1],-FyHazeHA*micro*[1,1],'Color','w','LineWidth',1.2);
    % text(micro*sind(90)/in.lambda0Red,0,'90^o','color','w','HorizontalAlignment','right','FontSize',12);
    % text(sind(out.HazeStopDeg)/in.lambda0Red*micro,0,'Stop','color','w','FontSize',12);
    LabelPlot('f_x [1/{\mu}m]','f_y [1/{\mu}m]',['Haze test: log_{10}(I), Haze = ',...
            num2str(round(100*(out.Haze(out.NMux*in.Nt+3,end)-out.HazeX(out.NMux*in.Nt+3,end)),2)),'%']);
    
    shading('interp'); colormap('bone');
    caxis(log10(CAxisMax)+[-6, 0]);
    daspect([1 1 1]);
    colorbar;
    xlim(1/in.lambda0Red*micro*[-1.1,1.1]);
    ylim(1/in.lambda0Red*micro*[-1.1,1.1]);
end % if writing

%--------------------------------------------------------------------------
% If multiple z slices for haze, plot evolution of haze vs tau and z
%--------------------------------------------------------------------------
if length(out.HazeZIndxs) > 1
    figure(1000)
    for it = 2:2:out.NMux*in.Nt
        ExcessHazePct = 100*(out.Haze(it,:)'-out.HazeX(it,:)' - (out.Haze(1,:)'-out.HazeX(1,:)'));
        f = fit((out.HazeZIndxs-1)'*out.dz/micro,ExcessHazePct,'poly3');
        plot(f,(out.HazeZIndxs-1)*out.dz/micro,ExcessHazePct);
        hold on;
    end
    legend('hide');
    LabelPlot('z [{\mu}m]','Haze - Initial Haze [%]','Haze vs exposure time and thickness');
end

%--------------------------------------------------------------------------
% Holographic replay xz slice
%--------------------------------------------------------------------------
figure (203);
Eref_fx0z = abs(fftshift(fft(out.Eref_x0z,in.Nx,1),1)).^2;
Eref_fx0z = Eref_fx0z/max(max(Eref_fx0z));

h = pcolor(out.z_z/micro,out.fx_fx*micro,log10(Eref_fx0z));
LabelPlot('z [{\mu}m]','f_x [1/{\mu}m]','log_{10}[ I(f_x,0,z) ]');
colormap(jet);   caxis([-3,0]);  colorbar;  

%--------------------------------------------------------------------------
% Efficiency and haze vs time
% Compare to Rayleigh prediction for TIS = n L sigma_theta 
% See "Holographic Model of Optical Haze.ppt"
%--------------------------------------------------------------------------
if out.Writing
    HazeStopDegInMat = asind(sind(out.HazeStopDeg)/out.n0);
    
    LineColors = [1 0 0
                  0 1 0
                  0 0 1
                  .5 .5 0
                  .5 0 .5
                  0 .5 .5
                  .5 .5 .5];         
    figure(204); hold off;
    colororder(LineColors);
    
    for iMux = 1:out.NMux           % Plot each multiplexed diffraction growth
        if in.TrackAllMux
            indices = (2+(iMux-1)*in.Nt):(out.NMux*in.Nt+1);      % Plot to end of illumination
        else
            indices = (2+(iMux-1)*in.Nt):(iMux*in.Nt+1);          % Plot only for this illumination
        end
        semilogy(tau_tau(indices),out.Eta(indices,iMux));
        hold on;
    end    
    semilogy(tau_tau,out.Haze(1:(out.NMux*in.Nt+1),end)-out.HazeX(1:(out.NMux*in.Nt+1),end),'r', ...
             tau_tau,out.HazeRef(1:(out.NMux*in.Nt+1))-out.HazeXRef(1:(out.NMux*in.Nt+1)),'b');
             
    % hold on;
    % % After diffusion
    % semilogy(tau_tau(end), out.Eta(out.NMux*in.Nt+2),'ko','MarkerFaceColor','k');
    % semilogy(tau_tau(end), out.Haze(out.NMux*in.Nt+2,end)-out.HazeX(out.NMux*in.Nt+2,end),'ro','MarkerFaceColor','r');
    % semilogy(tau_tau(end), out.HazeRef(out.NMux*in.Nt+2)-out.HazeXRef(out.NMux*in.Nt+2),'bo','MarkerFaceColor','b');
    % 
    % % After flood cure
    % semilogy(tau_tau(end), out.Eta(out.NMux*in.Nt+3),'ks','MarkerFaceColor','k');
    % semilogy(tau_tau(end), out.Haze(out.NMux*in.Nt+3,end)-out.HazeX(out.NMux*in.Nt+3,end),'rs','MarkerFaceColor','r');
    semilogy(tau_tau(end), out.HazeRef(out.NMux*in.Nt+3)-out.HazeXRef(out.NMux*in.Nt+3),'bs','MarkerFaceColor','b');

    ylim([10^(-5),1]);
    LabelPlot('\tau','\eta',...
        'Holographic efficiency and haze verus exposure time');
    % legend('Efficiency','Haze out of Bragg plane','Haze from reference beam',...
    %    'Location','southeast');
end % if writing
            
%--------------------------------------------------------------------------
% Bragg selectivity
%--------------------------------------------------------------------------
if in.CalcBragg
    
    % Kogelnik theory for mux'ed holograms
    FineTheta_theta = linspace(min(out.HoloThetaRed_theta),max(out.HoloThetaRed_theta),100);
    etaTheory_theta = zeros(1,100);        

    % Sum up expected DE from all holograms.  This isn't an exact solution
    % if there are multiple holograms, but it's simple.
    for iMux = 1:out.NMux

        if out.Writing
            dnTheoretical = out.DnThy(iMux * in.Nt +1) - out.DnThy((iMux-1) * in.Nt  + 2);
        else
            dnTheoretical = in.deltan;
        end

        % Calcualte expected plane wave diffraction efficiency for this Dn        
         etaTheory_theta = etaTheory_theta + ...
            Kogelnik_Transmission(in.Z, out.n0, dnTheoretical, ...
            out.Lambda(iMux), in.lambda0Red, out.phi(iMux), FineTheta_theta);
    end

    if out.Writing  % Only if Gaussian writing beams used
        % Correct for Gaussian beams, see Proc. of SPIE Vol. 10233 102330B-2
        % Figure 1.a can be derived analytyicaly (see Gaussian holographic
        % efficiency.nb) as 1/(1 + 2 (w0read/w0write)^2))    
        etaTheory_theta = etaTheory_theta / (1 + 2*(in.w0red/in.w0wrt)^2);          
    end
        
    for iBraggT = 1:out.NBraggT         % For all measurement times
        figure(300+iBraggT);

        plot(180/pi*FineTheta_theta,etaTheory_theta,'k');
        hold on;
        % This is the direct integration efficiency.  The following line is the
        % more sophisticated readout which corrects for finite read beam size
       plot(180/pi*out.HoloThetaRed_theta,out.Eta_theta,'b'); 

       plot(180/pi*out.HoloThetaRed_theta,out.Eta_Fz(:,iBraggT),'ro');
        legend('Theory','Direct numerical','Corrected numerical');
        LabelPlot('{\theta} [degrees]','\eta',['Bragg selectivity',]);  
    end

end

%--------------------------------------------------------------------------
% Compare simple predicted to measured polymer zero and first harmonic
%--------------------------------------------------------------------------
if out.Writing & (in.ThetaIncWrt ~= in.ThetaDifWrt)             % Skip if beams co-aligned
    figure(206);
    plot(tau_tau,1-exp(-tau_tau),'k.');hold on;
    plot(tau_tau,out.P0Num(1:out.NMux*in.Nt+1),'ko'); 
    plot(tau_tau,out.P0Thy(1:out.NMux*in.Nt+1),'k'); 
    plot(tau_tau,out.P1Num(1:out.NMux*in.Nt+1),'bo'); 
    plot(tau_tau,out.P1Thy(1:out.NMux*in.Nt+1),'b'); 
    plot(tau_tau,out.P2Num(1:out.NMux*in.Nt+1),'go'); 
    plot(tau_tau,out.P2Thy(1:out.NMux*in.Nt+1),'g'); 
    plot(tau_tau,out.P3Num(1:out.NMux*in.Nt+1),'ro'); 
    plot(tau_tau,out.P3Thy(1:out.NMux*in.Nt+1),'r'); 

    LabelPlot('\tau','P_i','Polymer harmonics');  
     legend('R_m \Rightarrow \infty','P_0 numerical','P_0 theory',...
         'P_1 numerical','P_1 theory','P_2 numerical','P_2 theory',...
         'P_3 numerical','P_3 theory','location','northwest');

    %----------------------------------------------------------------------
    % Compare theory and numerical index contrast
    %----------------------------------------------------------------------
    figure(207);
    plot(tau_tau,out.DnThy,'r'); hold on;
    plot(tau_tau,out.DnNum(1:out.NMux*in.Nt+1),'bo'); 
    LabelPlot('\tau','{\delta}n','Theory vs. numerical index');  
    legend('Theory','Numerical','location','northwest');

end

end % RenderHoloPolymer

%==========================================================================
% Common plot options + labels
%==========================================================================
function LabelPlot(xtitle,ytitle,plottitle)

shading('interp'); colormap('bone');set(gcf,'color','w');
set(gca,'FontSize',16,'FontName','Arial');
xlabel(xtitle); 
ylabel(ytitle);
title(plottitle);

end % LabelPlot

