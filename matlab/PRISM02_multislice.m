function [emdSTEM] = PRISM02_multislice(emdSTEM)

% Colin Ophus - 2020 Sept
% 02 - conventional multislice STEM for comparison, from same potential input

% Inputs
flagProgress = true;  % Display progress on console
% Output detector settings
if ~isfield(emdSTEM,'flagOutput3D'); emdSTEM.flagOutput3D = true; end
if ~isfield(emdSTEM,'flagOutput4D'); emdSTEM.flagOutput4D = false; end
% Microscope voltage [Volts]
if ~isfield(emdSTEM,'E0'); emdSTEM.E0 = 80e3; end
% probe convergence semiangle [rads]
if ~isfield(emdSTEM,'probeSemiangleArray'); emdSTEM.probeSemiangleArray = 20 / 1000; end
 % relative to the entrance surface of the potentials [Angstroms]
if ~isfield(emdSTEM,'probeDefocusDF'); emdSTEM.probeDefocusDF = 0; end
% spacing of bins in 3D output [rads]
if ~isfield(emdSTEM,'drBins3D'); emdSTEM.drBins3D = 1 / 1000; end
% single thickness output
if ~isfield(emdSTEM,'thicknessOutput'); emdSTEM.thicknessOutput = emdSTEM.cellDim(3); end

%  Probe positions
if ~isfield(emdSTEM,'xp')
    dxy = emdSTEM.cellDim(1:2) / 512;
    xR = [0 1]*emdSTEM.cellDim(1);
    yR = [0 1]*emdSTEM.cellDim(2);
    emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
    emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);

    % Single probe position
%     emdSTEM.xp = emdSTEM.cellDim(1) * 0.5;
%     emdSTEM.yp = emdSTEM.cellDim(2) * 0.5;
    emdSTEM.xp = emdSTEM.cellDim(1) * 0.52;
    emdSTEM.yp = emdSTEM.cellDim(2) * 0.54;
end

% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
emdSTEM.lambda = h/sqrt(2*m*e*emdSTEM.E0) ...
    /sqrt(1 + e*emdSTEM.E0/2/m/c^2) * 10^10; % wavelength in A
emdSTEM.sigma = (2*pi/emdSTEM.lambda/emdSTEM.E0) ...
    *(m*c^2+e*emdSTEM.E0)/(2*m*c^2+e*emdSTEM.E0);

% Fourier coordinates
qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
[emdSTEM.qya,emdSTEM.qxa] = meshgrid(qy,qx);
emdSTEM.q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;
emdSTEM.q1 = sqrt(emdSTEM.q2);

% Initial probe amplitude
qProbe = emdSTEM.probeSemiangleArray / emdSTEM.lambda;
if  emdSTEM.interpolationFactor(1) == 1 && ...
        emdSTEM.interpolationFactor(2) == 1
    dqx = qx(2) - qx(1);
    dqy = qy(2) - qy(1);
    emdSTEM.PsiInit = min(max( ...
        (qProbe*emdSTEM.q1 - emdSTEM.q2) ./ ...
        sqrt(dqx^2*emdSTEM.qxa.^2 + dqy^2*emdSTEM.qya.^2) ...
        + 0.5, 0), 1);
    emdSTEM.PsiInit(1,1) = 1;
    emdSTEM.PsiInit(:) = emdSTEM.PsiInit / sqrt(sum(abs(emdSTEM.PsiInit(:)).^2));

else
    %  Updated probes which take into account PRISM interpolation factors
    vxInter = 1:emdSTEM.interpolationFactor(1):emdSTEM.imageSize(1);
    vyInter = 1:emdSTEM.interpolationFactor(2):emdSTEM.imageSize(2);
    qxaInterp = emdSTEM.qxa(vxInter,vyInter);
    qyaInterp = emdSTEM.qya(vxInter,vyInter);
    q2Interp = qxaInterp.^2 + qyaInterp.^2;
    q1Interp = sqrt(q2Interp);
    dqx = qxaInterp(2,1) - qxaInterp(1,1);
    dqy = qyaInterp(1,2) - qyaInterp(1,1);
    
    PsiInitInterp = min(max( ...
        (qProbe*q1Interp - q2Interp) ./ ...
        sqrt(dqx^2*qxaInterp.^2 + dqy^2*qyaInterp.^2) ...
        + 0.5, 0), 1);
    PsiInitInterp(1,1) = 1;

    PsiInitInterp(:) = ifft2(PsiInitInterp);
    vx = [(1:(size(PsiInitInterp,1)/2)) ...
        (1-(size(PsiInitInterp,1)/2):0) + emdSTEM.imageSize(1)];
    vy = [(1:(size(PsiInitInterp,2)/2)) ...
        (1-(size(PsiInitInterp,2)/2):0) + emdSTEM.imageSize(2)];
    emdSTEM.PsiInit = zeros(emdSTEM.imageSize);
    emdSTEM.PsiInit(vx,vy) = PsiInitInterp;
    %     emdSTEM.PsiInit(:) = repmat(PsiInitInterp,emdSTEM.interpolationFactor);
    emdSTEM.PsiInit(:) = fft2(emdSTEM.PsiInit);
    emdSTEM.PsiInit(:) = emdSTEM.PsiInit / sqrt(sum(abs(emdSTEM.PsiInit(:)).^2));
end

% Probe defocus
emdSTEM.probeDefocusC1 = -1 * emdSTEM.probeDefocusDF;
chi = (pi*emdSTEM.lambda*emdSTEM.probeDefocusC1)*emdSTEM.q2;
PsiInit = emdSTEM.PsiInit .* exp(-1i*chi);

% Probe shift basis
qxShift = -2i*pi*emdSTEM.qxa;
qyShift = -2i*pi*emdSTEM.qya;
sub = abs(PsiInit(:)) > 0;
qxShiftSub = qxShift(sub);
qyShiftSub = qyShift(sub);

% AA aperture reduced coordinate indices and mask
emdSTEM.xAA = [(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-(emdSTEM.imageSize(1)/4)):0)+emdSTEM.imageSize(1)];
emdSTEM.yAA = [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-(emdSTEM.imageSize(2)/4)):0)+emdSTEM.imageSize(2)];

emdSTEM.maskAA = false(emdSTEM.imageSize);
emdSTEM.maskAA(emdSTEM.xAA,emdSTEM.yAA) = true;
emdSTEM.maskAAinv = ~emdSTEM.maskAA;

if emdSTEM.interpolationFactor(1) > 1 || ...
        emdSTEM.interpolationFactor(2) > 1
        emdSTEM.xAAinterp = [(1:emdSTEM.interpolationFactor(1):(emdSTEM.imageSize(1)/4)) ...
        ((1-(emdSTEM.imageSize(1)/4)):emdSTEM.interpolationFactor(1):0) ...
        + emdSTEM.imageSize(1)];
    emdSTEM.yAAinterp = [(1:emdSTEM.interpolationFactor(2):(emdSTEM.imageSize(2)/4)) ...
        ((1-(emdSTEM.imageSize(2)/4)):emdSTEM.interpolationFactor(2):0) ....
        + emdSTEM.imageSize(2)];
end

% propagator
emdSTEM.prop = exp( emdSTEM.q2(emdSTEM.xAA,emdSTEM.yAA) * ...
    (-1i*pi*emdSTEM.lambda*emdSTEM.sliceThickness));

% init output
PsiOutput = zeros(emdSTEM.imageSize / 2);

% Detector coordinates of 3D output
if emdSTEM.flagOutput3D == true
    % Generate the AA coordinate system
    N = emdSTEM.imageSize  / 2;
    qx = makeFourierCoords(N(1), 2*emdSTEM.pixelSize(1));
    qy = makeFourierCoords(N(2), 2*emdSTEM.pixelSize(2));
    [emdSTEM.qyaInterp,emdSTEM.qxaInterp] = meshgrid(qy,qx);
    emdSTEM.q2Interp = emdSTEM.qxaInterp.^2 + emdSTEM.qyaInterp.^2;
    emdSTEM.q1Interp = sqrt(emdSTEM.q2Interp);

    % total detector bins
    emdSTEM.qMax = min(max(emdSTEM.qxaInterp(:,1)),max(emdSTEM.qyaInterp(1,:)));
    alphaMax = emdSTEM.qMax * emdSTEM.lambda;
    emdSTEM.detectorAngles = (emdSTEM.drBins3D/2):emdSTEM.drBins3D:(alphaMax-emdSTEM.drBins3D/2);
    numDetBins = length(emdSTEM.detectorAngles);
    
    % detector bin indices
    dqDet = emdSTEM.drBins3D / emdSTEM.lambda;
    qDetInds = max(emdSTEM.q1Interp / dqDet,1);
    qF = floor(qDetInds);
    dq = qDetInds - qF;
    
    % lower index detector bins
    qDet1sub = qF <= numDetBins;
    qDetBins1 = qF(qDet1sub);
    qDetWeights1 = 1 - dq(qDet1sub);
    
    qDet2sub = qF <= numDetBins-1;
    qDetBins2 = qF(qDet2sub)+1;
    qDetWeights2 = dq(qDet2sub);
    
    % init 3D output array
    emdSTEM.output3D = zeros(length(emdSTEM.xp),length(emdSTEM.yp),numDetBins);
end


if emdSTEM.flagOutput4D == true
    % planes to output thickness
    z = (1:emdSTEM.numPlanes)*emdSTEM.sliceThickness;
    indOutput4D = zeros(length(emdSTEM.thicknessOutput),1);
    for a0 = 1:length(emdSTEM.thicknessOutput)
        [~,indOutput4D(a0)] = min(abs(z - emdSTEM.thicknessOutput(a0)));
    end
    
    emdSTEM.output4D = zeros( ...
        length(emdSTEM.xAA) / emdSTEM.interpolationFactor(1),...
        length(emdSTEM.yAA) / emdSTEM.interpolationFactor(2),...
        length(emdSTEM.xp),...
        length(emdSTEM.yp),...
        length(emdSTEM.thicknessOutput),...
        'single');
end

% Extra probe shift to align coordinate system with PRISM02_SMatrix
dxyAlign = -emdSTEM.pixelSize;

% Main loop
if flagProgress == true
    reverseStr = ''; % initialize console message piece
end
tic
for a0 = 1:emdSTEM.numFP
    trans = exp((1i*emdSTEM.sigma) * emdSTEM.pot(:,:,:,a0));
    
    for ax = 1:length(emdSTEM.xp)
        
        
        for ay = 1:length(emdSTEM.yp)
            
            % Shift probe - note 1 pixel shift added to agree with the 
            % coordinate convention used in the partitioning / PRISM sims
            Psi = PsiInit;
            Psi(sub) = Psi(sub) .* exp( ....
                qxShiftSub*(emdSTEM.xp(ax) + dxyAlign(1)) + ...
                qyShiftSub*(emdSTEM.yp(ay) + dxyAlign(2)));
            %             Psi(sub) = Psi(sub) .* exp( ....
            %                 qxShiftSub*emdSTEM.xp(ax) + qyShiftSub*emdSTEM.yp(ay));
            
            
            % Propagate through foil
            for a2 = 1:emdSTEM.numPlanes
                Psi(:) = fft2(ifft2(Psi) .* trans(:,:,a2));
                if a2 < emdSTEM.numPlanes
                    Psi(emdSTEM.xAA,emdSTEM.yAA) = ...
                        Psi(emdSTEM.xAA,emdSTEM.yAA) .* emdSTEM.prop;
                    Psi(emdSTEM.maskAAinv) = 0;
                end
                
                if emdSTEM.flagOutput4D == true
                    [val,ind] = min(abs(indOutput4D - a2));
                    if val == 0
                        if emdSTEM.interpolationFactor(1) > 1 || ...
                                emdSTEM.interpolationFactor(2) > 1
                            emdSTEM.output4D(:,:,ax,ay,ind) = ...
                                emdSTEM.output4D(:,:,ax,ay,ind) + ...
                                abs(Psi(emdSTEM.xAAinterp,emdSTEM.yAAinterp)).^2;
                        else
                            emdSTEM.output4D(:,:,ax,ay,ind) = ...
                                emdSTEM.output4D(:,:,ax,ay,ind) + ...
                                abs(Psi(emdSTEM.xAA,emdSTEM.yAA)).^2;
                        end
                    end
                end
            end
            PsiOutput(:) = abs(Psi(emdSTEM.xAA,emdSTEM.yAA)).^2;
            
            if emdSTEM.flagOutput3D == true
                emdSTEM.output3D(ax,ay,:) = squeeze(emdSTEM.output3D(ax,ay,:)) + ...
                    accumarray(qDetBins1, ...
                    PsiOutput(qDet1sub) .* qDetWeights1,[numDetBins 1]) + ...
                    accumarray(qDetBins2, ...
                    PsiOutput(qDet2sub) .* qDetWeights2,[numDetBins 1]);
            end
            
            if flagProgress == true
                comp = ((ay / length(emdSTEM.yp) ...
                    + ax - 1) / length(emdSTEM.xp) ...
                    + a0 - 1) / emdSTEM.numFP;
                msg = sprintf(['Multislice simulation is ' ...
                    sprintf('%0.2f',100*comp) ' percent complete']);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'),1,length(msg));
            end
        end
    end
end
if flagProgress == true 
    fprintf([reverseStr '']);
end
if emdSTEM.flagOutput3D == true
    emdSTEM.output3D(:) = emdSTEM.output3D / emdSTEM.numFP;
end
if emdSTEM.flagOutput4D == true
    emdSTEM.output4D(:) = emdSTEM.output4D ...
        * (prod(emdSTEM.interpolationFactor) / emdSTEM.numFP);
end



% figure(1002)
% clf
% plot(emdSTEM.detectorAngles * 1e3,...
%     squeeze(emdSTEM.output3D(ax,ay,:)),...
%     'linewidth',2,'color','r')
% % plot(emdSTEM.detectorAngles * 1e3,...
% %     squeeze(emdSTEM.output3D(ax,ay,:)) ./ emdSTEM.detectorAngles(:),...
% %     'linewidth',2,'color','r')



% figure(1001)
% clf
% EWamp = sqrt(PsiOutput);
% imagesc(fftshift(EWamp))
% % imagesc(fftshift(abs(Psi)))
% % imagesc(fftshift(abs(emdSTEM.PsiInit)))
% % imagesc((abs(ifft2(Psi))))
% axis equal off
% colormap(jetBlack)
% colorbar
% set(gca,'position',[0 0 0.85 1])

emdSTEM.time02multislice = toc;
end