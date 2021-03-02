function [emdSTEM] = PRISM03_probes(emdSTEM)
tic

% Colin Ophus - 2021 Feb
% 03 - Simulate STEM probes, save 3D and / or 4D output intensities

% Inputs
flagProgress = true;  % Display progress on console
% Output detector settings
if ~isfield(emdSTEM,'flagOutput3D'); emdSTEM.flagOutput3D = true; end
if ~isfield(emdSTEM,'flagOutput4D'); emdSTEM.flagOutput4D = false; end
% spacing of bins in 3D output (rads)
if ~isfield(emdSTEM,'drBins3D'); emdSTEM.drBins3D = 1 / 1000; end
% Probe positions at nearest wavefunction pixels (faster)
if ~isfield(emdSTEM,'flagProbePositionsNearestPixel'); emdSTEM.flagProbePositionsNearestPixel = true; end

% Probe positions
if ~isfield(emdSTEM,'xp')
    dxy = emdSTEM.cellDim(1:2) / 100;
    xR = [0 1]*emdSTEM.cellDim(1);
    yR = [0 1]*emdSTEM.cellDim(2);
    emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
    emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
    % Single probe position
    emdSTEM.xp = emdSTEM.cellDim(1) * 0.5;
    emdSTEM.yp = emdSTEM.cellDim(2) * 0.5;
end

if isempty(emdSTEM.partitionNumberRings)
    emdSTEM.xp = emdSTEM.xp - emdSTEM.pixelSize(1);
    emdSTEM.yp = emdSTEM.yp - emdSTEM.pixelSize(2);
end

% Output init
if emdSTEM.flagOutput3D == true
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
    emdSTEM.output4D = zeros( ...
        length(emdSTEM.xAA) / emdSTEM.interpolationFactor(1),...
        length(emdSTEM.yAA) / emdSTEM.interpolationFactor(2),...
        length(emdSTEM.xp),...
        length(emdSTEM.yp),...
        'single');
end

% Cropping box coordinate system
Nout = emdSTEM.imageSize ./ emdSTEM.interpolationFactor / 2;
emdSTEM.xCrop = (1:Nout(1)) - floor(Nout(1)/2) - 1;
emdSTEM.yCrop = (1:Nout(2)) - floor(Nout(2)/2) - 1;

% init outputs STEM probes, amplitude and masking
PsiOutput = zeros(Nout);
% Check for legacy PRISM or partitioned PRISM algorithm
if isempty(emdSTEM.partitionNumberRings)
    sub = emdSTEM.beamIndexInterp > 0;
    PsiCoefs = zeros(emdSTEM.beamNum,1);
    PsiCoefs(emdSTEM.beamIndexInterp(sub)) = emdSTEM.PsiInit(sub);
    scalePsiOutput = 1/prod(Nout)^2;
else
    Ntile = [prod(Nout) emdSTEM.beamNum];
    beamInds = cell2mat(emdSTEM.beamInds);
    
    if emdSTEM.flagProbePositionsNearestPixel == true
        % NOTE - should perhaps force this to be even after AA for centering
        %         beamPartition = reshape(circshift(ifft2( ...
        %             emdSTEM.beamPartition),[Nout/2 0]),Ntile);
        beamPartition = reshape(circshift(ifft2( ...
            emdSTEM.beamPartition),[floor(Nout/2) 0]),Ntile);
    else
        beamPartition = zeros(size(emdSTEM.beamPartition));
        
        % Pixel coordinate basis for subpixel shifting of the probe
        [qxa,qya] = makeFourierCoords(Nout,1);
        qxShiftBasis = -2i*pi*qxa(beamInds(:,1));
        qyShiftBasis = -2i*pi*qya(beamInds(:,1));
    end
    

end
probeCoefs = zeros(emdSTEM.beamNum,1,'single');

% pixel size and image, after AA cropping
emdSTEM.pixelSizeAA = 2 * emdSTEM.pixelSize;
emdSTEM.imageSizeAA = emdSTEM.imageSize / 2;


% Loop over all STEM probes
if flagProgress == true
    reverseStr = ''; % initialize console message piece
end
for ax = 1:length(emdSTEM.xp)
    
    % x dim cropping box, subpixel probe shift along x
    xP = emdSTEM.xp(ax) / emdSTEM.pixelSizeAA(1);
    xF = floor(xP);
    x = mod(emdSTEM.xCrop + xF, emdSTEM.imageSizeAA(1)) + 1;
    xa = repmat(x(:),[Nout(2) 1]);
    if ~isempty(emdSTEM.partitionNumberRings) && emdSTEM.flagProbePositionsNearestPixel == false
        dx = xP - xF + Nout(1)/2;
    end
    
    for ay = 1:length(emdSTEM.yp)
        % y dim cropping box, subpixel probe shift along x
        yP = emdSTEM.yp(ay) / emdSTEM.pixelSizeAA(2);
        yF = floor(yP);
        y = mod(emdSTEM.yCrop + yF, emdSTEM.imageSizeAA(2)) + 1;
        ya = repelem(y(:),Nout(1));
        if ~isempty(emdSTEM.partitionNumberRings) && emdSTEM.flagProbePositionsNearestPixel == false
            dy = yP - yF + Nout(2)/2;
        end
        
        % Realspace cropping coords
        inds = sub2ind(emdSTEM.imageSizeAA,xa,ya);
        
        % Check for legacy PRISM or partitioned PRISM algorithm
        if isempty(emdSTEM.partitionNumberRings)
            probeCoefs(:) = PsiCoefs .* exp( -2i*pi*( ...
                emdSTEM.beamList(:,3)*emdSTEM.xp(ax) + ...
                emdSTEM.beamList(:,4)*emdSTEM.yp(ay) ));
            
            for a0 = 1:emdSTEM.numFP
                
                if emdSTEM.interpolationFactor(1) == 1 && ...
                        emdSTEM.interpolationFactor(2) == 1 && ...
                        emdSTEM.numFP == 1
                    PsiOutput(:) = scalePsiOutput * abs(fft2(reshape( ...
                        emdSTEM.Scompact * probeCoefs,Nout))).^2;
                else
                    PsiOutput(:) = abs(fft2(reshape( ...
                        emdSTEM.Scompact(inds,:,a0) * probeCoefs,Nout))).^2;
                end
                
                if emdSTEM.flagOutput3D == true
                    emdSTEM.output3D(ax,ay,:) = squeeze(emdSTEM.output3D(ax,ay,:)) + ...
                        accumarray(qDetBins1, ...
                        PsiOutput(qDet1sub) .* qDetWeights1,[numDetBins 1]) + ...
                        accumarray(qDetBins2, ...
                        PsiOutput(qDet2sub) .* qDetWeights2,[numDetBins 1]);
                end
                
                if emdSTEM.flagOutput4D == true
                    emdSTEM.output4D(:,:,ax,ay) = ...
                        emdSTEM.output4D(:,:,ax,ay) + ...
                        PsiOutput;
                end
            end
            
            
        else
            if emdSTEM.flagProbePositionsNearestPixel == true
                for a0 = 1:emdSTEM.numFP
                    PsiOutput(:) = abs(fft2(reshape(sum(...
                        emdSTEM.Scompact(inds,:,a0) ...
                        .* beamPartition,2),Nout))).^2;
                    
                    if emdSTEM.flagOutput3D == true
                        emdSTEM.output3D(ax,ay,:) = squeeze(emdSTEM.output3D(ax,ay,:)) + ...
                            accumarray(qDetBins1, ...
                            PsiOutput(qDet1sub) .* qDetWeights1,[numDetBins 1]) + ...
                            accumarray(qDetBins2, ...
                            PsiOutput(qDet2sub) .* qDetWeights2,[numDetBins 1]);
                    end
                    
                    if emdSTEM.flagOutput4D == true
                        emdSTEM.output4D(:,:,ax,ay) = ...
                            emdSTEM.output4D(:,:,ax,ay) + ...
                            PsiOutput;
                    end
                end
                
            else
                % Generate initial STEM probe
                beamPartition(:) = emdSTEM.beamPartition;
                
                beamPartition(beamInds(:,2)) = ...
                    beamPartition(beamInds(:,2)) .* exp( ...
                    dx*qxShiftBasis + dy*qyShiftBasis);
                beamPartition(:) = ifft2(beamPartition);
                
                % Calculate the STEM probe using the masked, compact S-matrix, for
                % all FP configs
                for a0 = 1:emdSTEM.numFP
                    PsiOutput(:) = abs(fft2( ...
                        reshape(sum(emdSTEM.Scompact(inds,:,a0) ...
                        .* reshape(beamPartition,Ntile),2),Nout))).^2;
                    
                    if emdSTEM.flagOutput3D == true
                        emdSTEM.output3D(ax,ay,:) = squeeze(emdSTEM.output3D(ax,ay,:)) + ...
                            accumarray(qDetBins1, ...
                            PsiOutput(qDet1sub) .* qDetWeights1,[numDetBins 1]) + ...
                            accumarray(qDetBins2, ...
                            PsiOutput(qDet2sub) .* qDetWeights2,[numDetBins 1]);
                    end
                    
                    if emdSTEM.flagOutput4D == true
                        emdSTEM.output4D(:,:,ax,ay) = ...
                            emdSTEM.output4D(:,:,ax,ay) + ...
                            PsiOutput;
                    end
                end
            end
        end
    end
    
    if flagProgress == true
        comp =  ax / length(emdSTEM.xp);
        msg = sprintf(['Probe calculations are ' ...
            sprintf('%0.2f',100*comp) ' percent complete']);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
    end
end
if flagProgress == true 
    fprintf([reverseStr '']);
end
if emdSTEM.flagOutput3D == true
    emdSTEM.output3D(:) = emdSTEM.output3D / emdSTEM.numFP;
end
if emdSTEM.flagOutput4D == true
    emdSTEM.output4D(:) = emdSTEM.output4D / emdSTEM.numFP;
end

% % testing
% figure(12)
% clf
% EWamp = sqrt(PsiOutput);
% imagesc(fftshift(EWamp))
% axis equal off
% colorbar
% colormap(jetBlack(256))
% xlim([0 1]*Nout(1) + [1 0])
% ylim([0 1]*Nout(2) + [1 0])
% % xlim([0.3 0.7]*Nout(1) + [1 0])
% % ylim([0.3 0.7]*Nout(2) + [1 0])
% % % xlim([76 425])
% % % ylim([76 425])

emdSTEM.time03probes = toc;
end