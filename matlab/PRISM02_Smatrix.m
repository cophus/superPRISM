function [emdSTEM] = PRISM02_Smatrix(emdSTEM)
tic

% Colin Ophus - 2021 Feb
% 02 - Compute the interpolated "Compact S-Matrix", either for STEM probe
%      partitioning + PRISM interpolation, or legacy PRISM.

% Inputs
flagProgress = true;  % Display progress on console
% Microscope voltage [Volts]
if ~isfield(emdSTEM,'E0'); emdSTEM.E0 = 80e3; end    
% probe convergence semiangle [rads]
if ~isfield(emdSTEM,'probeSemiangleArray'); emdSTEM.probeSemiangleArray = 20 / 1000; end
% relative to the entrance surface of the potentials [Angstroms]
if ~isfield(emdSTEM,'probeDefocusDF'); emdSTEM.probeDefocusDF = 0; end
% relative to the entrance surface of the potentials [Angstroms]
if ~isfield(emdSTEM,'sMatrixRefocus'); emdSTEM.sMatrixRefocus = false; end
% Number of radial rings for the partitioning
% Set this value to empty to use original PRISM
if ~isfield(emdSTEM,'partitionNumberRings'); emdSTEM.partitionNumberRings = 4; end
% emdSTEM.partitionNumberRings = [];
% If this option is not empty, select N brightest beams from FFT of proj pot.
% if ~isfield(emdSTEM,'selectBrightestBeamsNum'); emdSTEM.selectBrightestBeamsNum = 6; end
if ~isfield(emdSTEM,'selectBrightestBeamsNum'); emdSTEM.selectBrightestBeamsNum = []; end

% Set this flag to true to use a sigmoidal basis
if ~isfield(emdSTEM,'partitionSigmoidal'); emdSTEM.partitionSigmoidal = false; end

% Probe C1 defocus, set to -DF value, to make defocusing more intuitive
emdSTEM.probeDefocusC1 = -1 * emdSTEM.probeDefocusDF;

% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
emdSTEM.lambda = h/sqrt(2*m*e*emdSTEM.E0) ...
    /sqrt(1 + e*emdSTEM.E0/2/m/c^2) * 10^10; % wavelength in A
emdSTEM.sigma = (2*pi/emdSTEM.lambda/emdSTEM.E0) ...
    *(m*c^2+e*emdSTEM.E0)/(2*m*c^2+e*emdSTEM.E0);


% Interpolated Fourier coordinates, with 0.5 AA aperture
N = emdSTEM.imageSize ./ emdSTEM.interpolationFactor / 2;
qx = makeFourierCoords(N(1),2 * emdSTEM.pixelSize(1));
qy = makeFourierCoords(N(2),2 * emdSTEM.pixelSize(2));
[emdSTEM.qyaInterp,emdSTEM.qxaInterp] = meshgrid(qy,qx);
emdSTEM.q2Interp = emdSTEM.qxaInterp.^2 + emdSTEM.qyaInterp.^2;
emdSTEM.q1Interp = sqrt(emdSTEM.q2Interp);
emdSTEM.qInterpTheta = atan2(emdSTEM.qyaInterp,emdSTEM.qxaInterp);

% Initial probe in the interpolated coordinate space
qProbe = emdSTEM.probeSemiangleArray / emdSTEM.lambda;
dqx = qx(2) - qx(1);
dqy = qy(2) - qy(1);
emdSTEM.PsiInit = min(max( ...
    (qProbe*emdSTEM.q1Interp - emdSTEM.q2Interp) ./ ...
    sqrt(dqx^2*emdSTEM.qxaInterp.^2 + dqy^2*emdSTEM.qyaInterp.^2) ...
    + 0.5, 0), 1);
emdSTEM.PsiInit(1,1) = 1;
PsiInitMask = abs(emdSTEM.PsiInit) > 0;
% normalize
emdSTEM.PsiInit(:) = emdSTEM.PsiInit / sqrt(sum(abs(emdSTEM.PsiInit(:)).^2));


% Select either regular PRISM or the new beam-partitioned PRISM
if isempty(emdSTEM.partitionNumberRings)
    % Conventional PRISM algorithm
    sub = PsiInitMask(:)>0;
    emdSTEM.beamNum = sum(sub);
    inds = find(sub);
    qxy = [emdSTEM.qxaInterp(inds) emdSTEM.qyaInterp(inds)];
    [~,indsOrder] = sort(sum(qxy.^2,2));
    inds = inds(indsOrder);
    
    emdSTEM.beamIndexInterp = zeros(N);
    emdSTEM.beamIndexInterp(inds) = 1:emdSTEM.beamNum;
    
else
    qxProbe = qProbe / dqx;
    qyProbe = qProbe / dqy;
    
    if isempty(emdSTEM.selectBrightestBeamsNum)
        % partition using equally spaced tiles
    qRing = (0:emdSTEM.partitionNumberRings)';
    numRingBeams = qRing*6;
    numRingBeams(1) = 1;
    qRadius = (qRing / emdSTEM.partitionNumberRings);

    
    % Loop over rings of beams
    p = zeros(sum(numRingBeams),2);
    if emdSTEM.partitionNumberRings > 0
        rgbPlot = zeros(sum(numRingBeams),3);
        for a0 = 1:length(qRadius)
            inds = (1:numRingBeams(a0)) + sum(numRingBeams(1:(a0-1)));
            theta = linspace(0,2*pi,numRingBeams(a0)+1);
            theta(end) = [];
            
            p(inds,1) = round(qxProbe * qRadius(a0) * cos(theta));
            p(inds,2) = round(qyProbe * qRadius(a0) * sin(theta));
            
            if a0 == 1
                rgbPlot(a0,:) = [1 1 1];
            else
                %             if mod(a0,2) == 0
                %                 sub = mod(inds,2) == 0;
                %                 rgbPlot(inds(sub),:) = repmat([1 0 0],[sum(sub) 1]);
                %                 rgbPlot(inds(~sub),:) = repmat([0 1 1],[sum(sub) 1]);
                %             else
                %                 sub = mod(inds,2) == 0;
                %                 rgbPlot(inds(sub),:) = repmat([1 0 1],[sum(sub) 1]);
                %                 rgbPlot(inds(~sub),:) = repmat([0 1 0],[sum(sub) 1]);
                %             end
                if mod(a0,3) == 2
                    sub = mod(inds,2) == 0;
                    rgbPlot(inds(sub),:) = repmat([1 0 0],[sum(sub) 1]);
                    rgbPlot(inds(~sub),:) = repmat([0 1 1],[sum(sub) 1]);
                elseif mod(a0,3) == 0
                    sub = mod(inds,2) == 0;
                    rgbPlot(inds(sub),:) = repmat([0 1 0],[sum(sub) 1]);
                    rgbPlot(inds(~sub),:) = repmat([1 0 1],[sum(sub) 1]);
                elseif mod(a0,3) == 1
                    sub = mod(inds,2) == 0;
                    rgbPlot(inds(sub),:) = repmat([0 0 1],[sum(sub) 1]);
                    rgbPlot(inds(~sub),:) = repmat([1 1 0],[sum(sub) 1]);
                end
            end
        end
    else
        p = [0 0];
        rgbPlot = [1 1 1];
    end
    
    else
        % partition from projected potential FFT spots
        potSumFFT = abs(fft2(mean(sum(emdSTEM.pot,3),4)));
        
        % find N+1 brightest FFT spots (+1 is for the center beam)
        [~,inds] = sort(potSumFFT(:),'descend');
        indsBeam = inds(1:(emdSTEM.selectBrightestBeamsNum+1));
        [xp,yp] = ind2sub(emdSTEM.imageSize,indsBeam);
        xp = mod(xp - 1 + emdSTEM.imageSize(1)/2, emdSTEM.imageSize(1)) ...
            - emdSTEM.imageSize(1)/2;
        yp = mod(yp - 1 + emdSTEM.imageSize(2)/2, emdSTEM.imageSize(2)) ...
            - emdSTEM.imageSize(2)/2;
        p = [xp yp];
        
        % colours
        rgbPlot = rand(emdSTEM.selectBrightestBeamsNum+1,3);
        rgbPlot(1,:) = [1 1 1];
    end
    
    % Protect interpolation against duplicate points:
    [p,indsKeep] = unique(p,'rows');
    rgbPlot = rgbPlot(indsKeep,:);
    % Output beams
    qxyBeams = p.*[dqx dqy];
    emdSTEM.beamVectors = qxyBeams;
    emdSTEM.beamNum = size(qxyBeams,1);
    
    % pixel basis
    [ya,xa] = meshgrid(1:N(2),1:N(1));
    
    % Coordinate system - spec up by cropping subset
    if emdSTEM.partitionNumberRings > 0
        if isempty(emdSTEM.selectBrightestBeamsNum)
            distMax = round(max(qxProbe,qyProbe) / emdSTEM.partitionNumberRings * 1.2);
        else
            distMax = round(max(qxProbe,qyProbe) * 1.2);
        end
        v = -distMax:distMax;
        
        % Calculate the partitions via level sets of each beam
        emdSTEM.beamPartition = zeros(N(1),N(2),emdSTEM.beamNum,'single');
        imageTemp = zeros(N(1:2));
        for a0 = 1:emdSTEM.beamNum
            pData = [p+N/2+1 zeros(size(p,1),1)];
            pData(a0,3) = 1;
            
            F = scatteredInterpolant(pData(:,1:2),pData(:,3));
            F.Method = 'natural';
            
            vx = v + pData(a0,1);
            vy = v + pData(a0,2);
            vx(vx<1) = [];
            vy(vy<1) = [];
            vx(vx>N(1)) = [];
            vy(vy>N(2)) = [];
            xSub = xa(vx,vy);
            ySub = ya(vx,vy);
            
            imageTemp(:) = 0;
            imageTemp(vx,vy) = reshape(max( ...
                F([xSub(:) ySub(:)]),0),...
                [length(vx) length(vy)]);
            
            emdSTEM.beamPartition(:,:,a0) = ifftshift(imageTemp);
            %         emdSTEM.beamPartition(:,:,a0) = ...
            %             ifftshift(max(reshape(F(xyAll),N),0));
        end
        
        if emdSTEM.partitionSigmoidal == true
            emdSTEM.beamPartition(:) = sin(emdSTEM.beamPartition(:)*(pi/2)).^2;
        end
        
        beamSum = sum(emdSTEM.beamPartition,3);
        emdSTEM.beamPartition(:) = ...
            emdSTEM.beamPartition ./ beamSum;
        emdSTEM.beamPartition(isnan(emdSTEM.beamPartition(:))) = 0;
        imagePlotPsi = emdSTEM.beamPartition .* ...
            (emdSTEM.PsiInit / max(abs(emdSTEM.PsiInit(:))));
        emdSTEM.beamPartition(:) = emdSTEM.beamPartition .* emdSTEM.PsiInit;
    else
        imagePlotPsi = repmat(emdSTEM.PsiInit / max(emdSTEM.PsiInit(:)),[1 1 3]);
        emdSTEM.beamPartition = emdSTEM.PsiInit;
    end
    

    
%     figure(41)
%     clf
%     Ip = fftshift(sum(emdSTEM.beamPartition,3));
%     sum(Ip(:)>0.5*max(Ip(:)))
%     imagesc(Ip)
%     axis equal off
%     colormap(gray(256))
    
%     % Subtract max from all beams to partition
%     beamMax = max(emdSTEM.beamPartition,[],3);
%     emdSTEM.beamPartition(:) = max(emdSTEM.beamPartition - beamMax + 0.5,0);
%     beamSum = sum(emdSTEM.beamPartition,3);
%     emdSTEM.beamPartition(:) = emdSTEM.beamPartition ./ beamSum;
%     emdSTEM.beamPartition(:) = emdSTEM.beamPartition .* emdSTEM.PsiInit;
    
    % Generate parent beams index array
    emdSTEM.beamIndexInterp = zeros(N);
    emdSTEM.beamIndexInterp(sub2ind(N,...
        mod(p(:,1),N(1))+1,...
        mod(p(:,2),N(2))+1 ...
        )) = 1:emdSTEM.beamNum;
    
    % find all nonzero elements in partitioned beam, size of each tile
    emdSTEM.beamInds = cell(emdSTEM.beamNum,2);
    emdSTEM.beamTileSize = zeros(emdSTEM.beamNum,2);
    scale = 1./max(emdSTEM.beamPartition(:));
    for a0 = 1:emdSTEM.beamNum
        sub = emdSTEM.beamPartition(:,:,a0) > 0;
        emdSTEM.beamInds{a0,1} = find(sub);
        emdSTEM.beamInds{a0,2} = emdSTEM.beamInds{a0,1} + (a0-1)*prod(N);
        
        emdSTEM.beamTileSize(a0,1) = sum(sub(:));
        emdSTEM.beamTileSize(a0,2) = sum(sum(emdSTEM.beamPartition(:,:,a0)))*scale;
    end
end



if nargout == 0
    
    figure(567)
    clf
    set(gcf,'outerposition',[1300 80 576 512+128])
    set(gcf,'color','w')
    
    if isempty(emdSTEM.partitionNumberRings )
        Ip = fftshift(emdSTEM.beamIndexInterp);
        Ip = repelem(Ip, ...
            emdSTEM.interpolationFactor(1),...
            emdSTEM.interpolationFactor(2));
         imagesc(Ip,...
            'xdata',fftshift(qy),...
            'ydata',fftshift(qx))
        colormap(jetBlack)
    else        
        imagePlot = zeros(N(1),N(2),3);
        for a0 = 1:size(rgbPlot,1)
            for a1 = 1:3
                imagePlot(:,:,a1) = imagePlot(:,:,a1) ...
                    + rgbPlot(a0,a1)*imagePlotPsi(:,:,a0);
            end
        end
        
        imagesc(circshift(imagePlot,N(1:2)/2),...
            'xdata',fftshift(qy),...
            'ydata',fftshift(qx))
%         imagesc(fftshift(sum( ...
%             emdSTEM.beamPartition .* ...
%             reshape(1:emdSTEM.beamNum,[1 1 emdSTEM.beamNum]),3)),...
%             'xdata',fftshift(qy),...
%             'ydata',fftshift(qx))
    end
    hold on
    t = linspace(0,2*pi,180+1);
    plot(cos(t)*qProbe,sin(t)*qProbe,'linewidth',2,'color','w')
    if ~isempty(emdSTEM.partitionNumberRings)
        scatter(emdSTEM.beamVectors(:,2),emdSTEM.beamVectors(:,1),...
            'k.','sizedata',500)
        scatter(emdSTEM.beamVectors(:,2),emdSTEM.beamVectors(:,1),'w.',...
            'sizedata',125)
    end
    hold off
    xlim([-1 1]*1)
    ylim([-1 1]*1)
    axis equal off
%     colormap(jetBlack)
%     caxis([0 1])
    title(['Number beams = ' num2str(emdSTEM.beamNum)])
    set(gca,'position',[0 0 1 0.95])
    drawnow;
    
    % xlim([0.4 0.6]*N(2)+[1 0])
    % ylim([0.4 0.6]*N(1)+[1 0])
    
else
    
    % Fourier coordinates for the full FOV (include the AA aperture regions)
    qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
    qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
    [qya,qxa] = meshgrid(qy,qx);
    q2 = qxa.^2 + qya.^2;
    
    
    % AA aperture reduced coordinate indices and mask
    emdSTEM.xAA = [(1:(emdSTEM.imageSize(1)/4)) ...
        ((1-(emdSTEM.imageSize(1)/4)):0)+emdSTEM.imageSize(1)];
    emdSTEM.yAA = [(1:(emdSTEM.imageSize(2)/4)) ...
        ((1-(emdSTEM.imageSize(2)/4)):0)+emdSTEM.imageSize(2)];
    emdSTEM.maskAA = false(emdSTEM.imageSize);
    emdSTEM.maskAA(emdSTEM.xAA,emdSTEM.yAA) = true;
    emdSTEM.maskAAinv = ~emdSTEM.maskAA;
    
    
    % Propagation operators
    emdSTEM.prop = exp( q2(emdSTEM.xAA,emdSTEM.yAA) * ...
        (-1i*pi*emdSTEM.lambda*emdSTEM.sliceThickness));
    %     if emdSTEM.sMatrixRefocus == true
    emdSTEM.propBack = exp( q2(emdSTEM.xAA,emdSTEM.yAA) * ...
        (-1i*pi*emdSTEM.lambda* ...
        (-emdSTEM.probeDefocusC1 - emdSTEM.cellDim(3) + emdSTEM.sliceThickness)));
    %     else
    %         emdSTEM.propBackFull = exp( q2 * ...
    %             (-1i*pi*emdSTEM.lambda* ...
    %             (emdSTEM.probeDefocusC1 - emdSTEM.cellDim(3) + emdSTEM.sliceThickness)));
    %     end
    emdSTEM.xTiltBasis = qxa(emdSTEM.xAA,emdSTEM.yAA) * ...
        (2i*pi*emdSTEM.lambda*( emdSTEM.cellDim(3) - emdSTEM.sliceThickness));
    emdSTEM.yTiltBasis = qya(emdSTEM.xAA,emdSTEM.yAA) * ...
        (2i*pi*emdSTEM.lambda*( emdSTEM.cellDim(3) - emdSTEM.sliceThickness));
    
    % expand from interpolated beam coordinates to the full field-of-view (FOV)
    beamIndexAA = zeros(emdSTEM.imageSize/2);
    vx = 1:emdSTEM.interpolationFactor(1):emdSTEM.imageSize(1)/2;
    vy = 1:emdSTEM.interpolationFactor(2):emdSTEM.imageSize(2)/2;
    beamIndexAA(vx,vy) = emdSTEM.beamIndexInterp;
    emdSTEM.beamIndex = zeros(emdSTEM.imageSize);
    emdSTEM.beamIndex(emdSTEM.xAA,emdSTEM.yAA) = beamIndexAA;
    
    % list of all beams, their indices and locations in the full size coordinate system
    emdSTEM.beamList = zeros(emdSTEM.beamNum,5);
    for a0 = 1:emdSTEM.beamNum
        [xp,yp] = find(emdSTEM.beamIndex == a0);
        inds = sub2ind(emdSTEM.imageSize,xp,yp);
        emdSTEM.beamList(a0,:) = [xp yp qxa(inds) qya(inds) inds];
    end
    
    % initialize scattering matrix in compact form, plane wave Psi
    emdSTEM.Scompact = zeros( ...
        prod(emdSTEM.imageSize)/4,...
        emdSTEM.beamNum,emdSTEM.numFP,...
        'single');
    Psi = zeros(emdSTEM.imageSize,'single');
    PsiAA = zeros(emdSTEM.imageSize / 2,'single');
    
    % Calculate all columns of the scattering matrix, for all FP configs
    if flagProgress == true
        reverseStr = ''; % initialize console message piece
    end
    for a0 = 1:emdSTEM.numFP
        trans = exp((1i*emdSTEM.sigma) * emdSTEM.pot(:,:,:,a0));
        for a1 = 1:emdSTEM.beamNum
            Psi(:) = 0;
            Psi(emdSTEM.beamList(a1,5)) = prod(emdSTEM.imageSize)/4;
            Psi(:) = ifft2(Psi);
            
            % Propagate through foil
            for a2 = 1:emdSTEM.numPlanes
                Psi(:) = fft2(Psi .* trans(:,:,a2));
                if a2 < emdSTEM.numPlanes
                    Psi(emdSTEM.xAA,emdSTEM.yAA) = ...
                        Psi(emdSTEM.xAA,emdSTEM.yAA) .* emdSTEM.prop;
                    Psi(emdSTEM.maskAAinv) = 0;
                    Psi(:) = ifft2(Psi);
                else
%                     if emdSTEM.sMatrixRefocus == true
                        Psi(emdSTEM.xAA,emdSTEM.yAA) = ...
                            Psi(emdSTEM.xAA,emdSTEM.yAA) .* emdSTEM.propBack;
%                     else
%                         Psi(emdSTEM.xAA,emdSTEM.yAA) = ...
%                             Psi(emdSTEM.xAA,emdSTEM.yAA) ...
%                             * emdSTEM.propBackFull(emdSTEM.beamList(a1,5));
%                     end
                end
            end
            PsiAA(:) = Psi(emdSTEM.maskAA);
            
            % Remove plane wave tilt for partitioned PRISM, inverse FFT2
            if ~isempty(emdSTEM.partitionNumberRings)
                PsiAA(:) = circshift(PsiAA, [1 1] - emdSTEM.beamList(a1,1:2));
            end
            PsiAA(:) = ifft2(PsiAA);
            
            % Save AA reduced plane wave into compact S-matrix
            emdSTEM.Scompact(:,a1,a0) = PsiAA(:);
            
            % Progress
            if flagProgress == true
            comp = (a1 / emdSTEM.beamNum ...
                + a0 - 1) / emdSTEM.numFP;
                msg = sprintf(['S-Matrix calculation is ' ...
                    sprintf('%0.2f',100*comp) ' percent complete']);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'),1,length(msg));
            end
        end
    end
end
if nargout > 0 && flagProgress == true 
    fprintf([reverseStr '']);
end

emdSTEM.time02Smatrix = toc;
end