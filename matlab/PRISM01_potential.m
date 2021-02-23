function [emdSTEM] = PRISM01_potential(atoms,cellDim,emdSTEM)
rng(1);  % for repeatable testing of thermal displacements
tic
% Colin Ophus - 2020 Sept
% 01 - calculate the projected potential from a collection of an orthogonal
%      unit cell, sliced up for multislice calculations.

% Inputs
% atoms -   N x 4 array containing cartesian atomic coordinates, with rows
%           given by [x y z atomic_number], for N atoms [positions in Angstroms]
% cellDim - Dimensions of the unit cell [x_cell y_cell z_cell] in Angstroms

flagAvgFP = true;   % Average all frozen phonon configurations
flagFlipCellDimZ = true;
flagProgress = 1;  % Display progress on console
% RMS atomic displacements (from the Debye waller coefficients).
u = ones(118,1) * 0.10;
if nargin < 3; emdSTEM = struct; end

% Realspace pixel size.
if ~isfield(emdSTEM,'pixelSize'); emdSTEM.pixelSize = 0.1; end   
% Radial distance to integrate atomic potentials. [Ang]
if ~isfield(emdSTEM,'potBound'); emdSTEM.potBound = 3; end
% Z sampling density for integrating atomic potentials. [Ang]
if ~isfield(emdSTEM,'potSamplingZ'); emdSTEM.potSamplingZ = 0.1; end
% [min max] units of Nyquist, so 0 to 1
if ~isfield(emdSTEM,'potBandLimit'); emdSTEM.potBandLimit = [0.75 0.95]; end
 % Number of frozen phonon configurations.
 if ~isfield(emdSTEM,'numFP'); emdSTEM.numFP = 1; end
 % Thickness of each potential slice.
if ~isfield(emdSTEM,'sliceThickness'); emdSTEM.sliceThickness = 2; end
% PRISM interpolation factor
if ~isfield(emdSTEM,'interpolationFactor'); emdSTEM.interpolationFactor = [1 1]*1; end

% Save atomic positions in struct
emdSTEM.atoms = atoms;
emdSTEM.cellDim = cellDim;

% Simulation size
f = 2*2*emdSTEM.interpolationFactor;  % Force size to be divisible by 4*f
emdSTEM.imageSize = round(cellDim(1:2)/emdSTEM.pixelSize./f).*f;
emdSTEM.pixelSize = cellDim(1:2) ./ emdSTEM.imageSize;

% Construct projected potentials
xyLeng = ceil(emdSTEM.potBound./emdSTEM.pixelSize);
xvec = -xyLeng(1):xyLeng(1);
yvec = -xyLeng(2):xyLeng(2);
xr = xvec*emdSTEM.pixelSize(1);
yr = yvec*emdSTEM.pixelSize(2);
zLeng = emdSTEM.potBound / emdSTEM.potSamplingZ;
zvec = (0.5-zLeng):(zLeng-0.5);
zr = zvec * emdSTEM.potSamplingZ;

% dzPot = emdSTEM.sliceThickness / emdSTEM.sliceThicknessFactor;
% zLeng = ceil(emdSTEM.potBound./dzPot);
% zvecCalc = -zLeng:zLeng;
% zr = zvecCalc*dzPot;
% zvec = (-zLeng-1):(zLeng+1);

% Lookup table for atom types
atomTypes = unique(atoms(:,4));
emdSTEM.potLookup = zeros( ...
    length(xvec),...
    length(yvec),...
    length(zvec),...
    length(atomTypes));
emdSTEM.uLookup = u(atomTypes);
for a0 = 1:length(atomTypes)
    % store potentials in FFT space in x and y
    emdSTEM.potLookup(:,:,:,a0) = ...
        fft2(projPotLookup(atomTypes(a0),xr,yr,zr));
end

% shift operators, bandwidth limit, circular aperture in realspace?
[qxa,qya] = makeFourierCoords([length(xvec) length(yvec)],1);
qxShift = -2i*pi*qxa;
qyShift = -2i*pi*qya;
% qBandLimit = (qxa.^2 + qya.^2) < (emdSTEM.potBandLimit/2)^2;
q1 = sqrt(qxa.^2 + qya.^2);
qBandLimit = (emdSTEM.potBandLimit(2) - 2*q1) ...
    / (emdSTEM.potBandLimit(2) - emdSTEM.potBandLimit(1));
qBandLimit(:) = min(max(qBandLimit,0),1);
qBandLimit(:) = sin(qBandLimit*(pi/2)).^2;
% realspace band limit (circle function)
rBandLimit = (xvec(:)/(xyLeng(1)+0.5)).^2 ...
    + (yvec(:)'/(xyLeng(2)+0.5)).^2 <= 1;

% Generate projected potentials for all atoms /  frozen phonon configs
emdSTEM.numPlanes = round(cellDim(3) / emdSTEM.sliceThickness);
Npot = [emdSTEM.imageSize(1:2)  emdSTEM.numPlanes emdSTEM.numFP];
emdSTEM.pot = zeros(Npot,'single');

% Main loop over atoms
if flagProgress == true
%     progressbar(0,2);
    comp = 0;
end
for a0 = 1:size(atoms,1)
    xyzID = atoms(a0,1:4);
    [~,ind] = min(abs(atomTypes-xyzID(4)));
    u = emdSTEM.uLookup(ind);
    
    if flagFlipCellDimZ == true
        xyzID(3) = emdSTEM.cellDim(3) - xyzID(3);
    end
    
    for a1 = 1:emdSTEM.numFP
        
        xyz = xyzID(1:3) + randn(1,3)*u;
        xPx = xyz(1) / emdSTEM.pixelSize(1);
        yPx = xyz(2) / emdSTEM.pixelSize(2);
        x0 = round(xPx);
        y0 = round(yPx);
        dxPx = xPx - x0;
        dyPx = yPx - y0;
        
        % pixel coordinates
        x = mod(x0 + xvec,emdSTEM.imageSize(1)) + 1;
        y = mod(y0 + yvec,emdSTEM.imageSize(2)) + 1;
        z = round((xyz(3) + zr) / emdSTEM.sliceThickness) + 1;
        z(:) = min(max(z,1),emdSTEM.numPlanes);
        
        % Collapse along z planes, apply subpixel shifts, output into potentials
        zVals = unique(z);
        for a2 = 1:length(zVals)
            emdSTEM.pot(x,y,zVals(a2),a1) = ...
                emdSTEM.pot(x,y,zVals(a2),a1) ...
                + rBandLimit .* ifft2( ...
                sum(emdSTEM.potLookup(:,:,z==zVals(a2),ind),3) ...
                .* exp(qxShift*dxPx + qyShift*dyPx) ...
                .* qBandLimit,...
                'symmetric');
        end
    end
    
    if flagProgress == true
        if mod(a0,100) == 0
            comp = a0 / size(atoms,1);
%             progressbar(comp,2);
        end
    end
end
if flagProgress == true && comp < 1
%     progressbar(1,2);
end


if flagAvgFP == true
    emdSTEM.pot = mean(emdSTEM.pot,4);
    emdSTEM.numFP = 1;
end

emdSTEM.time01potential = toc;
end