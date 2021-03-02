function [emd] = emdOutput(emdSTEM)

% Colin Ophus - 2021 Mar
% Function for saving on the 3D outputs and most important parameters of a
% PRISM or partioned PRISM simulation.

emd.output3D = emdSTEM.output3D;
emd.xp = emdSTEM.xp;
emd.yp = emdSTEM.yp;
emd.probeSemiangleArray = emdSTEM.probeSemiangleArray;
emd.E0 = emdSTEM.E0;
emd.detectorAngles = emdSTEM.detectorAngles;
emd.interpolationFactor = emdSTEM.interpolationFactor;
emd.cellDim = emdSTEM.cellDim;
emd.pixelSize = emdSTEM.pixelSize;
emd.pixelSizeAA = emdSTEM.pixelSizeAA;

% size of the potential array and scattering matrix
if isfield(emdSTEM,'pot')
    emd.sizePot = size(emdSTEM.pot);
end
if isfield(emdSTEM,'Scompact')
    emd.sizeScompact = size(emdSTEM.Scompact);
end
% Partioning
if isfield(emdSTEM,'partitionNumberRings')
    emd.partitionNumberRings = emdSTEM.partitionNumberRings;
end



% timing
if isfield(emdSTEM,'time01potential')
    emd.time01potential = emdSTEM.time01potential;
end
if isfield(emdSTEM,'time02multislice')
    emd.time02multislice = emdSTEM.time02multislice;
end
if isfield(emdSTEM,'time02Smatrix')
    emd.time02Smatrix = emdSTEM.time02Smatrix;
end
if isfield(emdSTEM,'time03probes')
    emd.time03probes = emdSTEM.time03probes;
end


end

