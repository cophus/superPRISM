function [emd] = emdOutput(emdSTEM)

emd.output3D = emdSTEM.output3D;
emd.xp = emdSTEM.xp;
emd.yp = emdSTEM.yp;
emd.probeSemiangleArray = emdSTEM.probeSemiangleArray;
emd.E0 = emdSTEM.E0;
emd.detectorAngles = emdSTEM.detectorAngles;
emd.interpolationFactor = emdSTEM.interpolationFactor;

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

