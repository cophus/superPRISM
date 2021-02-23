function [emd] = emdOutput(emdSTEM)


emd.output3D = emdSTEM.output3D;
emd.xp = emdSTEM.xp;
emd.yp = emdSTEM.yp;
emd.probeSemiangleArray = emdSTEM.probeSemiangleArray;
emd.E0 = emdSTEM.E0;
emd.detectorAngles = emdSTEM.detectorAngles;
emd.interpolationFactor = emdSTEM.interpolationFactor;


end

