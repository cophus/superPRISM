function [] = testPRISM04partition()

% Colin Ophus - 2021 Feb
% 04 - testing of partitioned PRISM simulations

% input .mat file for potentials, output of final results
test_file_name = 'Au_deca_pot.mat';
output_file_base = 'Au_deca_partition';

% import potentials
load(test_file_name);

% interp factor testing
emdSTEM.interpolationFactor = [1 1]*1; 
% emdSTEM.interpolationFactor = [1 1]*2; 
% emdSTEM.interpolationFactor = [1 1]*5; 
emdSTEM.partitionNumberRings = 0; 
% emdSTEM.partitionNumberRings = 1; 
% emdSTEM.partitionNumberRings = 2; 
% emdSTEM.partitionNumberRings = 4; 
% emdSTEM.partitionNumberRings = 8; 


% Probe positions
dxy = emdSTEM.cellDim(1:2) / 40;
xR = [0 1]*emdSTEM.cellDim(1);
yR = [0 1]*emdSTEM.cellDim(2);
emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);

% Other inputs
emdSTEM.E0 = 80e3;     
emdSTEM.probeSemiangleArray = 20 / 1000; 
emdSTEM.probeDefocusDF = 0; 
emdSTEM.partitionSigmoidal = false; 
emdSTEM.flagOutput3D = true; 
emdSTEM.flagOutput4D = false; 
emdSTEM.drBins3D = 1 / 1000; 
emdSTEM.flagProbePositionsNearestPixel = true;


% Run PRISM simulation
emdSTEM = PRISM02_Smatrix(emdSTEM);
emdSTEM = PRISM03_probes(emdSTEM);

% Output struct
emd = emdOutput(emdSTEM);

% save file
filename = [output_file_base '_interp' ...
    '_' num2str(emdSTEM.interpolationFactor(1)) ...
    '_' num2str(emdSTEM.interpolationFactor(2)) ...
    '_rings' ...
    '_' num2str(emdSTEM.partitionNumberRings) ...
    '.mat'];
save(filename,'emd');



end