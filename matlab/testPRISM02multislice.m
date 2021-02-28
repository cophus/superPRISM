function [] = testPRISM02multislice()

% Colin Ophus - 2021 Feb
% 02 - testing of multislice simulations

% input .mat file for potentials, output of final results
test_file_name = 'Au_deca_pot.mat';
output_file_base = 'Au_deca_multislice';

% import potentials
load(test_file_name);

% interp factor testing
emdSTEM.interpolationFactor = [1 1]*1; 
% emdSTEM.interpolationFactor = [1 1]*2; 
% emdSTEM.interpolationFactor = [1 1]*5; 

% Probe positions
dxy = emdSTEM.cellDim(1:2) / 400;
xR = [0 1]*emdSTEM.cellDim(1);
yR = [0 1]*emdSTEM.cellDim(2);
emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);

% Other inputs
emdSTEM.flagOutput3D = true; 
emdSTEM.flagOutput4D = false; 
emdSTEM.E0 = 80e3; 
emdSTEM.probeSemiangleArray = 20 / 1000; 
emdSTEM.probeDefocusDF = 0; 
emdSTEM.drBins3D = 1 / 1000; 


% Run multislice simulation
emdSTEM = PRISM02_multislice(emdSTEM);

% Output struct
emd = emdOutput(emdSTEM);

% save file
filename = [output_file_base '_interp' ...
    '_' num2str(emdSTEM.interpolationFactor(1)) ...
    '_' num2str(emdSTEM.interpolationFactor(2)) ...
    '.mat'];
save(filename,'emd');

end