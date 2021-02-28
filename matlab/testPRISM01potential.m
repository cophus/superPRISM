function [] = testPRISM01potential()

% Colin Ophus - 2021 Feb
% 01 - This testing script generates and saves the projected potentials,
% intended to be reused for testing. 

% emd struct .mat output for testing
test_file_name = 'Au_deca_pot.mat';

% load the decahedral nanoparticle testing
load('inputData01.mat')

% Inputs:
emdSTEM.pixelSize = 0.1;
emdSTEM.potBound = 3;
emdSTEM.potSamplingZ = 0.1; 
emdSTEM.potBandLimit = [0.75 0.95]; 
emdSTEM.numFP = 1; 
emdSTEM.sliceThickness = 2; 
emdSTEM.interpolationFactor = [1 1]*1; 

% Compute potentials
emdSTEM = PRISM01_potential(atoms,cellDim,emdSTEM);

% Save .mat output
save(test_file_name,'emdSTEM');

end