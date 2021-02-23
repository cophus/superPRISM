function [cmap] = jetBlackSqrt(Nc)

% Colin Ophus - 2020 Feb
% Square root scaling of jetblack

scale = 16;
powerLaw = 0.5;

if nargin < 1
    Nc = 1024;
end

NcScale = Nc * scale;
uInput = linspace(0,1,NcScale)';
uOutput = linspace(0,1,Nc)' .^ powerLaw;


cmapScale = jet(NcScale);
b = round(NcScale/8);
c = linspace(0,1,b)';
cmapScale(1:b,:) = cmapScale(1:b,:) .* c;

cmap = zeros(Nc,3);
for a0 = 1:3
    cmap(:,a0) = interp1(uInput,cmapScale(:,a0),uOutput,'linear');
end


end