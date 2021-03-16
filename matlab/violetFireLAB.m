function [cmap] = violetFireLAB(Nc)

if nargin == 0
    Nc = 256;
end

% Make custom colormap
% [R G B fraction]
c = [ 
    0       0       0       0;
    0.2     0       0.4     0.3;
    1       0       0       0.5;
    1       0.6     0       0.7;
    1       0.8       0     0.85;
    1       1       1       1];
cmap = zeros(Nc,3);
for a0 = 1:(size(c,1)-1)
    f1 = round(c(a0,4)*Nc+1);
    f2 = round(c(a0+1,4)*Nc);
    inds = f1:f2;
    cnew = [linspace(c(a0,1),c(a0+1,1),length(inds))' ...
        linspace(c(a0,2),c(a0+1,2),length(inds))' ...
        linspace(c(a0,3),c(a0+1,3),length(inds))'];
    cmap(inds,:) = cnew;
end



% Smoothing
sigma = 11 * Nc / 256;
cmap(:,1) = smooth(cmap(:,1),sigma,'moving');
cmap(:,2) = smooth(cmap(:,2),sigma,'moving');
cmap(:,3) = smooth(cmap(:,3),sigma,'moving');
cmap(cmap<0) = 0;
cmap(cmap>1) = 1;

% LAB lightness
% L = linspace(0,100,Nc);
L = 100 * linspace(0,1,Nc).^1.2;
cmap3 = reshape(cmap,[Nc 1 3]);
cmap3(:) = rgb2lab(cmap3);
cmap3(:,:,1) = L;
cmap3(:) = lab2rgb(cmap3);
cmap(:) = cmap3;

end