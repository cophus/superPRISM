function [cmap] = cyanToRed(Nc)

% Colin Ophus - 2021 Feb
% Returns a color map with these entries:

% If user doesn't specify number of colours, assume 256 (8 bit)
if nargin == 0
    Nc = 1024;
end

% Make custom colormap
% [R G B fraction]
c = [ 
    0       1       1       0.0
    0       0.7     1       0.2
    1       1       1       0.5
    1       0       0       0.8
    0.5       0       0       1.0];
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
% L = [linspace(0,100,Nc/2) fliplr(linspace(0,100,Nc/2))];
c = linspace(-1,1,Nc);
% L = 100 * (1 - abs(c).^1.5);
% L = 100 * (1 - abs(c).^1.5);
L = 100 * cos(c * pi/2);

cmap3 = reshape(cmap,[Nc 1 3]);
cmap3(:) = rgb2lab(cmap3);
cmap3(:,:,1) = L;
cmap3(:) = lab2rgb(cmap3);
cmap(:) = cmap3;


if nargout == 0
    figure(1)
    clf
    imagesc(cmap3)
    axis off
end

end