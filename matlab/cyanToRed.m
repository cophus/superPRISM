function [cmap] = cyanToRed(Nvalues)

% Colin Ophus - 2021 Feb
% Returns a color map with these entries:

% If user doesn't specify number of colours, assume 256 (8 bit)
if nargin == 0
    Nvalues = 256;
end

% coordinate system
x = linspace(0,1,Nvalues)';
% Red
r = -4*(x-1/2).^2 + 1;
r(x>=1/2) = 1;
% % Green
g = -4*(x-1/2).^2 + 1;
% Blue
b = -4*(1/2 - x).^2 + 1;
b(x<=1/2) = 1;
% Cyan map - set green to blue, or maybe average?
g = (b + g)/2;

% Assemble colormap
cmap = [r g b];

% % Testing
% % r = -16*x.^3 + 12*x.^2;
% % g = 16*(x-1/2).^4 - 8*(x-1/2).^2 + 1;
% % b = -16*(1-x).^3 + 12*(1-x).^2;
% figure(1)
% clf
% imagesc(reshape(cmap,[Nvalues 1 3]));
% % hold on
% % plot(x,r,'linewidth',2,'color','r')
% % plot(x,g,'linewidth',2,'color','g')
% % plot(x,b,'linewidth',2,'color','b')
% % hold off


end