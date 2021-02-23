function [imageRGB] = colorComplex(EW,ampRange)

imageRGB = ones(size(EW,1),size(EW,2),3);
imageRGB(:,:,1) = mod(angle(EW)/(2*pi),1);
imageRGB(:,:,3) = min(max( ....
     (abs(EW) - ampRange(1)) / (ampRange(2) - ampRange(1)),0),1);
imageRGB(:) = hsv2rgb(imageRGB);

end