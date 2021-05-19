
% generate dummy surface
nX = 500; nY = 600;
surfaceMin = 4e-3; surfaceMax = 7e-3;
x = 0:(size(surface, 1) - 1);
x = x - mean(x); % center x vector around 0
y = 0:(size(surface, 2) - 1);
y = y - mean(y); % center y vector around 0
radMap = sqrt(x.^2 * (y').^2); % generate radial map 
size(radMap)
surface = -radMap.^2;


dx = single(25e-6); % resolution in x direction
dy = single(25e-6); % resolution in y direction

fluenceField = zeros(500, 200, 500, 'single'); % iz, ir, iSurface
dr = single(5e-6);


res = single([dr, dx, dy]);
fluenceVol = buildFluenceVolume(surface, fluenceField, res);
