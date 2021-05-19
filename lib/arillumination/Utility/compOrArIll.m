% File: compOrArIll.m
% Author: Urs Hofmann
% Mail: hofmannu@bomed.ee.ethz.ch
% Date: 04.05.2020

% Description: simulates different fiber types and penetration depths

clear all;
close all;

[mua532, mus532] = getAbsorptionCoefficient(532);
[mua578, mus578] = getAbsorptionCoefficient(578);

tissue532 = single([mua532, mus532, 0.8, 1.41]);
tissue578 = single([mua578, mus578, 0.8, 1.41]);
% [0] absorption coefficient [1/m]
% [1] scattering coefficient [1/m]
% [2] anisotroy factor
% [3] optical index of tissue

fiberAR = single([200e-6 / 2, 0.22]);
fiberOR = single([10e-6 / 2, 0]);
% [0] core radius of multimode fiber
% [1] numerical aperture

field = single([3e-3, 15e-3, 7e-3, 5e-6, 5e-6]);
% [0] maximum simulated radial distance [m]
% [1] distance between lower field end and fiber output [m]
% [2] distance between fiber output and skin surface [m]
% [3] spatial resolution of output field in radial direction [m]
% [4] spatial resolution of output field in axial direction [m]

simOR = single([500, 1e9]);
simAR = single([500, 50e9]);
% [0] number of photons per thread
% [1] overall number of photons to simulate

tic
outputOR532 = m_mc(tissue532, fiberOR, field, simOR);
outputAR532 = m_mc(tissue532, fiberAR, field, simAR);
outputOR578 = m_mc(tissue578, fiberOR, field, simOR);
outputAR578 = m_mc(tissue578, fiberAR, field, simAR);
toc

outputOR532(1, 1) = outputOR532(2, 1);
outputAR532(1, 1) = outputAR532(2, 1);
outputOR578(1, 1) = outputOR578(2, 1);
outputAR578(1, 1) = outputAR578(2, 1);

[nR, nZ] = size(outputOR532);
r = ((0:(nR - 1)) + 0.5)* field(4); % radial vector
z = ((0:(nZ - 1)) + 0.5)* field(5); % axial vector

figure()

ax1 = subplot(2, 2, 1)
imagesc(r, z, log(outputAR532'));
xlabel('Radial distance [m]');
ylabel('Axial distance [m]');
title('Fluence field AR illumination');
colormap(ax1, bone(1024));
colorbar;
axis image;

ax2 = subplot(2, 2, 2)
imagesc(r, z, log(outputOR532'));
xlabel('Radial distance [m]');
ylabel('Axial distance [m]');
title('Fluence field OR illumination');
colormap(ax2, bone(1024));
colorbar;
axis image;

subplot(2, 2, 3)
yyaxis left;
plot(z, outputAR532(1, :));
hold on
plot(z, outputAR578(1, :));
hold off;
axis tight;
grid on;

xlabel('Axial distance [m]');
ylabel('Fluence');
title('Central axis AR illumination');
yyaxis right
deltaAR = abs(outputAR532(1, :) - outputAR578(1, :)) ./ ...
	((outputAR532(1, :) + outputAR578(1, :)) / 2);
plot(z, deltaAR);
ylabel('Delta / meanFluence');
legend('532 nm', '578 nm', 'relDelta');

subplot(2, 2, 4)
yyaxis left
plot(z, outputOR532(1, :));
hold on
plot(z, outputOR578(1, :));
hold off;
axis tight;
grid on;

xlabel('Axial distance [m]');
ylabel('Fluence');

yyaxis right
deltaOR = abs(outputOR532(1, :) - outputOR578(1, :)) ./ ...
	((outputOR532(1, :) + outputOR578(1, :)) / 2);
plot(z, deltaOR);
ylabel('Delta / meanFluence');
title('Central axis OR illumination');
legend('532 nm', '578 nm', 'relDelta');
