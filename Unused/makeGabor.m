function [] = makeGabor(lambda,imSize,thetaStep,sigma,phase,trim)

% INTERACTIVE
resp = input('Condition? 1-9 -> ');
theta = thetaStep(resp);

% DEFAULT PARAMETERS
if nargin < 6
    trim = 0.005; end
if nargin < 5
    phase = 0.5; end
if nargin < 4
    sigma = 15; end
if nargin < 3
    theta = 0; end
if nargin < 2
    imSize = 150; end

% if nargin < 2
%     imSize = 150;                           % image size: n X n
% elseif nargin < 3
%     theta = 0;                             % grating orientation
% elseif nargin < 4
%     sigma = 15; phase = .5; trim = .005;     % gaussian standard deviation in pixels / size of gabor patch
% elseif nargin < 5
%     phase = .5; trim = .005;                           % phase (0 -> 1)
% elseif nargin < 6
%     trim = .005;                            % trim off gaussian values smaller than this
% end

% debugging only
% trim = 0.005; phase = .5; sigma = 15; theta = 0; imSize = 150; theta = thetaStep(1);

% make linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5

% manipulate wavelength/phase
freq = imSize/lambda;                   % compute frequency from wavelength
Xf = X0 * freq * 2*pi;                  % convert X to radians: 0 -> ( 2*pi * frequency)
sinX = sin(Xf) ;                        % make new sinewave
% plot(sinX, 'r-');                       % plot in red
phaseRad = phase * 2* pi;               % convert to radians: 0 -> 2*pi
sinX = sin( Xf + phaseRad) ;            % make phase-shifted sinewave
% hold on;                                % superimpose next plot on last
% plot(sinX, 'g-');                       % plot in green
% hold off;                               % next plot overwrites this one

% 2D grating
[Xm Ym] = meshgrid(X0, X0);             % 2D matrices
% imagesc([Xm Ym]);                       % display Xm and Ym
% colorbar; axis off                      % add colour bar to see values

% put colormap thru sine
Xf = Xm * freq * 2*pi;
grating = sin( Xf + phaseRad);          % make 2D sinewave
% imagesc( grating, [-1 1] );             % display
% colormap gray(256);                     % use gray colormap (0: black, 1: white)
% axis off; axis image;                   % use gray colormap

% change orientation
thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt = [ Xt + Yt ];                      % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phaseRad);                   % make 2D sinewave
% imagesc( grating, [-1 1] );                     % display
% axis off; axis image;    % use gray colormap

% make a gaussian mask
s = sigma/imSize;                                   % gaussian width as fraction of imageSize
Xg = exp(-(((X0.^2))./(2*s^2)));                    % formula for 1D gaussian
Xg = normpdf(X0,0,(20/imSize)); Xg = Xg/max(Xg);    % alternative using normalized probability function (stats toolbox)
% plot(Xg)

% make 2D gaussian blob
gauss = exp(-(((Xm.^2)+(Ym.^2))./(2*s^2)));         % formula for 2D gaussian
% imagesc( gauss,[-1 1]);                             % display
% axis off; axis image;                               % use gray colormap

% multiply grating & gaussian to get a gabor
gauss(gauss < trim) = 0;                            % trim around edges (for 8-bit colour displays)
gabor = grating.* gauss;                            % use .* dot-product
imagesc( gabor,[-1 1]);                             % display
colormap gray(256);
axis off; axis image;                               % use gray colormap
axis image; axis off; colormap gray(256);
set(gca,'pos', [0 0 1 1]);                          % display nicely without borders
set(gcf, 'menu', 'none', 'Color',[.5 .5 .5]);       % without background

end