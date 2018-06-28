% function FastMaskedNoiseDemo(numRects, rectSize, scale)
% FastMaskedNoiseDemo([numRects=1][, rectSize=128][, scale=1])
%
% Demonstrates how to generate and draw noise patches on-the-fly in a fast
% way. The patches are shown through circular apertures by the use of
% alpha-blending.
%
% numRects = Number of random patches to generate and draw per frame.
%
% rectSize = Size of the generated random noise image: rectSize by rectSize
%            pixels. This is also the size of the Psychtoolbox noise
%            texture.
%
% scale = Scalefactor to apply to texture during drawing: E.g. if you'd set
% scale = 2, then each noise pixel would be replicated to draw an image
% that is twice the width and height of the input noise image. In this
% demo, a nearest neighbour filter is applied, i.e., pixels are just
% replicated, not bilinearly filtered -- Important to preserve statistical
% independence of the random pixel values!

% History:
% 4.11.2006 Written (MK).

% Abort script if it isn't executed on Psychtoolbox-3:
AssertOpenGL;

numRects = 3; 
rectSize = barLenPx; 
scale = 1; 

PsychDebugWindowConfiguration(0, 1);

% Open fullscreen onscreen window on that screen. Background color is
% gray, double buffering is enabled. Return a 'win'dowhandle and a
% rectangle 'winRect' which defines the size of the window.
[win, winRect] = Screen('OpenWindow', screenid, 128);

% Compute destination rectangle locations for the random noise patches:

% 'objRect' is a rectangle of the size 'rectSize' by 'rectSize' pixels of
% our Matlab noise image matrix:
objRect = SetRect(0,0, rectSize, rectSize);

% ArrangeRects creates 'numRects' copies of 'objRect', all nicely
% arranged / distributed in our window of size 'winRect':
dstRect = ArrangeRects(numRects, objRect, winRect);

% Now we rescale all rects
for i=1:numRects
    % Compute center position [xc,yc] of the i'th rectangle:
    [xc, yc] = RectCenter(dstRect(i,:));
    % Create a new rectange, centered at the same position, but 'scale'
    % times the size of our pixel noise matrix 'objRect':
    dstRect(i,:)=CenterRectOnPoint(objRect, xc, yc);
end

% Build a nice aperture texture: Offscreen windows can be used as
% textures as well, so we open an Offscreen window of exactly the same
% size 'objRect' as our noise textures, with a gray default background.
% This way, we can use the standard Screen drawing commands to 'draw'
% our aperture:
aperture=Screen('OpenOffscreenwindow', win, 128, objRect);

% First we clear out the alpha channel of the aperture disk to zero -
% In this area the noise stimulus will shine through:
Screen('FillOval', aperture, [255 255 255 0], objRect);

% Enable alpha blending
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Perform initial Flip and sync us to the retrace:
vbl = Screen('Flip', win);

% Generate and draw 'numRects' noise images:
for i=1:numRects
    % Compute 'noiseimg' noise image matrix with Matlab:
    % Normally distributed noise with mean 128 and stddev. 50, each
    % pixel computed independently with a size of rectSize x
    % rectSize noise pixels:
    noiseimg=(50*randn(rectSize, rectSize) + 128);
    
    % Convert it to a texture 'tex':
    tex=Screen('MakeTexture', win, noiseimg);
    
    % Draw the texture into the screen location defined by the
    % destination rectangle 'dstRect(i,:)'. If dstRect is bigger
    % than our noise image 'noiseimg', PTB will automatically
    % up-scale the noise image. We set the 'filterMode' flag for
    % drawing of the noise image to zero: This way the bilinear
    % filter gets disabled and replaced by standard nearest
    % neighbour filtering. This is important to preserve the
    % statistical independence of the noise pixels in the noise
    % texture! The default bilinear filtering would introduce local
    % correlations when scaling is applied:
    Screen('DrawTexture', win, tex, [], dstRect(i,:), [], 0);
    
    % Overdraw the rectangular noise image with our special
    % aperture image. The noise image will shine through in areas
    % of the aperture image where its alpha value is zero (i.e.
    % transparent):
    Screen('DrawTexture', win, aperture, [], dstRect(i,:), [], 0);
    Screen('Flip',win)
end