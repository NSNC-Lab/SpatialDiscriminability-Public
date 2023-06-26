function figuresize(width, height, h, units)
% FIGURESIZE(WIDTH, HEIGHT, H, UNITS)
%    WIDTH and HEIGHT are the size in UNITS. H is handle to figure (defaults to
%    gcf). For proper size display onscreen, you need to use the command,
%    SET(0, 'ScreenPixelsPerInch', DPI), where DPI is the number of pixels per
%    inch of your monitor (which you use a ruler to calculate).  For a 19-inch
%    monitor with 1280x1024 resolution, this value is approximately 86.23.

if nargin < 2
	error('At least two arguments are needed.')
end
if nargin < 3
    h = gcf;
end
if nargin < 4
	units = 'inches';
	if isempty(h)
		h = gcf;		
	end
end

% set(h, 'WindowStyle', 'normal');
drawnow
set(h, 'units', units, 'PaperUnits', units, 'resize', 'off', 'color', [1 1 1])
set(h, 'position', [0.1 0.1 width height], 'papersize', [width height], 'paperposition', [0 0 width height])
drawnow
