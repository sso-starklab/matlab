% plotranges    of a 2-column matrix
%
% call          lh = plotranges( mat, ypos, Fs, lineargs )
%
% mat           a 2-column matrix, each row is set of two numbers
%               indicating a range. does not have to be sorted
% ypos          vertical position to plot
% Fs            conversion for horizontal scale 
% lineargs      passed to line; e.g. ..., 'color', [ 1 0 0 ] ) will plot in red
% 
% see also: dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, setdiffranges, sortranges, uniteranges

% 27-jan-12

% revisions
% 03-mar-13 efficiency

% note: consider changing plot mode to patches

function [ lh ] = plotranges( mat, ypos, Fs, varargin )

if nargin < 2 || isempty( ypos ), ypos = 1; end
if nargin < 3 || isempty( Fs ), Fs = 1; end

if size( mat, 2 ) ~= 2
    error( 'format should be 2-column matrices' )
end
if isempty( mat )
    lh = NaN;
    return
end
mat = sortranges( mat );
mat = round( mat * Fs );

m = size( mat, 1 );
lh = zeros( m, 1);
for i = 1 : m
    t = mat( i, 1 ) : mat( i, 2 );
    lh( i ) = line( t, ypos * ones( size( t ) ), varargin{ : } );
end
xlim( [ 1 mat( end, 2 ) ] / Fs );

return

% EOF
