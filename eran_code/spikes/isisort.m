% isisort               sort spikes into ISI groups
%
% call                  [ xc bedges cpb ] = isisort( x, bedges, x2 )
% 
% gets                  x           spike train/s. can be either a sparse matrix or a vector of spike times
%                       n           number of ISI bins (scalar) or an n+1 vector of ISI edges
%                       x2          (optional) reference spike train (cross-cell model, e.g. xsta or xwiener)
%
% returns               xc          classified spike trains, either a sparse matrix of a vector
%                                       (depending on the input format)
%                       bedges      edges; (n+1 x 1 vector)
%                       cpb         summary of counts per bin (n x 1 vector)
%
% calls                 make_equal_bins             (general)
%                       calc_isis                   (spikes) 
%
% see also              streconstruct

% 27-oct-14 ES

% revisions
% 30-oct-14 accounted for bedges empty/<2 case
% 14-oct-19 cleaned up, documented

function [ xc, bedges, cpb ] = isisort( x, bedges, x2 )

nargs                   = nargin;
if nargs < 1 || isempty( x ) || issparse( x ) && sum( x( : ) ) == 0
    xc                  = x;
    bedges              = [];
    cpb                 = [];
    return
end
if nargs < 2 || isempty( bedges )
    bedges              = 1;
end
if nargs < 3
    x2                  = [];
end

ftype                   = '';
if isempty( x2 )
    if issparse( x ) || isvector( x )
        ftype           = 'cmodel';
    end
else
    if issparse( x ) && issparse( x2 ) && isequal( size( x ), size( x2 ) ) ...
            || ~issparse( x ) && ~issparse( x2 ) && isvector( x ) && isvector( x2 )
        ftype           = 'xmodel';
    end
end

% 1. calc the ISIs
if isequal( ftype, 'cmodel' )
    isis                = calc_isis( x );
else
    isis                = calc_isis( x, x2 );
end

% 2. get the ISI edges
if isempty( bedges ) || length( bedges ) == 1 && bedges( 1 ) < 2
    xc                  = x;
    cpb                 = [];
    return
elseif length( bedges ) == 1
    nbins               = bedges;
    [ sval, eval ]      = make_equal_bins( isis, nbins );
    nbins               = length( sval );
    eval( nbins )       = inf;
    bedges              = [ sval; eval( nbins ) ];
else
    if any( diff( round( bedges ) ) <= 0 )
        error( 'bedges must be a vector of increasing integers' )
    end
    nbins               = length( bedges ) - 1;
end
isis( isnan( isis ) )   = inf;

% 3. sort the spikes (first spike in trial gets inf)
[ cpb, sidx ]           = histc( isis, bedges );
cpb( nbins )            = sum( cpb( end + [ -1 0 ] ) );
cpb( nbins + 1 )        = [];
sidx( sidx == nbins + 1 ) = nbins;

if issparse( x )
    xc = x;
    xc( x( : ) )        = sidx; % keep the classified trains
else
    xc                  = sidx;
end

return

% EOF
