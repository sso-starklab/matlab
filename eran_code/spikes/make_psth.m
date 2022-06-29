% MAKE_PSTH             from spike and trigger times
%
% call                  [ psth, bins, rast, st ] = make_psth( spike_times, trigs, win, bs, normf, tf )
%
% gets                  spike_times         spike times
%                       trigs               trigger times (same sampling rate as ST)
%                       win                 relevant window (in samples; same sampling rate as ST)
%                       bs                  new bin size (in samples; same sampling rate as ST)
%                       normf               {1}; normalizes to counts/bin/trigger; otherwise counts/bin
%                       tf                  {1}; truncates to max 1 count/bin/trial; otherwise may be as many as BS
%
% returns               psth                counts/bin/trigger or counts/bin
%                       bins                relative to 0 time (in samples)
%                       rast                sparse matrix, bins x trials
%                       st                  in samples relative to trigger
%
% calls                 uhist
%
% see also              multipeth

% 09-dec-09 ES

% revisions
% 05-may-10 rast and st outputs added
% 07-nov-10 (1) empty spike_times supported
%           (2) empty psth changed to psth of zeros 
% 18-mar-11 (1) bins and edges modified to match those in MAKE_PSTH_RASTER
%           (2) edges of window included as well
% 14-oct-19 cleaned up and documented

function [ psth, bins, mat, st ] = make_psth( spike_times, trigs, win, bs, normf, tf )

nargs                   = nargin;
if nargs < 4
    error( '4 arguments necesary' );
end
if nargs < 5 || isempty( normf )
    normf               = 1;
end
if nargs < 6 || isempty( tf )
    tf                  = 0; 
end

% allocate memory
spike_times             = spike_times( : );
trigs                   = trigs( : );
edges                   = [ fliplr( -bs/2 : -bs : ( win( 1 ) - bs ) )  bs/2 : bs : ( win( 2 ) + bs ) ]';
bins                    = mean( [ edges( 1 : end - 1 ) edges( 2 : end ) ], 2 );
mat                     = sparse( length( bins ), length( trigs ) );

% get the spike times relative to triggers
st                      = [];
for i                   = 1 : length( trigs )
    tmp                 = spike_times( spike_times >= ( trigs( i ) + edges( 1 ) ) & spike_times <= ( trigs( i ) + edges( end ) ) ) - trigs( i );
    st                  = [ st; tmp ];
    idx                 = min( ceil( ( tmp - edges( 1 ) + 1 ) / bs ), length( bins ) );
    if tf
        mat( idx, i )   = 1;
    elseif ~isempty( idx )
        [ aa, bb ]      = uhist( idx );
        mat( bb, i )    = aa;
    end
end

% build PSTH
psth                    = histc( st, edges );
if isempty( psth )
    psth                = zeros( length( bins ), 1 );
else
    psth( end )         = [];
    psth                = psth( : );
end
if normf
    psth                = psth / length( trigs );
end

return

% EOF
