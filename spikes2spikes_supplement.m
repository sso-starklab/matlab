% spikes2spikes_supplement  add ASG and STC to s2s structure
%
% call                      s2s = spikes2spikes_supplement( s2s )
%
% gets and returns          s2s structure
%
% optional arguments (given as name/value pairs):
%                           BinSizeMS           { 1 }, see spikes2spikes
%                           halfWidthMS         { 50 }, see spikes2spikes
%                           jitWindowMS         { 5 }, see spikes2spikes
%                           supportEdges        { [ -1 3 ] }, i.e. keep
%                                                   h(t) only for the bins -5 : 15 
%
% calls                     calc_asg, ParseArgPairs
%
% see also                  spikes2spikes.m

% 31-may-20 ES & SSo

function s2s = spikes2spikes_supplement( s2s, varargin )

% argument handling
nargs                               = nargin;
if nargs < 1 || isempty( s2s )
    return
end
[ BinSizeMS, halfWidthMS, jitWindowMS ...
    , supportEdges ]                = ParseArgPairs( ...
    { 'BinSizeMS', 'halfWidthMS', 'jitWindowMS' ...
    , 'supportEdges' } ...
    , { 1, 50, 5 ...
    , [ -1 3 ] }...
    , varargin{ : } );

% compute CCH time lag vector
nBins                               = halfWidthMS / BinSizeMS;                      % number of bins on each side of the CCH
t                                   = ( -nBins : nBins )' * BinSizeMS;              % [ms]
idx                                 = ( supportEdges( 1 ) * jitWindowMS ) : ( supportEdges( 2 ) * jitWindowMS );
sidx                                = find( ismember( t, idx ) );
nsupport                            = length( sidx );

% allocate space
nunits                              = size( s2s.shankclu, 1 );
g1mat                               = NaN( nunits, nunits );
g2mat                               = NaN( nunits, nunits );
g1tidx                              = NaN( nsupport, nunits, nunits );
g1gcch                              = NaN( nsupport, nunits, nunits );
g2tidx                              = NaN( nsupport, nunits, nunits );
g2gcch                              = NaN( nsupport, nunits, nunits );

% go over pairs abd compute ASGs
for n1                              = 1 : nunits
    for n2                          = 1 : nunits
        if n1 == n2
            continue
        end
        n12                         = s2s.shankclu( [ n1 n2 ], : );
        [ g1, g2, ~, ~, s ]         = calc_asg( s2s, n12, 'jitWindowMS', jitWindowMS );
        % assign the ASG
        g1mat( n1, n2 )             = g1;
        g2mat( n1, n2 )             = g2;
        % determine the support of the ASGs
        idx1t                       = ismember( sidx, s.g1base );
        idx2t                       = ismember( sidx, s.g2base );
        [ ~, idx1s ]                = intersect( s.g1base, sidx );
        [ ~, idx2s ]                = intersect( s.g2base, sidx );
        % assign the STC for ASGe and ASGi
        g1tidx( idx1t, n1, n2 )     = s.g1base( idx1s );
        g1gcch( idx1t, n1, n2 )     = s.gcch( s.g1base( idx1s ) );
        g2tidx( idx2t, n1, n2 )     = s.g2base( idx2s );
        g2gcch( idx2t, n1, n2 )     = s.gcch( s.g2base( idx2s ) );
        
    end
end

% update the structure
s2s.g1mat                           = g1mat;
s2s.g2mat                           = g2mat;
s2s.g1tidx                          = g1tidx;
s2s.g1gcch                          = g1gcch;
s2s.g2tidx                          = g2tidx;
s2s.g2gcch                          = g2gcch;

return

% EOF

% use example (update the s2s on the disk):
s2sfname            = [ filebase '.s2s' ];
load( s2sfname, '-mat' )
s2s                 = spikes2spikes_supplement( s2s );
save( s2sfname, 's2s' )

x = [ 16:19];
for i=1:length(x)   
    
    filebase = sprintf('//media//shirly//C22865A128659567//mice//EPS//s2s//mV99_%d', x(i));
    s2sfname            = [ filebase '.s2s' ];
    load( s2sfname, '-mat' )
    s2s                 = spikes2spikes_supplement( s2s );
    save( s2sfname, 's2s' )

end
