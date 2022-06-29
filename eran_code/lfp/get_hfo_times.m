% get_hfo_times         from existing *spw* files
%
% [ periods trigs hchannels ] = get_hfo_times( filebase, channels );
%
% filebase          full path
% channels          list of channels. 
%                       -defaults to the peak channels listed in *sps
%                           or if doesn't exit to all *spw* files
%                       -otherwise uses the intersection between the
%                           existing *spw* files and the requested channels
%
% periods           ranges
% trigs             peak power, aligned to nearest cycle trough
% hchannels         source of 
%
% see also          find_hfos

% 19-may-13 ES

% revisions
% 17-aug-19 cleaned up

function [ periods, trigs, hchannels ] = get_hfo_times( filebase, channels )

% initialize output
periods                 = [];
trigs                   = [];
hchannels               = [];

% arguments
nargs                   = nargin;
if nargs < 2 || isempty( channels )
    channels            = [];
end
pathname                = fileparts( filebase );

% get the relevant spw file names
spsfname                = [ filebase '.sps' ];
if isempty( channels ) && exist( spsfname, 'file' )
    L                   = load( spsfname, '-mat' );
    channels            = L.vote;
end
spwfnames               = dir( [ filebase '.spw*' ] );
lchans                  = zeros( 1, length( spwfnames ) );
rmv                     = [];
for i                   = 1 : length( spwfnames )
    if isempty( str2num( spwfnames( i ).name( end - 2 : end ) ) )
        rmv             = [ rmv; i ];
    end
end
if ~isempty( rmv )
    spwfnames( rmv )    = [];
    lchans( rmv )       = [];
end
for i                   = 1 : length( spwfnames )
    [ ~, ~, extname ]   = fileparts( spwfnames( i ).name );
    lchans( i )         = str2num( extname( 2 : end ) );
end
if ~isempty( channels )
    [ lchans, i1 ]      = intersect( lchans, channels );
    spwfnames           = spwfnames( i1 );
end
if isempty( lchans )
    return
end

% load and get the parameters
for i                   = 1 : length( spwfnames )
   	spwfname            = [ pathname '/' spwfnames( i ).name ];
    L                   = load( spwfname, '-mat' );
    if length( L.rips.trigs ) ~= size( L.rips.edges )
        fprintf( '%s: Mismatch at %s; recomputing...', upper( mfilename ), spwfname )
        diffOrd             = 0;
        behaviorDurSEC      = 2;
        padBuffer           = [ -0.01 0.01 ];
        filtMode            = 'dog';
        powerMode           = 'LP';
        clipBase            = 1;
        L.rips = find_hfos( filebase, chan...
            , 'Overwrite', 1, 'ripBP', L.rips.bp, 'ripTH', L.rips.TH...
            , 'mindurSECrip',  L.rips.mindur, 'minisiSECrip', L.rips.minisi...
            , 'diffOrd', diffOrd, 'behaviorDurSEC', behaviorDurSEC...
            , 'padBuffer', padBuffer, 'filtMode', filtMode, 'powerMode', powerMode, 'clipBase', clipBase );
        fprintf( 'done!\n' )
    end
    trigs                   = [ trigs; L.rips.trigs ]; 
    periods                 = [ periods; L.rips.edges ]; 
    hchannels               = [ hchannels; ones( length( L.rips.trigs ), 1 ) * L.rips.chans( 1 ) ];
end

% sort
[ trigs, sidx ]             = sort( trigs );
periods                     = periods( sidx, : ); % i.e. not neccesarily sorted by onset times
hchannels                   = hchannels( sidx, : );

return

% EOF
