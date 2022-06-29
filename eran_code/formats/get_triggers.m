% get_triggers          from stim structures
%
% [ trigs tims durs vals stims ] = get_triggers( filebase, channels, uflag, multi, simOnly, params )
% 
% ARGUMENTS:
% filebase          full base + base or par structure
% channels          channel numbers
%
% the two modifiers apply to single channel and simultaneous channel selection:
% uflag         single-channel event selection:
%                   {1} forces single-channel events to be unique
%                   0   allows other channels to be active at the same time
%                   -1  forces events to occur sim. on multiple channels
% multi         multi-channel event selection:
%                   -'le' indicates any subset of channels (e.g. if [ 1 2 3 ]
%                   are specified, 1,2,12,13,23, and 123 events will all be
%                   included)
%                   -{'eq'} requires an exact match (out of the available
%                   channels). thus if 123 are specified, events 1234 will
%                   be ignored
%                   -'ge' allows additional channels to be active at the
%                   same time (e.g. if [ 1 2 3 ] are specified and [ 4 ] is
%                   coactive, the 1234 events will not be rejected but
%                   rather labelled as 123 events)
% simOnly          {0}: don't care
%                    1: use only simultaneous structures
%                   -1: use only non-simultaneous structures
%
% the rest of the arguments are parameter pairs refering to fields of stim,
% listed as 'fieldname1', fieldvals1, 'fieldname2', fieldvals2 (see
% stim_select for details)
%
% OUTPUT
% trigs         integers; a unique number for each trigger type:
%                   the category field of the stim structure (see stim_categorize.m)
%                   is converted from a vector to a scalar by multiplying
%                       [ channel type_category val_category dur_category ]
%                   by [ 1e9 1e6 1e3 1 ]. for simultaneous structures, the
%                   4th column of the category (channel combination category) replaces the channel
% tims          [samples]; trigger onset times
% durs          [sec]
% vals          [V]
% stims         the structures (filtered by channels, parameters, and logic);
%                   (can be used to bypass the stim_categorize results)
%
% EXAMPLE:
% to get the triggers which were pulses 40-80 ms long, <80 mA, on channels 33, 35, or both, type:
%   >> channels = [ 33 35 ];
%   >> params = { 'types', 'PULSE', 'durs',[ 0.04 0.08 ], 'vals', [ 0 0.08 ] };
%   >> [ trigs tims durs vals stims ] = get_triggers( filebase, channels, 1, 1, params )
% to get all 0:40 ascending chirps:
%   >> [ trigs tims durs vals stims ] = get_triggers( filebase, 36, 1, 'eq', 'types', 'ZAP', 'franges', [ 0 40 ] ); 
%
% CALLS                 LoadStims, LoadXml, get_stimchans                   (formats)
%                       stim_get, stim_get_channels, stim_make, stim_select (stim)
%                       uhist, verb                                         (general)
%
% See also              parseNchannels, get_spikes

% 25-feb-13 ES

% revisions
% 25-mar-13 modified uflag to enable -1, non-simultaneous only
% 18-aug-19 cleaned up

function [ trigs, tims, durs, vals, stims ] = get_triggers( filebase, channels, uflag, multi, simOnly, varargin )

%------------------------------------------------------------------
% constants
verbose                     = 1;

% initialize
trigs                       = [];                                       % labels
tims                        = [];                                       % [samples]
durs                        = [];                                       % [s]
vals                        = [];                                       % [V]
stims                       = [];

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    fprintf( '%s: input mismatch\n', upper( mfilename ) );
    return
end
if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    par                 = LoadXml( [ filebase '.prm.xml' ] );
elseif isa( filebase, 'struct' )
    if isfield( filebase, 'FileName' )
        par                 = filebase;
        filebase            = filebase.FileName;
    end
end
if ~exist( 'par', 'var' )
    fprintf( '%s: missing filebase\n', upper( mfilename ) );
    return
end
schannels                   = get_stimchans( par );
if nargs < 2 || isempty( channels )
    channels                = schannels;
else
    channels                = intersect( channels, schannels );
end
if isempty( channels )
    fprintf( '%s: missing channels %s\n', upper( mfilename ), num2str( channels( : ).' ) );
    return
end
if nargs < 3 || isempty( uflag )                                        % which events to choose (from single-channel structures)
    uflag                   = 1;
    %  0: don't care
    %  1: use only unique (non-simultaneous) events
    % -1: use only non-unique (simultaneous) events
end
if nargs < 4 || isempty( multi )
    multi                   = 'eq';
end
if nargs < 5 || isempty( simOnly )                                      % which structures to use (single/multi/any -channel )
    simOnly                 = 0;
    %  0: don't care
    %  1: use only simultaneous structures
    % -1: use only non-simultaneous structures
end

%------------------------------------------------------------------
% get the stims structures
[ ~, filename ]             = fileparts( filebase );
verb( sprintf( '%s: Loading %s trigger data... ', upper( mfilename ), filename ), -verbose )
[ ~, ~, stims ]             = LoadStims( filebase );

% choose the relevant channels
nstims                      = length( stims );
j                           = 0;
n0                          = zeros( 1, nstims );
for i = 1 : nstims
    n0( i )                 = length( stims( i ).types );
    out                     = stim_get_channels( stims( i ), channels, multi );   
    if ~isempty( out )
        if simOnly == 0 ...
                || simOnly == 1 && size( out.category, 2 ) == 4 ...
                || simOnly == -1 && size( out.category, 2 ) == 3
            j               = j + 1;
            s( j )          = out;
        end
    end
end
if ~exist( 's', 'var' )
    nstims                  = 0;
else
    nstims                  = length( s );
end
if nstims == 0
    fprintf( '\n%s: No valid stimuli for %s\n', upper( mfilename ), filebase )
    return
end

%------------------------------------------------------------------
% extract the relevant stimuli
stims                       = stim_make;
rmvi                        = [];
for i                       = 1 : nstims                                % get the proper stimuli
    si                      = stim_select( s( i ), varargin{ : } );
    if size( si.category, 2 ) == 3                                      % single-channel structure
        if uflag == 1                                                   % unique events
            si              = stim_get( si, si.index == 1 );
        elseif uflag == -1                                              % non-unique (simultaneous)
            si              = stim_get( si, si.index == 0 );
        end
    end
    % convert to a trigs/tims/durs convention
    n                       = length( si.types );
    if n == 0
        rmvi                = [ rmvi; i ]; 
        continue; 
    end
    if size( si.category, 2 ) == 3
        cmat                = [ si.chan * ones( n, 1 ) si.category( :, 1 : 3 ) ];
    else
        cmat                = si.category( :, [ 4 1 : 3 ] );       
    end
    itrigs                  = round( sum( cmat .* ( ones( n, 1 ) * [ 1e9 1e6 1e3 1 ] ), 2 ) );
    % [ channel/combination type val dur/freq ]
    trigs                   = [ trigs; itrigs ];
    tims                    = [ tims; si.times( :, 1 ) ];
    durs                    = [ durs; ( diff( si.times, [], 2 ) + 1 ) / par.SampleRate ]; % temporary fix
    vals                    = [ vals; si.vals ];
    stims( i )              = si;
end

rmvi( rmvi > length( stims ) ) = [];
if ~isempty( rmvi )
    stims( rmvi )           = [];
end

%------------------------------------------------------------------
% summarize
ntrigs                      = uhist( trigs );
verb( sprintf( '%d/%d triggers (%d-%d/cat %d cats) from %d/%d chan-combs were loaded; '...
    , length( trigs ), sum( n0 ), min( ntrigs ), max( ntrigs ), length( ntrigs )...
    , length( stims ), length( n0 ) ), -verbose )
verb( sprintf( 'duration, %6.3f (SD=%6.3f) ms'...
    , 1000 * mean( durs ), 1000 * std( durs ) ), verbose )

return

% EOF
