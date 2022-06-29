% LoadStims         load data from *stm* files
%
% CALL              [ vals, chans, stims ] = LoadStims( filebase )
%
% conventions:
% backwards compatibility is kept with 
%       [ Vals, Trigs ] = LoadVals( filebase )
% so vals and chans are sorted by onset time
% 
% forward compatibility is kept with stim_make, 
% so chans have negative values for simultaneous events
%
% thus vals( :, 1 : 2 ) is not a unique set 
% (neither is Vals( :, 1 : 2 ))
%
% calls             stim_check
%                   LoadVals, LoadXml, makesrslen, stim_make
%                   get_stimchans, parseNchannels

% 06-feb-13 ES

% revisions
% 25-feb-13 extended to support call with par structure (must have 'FileName' field)
% 13-nov-13 if no *stm* files, try to load *val* files
% 21-jul-15 added try-catch for *stm* loading (in case not confirms to format)
% 17-aug-19 cleaned up

function [ vals, chans, stims ] = LoadStims( filebase, writeOption )

vals                = [];
chans               = [];

ok                  = 0;
if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    ok              = 1;
elseif isa( filebase, 'struct' )
    if isfield( filebase, 'FileName' )
        filebase    = filebase.FileName;
        ok          = 1;
    end
end
if ~ok
    fprintf( '%s: check argument filebase\n', upper( mfilename ) );
    return
end

pathname                                = fileparts( filebase );
fnames                                  = dir( [ filebase '*stm*' ] );
if isempty( fnames )
    [ vals, chans ]                     = LoadVals( filebase );
    par                                 = LoadXml( filebase );
    srslen                              = makesrslen( filebase );
    % make de-novo...
    uchans                              = unique( chans );
    for i = 1 : length( uchans )
        stim                            = stim_make;
        stim.filebase                   = filebase;
        stim.chan                       = uchans( i );
        idx                             = chans == uchans( i );
        stim.times                      = vals( idx, 1 : 2 );
        stim.slopes                     = vals( idx, 3 : 4 );
        stim.durs                       = ( diff( vals( idx, 1 : 2 ), 1, 2 ) + 1 ) / par.SampleRate;
        stim.vals                       = vals( idx, 5 );
        [ ~, ~, vrange, srs ]           = get_stimchans( filebase, uchans( i ) );
        stim.voltagerange               = vrange;
        if isa( srs, 'cell' )
            stim.source                 = srs{ 1 };
        end
        stim.suffix                     = 'eeg';
        stim.duration                   = srslen;
        stim.generator                  = { computer, datestr( now, 'ddmmmyy' ) };
        % assume:
        stim.index                      = true( sum( idx ), 1 ); % assume unique
        stim.types                      = repmat( { 'PULSE' }, [ sum( idx ) 1 ] );
        stim.category                   = ones( sum( idx ), 3 );
        stim.franges                    = zeros( sum( idx ), 2 );
        stim.plateaus                   = zeros( sum( idx ), 1 );
        stim.median                     = 0;
        stim.stats                      = [ stim.vals zeros( sum( idx ), 1 ) ];
        stims( i )                      = stim;
    end
    if nargin <= 2 || writeOption == 0
        return
    end
    try
        fprintf( '\n%s\n%s: No stim structures for %s; creating w/ default arguments\n%s\n'...
            , repmat( '*', [ 1 160 ] ) , upper( mfilename ), filebase, repmat( '*', [ 1 160 ] ) )
        parseNchannels( filebase );%, parseOW, parseSuffix, parseMflag, parseGraphics, minAmp, minDur, minDC, sdGauss ); % OK
        fnames                          = dir( [ filebase '*stm*' ] );
    catch
        return
    end
end

vals                        = [];
times                       = [];
slopes                      = [];
v                           = [];
chans                       = [];
mismatched                  = 0;
for i                       = 1 : length( fnames )
    fname                   = [ pathname '/' fnames( i ).name ];
    try
        L                   = load( fname, '-mat' );
    catch
        continue
    end
    if ~isfield( L, 'stim' )
        continue
    end
    stim                    = L.stim; 
    if ~stim_check( stim )
        mismatched          = 1;
        continue
    end
    stims( i )              = stim;
    
    % now fill the other fields
    n                       = size( stim.category, 1 );
    chan                    = stim.chan;
    if length( chan ) == 1
        chans               = [ chans; ones( n, 1 ) * chan ];
    elseif size( stim.category, 2 ) == 4
        chans               = [ chans; -stim.category( :, 4 ) ];
    else
        chans               = [ chans; -100 * ones( n, 1 ) ];
    end
    times                   = [ times; stim.times ];
    slopes                  = [ slopes; stim.slopes ];
    v                       = [ v; stim.vals ];
end
vals                        = [ times slopes v ];

% sort
if ~isempty( vals )
    [ vals, sidx ]          = sortrows( vals, 1 );
    chans                   = chans( sidx );
end

% report
if mismatched
    fprintf( '%s\n%s: Mismatching stim structures; consider running\n\t>> parseAllChannels( ''%s'' );\n%s\n'...
        , repmat( '*', [ 1 160 ] ) , upper( mfilename ), filebase, repmat( '*', [ 1 160 ] ) )
end

return

% EOF

sum( diff( Vals( Trigs > 0, 1 : 2 ), [], 2 ) + 1  ) / 20000