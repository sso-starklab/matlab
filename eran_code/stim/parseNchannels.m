% parseNchannels                parse stimulus events from continous data (multiple channels)
%
% call                          stims = parseNchannels( filebase, channels )
%
% gets                  filebase
%                       channels            which channels to use
%
% CALLS: 
%                       get_stimchans, ParseArgPairs, LoadXml
%                       parse1channel, stim_make
%                       parseSimEvents
%                       stim_plot
%                       stim_categorize

% 01-aug-19 ES based on parseMultipleChannels and parseAllChannels and stim_make_files

% revisions
% 04-aug-19 added support for transformation matrices in formatted file (see details in interp_mat_example.m)
% 05-aug-19 added stim_categorize 
% 17-aug-19 list of dependencies added

function stims = parseNchannels( filebase, channels, varargin )

vflag                       = 1;
V2I                         = 0.01;             % CS conversion: from V command [V] to applied I [A]

%------------------------------------------------------------------%
% initialization
%------------------------------------------------------------------%
stim                        = stim_make;

% arguments
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( channels )
    channels                = get_stimchans( filebase );
end
nchans                      = length( channels );
[ suffix, graphics, Overwrite...
    , minAmpRelative, minDurationSEC, minDutyCycle, sdGaussSEC...
    , savetype, tbaseline...
    , minCC, sineCC ...
    ] = ParseArgPairs(...
    { 'suffix', 'graphics', 'Overwrite'...
    , 'minAmpRelative', 'minDurationSEC', 'minDutyCycle', 'sdGaussSEC'...
    , 'savetype', 'tbaseline'...
    , 'minCC', 'sineCC' ...
    }...
    , { 'eeg', [ 1 0 ], -2 ...
    , 0.001, [], 0.5, [] ...
    , 'png', 0 ...
    , 0.45, [] ...
    }...
    , varargin{ : } );
if length( graphics ) < 2
    graphics                = [ graphics( : ).' zeros( 1, 2 - length( graphics ) ) ];
end

% *.prm.xml file
prmfile                     = [ filebase '.prm.xml' ];
if ~exist( prmfile, 'file' )
    fprintf( 1, '%s: Missing parameter file %s\n', mfname, prmfile )
    return
end
par                         = LoadXml( prmfile );

% internal parameters
switch suffix
    case { 'eeg', 'lfp' }
        if isempty( minDurationSEC )
            minDurationSEC  = 0.002;
        end
        if isempty( sdGaussSEC )
            sdGaussSEC      = 0;
        end
    case 'dat'
        if isempty( minDurationSEC )
            minDurationSEC  = 0.0005;
        end
        if isempty( sdGaussSEC )
            sdGaussSEC      = 0.00025;
        end
end
spkFs                       = par.SampleRate;

% directory for figures
[ pathname, filename, extname ] = fileparts( filebase );
homedir                     = fileparts( fileparts( pathname ) );
figdir                      = sprintf( '%s/figs', homedir );
if ~exist( figdir, 'dir' )
    mkdir( homedir, 'figs' )
end

% load transformation table
ichans                      = [];
tabfile                     = sprintf( '%s/prm/%s.tab', homedir, filename );
if exist( tabfile, 'file' )
    L                       = load( tabfile, '-mat' );
else
    fprintf( 1, '%s: No table found for %s in %s...\n', mfname, filename, tabfile )
end
if exist( 'L', 'var' ) && all( isfield( L, { 'tab', 'chans' } ) )
    % first, check the data in the file
    [ m, n ]                 = size( L.tab );
    if m == 0 || n == 0 || ( n - 1 ) ~= length( L.chans )
        fprintf( 1, '%s: Table in %s not constructed properly - ignored!!\n', mfname, tabfile )
    else
        % second, transform the units:
        % V command [V] -> I out [A]    divide by 100
        % P [microW] -> P [mW]          divide by 1000
        L.tab( :, 1 )           = L.tab( :, 1 ) * V2I;
        L.tab( :, 2 : n )       = L.tab( :, 2 : n ) * 0.001;
        % third, add a row of zeros at the beginning
        if sum( L.tab( 1, : ) ) ~= 0
            L.tab               = [ zeros( 1, n ); L.tab ];
            m                   = m + 1;
        end
        % fourth, populate individual matrices
        [ ichans, ix ]          = intersect( L.chans, channels );
        imats                   = zeros( m, 2, length( ix ) );
        for i                   = 1 : length( ix )
            imats( :, :, i )    = [ L.tab( :, 1 ) L.tab( :, ix( i ) + 1 ) ];
        end
    end
end


%------------------------------------------------------------------%
% parse individual channels
%------------------------------------------------------------------%
fprintf( 1, '%s: Parsing individual channels...', mfname )
stims               = stim_make( nchans, 1 );
for i               = 1 : nchans
    chan            = channels( i );
    if length( channels ) == length( minAmpRelative )
        ma          = minAmpRelative( i );
    else
        ma          = minAmpRelative;
    end
    if isreal( minCC )
        minCCparse  = minCC;
    else
        minCCparse  = imag( minCC );
    end
    if isempty( ichans )
        ichan       = [];
    else
        ichan       = ismember( ichans, chan );
    end
    if isempty( ichan ) 
        imat        = [];
    else
        imat        = imats( :, :, ichan );
    end
    stims( i ) = parse1channel( filebase, chan, 'suffix', suffix...
        , 'Overwrite', Overwrite...
        , 'minAmpRelative', ma, 'minDurationSEC', minDurationSEC...
        , 'minDutyCycle', minDutyCycle, 'sdGaussSEC', sdGaussSEC...
        , 'tbaseline', tbaseline , 'imat', imat ...
        , 'wnCC', minCCparse, 'sineCC', sineCC );
end

%------------------------------------------------------------------%
% determine whether to recompute sim structure
%------------------------------------------------------------------%
savenameSim                 = sprintf( '%s.stm.sim', filebase );
figname                     = sprintf( '%s/%s.stm.sim', figdir, [ filename extname ] );
if Overwrite >= 0 || ~exist( savenameSim, 'file' )
    [ stim, stims ]         = parseSimEvents( stims, channels, spkFs );
else
    fprintf( 1, '%s: loading %s...\n', mfname, savenameSim )
    load( savenameSim, '-mat' )
    if stim_check( stim )
        [ nums, types ]     = uhist( stim.types );
        str = '';
        for i               = 1 : length( nums )
            str             = sprintf( '%s%d %s; ', str, nums( i ), types{ i } );
        end
        if i > 0
            verb( sprintf( '%ss.', str( 1 : end - 2 ) ), vflag )
        else
            verb( 'No stimuli detected.', vflag )
        end
    else
        fprintf( 1, 'format mismatch!\n' )
        Overwrite = abs( Overwrite );
        clear stim
    end
end

%------------------------------------------------------------------%
% save to disk, append merged structure, and summarize graphically
%------------------------------------------------------------------%
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( savenameSim, 'file' ) )
    fprintf( 1, '%s: saving %s...\n', mfname, savenameSim ) 
    save( savenameSim, 'stim', '-v6' )
    simstim                 = stim;
    for i                   = 1 : nchans
        stim                = stims( i );
        savename1           = sprintf( '%s.stm.%s', filebase, num3str( channels( i ) ) );
        fprintf( 1, '%s: saving %s...\n', mfname, savename1 )
        save( savename1, 'stim', '-v6' )
    end
    stim                    = simstim;
    graphics( 1 )           = abs( graphics( 1 ) );
end

% append sim and plot
stims( end + 1 )            = stim;
if graphics( 1 ) > 0
    stim_plot( stims, figname, savetype, abs( Overwrite ) );
end

%------------------------------------------------------------------%
% categorize stimuli
%------------------------------------------------------------------%
stims                       = stim_categorize( stims, 0, 0, Overwrite, filebase );

return

% EOF
