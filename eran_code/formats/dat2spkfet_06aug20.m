% dat2spkfet                    extract spk and fet files from dat and res
%
% call                          rc = dat2spkfet( filebase )
%
% gets                          filebase
%
% optional arguments (name/value pairs)
%
%                               nrank           {3}     number of PC coefficients
%                               nsamples        {32}    number of samples to include in every spike
%                               peaksample      {16}    number of samples before and including the peak
%                               useres          {1}     add a column of res (spike times in samples) to each entry in the *fet*
%
% returns                       rc              return code - matrix of ngroups by 2 columns (spk, fet success)          
% 
% does                          goes over electrode groups according to *xml file:
%                               -extracts *spk* file from *dat according to *res* file
%                               -extracts *fet* file from *spk* by PCA of the data 
%
% requirements                  filebase.xml
%                               filebase.dat
%                               filebase.res.#
%
% calls                         ParseArgPairs, LoadXml, dat2spk, spk2fet
%
% see also                      ks2ndm

% 05-aug-20 ES+AL

% revisions
% 06-aug-20 (1) typo correction
%           (2) help expanded

function rc = dat2spkfet( filebase, varargin )


%------------------------------------------------------------
% constants
%------------------------------------------------------------
% constants
fetBits                         = 16;       % number of bits in fet file (scaling factor, it is an ASCII file anyhow)
NSAMPLES                        = 32;       % number of samples to include in every spike
PEAKSAMPLE                      = 16;       % number of samples before and including the peak

% default values
NRANK                           = 3;            % relevant for 'spk' only
USERES                          = 1;

%------------------------------------------------------------
% arguments
%------------------------------------------------------------
nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ nrank ...
    , nsamples, peaksample ...
    , useres ] = ParseArgPairs(...
    { 'nrank' ...
    , 'nsamples', 'peaksample' ...
    , 'useres' } ...
    , {  NRANK ...
    , NSAMPLES, PEAKSAMPLE ...
    ,  USERES } ...
    , varargin{ : } );

%------------------------------------------------------------
% collect information
%------------------------------------------------------------
mfname                          = upper( mfilename );
fprintf( 1, '%s: working on %s. Checking for *dat file...', mfname, filebase )
datfname                        = [ filebase '.dat' ];
if ~exist( datfname, 'file' )
    fprintf( 1, '%s: missing file %s, exiting\n', mfname, datfname )
    return
else
    fprintf( 1, 'OK\n' )
end

fprintf( 1, '%s: Loading *xml file...', mfname, filebase )
par                             = LoadXml( filebase );
fprintf( 1, 'done.\n' )

nallchans                       = par.nChannels;
nshanks                         = par.nElecGps;

%------------------------------------------------------------
% go over shanks
%------------------------------------------------------------
rc                              = NaN * ones( nshanks, 2 );
ushanks                         = 1 : nshanks;
for sidx                        = 1 : nshanks
    
    % determine the shank-specific spikes
    ashank                      = ushanks( sidx );
    fetfname                    = sprintf( '%s.fet.%d', filebase, ashank );
    resfname                    = sprintf( '%s.res.%d', filebase, ashank );
    spkfname                    = sprintf( '%s.spk.%d', filebase, ashank );
    chans                       = sort( par.SpkGrps( sidx ).Channels ) + 1;
    nchans                      = length( chans );
    
    %------------------------------------------------------------
    % res file
    fprintf( 1, '%s: Loading shank#%d res... ', mfname, ashank )
    res                         = load( resfname );
    nspks                       = length( res );
    fprintf( 1, '%d spikes. ', nspks )
    
    %------------------------------------------------------------
    % spk file
    % structure: binary file, nsamples around each spike
    fprintf( 1, 'spk... ' )
    rc( sidx, 1 )               = dat2spk( datfname, spkfname, double( res ), nallchans, chans ...
        , 'nsamples', nsamples, 'peaksample', peaksample );
    if rc( sidx, 1 ) == 0
        fprintf( 1, 'written spk file. fet...' )
    else
        fprintf( 1, 'failed writing file: %s!!', spkfname  )
    end
    
    %------------------------------------------------------------
    % fet file
    % structure: ASCII file. first, row with number of features
    % then, 3 PCA features per channel, then some additional features (nmissing),
    % then spike times
    % all are integers, positive or negative
    ngchans                     = nchans;                                   % see ks2ndm for other options
    nmissing                    = ngchans;
    nfets                       = ngchans * nrank + nmissing + useres;

    rc( sidx, 2 )               = spk2fet( spkfname, resfname, fetfname, ngchans ...
        , 'fetBits', fetBits, 'verbose', 1, 'npcs', nrank );

    if rc( sidx, 2 ) == 0
        fprintf( 1, '%d features. done!\n', nfets )
    else
        fprintf( 1, 'failed writing file: %s!!\n', fetfname  )
    end
    
end

return

% EOF