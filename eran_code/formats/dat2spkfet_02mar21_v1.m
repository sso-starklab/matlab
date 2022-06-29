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
%                               Overwrite       1       compute and overwrite
%                                               {-2}    compute and write only if does not exist
%                               filebaseTarget  {''}    defaults to filebase
%                                                       if distinct, will save *spk.#, *fet.#, and *alg.# to the target directory
%                               filebaseDat     {''}    defaults to filebase
%                                                       if distinct, will load *xml and *dat from the filebaseDat directory
%                               shanknums       {[]}    work on specific shanks
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
% flow control                  is implemented using an *.alg.# file, which is
%                               a *mat file with a 4-element integer flag
%                               each element represents the status of realignment for a given file type: 
%                                   [ res clu spk fet ]
%                               initial status is 0
%                               realignres updates flag( 1 ) -> 1
%                               reunique updates flag( 1 : 2 ) -> 2
%                               dat2spkfet updates flag( 3 : 4 ) -> 2
%
% calls                         ParseArgPairs, LoadXml, dat2spk, spk2fet
%
% see also                      ks2ndm, realignres, resunique

% 05-aug-20 ES+AL

% revisions
% 06-aug-20 (1) typo correction
%           (2) help expanded
% 16-feb-21 (1) added *.alg.# file with support for overwriting
%           (2) added filebaseTarget
% 02-mar-21 (1) added support for a subset of shanks 
%           (2) fixed bug in ashank/sidx

function rc = dat2spkfet( filebase, varargin )

%------------------------------------------------------------
% constants
%------------------------------------------------------------
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
    , shanknums ...
    , useres ...
    , filebaseTarget, filebaseDat, Overwrite ] = ParseArgPairs(...
    { 'nrank' ...
    , 'nsamples', 'peaksample' ...
    , 'shanknums' ...
    , 'useres' ...
    , 'filebaseTarget', 'filebaseDat', 'Overwrite' } ...
    , {  NRANK ...
    , NSAMPLES, PEAKSAMPLE ...
    , [] ...
    ,  USERES ...
    , '', '', -2 } ...
    , varargin{ : } );
if isempty( filebaseTarget )
    filebaseTarget              = filebase;
end
if isempty( filebaseDat )
    filebaseDat                 = filebase;
end

%------------------------------------------------------------
% collect information
%------------------------------------------------------------
mfname                          = upper( mfilename );
fprintf( 1, '%s: working on %s. Checking for *dat file...', mfname, filebaseDat )
datfname                        = [ filebaseDat '.dat' ];
if ~exist( datfname, 'file' )
    fprintf( 1, '%s: missing file %s, exiting\n', mfname, datfname )
    return
else
    fprintf( 1, 'OK\n' )
end

fprintf( 1, '%s: Loading *xml file...', mfname, filebaseDat )
par                             = LoadXml( filebaseDat );
fprintf( 1, 'done.\n' )

nallchans                       = par.nChannels;
nshanks                         = par.nElecGps;
allshanks                       = 1 : nshanks;
if isempty( shanknums )
    shanknums                   = allshanks;
else
    shanknums                   = intersect( allshanks, shanknums );
end
[ ~, i1 ]                       = intersect( allshanks, shanknums );
ushanks                         = allshanks( i1 );
nshanks                         = length( i1 );

%------------------------------------------------------------
% go over shanks
%------------------------------------------------------------
rc                              = NaN * ones( nshanks, 2 );
for sidx                        = 1 : nshanks
    
    % determine the shank-specific spikes
    ashank                      = ushanks( sidx );
    algfname                    = sprintf( '%s.alg.%d', filebaseTarget, ashank );
    resfname                    = sprintf( '%s.res.%d', filebase, ashank );
    spkfname                    = sprintf( '%s.spk.%d', filebaseTarget, ashank );
    fetfname                    = sprintf( '%s.fet.%d', filebaseTarget, ashank );
    chans                       = sort( par.SpkGrps( ashank ).Channels ) + 1;
    nchans                      = length( chans );
    
    % check if spk and fet have already been "uniquified"
    if ~exist( algfname, 'file' )
        fprintf( 1, '%s: Shank#%d res has not been even realigned yet - run realignres.m!!\n', mfname, ashank )
        continue
    end
    if Overwrite == -2
        load( algfname, 'flag', '-mat' )
        if ~isequal( size( flag ), [ 1 4 ] ) || ~isa( flag, 'numeric' )
            error( 'format mismatch for %s', algfname )
        end
        if flag( 3 ) == 2 && flag( 4 ) == 2
            fprintf( 1, '%s: Shank#%d spk + fet have already been realigned and uniquified, skipping\n', mfname, ashank )
            continue
        end
    end
    
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
    
    %------------------------------------------------------------
    % update alg file
    if rc( sidx, 1 ) == 0 && rc( sidx, 2 ) == 0
        load( algfname, 'flag', '-mat' )
        if ~isequal( size( flag ), [ 1 4 ] ) || ~isa( flag, 'numeric' )
            error( 'format mismatch for %s', algfname )
        end
        flag( 3 )               = 2;
        flag( 4 )               = 2;
        save( algfname, 'flag', '-v6' )
    end
    
end

return

% EOF