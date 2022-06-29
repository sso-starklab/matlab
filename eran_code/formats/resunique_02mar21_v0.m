% resunique      	removes same-unit spikes which have similar time samples 
%
% call            	rc              = resunique( filebase )
%
% gets              filebase        path to data
%
% optional arguments        
%                   DeadTime            {6}     [samples]
%                   Overwrite           1       compute and overwrite
%                                       {-2}    compute and write only if does not exist
%                   filebaseTarget      {''}    defaults to filebase
%                                               if distinct, will save *clu.#, *res.# and *alg.# to the target directory 
%                   shanknums           {[]}    work on specific shanks
%
% does              (1) loads *.clu.# and *.res.# files from filebase
%                   (2) dilutes spike occurrences for each unit
%                   (3) saves *.clu.# and *.res.# files to filebaseTarget
%
% requires          files named filebase.xml filebase.clu.shankNumber and filebase.res.shankNumber
%
% calls             LoadXml, ParseArgPairs
%
% usage             After the proccess of spike sorting manual curation, run the realignres function.
%                   Then run this function
%
% note              the output of this routine CANNOT be opened in klusters
%                   without updating the spk and fet files
%
% flow control   	is implemented using an *.alg.# file, which is
%                 	a *mat file with a 4-element integer flag
%                  	each element represents the status of realignment for a given file type: 
%                   	[ res clu spk fet ]
%                	initial status is 0
%                	realignres updates flag( 1 ) -> 1
%                   reunique updates flag( 1 : 2 ) -> 2
%                   dat2spkfet updates flag( 3 : 4 ) -> 2
%
% see also          realignres, dat2spkfet

% 25-aug-20 LS

% revisions
% 04-nov-20 added dead time of 7 samples for each unit (similar to NDM)
% 16-feb-21 (1) modified deleted spikes to be non-first spikes in "burst"
%           (2) modified DeadTime to be an argument
%           (3) added varagin
%           (4) added *.alg.# file with support for overwriting
%           (5) added filebaseTarget
% 02-mar-21 (1) added support for a subset of shanks 

function rc = resunique( filebase, varargin )

DeadTime_DFLT                   = 7;                                        % [samples]

%------------------------------------------------------------
% arguments
%------------------------------------------------------------
nargs                           = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ DeadTime ...
    , filebaseTarget ...
    , shanknums ...
    , Overwrite ]        = ParseArgPairs(...
    { 'DeadTime' ...
    , 'filebaseTarget' ...
    , 'shanknums' ...
    , 'Overwrite' } ...
    , {  DeadTime_DFLT ...
    , '' ...
    , [] ...
    ,  -2 } ...
    , varargin{ : } );
if ~ismember( Overwrite, [ 1 -2 ] )
    error( 'unsupported Overwrite mode %d', Overwrite )
end
if isempty( filebaseTarget )
    filebaseTarget              = filebase;
end

%------------------------------------------------------------
% preparations
%------------------------------------------------------------
mfname                          = upper( mfilename );
fprintf( 1, '%s: Loading *xml file...', mfname, filebase )
par                             = LoadXml( sprintf( '%s.xml', filebase ) );
fprintf( 1, 'done.\n' )

nshanks                         = par.nElecGps;
%ushanks                         = 1 : nshanks;
allshanks                       = 1 : nshanks;
if isempty( shanknums )
    shanknums                   = allshanks;
else
    shanknums                   = intersect( allshanks, shanknums );
end
[ ~, i1 ]                       = intersect( allshanks, shanknums );
ushanks                         = allshanks( i1 );
nshanks                         = length( i1 );

rc                              = NaN * ones( nshanks, 2 );
for sidx                        = 1 : nshanks
    
    %------------------------------------------------------------
    % prepare
    
    % determine the shank-specific spikes
    shanknum                    = ushanks( sidx );
    algfname                    = sprintf( '%s.alg.%d', filebaseTarget, shanknum );
    clufnameSource              = sprintf( '%s.clu.%d', filebase, shanknum );
    clufnameTarget              = sprintf( '%s.clu.%d', filebaseTarget, shanknum );
    resfnameSource              = sprintf( '%s.res.%d', filebase, shanknum );
    resfnameTarget              = sprintf( '%s.res.%d', filebaseTarget, shanknum );
    
    % check if res and clu have already been "uniquified"
    if ~exist( algfname, 'file' )
        fprintf( 1, '%s: Shank#%d res has not been realigned yet - run realignres.m!!\n', mfname, shanknum )
        continue
    end
    if Overwrite == -2
        load( algfname, 'flag', '-mat' )
        if ~isequal( size( flag ), [ 1 4 ] ) || ~isa( flag, 'numeric' )
            error( 'format mismatch for %s', algfname )
        end
        if flag( 1 ) == 0
            fprintf( 1, '%s: Shank#%d res has not been realigned yet - run realignres.m!!\n', mfname, shanknum )
            continue
        end
        if flag( 1 ) == 2 && flag( 2 ) == 2
            fprintf( 1, '%s: Shank#%d res + clu have already been realigned and uniquified, skipping\n', mfname, shanknum )
            continue
        end
    end
    
    %------------------------------------------------------------
    % load res file
    fprintf( 1, '\n %s: Loading shank#%d res... ', mfname, shanknum )
    res                         = load( resfnameSource );
    nspks                       = length( res );
    fprintf( 1, '%d spikes. ', nspks )
    
    %------------------------------------------------------------
    % load clu file
    fprintf( 1, 'clu... ' )
    clu                         = load( clufnameSource );
    nclu0                       = clu( 1 );
    clu( 1 )                    = [];
    uclu                        = unique( clu );
    nclu                        = length( uclu );
    if nclu0 ~= nclu
        warning( 'weird clu file' )
    end
    fprintf( 1, '%d clu. ', nclu )
    
    %------------------------------------------------------------
    % dilute
    Vecunits                    = unique( clu );
    resU                        = nan( size( res ) );
    cluU                        = nan( size( clu ) );
    for i                       = 1 : length( Vecunits )
        U                       = Vecunits( i );
        idx                     = clu == U;
        clu1                    = clu( idx );
        res1                    = res( idx );
        res0                    = res1;
        clu0                    = clu1;
        Uidx                    = find( diff( res0 ) <= DeadTime ) + 1;  	% keep only the first spike in an intra-unit "burst" 
        res0( Uidx )            = nan;
        clu0( Uidx )            = nan;
        resU( idx )             = res0;
        cluU( idx )             = clu0;
    end
    idxNan                      = isnan( cluU );
    cluU( idxNan )              = [];
    resU( idxNan )              = [];
    
    %------------------------------------------------------------
    % save res file
    % structure: ASCII file, one row per spike time, in samples
    
    res                         = int32( resU );           % this supports 29 hours of recordings
    fid                         = fopen( resfnameTarget, 'w' );
    [ ~ ]                       = fprintf( fid, '%d\n', res );
    rc( sidx, 1 )               = fclose( fid );
    if rc( sidx, 1 ) == 0
        fprintf( 1, '%d spikes. clu...', length( res ) )
    else
        fprintf( 1, 'failed writing file: %s!!\n', resfnameTarget  )
    end
    
    %------------------------------------------------------------
    % save clu file
    % use shank-specific clu numbers starting from 2 (1 is noise)
    % structure: ASCII file. first row, number of clusters. then, one row per spike.
    
    ncluV                       = unique( cluU );
    ncluW                       = length( ncluV );
    
    fid                         = fopen( clufnameTarget, 'w' );
    [ ~ ]                       = fprintf( fid, '%d\n', ncluW );
    [ ~ ]                       = fprintf( fid, '%d\n', cluU );
    rc( sidx, 2 )               = fclose( fid );
    if rc( sidx, 2 ) == 0
        fprintf( 1, '%d clusters.', ncluW )
    else
        fprintf( 1, 'failed writing file: %s!!\n', clufnameTarget  )
    end
    
    %------------------------------------------------------------
    % update alg file
    if rc( sidx, 1 ) == 0 && rc( sidx, 2 ) == 0
        load( algfname, 'flag', '-mat' )
        if ~isequal( size( flag ), [ 1 4 ] ) || ~isa( flag, 'numeric' )
            error( 'format mismatch for %s', algfname )
        end
        flag( 1 )               = 2;
        flag( 2 )               = 2;
        save( algfname, 'flag', '-v6' )
    end
    
end

fprintf( 1, '\n' )

return

% EOF
