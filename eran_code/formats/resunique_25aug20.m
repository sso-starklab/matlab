% resunique          cleans the spikes which has the same time sample form each unit 
%
% call               rc              = resunique( filebase)
%
% gets              filebase        path to data
%
%
% 
% requires          files named filebase.xml filebase.clu.shankNumber and filebase.res.shankNumber
%
% calls             LoadXml
%
% usage             After the proccess of spike sorting manual curation, run the realignres function.
%                   Then run this function
%
%
% note              1. the output of this routine can be opened in klusters
%
%
%
% does
%   Rewrites  the *.clu.# and *.res.# files  
%
% 25-aug-20 LS




function    rc =    resunique(filebase)

    mfname                          = upper( mfilename ); 
    fprintf( 1, '%s: Loading *xml file...', mfname, filebase )
    par                             = LoadXml( sprintf( '%s.xml', filebase ) );
    fprintf( 1, 'done.\n' )

    nshanks                         = par.nElecGps;
    rc                              = NaN * ones( nshanks, 2 );
    ushanks                         = 1 : nshanks;
    
for sidx                        = 1 : nshanks
    
    % determine the shank-specific spikes
    shanknum                    = ushanks( sidx );
    clufname                    = sprintf( '%s.clu.%d', filebase, shanknum );
    resfname                    = sprintf( '%s.res.%d', filebase, shanknum );

    
    %------------------------------------------------------------
    % load res file
    fprintf( 1, '\n %s: Loading shank#%d res... ', mfname, shanknum )
    res                         = load( resfname );
    nspks                       = length( res );
    fprintf( 1, '%d spikes. ', nspks )


    %------------------------------------------------------------
    % load clu file
    fprintf( 1, 'clu... ', mfname, shanknum )
    clu                         = load( clufname );
    nclu0                       = clu( 1 );
    clu( 1 )                    = [];
    uclu                        = unique( clu );
    nclu                        = length( uclu );
    if nclu0 ~= nclu
        warning( 'weird clu file' )
    end
    fprintf( 1, '%d clu. ', nclu )

    Vecunits            = unique(clu);
    resU                = nan(size(res));
    cluU                = nan(size(clu));
    for i =1:length(Vecunits);
        U               = Vecunits(i);
        idx             = clu==U;
        clu1            = clu(idx);
        res0            = nan(sum(idx),1);
        clu0            = nan(sum(idx),1);
        [res1,Uidx]     = unique(res(idx));
        res0(Uidx)      = res1;
        clu0(Uidx)      = clu1(Uidx); 
        resU(idx)       = res0;
        cluU(idx)       = clu0;
    end

    idxNan              = isnan(cluU);
    cluU(idxNan)        = [];
    resU(idxNan)        = [];
    
    %------------------------------------------------------------
    % save res file
    % structure: ASCII file, one row per spike time, in samples
    
    res                         = int32( resU );           % this supports 29 hours of recordings
    fid                         = fopen( resfname, 'w' );
    [ ~ ]                       = fprintf( fid, '%d\n', res );
    rc( sidx, 1 )               = fclose( fid );
    if rc( sidx, 1 ) == 0
        fprintf( 1, '%d spikes. clu...', length( res ) )
    else
        fprintf( 1, 'failed writing file: %s!!\n', resfname  )
    end
    
        %------------------------------------------------------------
    % clu file
    % use shank-specific clu numbers starting from 2 (1 is noise)
    % structure: ASCII file. first row, number of clusters. then, one row per spike. 

    ncluV                       = unique( cluU );
    ncluW                       = length( ncluV );

    fid                         = fopen( clufname, 'w' );
    [ ~ ]                       = fprintf( fid, '%d\n', ncluW );
    [ ~ ]                       = fprintf( fid, '%d\n', cluU );
    rc( sidx, 2 )               = fclose( fid );
    if rc( sidx, 2 ) == 0
        fprintf( 1, '%d clusters. ind...', ncluW )
    else
        fprintf( 1, 'failed writing file: %s!!\n', clufname  )
    end
end
end
    

