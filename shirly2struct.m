% WRITE HELP


function s = shirly2struct( filebases )

nsess                       = length( filebases );
for i                       = 1 : nsess
    filebase                = filebases{ i };
    try
        mono                = check_mono( filebase );
    catch
        fprintf( 1, 'Cannot load %s\n', filebase )
        continue
    end
    s1                      = mono2struct( mono );
    if i == 1
        s                   = s1;
    else
        s                   = struct_cat( s, s1 );
    end
end

return

% take the excitatory pairs:
function s1 = mono2struct( mono )
    [~, filename] = fileparts(mono.filebase);
for i                       = 1 : 4
    switch i
        case 1
            tmp             = mono.pairsExc;
            whatI           = 1;
        case 2
            tmp             = mono.pairsInh;
            whatI           = -1;
        case 3
            tmp             = mono.pairsSync;
            whatI           = 0;
        case 4
            tmp             = mono.pairsDesync;
            whatI           = -2;
    end
    npairs                  = size( tmp, 1 );
    filenames               = repmat( { filename }, [ npairs 1 ] );
    shankclu1               = mono.shankclu( tmp( :, 1 ), : );
    shankclu2               = mono.shankclu( tmp( :, 2 ), : );
    pval                    = tmp( :, 3 );
    timelag                 = tmp( :, 4 );
    STG                     = tmp( :, 5 );
    what                    = whatI * ones( npairs, 1 ); % 1 == exc; -1 == inh; sync = 0; desync = -2;
    sI.filebase             = filenames;
    sI.shankclu1            = shankclu1;
    sI.shankclu2            = shankclu2;
    sI.pval                 = pval;
    sI.timelag              = timelag;
    sI.STG                  = STG;
    sI.what                 = what;
    if i == 1
        s1                  = sI;
    else
        s1 = struct_cat( s1, sI );
    end
    
end

return

%% create filebases
datadir     = '/probox1/mice/';
filenames   = unique(sst.filebase);
nfiles                  = length( filenames );
mouse       = extractBefore(filenames,"_")
for i                   = 1 : nfiles
    filepath{i,1}    = [datadir mouse{i} '/dat/' filenames{i} '/' filenames{i}  ];
        
end        
