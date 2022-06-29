

OverwriteRealign            = -2;
realignResMethod            = 'extremum';

x = [ 16];
str = 'mP23';
str_ = 'mP23_';

for i=1:length(x)
    try
        eval(sprintf('!rsync -avh --progress //odin2//Recordings//%s//dat//%s%d//%s%d.clu* /media/shirly/C22865A128659567/mice/%s/dat/%s%d/', str, str_, x(i), str_, x(i), str, str_, x(i)));
        eval(sprintf('!rsync -avh --progress //odin2//Recordings//%s//dat//%s%d//%s%d.res* /media/shirly/C22865A128659567/mice/%s/dat/%s%d/', str, str_, x(i), str_, x(i), str, str_, x(i)));
        eval(sprintf('!rsync -avh --progress //odin2//Recordings//%s//dat//%s%d//%s%d.spk* /media/shirly/C22865A128659567/mice/%s/dat/%s%d/', str, str_, x(i), str_, x(i), str, str_, x(i)));
%         eval(sprintf('!rsync -avh --progress //odin2//Recordings//%s//dat//%s%d//%s%d.xml /media/shirly/C22865A128659567/mice/%s/dat/%s%d/', str, str_, x(i), str_, x(i), str, str_,x(i)));
        filebaseTarget  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
%         filebaseDat     = sprintf('//odin2//Recordings//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i)); % only for *dat file
        filename        = sprintf('%s%d', str_, x(i));
        setup           = sprintf('%s', str);

        par                     = LoadXml( [ filebaseTarget '.xml' ] );
        rc1                     = realignres( filebaseTarget, 'method', realignResMethod, 'filebaseTarget', filebaseTarget, 'Overwrite', OverwriteRealign );
        rc2                     = resunique( filebaseTarget, 'filebaseTarget', filebaseTarget, 'Overwrite', OverwriteRealign  );
        nsamples                = par.SpkGrps( 1 ).nSamples;                                    % assume same for all spike groups
        peaksample              = par.SpkGrps( 1 ).PeakSample;                                  % assume same for all spike groups
        rc3                     = dat2spkfet( filebaseTarget, 'nsamples', nsamples, 'peaksample', peaksample, 'filebaseDat', filebaseTarget, 'Overwrite', OverwriteRealign );
    catch
        fprintf('session %s%d alignment failed\n', str_, x(i));
    end
end

toFlip                      = [ 1 0 ];
flipLFP                     = toFlip( 1 );
flipSPK                     = toFlip( 2 );
OverwriteSST                = 1;
byParSST                    = 0;
spkNotDetrended             = 0;                    % typically, spikes are detrended during extraction, and the cell type classifiers assume this. If no detrending was done, this flag should be 1
savetype                    = 'png';

x = [11];
str = 'mK01';
str_ = 'mK01_';

for i=1:length(x)
    try
        filebaseTarget  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
toFlip                      = [ 1 0 ];
flipLFP                     = toFlip( 1 );
flipSPK                     = toFlip( 2 );
OverwriteSST                = 1;
byParSST                    = 0;
spkNotDetrended             = 0;              [ ~, sst ]              = run_spikes_stats( filebaseTarget, 'Overwrite', OverwriteSST, 'byPar', byParSST );
%        [ ~, sst ]              = run_spikes_stats( filebaseTarget, 'Overwrite', OverwriteSST, 'byPar', byParSST, 's2s_only', 1 );
        spikes_stats_plot( sst, filebaseTarget, 'savetype', savetype, 'graphics', 1 ...
        , 'flipLFP', flipLFP, 'flipSPK', flipSPK, 'hpfUsed', spkNotDetrended );
     catch
        fprintf('session %s%d spike_stats failed\n', str_, x(i));
    end
end


datadir = '/media/shirly/C22865A128659567/mice/EPS/sst/new_realigned';
[ res, sst, ff ] = shirly_eps_analysis( datadir, 'onlygather', 1);
ispos=sst.extremum>0;
idx = find(ispos);

for pos = 241:250%1:length(idx)
    Pfilename = sst.filebase {idx(pos)};
    Pfilebase = filebase_lookup(Pfilename);
    Pshankclu = sst.shankclu(idx(pos),:);
    figure(pos);
    plot_ss(Pfilebase,Pshankclu);
end

Pfilename = sst.filebase {idx(162)};

Pfilebase = filebase_lookup(Pfilename);
Pshankclu = sst.shankclu(idx(162),:);
shankclu1  = Pshankclu;
shankclu2 = [2 14];
filebase = Pfilebase



ignoredChannels = [];
flipLFP = 1;

x = [ 7 9 ];
str = 'mV99';
str_ = 'mV99_0';

for i=1:length(x)
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
    plotProbeHFOs( filebase, flipLFP, ignoredChannels ); % flipLFP is used here to flip the probe
    sst                 = spikes_stats_depth( filebase, 'graphics', 1,  'Overwrite', 1, 'flipLFP', 1 );
    sstfname            = [ filebase '.sst' ];
    save( sstfname, 'sst' )
end

x = [10:19];
str = 'mV99';
str_ = 'mV99_';
OverwriteHFOspkCount = -2;
OverwriteHFOspk = -2;
OverwriteHFOphs = -2;
savetype = 'png';
ilevel = 'B';

for i=1:length(x)
    try
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
    hfoAnalysisSpikingCount( filebase , 'Overwrite', OverwriteHFOspkCount, 'graphics', 1 ...
        , 'ilevel', ilevel, 'savetype', savetype );
        hfoAnalysisSpiking( filebase, 'spontaneous' ...
        , 'Overwrite', OverwriteHFOspk, 'OverwritePHS', OverwriteHFOphs ...
        , 'graphics', [ 0 1 0 0 ], 'ilevel', ilevel, 'savetype', savetype );
        hfoAnalysisSpiking_channelinduced( filebase, 'induced' ...
        , 'Overwrite', OverwriteHFOspk, 'OverwritePHS', OverwriteHFOphs ...
        , 'graphics', [ 0 1 0 0 ], 'ilevel', ilevel, 'savetype', savetype );
    catch
    end
end