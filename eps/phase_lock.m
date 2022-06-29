%% analysis for theta and ripple locking
function phase_lock (mode, state, varargin)
%% ripple
% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ colors_NPB ]            = ParseArgPairs( ...
    {  } ...
    , {  }...
    , varargin{ : } );

nargs                   = nargin;
if nargs < 1
    mode = 'ripple';
end
if nargs < 2
    state = 'spontaneous';
end
    
colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR

% ripple lock
 % for one session mA234_17
phsS =load ('/media/shirly/C22865A128659567/mice/mA234/mat/hfo/mA234_17.hfo.spiking.spontaneous.details', '-mat');
idxS = ismember(sst.filebase,'mA234_17');
sstS = struct_select(sst,idxS);
% sstD = load ('/media/shirly/C22865A128659567/mice/EPS/sst/depth/sst/mA234_17.sst', '-mat')
% sstD = sstD.sst
% gidx                = check_cluster_quality( sstD, ilevel );
% sstD = struct_select (sstD, gidx);
phsS.stats.depth         = phsS.stats.depth* 20;        % [um]

isbipS           = ~isnan(sstS.bpi) & ~isinf( sstS.bpi );
isposS           = sstS.extremum > 0 & ~isbipS;
idxI = phsS.stats.pval<0.05;
figure,
scatter (phsS.stats.mPhase(idxI&isbipS),sstS.depth(idxI&isbipS))
hold on,
scatter (phsS.stats.mPhase(idxI&isposS),sstS.depth(idxI&isposS))
scatter (phsS.stats.mPhase(idxI&~isbipS&~isposS),phsS.stats.depth(idxI&~isbipS&~isposS))

% for the full dataset
sst.mPhase = NaN(length(sst.shankclu),1);
sst.rPhase = NaN(length(sst.shankclu),1);
sst.pPhase = NaN(length(sst.shankclu),1);

% x = [17:18, 20:33, 35:36, 38:42 ];
% str = 'mA234';
% str_ = 'mA234_';

% x = [20:21 23:26 30:45 47:48 52:58 60 ];
% str = 'mC41';
% str_ = 'mC41_';

x = [10:11 13:15 18:24 27:29 31 ];
str = 'mDS1';
str_ = 'mDS1_';

for i=1:length(x)
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//mat//hfo//%s%d.hfo.spiking.spontaneous.details', str, str_, x(i));
    phsS =load (filebase, '-mat');
    phsS = phsS.stats;
    str_f = sprintf ('%s%d', str_, x(i));
    idxS = ismember(sst.filebase,str_f);
    idxSf = find(idxS);
    sst.mPhase (idxSf) = phsS.mPhase;
    sst.rPhase (idxSf) = phsS.rPhase;
    sst.pPhase (idxSf) = phsS.pval;
    sst.depth (idxSf)  = phsS.depth;
end
%-----------------------------------%

idxI = sst.pPhase<0.05;
% depth         = sst.depth* 20;        % [um] - fix for mC41!!!!!!!!
   % compute depth in um based on probe resolution
    fnames              = sst.filebase;
    ufnames             = unique( sst.filebase );
    nufnames            = length( ufnames );
    tmpcell             = strfind( ufnames, 'mC41' ); % presently supporting only one animal
    % to support multiple animals (strings), must expand the code
    tmpvec              = NaN( nufnames, 1 );
    for i               = 1 : nufnames
        if isempty( tmpcell{ i } )
            tmpvec( i ) = 0;
        else
            tmpvec( i ) = tmpcell{ i };
        end
    end
    tmpvec              = logical( tmpvec );
    hridx               = ismember( fnames, ufnames( tmpvec ) );
    depth( hridx )          = sst.depth( hridx ) * 15;         % [um]
    depth( ~hridx )         = sst.depth( ~hridx ) * 20;        % [um]
    


figure,
ph( 1 ) = scatter (sst.mPhase(idxI&isbip),depth(idxI&isbip), '.' ) % Biphasic
hold on,
ph( 2 ) = scatter (sst.mPhase(idxI&ispos),depth(idxI&ispos), '.' ) % Puntis

ph( 3 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&sst.pyr),depth(idxI&~isbip&~ispos&sst.pyr), '.' ) % Pyr

ph( 4 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&~sst.pyr),depth(idxI&~isbip&~ispos&~sst.pyr), '.' ) % Int

% Divide the depth to bins of 10um
range = [-150 150];
bins_size = 10;
bin_vec = range(1):bins_size:range(2);
depth_binned = NaN(length(depth),1);
for i = 1:(length(bin_vec)-1)
    idx1 = depth<bin_vec(i+1) & depth>=bin_vec(i);
    idx2 = find(idx1);
    depth_binned(idx2) = (bin_vec(i+1)+bin_vec(i))/2;
end

figure;
nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));
% for each uniqe bin find all the phases
binned_struct = [];
for i = 1:length(uniq_depth)
    idx3 = depth_binned==uniq_depth(i);
    binned_struct.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
    binned_struct.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
    binned_struct.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
    binned_struct.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));
    
    binned_struct.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
    binned_struct.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   
    
    binned_struct.depth(i) = uniq_depth(i);
    
    binned_struct.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
    sumbpiS(i) = sum(idx3&isbip);
    sumbpiSL(i) = sum(idx3&isbip&idxI);
   
    sumPYRS(i) = sum(idx3&~isbip&~ispos&sst.pyr);
    sumPYRSL(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);

    sumPVS(i) = sum(idx3&~isbip&~ispos&~sst.pyr);
    sumPVSL(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);   
end


figure;
ph( 1 )=subplot(2,2,1);
errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*m', 'horizontal')
ph( 2 )=subplot(2,2,2);
errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*g', 'horizontal')
ph( 3 )=subplot(2,2,3);
errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*b', 'horizontal')
ph( 4 )=subplot(2,2,4);
errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*r', 'horizontal')
for i =1:4
    subplot(2,2,i)
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );


figure;
errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
hold on,
errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
%  errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'+', 'horizontal', 'Color', colors_NPB( 2, : ))

for i = 1:length(binned_struct.depth)
    if binned_struct.pval_PYRbpi(i)<0.05
        errorbar(6.3,binned_struct.depth(i),[],'*k')
        text(6.35,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRbpi(i)))   
    end
    text(6.7,binned_struct.depth(i),sprintf('%d,%d,%d',sumbpiSL(i),sumPYRSL(i),sumPVSL(i)))
end
for i =1:4
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end


set( ph( 1 ), 'CData', colors_NPB( 3, : ) )    
set( ph( 2 ), 'CData', colors_NPB( 2, : ) )    
set( ph( 3 ), 'CData', colors_IP_light( 2, : ) )    
set( ph( 4 ), 'CData', colors_IP_light( 1, : ) )    

set( ph( 1 ), 'Colormap', colors_NPB( 3, : ) )    
set( ph( 2 ), 'Colormap', colors_NPB( 2, : ) )    
set( ph( 3 ), 'Colormap', colors_IP_light( 2, : ) )    
set( ph( 4 ), 'Colormap', colors_IP_light( 1, : ) ) 
set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
%-----------------------------------------------%
% induced ripples


% for the full dataset
sst.ImPhas = NaN(length(sst.shankclu),1);
sst.IrPhas = NaN(length(sst.shankclu),1);
sst.IpPhas = NaN(length(sst.shankclu),1);

% x = [17:18, 20:33, 35:36, 38:42 ];
% str = 'mA234';
% str_ = 'mA234_';

x = [20:21 23:26 30:45 47:48 52:58 60 ];
str = 'mC41';
str_ = 'mC41_';

% x = [10:15 18:24 27:31 ];
% str = 'mDS1';
% str_ = 'mDS1_';

for i=1:length(x)
    try
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//mat//hfo//induced//%s%d.hfo_spiking_induced.details', str, str_, x(i));
    phsS =load (filebase, '-mat');
    phsS = phsS.stats;
    str_f = sprintf ('%s%d', str_, x(i));
    idxS = ismember(sst.filebase,str_f);
    idxSf = find(idxS);
    sst.ImPhas (idxSf) = phsS.mPhase;
    sst.IrPhas (idxSf) = phsS.rPhase;
    sst.IpPhas (idxSf) = phsS.pval;
    catch
    end
end

idxI = sst.IpPhas<0.05;

% Divide the depth to bins of 10um
range = [-150 150];
bins_size = 10;
bin_vec = range(1):bins_size:range(2);
depth_binned = NaN(length(depth),1);
for i = 1:(length(bin_vec)-1)
    idx1 = depth<bin_vec(i+1) & depth>=bin_vec(i);
    idx2 = find(idx1);
    depth_binned(idx2) = (bin_vec(i+1)+bin_vec(i))/2;
end

figure;
nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));
% for each uniqe bin find all the phases
binned_struct = [];
for i = 1:length(uniq_depth)
    idx3 = depth_binned==uniq_depth(i);
    binned_struct.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
    binned_struct.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
    binned_struct.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
    binned_struct.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));
    
    binned_struct.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
    binned_struct.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   
    
    binned_struct.depth(i) = uniq_depth(i);
    
    sumbpiI(i) = sum(idx3&isbip);
    sumbpiIL(i) = sum(idx3&isbip&idxI);
   
    sumPYRI(i) = sum(idx3&~isbip&~ispos&sst.pyr);
    sumPYRIL(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);

    sumPVI(i) = sum(idx3&~isbip&~ispos&~sst.pyr);
    sumPVIL(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);
end


figure;
ph( 1 )=subplot(2,2,1);
errorbar(binned_struct.mIPhasePyr, binned_struct.depth,binned_struct.SEMIPhasePyr,'*m', 'horizontal')
ph( 2 )=subplot(2,2,2);
errorbar(binned_struct.mIPhaseInt, binned_struct.depth,binned_struct.SEMIPhaseInt,'*g', 'horizontal')
ph( 3 )=subplot(2,2,3);
errorbar(binned_struct.mIPhaseBip, binned_struct.depth,binned_struct.SEMIPhaseBip,'*b', 'horizontal')
ph( 4 )=subplot(2,2,4);
errorbar(binned_struct.mIPhasePos, binned_struct.depth,binned_struct.SEMIPhasePos,'*r', 'horizontal')
for i =1:4
    subplot(2,2,i)
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );


figure;
errorbar(binned_struct.mIPhasePyr, binned_struct.depth,binned_struct.SEMIPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
hold on,
errorbar(binned_struct.mIPhaseInt, binned_struct.depth,binned_struct.SEMIPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
errorbar(binned_struct.mIPhaseBip, binned_struct.depth,binned_struct.SEMIPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
% errorbar(binned_struct.mIPhasePos, binned_struct.depth,binned_struct.SEMIPhasePos,'*', 'horizontal')
for i =1:4
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end


%plot sum/sum
figure;
bar([1,2,3,4,5,6], [sum(sumbpiSL)/sum(sumbpiS) sum(sumbpiIL)/sum(sumbpiI)...
    sum(sumPYRSL)/sum(sumPYRS) sum(sumPYRIL)/sum(sumPYRI) sum(sumPVSL)/sum(sumPVS) sum(sumPVIL)/sum(sumPVI)])
% subplot(1,2,2)
% bar([1,2,3], [sum(sumbpiIL)/sum(sumbpiI) sum(sumPYRIL)/sum(sumPYRI) sum(sumPVIL)/sum(sumPVI)])
%---------------------------------------------------%
%%
% calculating all units (not only significant)
% for the full dataset
sst.mPhase = NaN(length(sst.shankclu),1);
sst.rPhase = NaN(length(sst.shankclu),1);
sst.pPhase = NaN(length(sst.shankclu),1);

x = [17:18, 20:33, 35:36, 38:42 ];
str = 'mA234';
str_ = 'mA234_';


for i=1:length(x)
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//mat//hfo//%s%d.hfo.spiking.spontaneous.details', str, str_, x(i));
    phsS =load (filebase, '-mat');
    phsS = phsS.stats;
    str_f = sprintf ('%s%d', str_, x(i));
    idxS = ismember(sst.filebase,str_f);
    idxSf = find(idxS);
    sst.mPhase (idxSf) = phsS.mPhase;
    sst.rPhase (idxSf) = phsS.rPhase;
    sst.pPhase (idxSf) = phsS.pval;
    sst.depth (idxSf)  = phsS.depth;
end

idxI = sst.pPhase>0;
depth         = sst.depth* 20;        % [um]

figure,
ph( 1 ) = scatter (sst.mPhase(idxI&isbip),depth(idxI&isbip), '.' ) % Biphasic
hold on,
ph( 2 ) = scatter (sst.mPhase(idxI&ispos),depth(idxI&ispos), '.' ) % Puntis

ph( 3 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&sst.pyr),depth(idxI&~isbip&~ispos&sst.pyr), '.' ) % Pyr

ph( 4 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&~sst.pyr),depth(idxI&~isbip&~ispos&~sst.pyr), '.' ) % Int

% Divide the depth to bins of 10um
range = [-150 150];
bins_size = 10;
bin_vec = range(1):bins_size:range(2);
depth_binned = NaN(length(depth),1);
for i = 1:(length(bin_vec)-1)
    idx1 = depth<bin_vec(i+1) & depth>=bin_vec(i);
    idx2 = find(idx1);
    depth_binned(idx2) = (bin_vec(i+1)+bin_vec(i))/2;
end

figure;
nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));
% for each uniqe bin find all the phases
binned_struct = [];
for i = 1:length(uniq_depth)
    idx3 = depth_binned==uniq_depth(i);
    binned_struct.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
    binned_struct.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
    binned_struct.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
    binned_struct.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));
    
    binned_struct.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
    binned_struct.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   
    
    binned_struct.depth(i) = uniq_depth(i);
end


figure;
ph( 1 )=subplot(2,2,1);
errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*m', 'horizontal')
ph( 2 )=subplot(2,2,2);
errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*g', 'horizontal')
ph( 3 )=subplot(2,2,3);
errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*b', 'horizontal')
ph( 4 )=subplot(2,2,4);
errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*r', 'horizontal')
for i =1:4
    subplot(2,2,i)
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );


figure;
errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
hold on,
errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
% errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*', 'horizontal')
for i =1:4
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

%-----------------------------------------------%
% induced ripples


% for the full dataset
sst.ImPhas = NaN(length(sst.shankclu),1);
sst.IrPhas = NaN(length(sst.shankclu),1);
sst.IpPhas = NaN(length(sst.shankclu),1);

% x = [17:18, 20:33, 35:36, 38:42 ];
% str = 'mA234';
% str_ = 'mA234_';
% 
% x = [20:21 23:26 30:45 47:48 52:58 60 ];
% str = 'mC41';
% str_ = 'mC41_';

x = [10:11 13:15 18:24 27:29 31 ];
str = 'mDS1';
str_ = 'mDS1_';

for i=1:length(x)
    try
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//EPS//hfo_ind//%s%d.hfo_spiking_induced.details', str_, x(i));
    phsS =load (filebase, '-mat');
    phsS = phsS.stats;
    str_f = sprintf ('%s%d', str_, x(i));
    idxS = ismember(sst.filebase,str_f);
    idxSf = find(idxS);
    sst.ImPhas (idxSf) = phsS.mPhase;
    sst.IrPhas (idxSf) = phsS.rPhase;
    sst.IpPhas (idxSf) = phsS.pval;
    catch
    end
end

idxI = sst.IpPhas>0;

% Divide the depth to bins of 10um
range = [-150 150];
bins_size = 10;
bin_vec = range(1):bins_size:range(2);
depth_binned = NaN(length(depth),1);
for i = 1:(length(bin_vec)-1)
    idx1 = depth<bin_vec(i+1) & depth>=bin_vec(i);
    idx2 = find(idx1);
    depth_binned(idx2) = (bin_vec(i+1)+bin_vec(i))/2;
end

figure;
nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));
% for each uniqe bin find all the phases
binned_struct = [];
for i = 1:length(uniq_depth)
    idx3 = depth_binned==uniq_depth(i);
    binned_struct.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
    binned_struct.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
    binned_struct.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
    binned_struct.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));
    
    binned_struct.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
    binned_struct.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   
    
    binned_struct.depth(i) = uniq_depth(i);
    
    sumbpi(i) = sum(idx3&isbip);
end


figure;
ph( 1 )=subplot(2,2,1);
errorbar(binned_struct.mIPhasePyr, binned_struct.depth,binned_struct.SEMIPhasePyr,'*m', 'horizontal')
ph( 2 )=subplot(2,2,2);
errorbar(binned_struct.mIPhaseInt, binned_struct.depth,binned_struct.SEMIPhaseInt,'*g', 'horizontal')
ph( 3 )=subplot(2,2,3);
errorbar(binned_struct.mIPhaseBip, binned_struct.depth,binned_struct.SEMIPhaseBip,'*b', 'horizontal')
ph( 4 )=subplot(2,2,4);
errorbar(binned_struct.mIPhasePos, binned_struct.depth,binned_struct.SEMIPhasePos,'*r', 'horizontal')
for i =1:4
    subplot(2,2,i)
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );


figure;
errorbar(binned_struct.mIPhasePyr, binned_struct.depth,binned_struct.SEMIPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
hold on,
errorbar(binned_struct.mIPhaseInt, binned_struct.depth,binned_struct.SEMIPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
errorbar(binned_struct.mIPhaseBip, binned_struct.depth,binned_struct.SEMIPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
% errorbar(binned_struct.mIPhasePos, binned_struct.depth,binned_struct.SEMIPhasePos,'*', 'horizontal')
for i =1:4
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

%% calculate theta by waveform
BP                             = [1 60];
dtFactor                       = 14;
fir                            = [BP, dtFactor];

x = [18, 20:33, 35:36, 38:42 ];
str = 'mA234';
str_ = 'mA234_';
for i=1:length(x)
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
    temp                           = load( [ filebase '.phs' ], 'eegchan', '-mat' );                     % the 'best theta' channel
    eegchan                        = temp.eegchan;
    phs                            = eegPhase( filebase, eegchan, [], fir, [], [], [], [], 'extrema');
    save( [ filebase '.phs_twf_spi' ], 'phs', 'eegchan' )
end
% ------------------------------------------%
% once the phs is ran over all the sessions:

% for the full dataset
% sst.mThphs = NaN(length(sst.shankclu),1);
% sst.rThphs = NaN(length(sst.shankclu),1);
% sst.pThphs = NaN(length(sst.shankclu),1);

x = [17:18, 21:33, 35:36, 38:42 ];
str = 'mA234';
str_ = 'mA234_';

for i=1:length(x)
    filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
    % load the phases:
    load ([ filebase '.phs_twf_spi' ], '-mat');

    % get theta times
    THETAper                       = get_states( filebase, 'the' ); % Periods of theta, phsFs

    phsFs = 1250; % [Hz]
    % get histogram time per each side
    [~,~, phsEdges]             = circ_hist( [0, 2*pi], 100 );
    phsdt                       = ( phsEdges(2) - phsEdges(1) ) / 2;
    phsCenters                  = [ phsEdges(1:end-1) + phsdt ]';
    timeHist                    = histc( phs(enumerate(THETAper)), phsEdges ); % theta-phase distribution within theta periods
    timeHist                    = timeHist(1:end-1);
    timeHist                    = timeHist / phsFs; %[s]
    timeHist                    = timeHist(:)';

    % calc mean phase
    spkFs = 20000;
    [s, ~]                      = load_spikes( filebase, [], 'B' );
    nunits                      = size( s.shankclu, 1 );
    mThphs = NaN(nunits,1);
    rThphs = NaN(nunits,1);
    pThphs = NaN(nunits,1);
    for j = 1 : nunits
        % get unit data
        aclunum                 = s.map( j, 1 );
        idx                     = s.clu == aclunum;
        spk                     = s.res( idx );         % [spkFs]

        spk_phsFs               = ceil( spk / spkFs * phsFs );     % spk_phsFs - spike times in phsFs
        tspk_phsFs              = ismember(spk_phsFs,THETAper);
        spk_phs                 = phs(tspk_phsFs);

        phshistCOUNT            = histc(spk_phs,phsEdges);
        phshistCOUNT            = phshistCOUNT(1:end-1);
        phshistCOUNT            = phshistCOUNT(:)';
        phshistRate             = phshistCOUNT ./ timeHist;
        [phi,R,SD]              = circ_mean(phsCenters(:), phshistRate(:));
        lock_sig                = ray_test(phsCenters(:), phshistRate(:));
        mThphs(j)               = phi;
        rThphs(j)               = R;
        pThphs(j)               = lock_sig;
    end
    % populate the data in full sst fields:
    str_f = sprintf ('%s%d', str_, x(i));
    idxS = ismember(sst.filebase,str_f);
    idxSf = find(idxS);
    sst.mThphs (idxSf) = mThphs;
    sst.rThphs (idxSf) = rThphs;
    sst.pThphs (idxSf) = pThphs;
end

idxI = sst.pThphs<0.05;
depth         = sst.depth* 20;        % [um]

% Divide the depth to bins of 10um
range = [-150 150];
bins_size = 10;
bin_vec = range(1):bins_size:range(2);
depth_binned = NaN(length(depth),1);
for i = 1:(length(bin_vec)-1)
    idx1 = depth<bin_vec(i+1) & depth>=bin_vec(i);
    idx2 = find(idx1);
    depth_binned(idx2) = (bin_vec(i+1)+bin_vec(i))/2;
end

figure;
nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));
% for each uniqe bin find all the phases
binned_struct = [];
for i = 1:length(uniq_depth)
    idx3 = depth_binned==uniq_depth(i);
    binned_struct.mPhaseBip(i) = nanmean(sst.mThphs(idx3&idxI&isbip));
    binned_struct.SEMPhaseBip(i) = calc_sem(sst.mThphs(idx3&idxI&isbip));
    binned_struct.mPhasePos(i) = nanmean(sst.mThphs(idx3&idxI&ispos));
    binned_struct.SEMPhasePos(i) = calc_sem(sst.mThphs(idx3&idxI&ispos));
    
    binned_struct.mPhasePyr(i) = nanmean(sst.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.SEMPhasePyr(i) = calc_sem(sst.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
    binned_struct.mPhaseInt(i) = nanmean(sst.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));
    binned_struct.SEMPhaseInt(i) = calc_sem(sst.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));   
    
    binned_struct.depth(i) = uniq_depth(i);
end

figure;
errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
hold on,
errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
% errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*', 'horizontal')
for i =1:4
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end