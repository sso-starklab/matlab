%% analysis for theta and ripple locking
% 
% gets:                 mode {'ripple'}, type of event to analyze. Can also support 'theta'
%                       state {'spontaneous'}, in the case of ripples can alsp support 'induced'
%                       lock {1}, analysis of ripple/theta locked units
%
% arguments:            gather {0} to repopulate the sst
%                       str        animal name (e.g. mC41)
%                       x          a vector of session numbers (e.g. [15:16 20])
%                       save {0}   save sst structure after run
%                       savedate   string of suffix for sst save
%                       

function phase_lock_bpi (mode, state, lock, varargin)
%% ripple
% argument handling
[ gather, x, str, ...
    saveST, savedate, ...
    graphics]            = ParseArgPairs( ...
    { 'gather', 'x', 'str',  ...
    'saveST', 'savedate'...
    'graphics'} ...
    , { 0, [], [], ...
    0, [], ...
     [0 0 1 1]}...
    , varargin{ : } );

nargs                   = nargin;
if nargs < 1 || isempty( mode )
    mode = 'ripple';
end
if nargs < 2 || isempty( state )
    state = 'spontaneous';
end
if nargs < 2 || isempty( lock )
    lock = 1;
end
if ~isempty(str)
    str_ = sprintf('%s_',str );
end
if isunix
    dir = '/media/shirly/C22865A128659567/mice/EPS/';
end
load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_29mar22.mat', '-mat');
colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR
isbip           = ~isnan(sst.bpi) & ~isinf( sst.bpi );
ispos           = sst.extremum > 0 & ~isbip;

% compute depth in um based on probe resolution
fnames              = sst.filebase;
ufnames             = unique( sst.filebase );
% nufnames            = length( ufnames );
% tmpcell             = strfind( ufnames, 'mC41' ); % presently supporting only one animal - should add mP20, nF79, mP93, mF105, mF108
% tmpvec              = NaN( nufnames, 1 );
% for i               = 1 : nufnames
%     if isempty( tmpcell{ i } )
%         tmpvec( i ) = 0;
%     else
%         tmpvec( i ) = tmpcell{ i };
%     end
% end
% tmpvec              = logical( tmpvec );
idx15 = (contains(sst.filebase, 'mC41') | contains(sst.filebase, 'mP20')...
    | contains(sst.filebase, 'nF79')| contains(sst.filebase, 'mP93')...
    | contains(sst.filebase, 'mF105')| contains(sst.filebase, 'mF108'));
% hridx               = ismember( fnames, ufnames( tmpvec ) );
% depth( hridx )          = sst.depth( hridx ) * 15;         % [um]
% depth( ~hridx )         = sst.depth( ~hridx ) * 20;        % [um]
depth( idx15 ,1)          = sst.depth( idx15 ) * 15;         % [um]
depth(~idx15,1)         = sst.depth( ~idx15 ) * 20;        % [um]

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

nanidx = ~isnan(depth_binned);
uniq_depth = unique(depth_binned(nanidx));

% ripple lock
switch mode
    case 'ripple' 
          if lock && isequal(state,'spontaneous')
               idxI = sst.pPhase < 0.05;
           else
               idxI = sst.pPhase >= 0;
          end   
          if lock && isequal(state,'induced')
                idxI = sst.IpPhas<0.05;
          else
                idxI = sst.IpPhas >= 0;
          end 
            % for each uniqe bin find all the phases
            binned_struct = [];
            for i = 1:length(uniq_depth)
                idx3 = depth_binned==uniq_depth(i);
                % mean preffered phases
                binned_struct.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
                binned_struct.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
                binned_struct.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
                binned_struct.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));

                binned_struct.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                binned_struct.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   
                
                % mean R
                binned_struct.rPhaseBip(i) = nanmean(sst.rPhase(idx3&idxI&isbip));
                binned_struct.SEMrPhaseBip(i) = calc_sem(sst.rPhase(idx3&idxI&isbip));
                binned_struct.rPhasePos(i) = nanmean(sst.rPhase(idx3&idxI&ispos));
                binned_struct.SEMrPhasePos(i) = calc_sem(sst.rPhase(idx3&idxI&ispos));

                binned_struct.rPhasePyr(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.SEMrPhasePyr(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.rPhaseInt(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                binned_struct.SEMrPhaseInt(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 
                binned_struct.depth(i) = uniq_depth(i);

                binned_struct.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.rpval_PYRbpi(i) = utest( sst.rPhase(idx3&idxI&isbip),sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

%                 binned_struct.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
%                 binned_struct.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
%                 binned_struct.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
%                 binned_struct.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));
% 
%                 binned_struct.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
%                 binned_struct.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                binned_struct.depth(i) = uniq_depth(i);
                
%                 binned_struct.pval_IPYRbpi(i) = utest( sst.ImPhas(idx3&idxI&isbip),sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));

                
                sumbpiS(i) = sum(idx3&isbip);
                sumbpiSL(i) = sum(idx3&isbip&idxI);

                sumPYRS(i) = sum(idx3&~isbip&~ispos&sst.pyr);
                sumPYRSL(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);

                sumINTS(i) = sum(idx3&~isbip&~ispos&~sst.pyr);
                sumINTSL(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);   
%                 
%                 sumbpiI(i) = sum(idx3&isbip);
%                 sumbpiIL(i) = sum(idx3&isbip&idxI);
% 
%                 sumPYRI(i) = sum(idx3&~isbip&~ispos&sst.pyr);
%                 sumPYRIL(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);
% 
%                 sumINTI(i) = sum(idx3&~isbip&~ispos&~sst.pyr);
%                 sumINTIL(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);
            end
            %-----------------------------------%          
        switch state
            case 'spontaneous'
                if gather
                 % for the full dataset
%                         sst.mPhase = NaN(length(sst.shankclu),1);
%                         sst.rPhase = NaN(length(sst.shankclu),1);
%                         sst.pPhase = NaN(length(sst.shankclu),1);
                    for i=1:length(x)
                        filebase  = sprintf('//media//shirly//C22865A128659567//mice//EPS//hfo//%s%d.hfo.spiking.spontaneous.details', str_, x(i));
                        phsS =load (filebase, '-mat');
                        phsS = phsS.stats;
                        str_f = sprintf ('%s%d', str_, x(i));
                        idxS = ismember(sst.filebase,str_f);
                        idxSf = find(idxS);
                        sst.mPhase (idxSf) = phsS.mPhase;
                        sst.rPhase (idxSf) = phsS.rPhase;
                        sst.pPhase (idxSf) = phsS.pval;
                        sst.depth (idxSf)  = phsS.depth;
%                         sst.mPhase (idxSf) = NaN;
%                         sst.rPhase (idxSf) = NaN;
%                         sst.pPhase (idxSf) = NaN;
%                         sst.depth (idxSf)  = NaN;
                    end
                end
               %-----------------------------------%
         
        
            if graphics(1)
                figure,
                ph( 1 ) = scatter (sst.mPhase(idxI&isbip),depth(idxI&isbip), '.' ); % Biphasic
                hold on,
                ph( 2 ) = scatter (sst.mPhase(idxI&ispos),depth(idxI&ispos), '.' ); % Puntis
                ph( 3 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&sst.pyr),depth(idxI&~isbip&~ispos&sst.pyr), '.' ); % Pyr
                ph( 4 ) = scatter (sst.mPhase(idxI&~isbip&~ispos&~sst.pyr),depth(idxI&~isbip&~ispos&~sst.pyr), '.' ); % Int
            end

            % for each uniqe bin find all the phases
            for i = 1:length(uniq_depth)
            
            end

            if graphics(2)
                figure;
                ph( 1 )=subplot(2,2,1);
                errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*','horizontal', 'Color', colors_IP_light( 2, : ))
                ph( 2 )=subplot(2,2,2);
                errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
                ph( 3 )=subplot(2,2,3);
                errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*', 'horizontal', colors_NPB( 3, : ))
                ph( 4 )=subplot(2,2,4);
                errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*', 'horizontal', colors_NPB( 2, : ))
                for i =1:4
                    subplot(2,2,i)
                    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
                    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
                    set( gca, 'tickdir', 'out', 'box', 'off' )
                end
                alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
            end
            if graphics(3)           
                figure;
                errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
                hold on,
                errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
                errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
                 errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'+', 'horizontal', 'Color', colors_NPB( 2, : ))
                title ('Ripple-phase lock vs depth')
                xlabel ('Phase during ripple')
                ylabel ('Depth')
                for i = 1:length(binned_struct.depth)
                    if binned_struct.pval_PYRbpi(i)<0.05
                        errorbar(6.25,binned_struct.depth(i),[],'*k')
                        text(6.3,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRbpi(i)))   
                    end
                    text(6.75,binned_struct.depth(end)+7,sprintf('Bpi,Pyr,Int',binned_struct.pval_PYRbpi(i)))                     
                    text(6.85,binned_struct.depth(i),sprintf('%d,%d,%d',sumbpiSL(i),sumPYRSL(i),sumINTSL(i)))
                end
                for i =1:4
                    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
                    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
                    set( gca, 'tickdir', 'out', 'box', 'off' )
                end
            end
            if graphics(4)           
                figure;
                errorbar(binned_struct.rPhasePyr, binned_struct.depth,binned_struct.SEMrPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
                hold on,
                errorbar(binned_struct.rPhaseInt, binned_struct.depth,binned_struct.SEMrPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
                errorbar(binned_struct.rPhaseBip, binned_struct.depth,binned_struct.SEMrPhaseBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
%                  errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'+', 'horizontal', 'Color', colors_NPB( 2, : ))
                xlim ([0 1])
                title ('Ripple-R lock vs depth')
                xlabel ('R during ripple')
                ylabel ('Depth')               
                for i = 1:length(binned_struct.depth)
                    if binned_struct.rpval_PYRbpi(i)<0.05
                        errorbar(0.95,binned_struct.depth(i),[],'*k')
                        text(0.955,binned_struct.depth(i),sprintf('%0.2g',binned_struct.rpval_PYRbpi(i)))   
                    end
                    text(1.02,binned_struct.depth(end)+7,sprintf('Bpi,Pyr,Int',binned_struct.rpval_PYRbpi(i)))                     
                    text(1.03,binned_struct.depth(i),sprintf('%d,%d,%d',sumbpiSL(i),sumPYRSL(i),sumINTSL(i)))
                end
                set( gca, 'tickdir', 'out', 'box', 'off' )
            end
    %-----------------------------------------------%
    % induced ripples
            case 'induced'
                if gather
%                 % for the full dataset
%                 sst.ImPhas = NaN(length(sst.shankclu),1);
%                 sst.IrPhas = NaN(length(sst.shankclu),1);
%                 sst.IpPhas = NaN(length(sst.shankclu),1);
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
                end

            if graphics(2)           
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
            end
            
            if graphics(3)
                figure;
                errorbar(binned_struct.mIPhasePyr, binned_struct.depth,binned_struct.SEMIPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
                hold on
                errorbar(binned_struct.mIPhaseInt, binned_struct.depth,binned_struct.SEMIPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
                errorbar(binned_struct.mIPhaseBip, binned_struct.depth,binned_struct.SEMIPhaseBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
                errorbar(binned_struct.mIPhasePos, binned_struct.depth,binned_struct.SEMIPhasePos,'.', 'horizontal', 'Color', colors_NPB( 2, : ))
                    title ('Induced ripple-phase lock vs depth')
                    xlabel ('Phase during induced ripples')
                    ylabel ('Depth')
                for i =1:4
                    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
                    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
                    set( gca, 'tickdir', 'out', 'box', 'off' )
                end
                for i = 1:length(binned_struct.depth)
                    if binned_struct.pval_IPYRbpi(i)<0.05
                        errorbar(6.3,binned_struct.depth(i),[],'*k')
                        text(6.35,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_IPYRbpi(i)))   
                    end
                    text(6.65,binned_struct.depth(end)+5,sprintf('Bpi,Pyr,Int',binned_struct.pval_IPYRbpi(i)))                     
                    text(6.7,binned_struct.depth(i),sprintf('%d,%d,%d',sumbpiIL(i),sumPYRIL(i),sumINTIL(i)))
                end
            end
            
            if graphics(4)
                %plot sum/sum
                figure;
                bar([1,2,3,4,5,6], [sum(sumbpiSL)/sum(sumbpiS) sum(sumbpiIL)/sum(sumbpiI)...
                    sum(sumPYRSL)/sum(sumPYRS) sum(sumPYRIL)/sum(sumPYRI)...
                    sum(sumINTSL)/sum(sumINTS) sum(sumINTIL)/sum(sumINTI)])
                % subplot(1,2,2)
                % bar([1,2,3], [sum(sumbpiIL)/sum(sumbpiI) sum(sumPYRIL)/sum(sumPYRI) sum(sumINTIL)/sum(sumINTI)])
            end
        end
    %---------------------------------------------------%
%% theta
    case 'theta' 
        % calculate theta by waveform
        BP                             = [1 60];
        dtFactor                       = 14;
        fir                            = [BP, dtFactor];
        minspk                         = 5;
        
        if gather
%             for i=1:length(x)
%                 filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
%                 temp                           = load( [ filebase '.phs' ], 'eegchan', '-mat' );                     % the 'best theta' channel
%                 eegchan                        = temp.eegchan;
%                 phs                            = eegPhase( filebase, eegchan, [], fir, [], [], [], [], 'extrema');
%                 save( [ filebase '.phs_twf_spi' ], 'phs', 'eegchan' )
%             end
            for i=1:length(x)
                filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
                filephs   = [ filebase '.phs_ripp' ];
                if ~exist(filephs, 'file')
                    L                          = load([filebase '.sps'],'-mat');
                    mstats                     = L.stats;
                    f                          = mstats(:,5);
                    pf                         = mstats(:,6);
                    minRf                   = 130;  
                    maxRf                   = 190;
                    pf(f < minRf | f > maxRf ) = NaN;
                    [~,nrp]                    = max( pf ); % take the shank with the highest ripple power
                    central_ch                 = mstats(nrp,2); % take it's central channel
                    % calculate theta phase
                    phs                   = eegPhase( filebase, central_ch, [], fir, [], [], [], [], 'extrema');
                    save( [ filebase '.phs_ripp' ], 'phs', 'central_ch' )
                end
            end    
        
            % ------------------------------------------%
            % once the phs is ran over all the sessions:

            % for the full dataset
            % sst.mThphs = NaN(length(sst.shankclu),1);
            % sst.rThphs = NaN(length(sst.shankclu),1);
            % sst.pThphs = NaN(length(sst.shankclu),1);

            for i=1:length(x)
                filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
                % load the phases:
                load ([ filebase '.phs_ripp' ], '-mat');

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

%                     if length(spk_phs) < minspk
%                         continue;
%                     end
                    
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
        end
        
        if lock
            idxI = sst.pThphs<0.05;
        else 
            idxI = sst.pThphs>=0;
        end

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

        nanidx = ~isnan(depth_binned);
        uniq_depth = unique(depth_binned(nanidx));
        % for each uniqe bin find all the phases
        binned_struct = [];
        for i = 1:length(uniq_depth)
            idx3 = depth_binned==uniq_depth(i);
            % phase
            binned_struct.mPhaseBip(i) = nanmean(sst.mThphs(idx3&idxI&isbip));
            binned_struct.SEMPhaseBip(i) = calc_sem(sst.mThphs(idx3&idxI&isbip));
            binned_struct.mPhasePos(i) = nanmean(sst.mThphs(idx3&idxI&ispos));
            binned_struct.SEMPhasePos(i) = calc_sem(sst.mThphs(idx3&idxI&ispos));

            binned_struct.mPhasePyr(i) = nanmean(sst.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.SEMPhasePyr(i) = calc_sem(sst.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.mPhaseInt(i) = nanmean(sst.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));
            binned_struct.SEMPhaseInt(i) = calc_sem(sst.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));   
           
            % R
            binned_struct.rPhaseBip(i) = nanmean(sst.rThphs(idx3&idxI&isbip));
            binned_struct.SEMrPhaseBip(i) = calc_sem(sst.rThphs(idx3&idxI&isbip));
            binned_struct.rPhasePos(i) = nanmean(sst.rThphs(idx3&idxI&ispos));
            binned_struct.SEMrPhasePos(i) = calc_sem(sst.rThphs(idx3&idxI&ispos));

            binned_struct.rPhasePyr(i) = nanmean(sst.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.SEMrPhasePyr(i) = calc_sem(sst.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.rPhaseInt(i) = nanmean(sst.rThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));
            binned_struct.SEMrPhaseInt(i) = calc_sem(sst.rThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));   
            
            binned_struct.pval_TH_PYRbpi(i) = utest( sst.mThphs(idx3&idxI&isbip),sst.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.pval_THr_PYRbpi(i) = utest( sst.rThphs(idx3&idxI&isbip),sst.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));

            sumbpiTH(i) = sum(idx3&isbip);
            sumbpiTHL(i) = sum(idx3&isbip&idxI);

            sumPYRTH(i) = sum(idx3&~isbip&~ispos&sst.pyr);
            sumPYRTHL(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);

            sumINTTH(i) = sum(idx3&~isbip&~ispos&~sst.pyr);
            sumINTTHL(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr); 
            
            sumposTH(i) = sum(idx3&ispos);
            sumposTHL(i) = sum(idx3&ispos&idxI);
            
            binned_struct.depth(i) = uniq_depth(i);
        end

        figure;
        errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
        hold on,
        errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
        errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
        errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*', 'horizontal', 'Color', colors_NPB( 2, : ))
        title ('Theta-phase lock vs depth')
        xlabel ('Phase during theta')
        ylabel ('Depth')
        for i =1:4
            set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
            set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
            set( gca, 'tickdir', 'out', 'box', 'off' )
        end
        for i = 1:length(binned_struct.depth)
            if binned_struct.pval_TH_PYRbpi(i)<0.05
                errorbar(6.3,binned_struct.depth(i),[],'*k')
                text(6.35,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_TH_PYRbpi(i)))   
            end
            text(6.65,binned_struct.depth(end)+5,sprintf('Bpi,Pyr,Int, Pos',binned_struct.pval_TH_PYRbpi(i)))                     
            text(6.7,binned_struct.depth(i),sprintf('%d,%d,%d,%d',sumbpiTHL(i),sumPYRTHL(i),sumINTTHL(i),sumposTHL(i)))
        end
        ylim ([-150 150]);
        
        figure;
        errorbar(binned_struct.rPhasePyr, binned_struct.depth,binned_struct.SEMrPhasePyr,'*', 'horizontal', 'Color', colors_IP_light( 2, : ))
        hold on,
        errorbar(binned_struct.rPhaseInt, binned_struct.depth,binned_struct.SEMrPhaseInt,'*', 'horizontal', 'Color', colors_IP_light( 1, : ))
        errorbar(binned_struct.rPhaseBip, binned_struct.depth,binned_struct.SEMrPhaseBip,'*', 'horizontal', 'Color', colors_NPB( 3, : ))
        errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'*', 'horizontal', 'Color', colors_NPB( 2, : ))
        title ('Theta-R lock vs depth')
        xlabel ('R during theta')
        ylabel ('Depth')
        for i =1:4
            set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
            set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
            set( gca, 'tickdir', 'out', 'box', 'off' )
        end
        for i = 1:length(binned_struct.depth)
            if binned_struct.pval_TH_PYRbpi(i)<0.05
                errorbar(0.95,binned_struct.depth(i),[],'*k')
                text(0.955,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_TH_PYRbpi(i)))   
            end
            text(1.02,binned_struct.depth(end)+5,sprintf('Bpi,Pyr,Int',binned_struct.pval_TH_PYRbpi(i)))                     
            text(1.03,binned_struct.depth(i),sprintf('%d,%d,%d',sumbpiTHL(i),sumPYRTHL(i),sumINTTHL(i)))
        end
        xlim ([0 1]);
end
if saveST
    save ([dir 'struct_punits_biphasic_' savedate], 'sst')
end
return