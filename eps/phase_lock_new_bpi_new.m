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
% load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_29mar22.mat', '-mat');
load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_lean_28apr22.mat', '-mat')
load ('/media/shirly/C22865A128659567/mice/EPS/struct_phaseRippleTheta_linearDepth_12may22.mat', '-mat')
colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_PB_light        = [90, 219, 179; 106, 151, 232] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR
isbip           = ~isnan(sst.bpi) & ~isinf( sst.bpi );
ispos           = sst.extremum > 0 & ~isbip;

%add idx for multimodal biphasic
Nbip = sum(isbip);
wherebip = find(isbip);
idx_bpi_mult = false(length(isbip),1);
for i = 1:Nbip
    Vb = sst.vB(wherebip(i),:);
    Va = sst.vA(wherebip(i),:);
    nanidx = isnan(Vb);
    thidx = Va<-40;
    idx_multi = nanidx&thidx;
    if sum(idx_multi)>0
        idx_bpi_mult(wherebip(i)) = true;
    end
end

% add idx for bpi that are also place cellslslslsls
load ('/media/shirly/C22865A128659567/mice/EPS/bpi_field_struct.mat', '-mat')
for i = 1:length(bpi_field_struct.filename)
    afilename = bpi_field_struct.filename{i};
    idxname = contains(sst.filebase,afilename);
    idxclu = ismember(sst.shankclu,bpi_field_struct.shankclu(i,1:2),'rows');
    idxbpiplace(i) = find( idxname&idxclu);
end
bpiplace = unique(idxbpiplace);
idx_bpiplace = false(length(isbip),1);
idx_bpiplace(bpiplace)=true;

                if gather
                 % for the full dataset
%                         sst.mPhase = NaN(length(sst.shankclu),1);
%                         sst.rPhase = NaN(length(sst.shankclu),1);
%                         sst.pPhase = NaN(length(sst.shankclu),1);
%                          phasBP.ripGain = NaN(length(sst.shankclu),size(phsS.phsBins,2));

                    for i=1:length(x)
                        filebase  = sprintf('//media//shirly//C22865A128659567//mice//EPS//hfo//%s%d.hfo.spiking.spontaneous.details', str_, x(i));
                        phsS =load (filebase, '-mat');
                        phsS = phsS.stats;
                        str_f = sprintf ('%s%d', str_, x(i));
                        idxS = ismember(sst.filebase,str_f);
                        idxSf = find(idxS);
%                         sst.mPhase (idxSf) = phsS.mPhase;
%                         sst.rPhase (idxSf) = phsS.rPhase;
%                         sst.pPhase (idxSf) = phsS.pval;
%                         sst.depth (idxSf)  = phsS.depth;
                        phasBP.ripGain(idxSf,:)  = phsS.gain;

%                         sst.mPhase (idxSf) = NaN;
%                         sst.rPhase (idxSf) = NaN;
%                         sst.pPhase (idxSf) = NaN;
%                         sst.depth (idxSf)  = NaN;
                    end
                end
               %-----------------------------------%
% compute depth in um based on probe resolution
fnames              = phasBP.filebase;
ufnames             = unique( phasBP.filebase );
idx15 = (contains(phasBP.filebase, 'mC41') | contains(phasBP.filebase, 'mP20')...
    | contains(phasBP.filebase, 'nF79')| contains(phasBP.filebase, 'mP93')...
    | contains(phasBP.filebase, 'mF105')| contains(phasBP.filebase, 'mF108'));
depth( idx15 ,1)          = phasBP.depth( idx15 ) * 15;         % [um]
depth(~idx15,1)         = phasBP.depth( ~idx15 ) * 20;        % [um]

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
               idxI = phasBP.pPhase < 0.05;
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
                binned_struct.mPhaseBip(i) = nanmean(phasBP.mPhase(idx3&idxI&isbip));
                binned_struct.SEMPhaseBip(i) = calc_sem(phasBP.mPhase(idx3&idxI&isbip));
                binned_struct.mPhasePos(i) = nanmean(phasBP.mPhase(idx3&idxI&ispos));
                binned_struct.SEMPhasePos(i) = calc_sem(phasBP.mPhase(idx3&idxI&ispos));

                binned_struct.mPhasePyr(i) = nanmean(phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.SEMPhasePyr(i) = calc_sem(phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.mPhaseInt(i) = nanmean(phasBP.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                binned_struct.SEMPhaseInt(i) = calc_sem(phasBP.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   
                
%                     % mean preffered phases of multimodal biphasic
%                     binned_struct.mPhaseMBip(i) = nanmean(phasBP.mPhase(idx3&idxI&isbip&idx_bpi_mult));
%                     binned_struct.SEMPhaseMBip(i) = calc_sem(phasBP.mPhase(idx3&idxI&isbip&idx_bpi_mult));
%                     binned_struct.mPhaseNBip(i) = nanmean(phasBP.mPhase(idx3&idxI&isbip&~idx_bpi_mult));
%                     binned_struct.SEMPhaseNBip(i) = calc_sem(phasBP.mPhase(idx3&idxI&isbip&~idx_bpi_mult));
%                     % mean preffered phases of biphasic who are place cells
%                     binned_struct.mPhasePBip(i) = nanmean(phasBP.mPhase(idx3&idxI&isbip&idx_bpiplace));
%                     binned_struct.SEMPhasePBip(i) = calc_sem(phasBP.mPhase(idx3&idxI&isbip&idx_bpiplace));
%                     binned_struct.mPhaseNPBip(i) = nanmean(phasBP.mPhase(idx3&idxI&isbip&~idx_bpiplace));
%                     binned_struct.SEMPhaseNPBip(i) = calc_sem(phasBP.mPhase(idx3&idxI&isbip&~idx_bpiplace));                    
                
                % mean R
                binned_struct.rPhaseBip(i) = nanmean(phasBP.rPhase(idx3&idxI&isbip));
                binned_struct.SEMrPhaseBip(i) = calc_sem(phasBP.rPhase(idx3&idxI&isbip));
                binned_struct.rPhasePos(i) = nanmean(phasBP.rPhase(idx3&idxI&ispos));
                binned_struct.SEMrPhasePos(i) = calc_sem(phasBP.rPhase(idx3&idxI&ispos));

                binned_struct.rPhasePyr(i) = nanmean(phasBP.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.SEMrPhasePyr(i) = calc_sem(phasBP.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.rPhaseInt(i) = nanmean(phasBP.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                binned_struct.SEMrPhaseInt(i) = calc_sem(phasBP.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 
                binned_struct.depth(i) = uniq_depth(i);

                binned_struct.pval_PYRbpi(i) = utest( phasBP.mPhase(idx3&idxI&isbip),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.rpval_PYRbpi(i) = utest( phasBP.rPhase(idx3&idxI&isbip),phasBP.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.pval_PYRMbpi(i) = utest( phasBP.mPhase(idx3&idxI&isbip&idx_bpi_mult),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.pval_PYRNbpi(i) = utest( phasBP.mPhase(idx3&idxI&isbip&~idx_bpi_mult),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.pval_PYRPbpi(i) = utest( phasBP.mPhase(idx3&idxI&isbip&idx_bpiplace),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
%                 binned_struct.pval_PYRNPbpi(i) = utest( phasBP.mPhase(idx3&idxI&isbip&~idx_bpiplace),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

                binned_struct.pval_PYRpos(i) = utest( phasBP.mPhase(idx3&idxI&ispos),phasBP.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                binned_struct.rpval_PYRpos(i) = utest( phasBP.rPhase(idx3&idxI&ispos),phasBP.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                
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
                
%                 sumMbpiS(i) = sum(idx3&isbip&idx_bpi_mult);
%                 sumMbpiSL(i) = sum(idx3&isbip&idxI&idx_bpi_mult);
%                 sumNbpiS(i) = sum(idx3&isbip&~idx_bpi_mult);
%                 sumNbpiSL(i) = sum(idx3&isbip&idxI&~idx_bpi_mult);
%                 sumPbpiS(i) = sum(idx3&isbip&idx_bpiplace);
%                 sumPbpiSL(i) = sum(idx3&isbip&idxI&idx_bpiplace);                    
%                 sumNPbpiS(i) = sum(idx3&isbip&~idx_bpiplace);
%                 sumNPbpiSL(i) = sum(idx3&isbip&idxI&~idx_bpiplace);                
                
                sumposS(i) = sum(idx3&ispos);
                sumposSL(i) = sum(idx3&ispos&idxI);
                
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
%                  errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'.', 'horizontal', 'Color', colors_NPB( 2, : ))
                title ('Ripple-phase lock vs depth')
                xlabel ('Phase during ripple')
                ylabel ('Depth')
                for i = 1:length(binned_struct.depth)
                    if binned_struct.pval_PYRbpi(i)<0.05
                        errorbar(6.25,binned_struct.depth(i),[],'*k')
                        text(6.3,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRbpi(i)))   
                    end
                    if binned_struct.pval_PYRpos(i)<0.05
                        errorbar(5.8,binned_struct.depth(i),[],'*m')
                        text(5.85,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRpos(i)))   
                    end
                    text(6.75,binned_struct.depth(end)+7,sprintf('Pos,Bpi,Pyr,Int',binned_struct.pval_PYRbpi(i)))                     
                    text(6.85,binned_struct.depth(i),sprintf('%d,%d,%d,%d',sumposSL(i),sumbpiSL(i),sumPYRSL(i),sumINTSL(i)))
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
            
            if graphics(5)           
                figure;
                errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
                hold on,
%                 errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
                errorbar(binned_struct.mPhaseMBip, binned_struct.depth,binned_struct.SEMPhaseMBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
%                 errorbar(binned_struct.mPhaseNBip, binned_struct.depth,binned_struct.SEMPhaseNBip,'.', 'horizontal', 'Color', colors_PB_light( 2, : ))
%                  errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'.', 'horizontal', 'Color', colors_NPB( 2, : ))
                title ('Ripple-phase lock vs depth')
                xlabel ('Phase during ripple')
                ylabel ('Depth')
                for i = 1:length(binned_struct.depth)
                    if binned_struct.pval_PYRMbpi(i)<0.05
                        errorbar(6.25,binned_struct.depth(i),[],'*k')
                        text(6.3,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRMbpi(i)))   
                    end
                    if binned_struct.pval_PYRNbpi(i)<0.05
                        errorbar(5.8,binned_struct.depth(i),[],'*b')
                        text(5.85,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRNbpi(i)))   
                    end
                    text(6.75,binned_struct.depth(end)+7,sprintf('BpiN,BpiM,Pyr',binned_struct.pval_PYRbpi(i)))                     
                    text(6.85,binned_struct.depth(i),sprintf('%d,%d,%d,%d',sumNbpiSL(i),sumMbpiSL(i),sumPYRSL(i)))
                end
                for i =1:4
                    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
                    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
                    set( gca, 'tickdir', 'out', 'box', 'off' )
                end
            end
            if graphics(6)           
                figure;
                errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'.', 'horizontal', 'Color', colors_IP_light( 2, : ))
                hold on,
%                 errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'.', 'horizontal', 'Color', colors_IP_light( 1, : ))
                errorbar(binned_struct.mPhasePBip, binned_struct.depth,binned_struct.SEMPhasePBip,'.', 'horizontal', 'Color', colors_NPB( 3, : ))
                errorbar(binned_struct.mPhaseNPBip, binned_struct.depth,binned_struct.SEMPhaseNPBip,'.', 'horizontal', 'Color', colors_PB_light( 2, : ))
%                  errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'.', 'horizontal', 'Color', colors_NPB( 2, : ))
                title ('Ripple-phase lock vs depth')
                xlabel ('Phase during ripple')
                ylabel ('Depth')
                for i = 1:length(binned_struct.depth)
                    if binned_struct.pval_PYRPbpi(i)<0.05
                        errorbar(6.25,binned_struct.depth(i),[],'*k')
                        text(6.3,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRMbpi(i)))   
                    end
                    if binned_struct.pval_PYRNPbpi(i)<0.05
                        errorbar(5.8,binned_struct.depth(i),[],'*b')
                        text(5.85,binned_struct.depth(i),sprintf('%0.2g',binned_struct.pval_PYRNPbpi(i)))   
                    end
                    text(6.75,binned_struct.depth(end)+7,sprintf('BpiNP,BpiP,Pyr',binned_struct.pval_PYRbpi(i)))                     
                    text(6.85,binned_struct.depth(i),sprintf('%d,%d,%d,%d',sumNPbpiSL(i),sumPbpiSL(i),sumPYRSL(i)))
                end
                for i =1:4
                    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
                    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
                    set( gca, 'tickdir', 'out', 'box', 'off' )
                end
            end   
            
            if graphics(7)           
                figure;
%                 t_pyr = ~isbip&~ispos&sst.pyr;
%                 t_int = ~isbip&~ispos&~sst.pyr;
%                 t_pun = ispos;
%                 t_bip = isbip;
                t_pyr = ~isbip&~ispos&sst.pyr&idxI;
                t_int = ~isbip&~ispos&~sst.pyr&idxI;
                t_pun = ispos&idxI;
                t_bip = isbip&idxI;
                tstr                    = 'Phase lock ripple angle (sig)';
                id                     = [ t_pyr, t_int, t_pun, t_bip];
                param_type              = 'circ';
                prm                   = phasBP.mPhase;
                make_TPR_TPP_figures_histograms( prm, id, tstr, param_type, 'color_gro', 1:4, 'edges',  phsCenters)

                % test for equal medians using circular non-parametric "ANOVA" (CM test)
                vidx                    = [ t_pyr | t_int | t_pun| t_bip];
                xx                      = prm( vidx );
                [ ~, gidx ]             = find( id( vidx, : ) );
                pval_all                = circ_cmtest( xx, gidx );
                fprintf( 1, '%s, p=%0.3g\n', tstr, pval_all )
                pairs                   = make_ai( length( unique( gidx ) ) );
                npairs                  = size( pairs, 1 );
                pvals                   = NaN( npairs, 1 );
                for pidx                = 1 : npairs
                    g1                  = xx( gidx == pairs( pidx, 1 ) );
                    g2                  = xx( gidx == pairs( pidx, 2 ) );
                    pvals( pidx, 1 )    = circ_cmtest( g1, g2 );
                    %[ ~, pvals( pidx, 2 ) ]         = wheeler( g1, g2 );
                    fprintf( 1, '%s, g%d vs. g%d, p=%0.3g\n', tstr, pairs( pidx, 1 ), pairs( pidx, 2 ), pvals( pidx, 1 ) )
                end
            end

            if graphics(8)           
                figure;            
                t_pyr = ~isbip&~ispos&sst.pyr;
                t_int = ~isbip&~ispos&~sst.pyr;
                t_pun = ispos;
                t_bip = isbip;
                t_pyrL = ~isbip&~ispos&sst.pyr&idxI;
                t_intL = ~isbip&~ispos&~sst.pyr&idxI;
                t_punL = ispos&idxI;
                t_bipL = isbip&idxI;
                
                barwerror([1 2 3 4],[sum(t_pyrL) sum(t_intL)...
                sum(t_punL)  sum(t_bipL)], [], colors_NPB(3,:));
                hold on
                h = barwerror([1 2 3 4],[sum(t_pyr) sum(t_int)...
                sum(t_pun)  sum(t_bip)]);
                set(h,'FaceColor','none');
                set(h,'EdgeColor','k');
                set( gca, 'tickdir', 'out', 'box', 'off' );
                title('all units locked/all spontaneous');
                ylabel('count');
                set(gca, 'XTickLabel',{ 'PYR' 'INT' 'POS' 'BPI' })
                text(1, 7500 ,sprintf('%.2f',sum(t_pyrL)/sum(t_pyr)));
                text(2, 7500 ,sprintf('%.2f',sum(t_intL)/sum(t_int)));
                text(3, 7500 ,sprintf('%.2f',sum(t_punL)/sum(t_pun)));
                text(4, 7500 ,sprintf('%.2f',sum(t_bipL)/sum(t_bip)));
            end
            
            if graphics(9)           
                figure;            
                t_pyr = ~isbip&~ispos&sst.pyr;
                t_int = ~isbip&~ispos&~sst.pyr;
                t_pun = ispos;
                t_bip = isbip;
                t_pyrL = ~isbip&~ispos&sst.pyr&idxI;
                t_intL = ~isbip&~ispos&~sst.pyr&idxI;
                t_punL = ispos&idxI;
                t_bipL = isbip&idxI;
                
                barwerror([1 2 3 4],[sum(t_pyrL)./sum(t_pyr) sum(t_intL)./sum(t_int)...
                sum(t_punL)./sum(t_pun) sum(t_bipL)./sum(t_bip)], [], colors_NPB(3,:));
                set(h,'FaceColor','none');
                set(h,'EdgeColor','k');
                set( gca, 'tickdir', 'out', 'box', 'off' );
                title('fraction of locked/all units');
                ylabel('count');
                set(gca, 'XTickLabel',{ 'PYR' 'INT' 'POS' 'BPI' })
                text(1, 0.8 ,sprintf('%.2f',sum(t_pyrL)/sum(t_pyr)));
                text(2, 0.8 ,sprintf('%.2f',sum(t_intL)/sum(t_int)));
                text(3, 0.8 ,sprintf('%.2f',sum(t_punL)/sum(t_pun)));
                text(4, 0.8 ,sprintf('%.2f',sum(t_bipL)/sum(t_bip)));
                ylim([0 1])
            end            
            
            if graphics(10)
                gain_pyr = phasBP.ripGain(t_pyrL,:);
                gain_int = phasBP.ripGain(t_intL,:);
                gain_bip = phasBP.ripGain(t_bipL,:);
                gain_pun = phasBP.ripGain(t_punL,:);
%                 gain_pyr = gain_pyr(:,80:360);               
%                 gain_int = gain_int(:,80:360);                
%                 gain_bip = gain_bip(:,80:360);                
%                 gain_pun = gain_pun(:,80:360);


                figure,patch_band([1:421],nanmean(gain_pyr),calc_sem(gain_pyr),colors_IP_light( 2, : ))
                hold on
                patch_band([1:421],nanmean(gain_int),calc_sem(gain_int),colors_IP_light( 1, : ), colors_IP_light( 1, : ))
                patch_band([1:421],nanmean(gain_bip),calc_sem(gain_bip),colors_NPB( 3, : ), colors_NPB( 3, : ))
                patch_band([1:421],nanmean(gain_pun),calc_sem(gain_pun),colors_NPB( 2, : ), colors_NPB( 2, : ))

                figure,patch_band([1:421],nanmean(gain_pyr),[],colors_IP_light( 2, : ))
                hold on
                patch_band([1:421],nanmean(gain_int),[],colors_IP_light( 1, : ), colors_IP_light( 1, : ))
                patch_band([1:421],nanmean(gain_bip),[],colors_NPB( 3, : ), colors_NPB( 3, : ))
                patch_band([1:421],nanmean(gain_pun),[],colors_NPB( 2, : ), colors_NPB( 2, : ))                
            end
            
%             examples:
            [~,~, phsEdges]             = circ_hist( [0, 2*pi], 20 );
%             pyr 
            idxPYREX = ~isbip&~ispos&sst.pyr&phasBP.rPhase<0.81&phasBP.rPhase>0.78&phasBP.mPhase>0.97*pi&phasBP.mPhase<1*pi;
            phsS = load ('/media/shirly/C22865A128659567/mice/EPS/hfo/mA234_18.hfo.spiking.spontaneous.details', '-mat' )
            
            filename = 'mA234_18';filebase = filebase_lookup(filename,1);
            figure, plot_ss(filebase, [4 31]);
            binsize = phsS.stats.phsbins(2)-phsS.stats.phsbins(1);
            edges       = ( floor( min( phsS.stats.phsbins ) ) - binsize/2 ) : binsize : ( ceil( max( phsS.stats.phsbins ) ) + binsize/2 );
            figure;
            bh = polarhistogram('BinCounts',phsS.stats.phshists(50,:), 'BinEdges',edges);
            hold on;
            set( bh, 'EdgeColor', colors_IP_light( 2, : ), 'FaceColor', colors_IP_light( 2, : ) );
            phi = phsS.stats.mPhase(50,:);
            R   = phsS.stats.rPhase(50,:);
            maxcount = max(bh.Values);
            polarplot( [phi; phi], [0 R*maxcount], 'k', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            tit = sprintf( 'pyr, %s %d.%d med=%0.3g%c, R=%.2f',filename, phsS.stats.shankclu(50,1), phsS.stats.shankclu(50,2), phi/pi, char(960), R);            
            title( tit )
            %             other potentials
%             filename = 'mA234_29';filebase = filebase_lookup(filename,1);
%             figure, plot_ss(filebase, [1 55]);
%             filename = 'mA234_41';filebase = filebase_lookup(filename,1);
%             figure, plot_ss(filebase, [2 46]);
%             figure, plot_ss(filebase, [2 51]);

%             int 
            idxINTEX = ~isbip&~ispos&~sst.pyr&phasBP.rPhase<0.4&phasBP.rPhase>0.35&phasBP.mPhase>1.3*pi&phasBP.mPhase<1.4*pi;
            filename = 'mC41_36';filebase = filebase_lookup(filename,1);
            figure, plot_ss(filebase, [1 60]);
            phsS = load ('/media/shirly/C22865A128659567/mice/EPS/hfo/mC41_36.hfo.spiking.spontaneous.details', '-mat')
            figure;
            bh = polarhistogram('BinCounts',phsS.stats.phshists(38,:), 'BinEdges',edges);
            hold on;
            set( bh, 'EdgeColor', colors_IP_light( 1, : ), 'FaceColor', colors_IP_light( 1, : ) );
            phi = phsS.stats.mPhase(38,:);
            R   = phsS.stats.rPhase(38,:);
            maxcount = max(bh.Values);
            polarplot( [phi; phi], [0 R*maxcount], 'k', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            tit = sprintf( 'int, %s %d.%d med=%0.3g%c, R=%.2f',filename, phsS.stats.shankclu(38,1), phsS.stats.shankclu(38,2), phi/pi, char(960), R);            
            title( tit )
            
%             pun 
            idxposEX = ~isbip&ispos&phasBP.rPhase<0.42&phasBP.rPhase>0.25&phasBP.mPhase>0.28*pi&phasBP.mPhase<0.4*pi;
            filename = 'mC41_21';filebase = filebase_lookup(filename,1);
            figure, plot_ss(filebase, [2 9]);
            phsS = load ('/media/shirly/C22865A128659567/mice/EPS/hfo/mC41_21.hfo.spiking.spontaneous.details', '-mat')
            figure;
            bh = polarhistogram('BinCounts',phsS.stats.phshists(17,:), 'BinEdges',edges);
            hold on;
            set( bh, 'EdgeColor', colors_NPB( 2, : ), 'FaceColor', colors_NPB( 2, : ) );
            phi = phsS.stats.mPhase(17,:);
            R   = phsS.stats.rPhase(17,:);
            maxcount = max(bh.Values);
            polarplot( [phi; phi], [0 R*maxcount], 'k', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            tit = sprintf( 'pun, %s %d.%d med=%0.3g%c, R=%.2f',phsS.stats.filename{1}, phsS.stats.shankclu(17,1), phsS.stats.shankclu(17,2), phi/pi, char(960), R);            
            title( tit )
            
            %             bip 
            idxbipEX = isbip&~ispos&phasBP.rPhase<0.63&phasBP.rPhase>0.5&phasBP.mPhase>0.8*pi&phasBP.mPhase<1*pi;
            filename = 'mDS2_13';filebase = filebase_lookup(filename,1);
            figure, plot_ss(filebase, [1 16]);
            phsS = load ('/media/shirly/C22865A128659567/mice/EPS/hfo/mDS2_13.hfo.spiking.spontaneous.details', '-mat')
            figure;
            bh = polarhistogram('BinCounts',phsS.stats.phshists(14,:), 'BinEdges',edges);
            hold on;
            set( bh, 'EdgeColor', colors_NPB( 3, : ), 'FaceColor', colors_NPB( 3, : ) );
            phi = phsS.stats.mPhase(14,:);
            R   = phsS.stats.rPhase(14,:);
            maxcount = max(bh.Values);
            polarplot( [phi; phi], [0 R*maxcount], 'k', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            tit = sprintf( 'bip, %s %d.%d med=%0.3g%c, R=%.2f',phsS.stats.filename{1}, phsS.stats.shankclu(14,1), phsS.stats.shankclu(14,2), phi/pi, char(960), R);            
            title( tit )
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
%             phasBP.mThphs = NaN(length(phasBP.shankclu),1);
%             phasBP.rThphs = NaN(length(phasBP.shankclu),1);
%             phasBP.pThphs = NaN(length(phasBP.shankclu),1);
            for i=1:length(x)
                filebase  = sprintf('//media//shirly//C22865A128659567//mice//%s//dat//%s%d//%s%d', str, str_, x(i), str_, x(i));
                % load the phases:
                load ([ filebase '.phs_ripp' ], '-mat');
                phs                                = wrapTo2Pi(phs);

                % get theta times
                THETAper                       = get_states( filebase, 'the' ); % Periods of theta, phsFs

                phsFs = 1250; % [Hz]
                % get histogram time per each side
                [~,~, phsEdges]             = circ_hist( [0, 2*pi], 20 );
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
                idxS = ismember(phasBP.filebase,str_f);
                idxSf = find(idxS);
                phasBP.mThphs (idxSf) = mThphs;
                phasBP.rThphs (idxSf) = rThphs;
                phasBP.pThphs (idxSf) = pThphs;
            end
        end
        
        if lock
            idxI = phasBP.pThphs<0.05;
        else 
            idxI = phasBP.pThphs>=0;
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
            binned_struct.mPhaseBip(i) = nanmean(phasBP.mThphs(idx3&idxI&isbip));
            binned_struct.SEMPhaseBip(i) = calc_sem(phasBP.mThphs(idx3&idxI&isbip));
            binned_struct.mPhasePos(i) = nanmean(phasBP.mThphs(idx3&idxI&ispos));
            binned_struct.SEMPhasePos(i) = calc_sem(phasBP.mThphs(idx3&idxI&ispos));

            binned_struct.mPhasePyr(i) = nanmean(phasBP.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.SEMPhasePyr(i) = calc_sem(phasBP.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.mPhaseInt(i) = nanmean(phasBP.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));
            binned_struct.SEMPhaseInt(i) = calc_sem(phasBP.mThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));   
           
            % R
            binned_struct.rPhaseBip(i) = nanmean(phasBP.rThphs(idx3&idxI&isbip));
            binned_struct.SEMrPhaseBip(i) = calc_sem(phasBP.rThphs(idx3&idxI&isbip));
            binned_struct.rPhasePos(i) = nanmean(phasBP.rThphs(idx3&idxI&ispos));
            binned_struct.SEMrPhasePos(i) = calc_sem(phasBP.rThphs(idx3&idxI&ispos));

            binned_struct.rPhasePyr(i) = nanmean(phasBP.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.SEMrPhasePyr(i) = calc_sem(phasBP.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.rPhaseInt(i) = nanmean(phasBP.rThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));
            binned_struct.SEMrPhaseInt(i) = calc_sem(phasBP.rThphs(idx3&idxI&~isbip&~ispos&~sst.pyr));   
            
            binned_struct.pval_TH_PYRbpi(i) = utest( phasBP.mThphs(idx3&idxI&isbip),phasBP.mThphs(idx3&idxI&~isbip&~ispos&sst.pyr));
            binned_struct.pval_THr_PYRbpi(i) = utest( phasBP.rThphs(idx3&idxI&isbip),phasBP.rThphs(idx3&idxI&~isbip&~ispos&sst.pyr));

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
        errorbar(binned_struct.mPhasePyr, binned_struct.depth,binned_struct.SEMPhasePyr,'+', 'horizontal', 'Color', colors_IP_light( 2, : ))
        hold on,
        errorbar(binned_struct.mPhaseInt, binned_struct.depth,binned_struct.SEMPhaseInt,'+', 'horizontal', 'Color', colors_IP_light( 1, : ))
        errorbar(binned_struct.mPhaseBip, binned_struct.depth,binned_struct.SEMPhaseBip,'+', 'horizontal', 'Color', colors_NPB( 3, : ))
%         errorbar(binned_struct.mPhasePos, binned_struct.depth,binned_struct.SEMPhasePos,'+', 'horizontal', 'Color', colors_NPB( 2, : ))
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
        
        if graphics(7)           
            figure;
                t_pyr = ~isbip&~ispos&sst.pyr;
                t_int = ~isbip&~ispos&~sst.pyr;
                t_pun = ispos;
                t_bip = isbip;
%             t_pyr = ~isbip&~ispos&sst.pyr&idxI;
%             t_int = ~isbip&~ispos&~sst.pyr&idxI;
%             t_pun = ispos&idxI;
%             t_bip = isbip&idxI;
            tstr                    = 'Phase lock theta angle (sig)';
            id                     = [ t_pyr, t_int, t_pun, t_bip];
            param_type              = 'circ';
            prm                   = phasBP.mThphs;
            make_TPR_TPP_figures_histograms( prm, id, tstr, param_type, 'color_gro', 1:4 )

            % test for equal medians using circular non-parametric "ANOVA" (CM test)
            vidx                    = [ t_pyr | t_int | t_pun| t_bip];
            xx                      = prm( vidx );
            [ ~, gidx ]             = find( id( vidx, : ) );
            pval_all                = circ_cmtest( xx, gidx );
            fprintf( 1, '%s, p=%0.3g\n', tstr, pval_all )
            pairs                   = make_ai( length( unique( gidx ) ) );
            npairs                  = size( pairs, 1 );
            pvals                   = NaN( npairs, 1 );
            for pidx                = 1 : npairs
                g1                  = xx( gidx == pairs( pidx, 1 ) );
                g2                  = xx( gidx == pairs( pidx, 2 ) );
                pvals( pidx, 1 )    = circ_cmtest( g1, g2 );
                %[ ~, pvals( pidx, 2 ) ]         = wheeler( g1, g2 );
                fprintf( 1, '%s, g%d vs. g%d, p=%0.3g\n', tstr, pairs( pidx, 1 ), pairs( pidx, 2 ), pvals( pidx, 1 ) )
            end
        end
end
if saveST
    save ([dir 'struct_punits_biphasic_' savedate], 'sst')
end
return