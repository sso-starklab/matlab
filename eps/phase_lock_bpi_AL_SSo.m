%% analysis for theta and ripple locking
%
%
% arguments:            gather {0} to repopulate the sst
%                       str        animal name (e.g. mC41)
%                       x          a vector of session numbers (e.g. [15:16 20])
%                       save {0}   save sst structure after run
%                       savedate   string of suffix for sst save
%                       
% note - support for theta analysis not finished

function phase_lock_bpi (mode, state, lock, varargin)
%% ripple
% argument handling
[ gather, x, str, ...
    save, savedate, ...
    graphics]            = ParseArgPairs( ...
    { 'gather', 'x', 'str',  ...
    'save', 'savedate'...
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
load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_06apr22.mat', '-mat');
colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR
isbip           = ~isnan(sst.bpi) & ~isinf( sst.bpi );
ispos           = sst.extremum > 0 & ~isbip;

% compute depth in um based on probe resolution
fnames              = sst.filebase;
ufnames             = unique( sst.filebase );
nufnames            = length( ufnames );
tmpcell             = strfind( ufnames, 'mC41' ); % presently supporting only one animal
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
 binned_struct_LS = [];
 binned_struct_US = [];
 binned_struct_LI = [];
 binned_struct_UI = [];
 namesSI = {'spontaneous' , 'induced'};
        for lock = [0 1]
            for j = 1:2
                state = namesSI{j};
                for i = 1:length(uniq_depth)
                  if lock && isequal(state,'spontaneous')
                        idxI = sst.pPhase < 0.05;
                        idx3 = depth_binned==uniq_depth(i);
                        % mean preffered phases
                        binned_struct_LS.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_LS.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_LS.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
                        binned_struct_LS.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));

                        binned_struct_LS.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LS.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        % mean R
                        binned_struct_LS.rPhaseBip(i) = nanmean(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_LS.SEMrPhaseBip(i) = calc_sem(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_LS.rPhasePos(i) = nanmean(sst.rPhase(idx3&idxI&ispos));
                        binned_struct_LS.SEMrPhasePos(i) = calc_sem(sst.rPhase(idx3&idxI&ispos));

                        binned_struct_LS.rPhasePyr(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.SEMrPhasePyr(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.rPhaseInt(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LS.SEMrPhaseInt(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 

                        binned_struct_LS.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.rpval_PYRbpi(i) = utest( sst.rPhase(idx3&idxI&isbip),sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

                        binned_struct_LS.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_LS.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_LS.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
                        binned_struct_LS.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));

                        binned_struct_LS.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LS.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LS.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        binned_struct_LS.depth(i) = uniq_depth(i);

                        binned_struct_LS.pval_IPYRbpi(i) = utest( sst.ImPhas(idx3&idxI&isbip),sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));   
                         binned_struct_LS.sumbpi(i) = sum(idx3&isbip&idxI);
                         binned_struct_LS.sumpos(i) = sum(idx3&ispos&idxI);                         
                         binned_struct_LS.sumPYR(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);
                         binned_struct_LS.sumINT(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);                       
                   else
                       idxI = sst.pPhase > 0;
                        idx3 = depth_binned==uniq_depth(i);
                        % mean preffered phases
                        binned_struct_US.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_US.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_US.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
                        binned_struct_US.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));

                        binned_struct_US.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_US.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        % mean R
                        binned_struct_US.rPhaseBip(i) = nanmean(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_US.SEMrPhaseBip(i) = calc_sem(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_US.rPhasePos(i) = nanmean(sst.rPhase(idx3&idxI&ispos));
                        binned_struct_US.SEMrPhasePos(i) = calc_sem(sst.rPhase(idx3&idxI&ispos));

                        binned_struct_US.rPhasePyr(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.SEMrPhasePyr(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.rPhaseInt(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_US.SEMrPhaseInt(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 

                        binned_struct_US.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.rpval_PYRbpi(i) = utest( sst.rPhase(idx3&idxI&isbip),sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

                        binned_struct_US.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_US.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_US.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
                        binned_struct_US.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));

                        binned_struct_US.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_US.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_US.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        binned_struct_US.depth(i) = uniq_depth(i);

                         binned_struct_US.pval_IPYRbpi(i) = utest( sst.ImPhas(idx3&idxI&isbip),sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr)); 
                         binned_struct_US.sumbpi(i) = sum(idx3&isbip&idxI);
                         binned_struct_US.sumpos(i) = sum(idx3&ispos&idxI);                         
                         binned_struct_US.sumPYR(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);
                         binned_struct_US.sumINT(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);
                  end   
                  if lock && isequal(state,'induced')
                        idxI = sst.IpPhas<0.05;
                        idx3 = depth_binned==uniq_depth(i);
                        % mean preffered phases
                        binned_struct_LI.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_LI.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_LI.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
                        binned_struct_LI.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));

                        binned_struct_LI.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LI.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        % mean R
                        binned_struct_LI.rPhaseBip(i) = nanmean(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_LI.SEMrPhaseBip(i) = calc_sem(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_LI.rPhasePos(i) = nanmean(sst.rPhase(idx3&idxI&ispos));
                        binned_struct_LI.SEMrPhasePos(i) = calc_sem(sst.rPhase(idx3&idxI&ispos));

                        binned_struct_LI.rPhasePyr(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.SEMrPhasePyr(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.rPhaseInt(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LI.SEMrPhaseInt(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 

                        binned_struct_LI.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.rpval_PYRbpi(i) = utest( sst.rPhase(idx3&idxI&isbip),sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

                        binned_struct_LI.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_LI.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_LI.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
                        binned_struct_LI.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));

                        binned_struct_LI.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_LI.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_LI.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        binned_struct_LI.depth(i) = uniq_depth(i);

                        binned_struct_LI.pval_IPYRbpi(i) = utest( sst.ImPhas(idx3&idxI&isbip),sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr)); 
                         binned_struct_LI.sumbpi(i) = sum(idx3&isbip&idxI);
                         binned_struct_LI.sumpos(i) = sum(idx3&ispos&idxI);                         
                         binned_struct_LI.sumPYR(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);
                         binned_struct_LI.sumINT(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);                         
                  else
                        idxI = sst.IpPhas > 0;
 idx3 = depth_binned==uniq_depth(i);
                        % mean preffered phases
                        binned_struct_UI.mPhaseBip(i) = nanmean(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_UI.SEMPhaseBip(i) = calc_sem(sst.mPhase(idx3&idxI&isbip));
                        binned_struct_UI.mPhasePos(i) = nanmean(sst.mPhase(idx3&idxI&ispos));
                        binned_struct_UI.SEMPhasePos(i) = calc_sem(sst.mPhase(idx3&idxI&ispos));

                        binned_struct_UI.mPhasePyr(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.SEMPhasePyr(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.mPhaseInt(i) = nanmean(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_UI.SEMPhaseInt(i) = calc_sem(sst.mPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        % mean R
                        binned_struct_UI.rPhaseBip(i) = nanmean(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_UI.SEMrPhaseBip(i) = calc_sem(sst.rPhase(idx3&idxI&isbip));
                        binned_struct_UI.rPhasePos(i) = nanmean(sst.rPhase(idx3&idxI&ispos));
                        binned_struct_UI.SEMrPhasePos(i) = calc_sem(sst.rPhase(idx3&idxI&ispos));

                        binned_struct_UI.rPhasePyr(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.SEMrPhasePyr(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.rPhaseInt(i) = nanmean(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_UI.SEMrPhaseInt(i) = calc_sem(sst.rPhase(idx3&idxI&~isbip&~ispos&~sst.pyr)); 

                        binned_struct_UI.pval_PYRbpi(i) = utest( sst.mPhase(idx3&idxI&isbip),sst.mPhase(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.rpval_PYRbpi(i) = utest( sst.rPhase(idx3&idxI&isbip),sst.rPhase(idx3&idxI&~isbip&~ispos&sst.pyr));

                        binned_struct_UI.mIPhaseBip(i) = nanmean(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_UI.SEMIPhaseBip(i) = calc_sem(sst.ImPhas(idx3&idxI&isbip));
                        binned_struct_UI.mIPhasePos(i) = nanmean(sst.ImPhas(idx3&idxI&ispos));
                        binned_struct_UI.SEMIPhasePos(i) = calc_sem(sst.ImPhas(idx3&idxI&ispos));

                        binned_struct_UI.mIPhasePyr(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.SEMIPhasePyr(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr));
                        binned_struct_UI.mIPhaseInt(i) = nanmean(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));
                        binned_struct_UI.SEMIPhaseInt(i) = calc_sem(sst.ImPhas(idx3&idxI&~isbip&~ispos&~sst.pyr));   

                        binned_struct_UI.depth(i) = uniq_depth(i);

                        binned_struct_UI.pval_IPYRbpi(i) = utest( sst.ImPhas(idx3&idxI&isbip),sst.ImPhas(idx3&idxI&~isbip&~ispos&sst.pyr)); 
                         binned_struct_UI.sumbpi(i) = sum(idx3&isbip&idxI);
                         binned_struct_UI.sumpos(i) = sum(idx3&ispos&idxI);                         
                         binned_struct_UI.sumPYR(i) = sum(idx3&idxI&~isbip&~ispos&sst.pyr);
                         binned_struct_UI.sumINT(i) = sum(idx3&idxI&~isbip&~ispos&~sst.pyr);   
                  end 

                end
            end
        end
struct_final.LI = binned_struct_LI;
struct_final.LS = binned_struct_LS;
struct_final.UI = binned_struct_UI;
struct_final.US = binned_struct_US;

k = sum(binned_struct_LS.sumbpi);
n = sum(binned_struct_US.sumbpi+binned_struct_LS.sumbpi);
bino_se_norm_bpi = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LS.sumPYR);
n = sum(binned_struct_US.sumPYR+binned_struct_LS.sumPYR);
bino_se_norm_pyr = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LS.sumINT);
n = sum(binned_struct_US.sumINT+binned_struct_LS.sumINT);
bino_se_norm_int = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LS.sumpos);
n = sum(binned_struct_US.sumpos+binned_struct_LS.sumpos);
bino_se_norm_pos = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )

k = sum(binned_struct_LI.sumbpi);
n = sum(binned_struct_UI.sumbpi+binned_struct_LI.sumbpi);
bino_se_norm_bpiI = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LI.sumPYR);
n = sum(binned_struct_UI.sumPYR+binned_struct_LI.sumPYR);
bino_se_norm_pyrI = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LI.sumINT);
n = sum(binned_struct_UI.sumINT+binned_struct_LI.sumINT);
bino_se_norm_intI = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )
k = sum(binned_struct_LI.sumpos);
n = sum(binned_struct_UI.sumpos+binned_struct_LI.sumpos);
bino_se_norm_posI = sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )

%Can you compute the mean ripple phases of the four groups regardless of depth?
mPhase = sst.mPhase;
rPhase = sst.rPhase;
ispyr = sst.pyr&~isbip&~ispos;
isINT = ~sst.pyr&~isbip&~ispos;

mIPhase = sst.ImPhas;

fig1 = figure;
subplot(1,2,1)
barwerror([1 2 3 4],[nanmean(mPhase(isbip)) nanmean(mPhase(ispos)) nanmean(mPhase(ispyr)) nanmean(mPhase(isINT))],...
    [calc_sem(mPhase(isbip)) calc_sem(mPhase(ispos)) calc_sem(mPhase(ispyr)) calc_sem(mPhase(isINT))], colors_NPB(3,:));% colors_NPB(2,:) colors_IP_light(2,:) colors_IP_light(1,:)] );
set( gca, 'tickdir', 'out', 'box', 'off' )
title('all units mean phase spontaneous')
ylabel('mean phase');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})
text(0.5, 4.3 ,sprintf('%d',sum(struct_final.LS.sumbpi)));
text(1.5, 4.3 ,sprintf('%d',sum(struct_final.LS.sumpos)));
text(2.5, 4.3 ,sprintf('%d',sum(struct_final.LS.sumpyr)));
text(3.5, 4.3 ,sprintf('%d',sum(struct_final.LS.sumint)));

subplot(1,2,2)
barwerror([1 2 3 4],[nanmean(mIPhase(isbip)) nanmean(mIPhase(ispos)) nanmean(mIPhase(ispyr)) nanmean(mIPhase(isINT))],...
    [calc_sem(mIPhase(isbip)) calc_sem(mIPhase(ispos)) calc_sem(mIPhase(ispyr)) calc_sem(mIPhase(isINT))], colors_NPB(3,:));
set( gca, 'tickdir', 'out', 'box', 'off' )
title('all units mean phase Induced')
ylabel('mean phase');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})
ylim([0 4.5])
text(0.5, 4.3 ,sprintf('%d',sum(struct_final.LI.sumbpi)));
text(1.5, 4.3 ,sprintf('%d',sum(struct_final.LI.sumpos)));
text(2.5, 4.3 ,sprintf('%d',sum(struct_final.LI.sumpyr)));
text(3.5, 4.3 ,sprintf('%d',sum(struct_final.LI.sumint)));

figname = '/shirly1/amir1/colors/mean_phase_06apr22';
save_print(figname);

% pooled R
figure;
barwerror([1 2 3 4],[nanmean(rPhase(isbip)) nanmean(rPhase(ispos)) nanmean(rPhase(ispyr)) nanmean(rPhase(isINT))],...
    [calc_sem(rPhase(isbip)) calc_sem(rPhase(ispos)) calc_sem(rPhase(ispyr)) calc_sem(rPhase(isINT))], colors_NPB(3,:));% colors_NPB(2,:) colors_IP_light(2,:) colors_IP_light(1,:)] );
set( gca, 'tickdir', 'out', 'box', 'off' )
title('all units R spontaneous')
ylabel('R');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})
ylim([0 0.5])
text(1, 0.46 ,sprintf('%d',sum(struct_final.LS.sumbpi)));
text(2, 0.46 ,sprintf('%d',sum(struct_final.LS.sumpos)));
text(3, 0.46 ,sprintf('%d',sum(struct_final.LS.sumpyr)));
text(4, 0.46 ,sprintf('%d',sum(struct_final.LS.sumint)));

figname = '/shirly1/amir1/colors/pooled_R';
save_print(figname);


%Can you compute the fraction of units (PYR, INT, biphasic, punits)
%that exhibit phase locking - as a function of depth, and pooled?
%as function of depth:
fig1 = figure;
subplot(4,2,1)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumbpi./(binned_struct_US.sumbpi+binned_struct_LS.sumbpi),[],colors_NPB(3,:), [], colors_NPB(3,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,2)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumbpi./(binned_struct_UI.sumbpi+binned_struct_LI.sumbpi), [], colors_NPB(3,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all induced')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,3)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumpos./(binned_struct_US.sumpos+binned_struct_LS.sumpos), [], colors_NPB(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,4)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumpos./(binned_struct_UI.sumpos+binned_struct_LI.sumpos), [], colors_NPB(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all induced')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,5)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumPYR./(binned_struct_US.sumPYR+binned_struct_LS.sumPYR), [], colors_IP_light(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PYRs locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,6)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumPYR./(binned_struct_UI.sumPYR+binned_struct_LI.sumPYR), [], colors_IP_light(2,:));
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('PYRs locked/all induced');
ylabel('fraction');
xlabel('depth');
ylim ([0 1]);

subplot(4,2,7)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumINT./(binned_struct_US.sumINT+binned_struct_LS.sumINT), [], colors_IP_light(1,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('INT locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,8)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumINT./(binned_struct_UI.sumINT+binned_struct_LI.sumINT), [], colors_IP_light(1,:));
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('INT locked/all induced');
xlabel('depth');
ylabel('fraction');
ylim ([0 1]);

figname = '/shirly1/amir1/colors/fracUnit_bars';
save_print(figname);

%pooled
figure;
subplot(2,2,1)
barwerror([1 2 3 4],[sum(binned_struct_LS.sumbpi) sum(binned_struct_LS.sumpos)...
    sum(binned_struct_LS.sumPYR)  sum(binned_struct_LS.sumINT)], [], colors_NPB(3,:));
hold on
h = barwerror([1 2 3 4],[sum(binned_struct_US.sumbpi)+sum(binned_struct_LS.sumbpi) ...
sum(binned_struct_LS.sumpos)+sum(binned_struct_US.sumpos)... 
sum(binned_struct_LS.sumPYR)+sum(binned_struct_US.sumPYR) ...
 sum(binned_struct_LS.sumINT)+sum(binned_struct_US.sumINT)]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('all units locked/all spontaneous');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})

subplot(2,2,2)
barwerror([1 2 3 4],[sum(binned_struct_LI.sumbpi) sum(binned_struct_LI.sumpos)...
    sum(binned_struct_LI.sumPYR)  sum(binned_struct_LI.sumINT)], [], colors_NPB(3,:));
hold on
h = barwerror([1 2 3 4],[sum(binned_struct_UI.sumbpi)+sum(binned_struct_LI.sumbpi) ...
sum(binned_struct_LI.sumpos)+sum(binned_struct_UI.sumpos)... 
sum(binned_struct_LI.sumPYR)+sum(binned_struct_UI.sumPYR) ...
 sum(binned_struct_LI.sumINT)+sum(binned_struct_UI.sumINT)]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('all units locked/all Induced');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})

subplot(2,2,3)
barwerror([1 2 3 4],[sum(binned_struct_LS.sumbpi)./sum((binned_struct_US.sumbpi+binned_struct_LS.sumbpi))...
    sum(binned_struct_LS.sumpos)./sum((binned_struct_US.sumpos+binned_struct_LS.sumpos))...
    sum(binned_struct_LS.sumPYR)./sum((binned_struct_US.sumPYR+binned_struct_LS.sumPYR))...
    sum(binned_struct_LS.sumINT)./sum((binned_struct_US.sumINT+binned_struct_LS.sumINT)) ]...
    , [bino_se_norm_bpi,bino_se_norm_pos,bino_se_norm_pyr,bino_se_norm_int], colors_NPB(3,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('fraction of all units locked/all spontaneous');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})
ylim ([0 1]);

subplot(2,2,4)
barwerror([1 2 3 4],[sum(binned_struct_LI.sumbpi)./sum((binned_struct_UI.sumbpi+binned_struct_LI.sumbpi))...
    sum(binned_struct_LI.sumpos)./sum((binned_struct_UI.sumpos+binned_struct_LI.sumpos))...
    sum(binned_struct_LI.sumPYR)./sum((binned_struct_UI.sumPYR+binned_struct_LI.sumPYR))...
    sum(binned_struct_LI.sumINT)./sum((binned_struct_UI.sumINT+binned_struct_LI.sumINT))]...
    , [bino_se_norm_bpiI,bino_se_norm_posI,bino_se_norm_pyrI,bino_se_norm_intI], colors_NPB(3,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('fraction of all units locked/all induced');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})
ylim ([0 1]);

figname = '/shirly1/amir1/colors/numUnit_pooled_bars_06apr22';
save_print(figname);


% Can you add histograms of number of cells (empty bars) and phase locked cells (full bars)?

%as function of depth:
fig1 = figure;
subplot(4,2,1)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumbpi, [], colors_NPB(3,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumbpi+binned_struct_LS.sumbpi);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI: %d/%d locked/all spontaneous', binned_struct_LS.sumbpi, binned_struct_US.sumbpi+binned_struct_LS.sumbpi)
ylabel('count');
xlabel('depth')

subplot(4,2,2)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumbpi, [], colors_NPB(3,:))
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumbpi+binned_struct_LI.sumbpi);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI: %d/%d locked/all induced', binned_struct_LI.sumbpi, binned_struct_UI.sumbpi+binned_struct_LI.sumbpi)
ylabel('count');
xlabel('depth')

subplot(4,2,3)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumpos, [], colors_NPB(2,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumpos+binned_struct_LS.sumpos);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('Pos: %d/%d locked/all spontaneous', binned_struct_LS.sumpos, binned_struct_US.sumpos+binned_struct_LS.sumpos)
ylabel('count');
xlabel('depth')

subplot(4,2,4)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumpos, [], colors_NPB(2,:))
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumpos+binned_struct_LI.sumpos);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('Pos: %d/%d locked/all induced', binned_struct_LI.sumpos, binned_struct_UI.sumpos+binned_struct_LI.sumpos)
ylabel('count');
xlabel('depth')

subplot(4,2,5)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumPYR, [], colors_IP_light(2,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumPYR+binned_struct_LS.sumPYR);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PYR: %d/%d locked/all spontaneous', binned_struct_LS.sumPYR, binned_struct_US.sumPYR+binned_struct_LS.sumPYR)
ylabel('count');
xlabel('depth')

subplot(4,2,6)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumPYR, [], colors_IP_light(2,:));
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumPYR+binned_struct_LI.sumPYR);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('PYR: %d/%d locked/all induced', binned_struct_LI.sumPYR, binned_struct_UI.sumPYR+binned_struct_LI.sumPYR);
ylabel('count');
xlabel('depth');

subplot(4,2,7)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumINT, [], colors_IP_light(1,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumINT+binned_struct_LS.sumINT);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('INT: %d/%d locked/all spontaneous', binned_struct_LS.sumINT, binned_struct_US.sumINT+binned_struct_LS.sumINT)
ylabel('count');
xlabel('depth')

subplot(4,2,8)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumINT, [], colors_IP_light(1,:));
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumINT+binned_struct_LI.sumINT);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('INT: %d/%d locked/all induced', binned_struct_LI.sumINT, binned_struct_UI.sumINT+binned_struct_LI.sumINT);
xlabel('depth');
ylabel('count');

figname = '/shirly1/amir1/colors/numUnit_bars';
save_print(figname);

%pooled
figure;
subplot(1,2,1)
barwerror([1 2 3 4],[sum(binned_struct_LS.sumbpi) sum(binned_struct_LS.sumpos)...
    sum(binned_struct_LS.sumPYR)  sum(binned_struct_LS.sumINT)], [], colors_NPB(3,:));
hold on
h = barwerror([1 2 3 4],[sum(binned_struct_US.sumbpi)+sum(binned_struct_LS.sumbpi) ...
sum(binned_struct_LS.sumpos)+sum(binned_struct_US.sumpos)... 
sum(binned_struct_LS.sumPYR)+sum(binned_struct_US.sumPYR) ...
 sum(binned_struct_LS.sumINT)+sum(binned_struct_US.sumINT)]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('all units locked/all spontaneous');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})

subplot(1,2,2)
barwerror([1 2 3 4],[sum(binned_struct_LI.sumbpi) sum(binned_struct_LI.sumpos)...
    sum(binned_struct_LI.sumPYR)  sum(binned_struct_LI.sumINT)], [], colors_NPB(3,:));
hold on
h = barwerror([1 2 3 4],[sum(binned_struct_UI.sumbpi)+sum(binned_struct_LI.sumbpi) ...
sum(binned_struct_LI.sumpos)+sum(binned_struct_UI.sumpos)... 
sum(binned_struct_LI.sumPYR)+sum(binned_struct_UI.sumPYR) ...
 sum(binned_struct_LI.sumINT)+sum(binned_struct_UI.sumINT)]);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('all units locked/all Induced');
ylabel('count');
set(gca, 'XTickLabel',{'BPI' 'POS' 'PYR' 'INT'})

figname = '/shirly1/amir1/colors/numUnit_pooled_bars';
save_print(figname);


% combined fracUnit_bars and numUnit_bars
%spontaneous

fig1 = figure;
subplot(4,2,1)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumbpi, [], colors_NPB(3,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumbpi+binned_struct_LS.sumbpi);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all spontaneous')
ylabel('count');
xlabel('depth')

subplot(4,2,2)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumbpi./(binned_struct_US.sumbpi+binned_struct_LS.sumbpi),[],colors_NPB(3,:), [], colors_NPB(3,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);


subplot(4,2,3)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumpos, [], colors_NPB(2,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumpos+binned_struct_LS.sumpos);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all spontaneous')
ylabel('count');
xlabel('depth')

subplot(4,2,4)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumpos./(binned_struct_US.sumpos+binned_struct_LS.sumpos), [], colors_NPB(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,5)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumPYR, [], colors_IP_light(2,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumPYR+binned_struct_LS.sumPYR);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PYRs locked/all spontaneous')
ylabel('count');
xlabel('depth')

subplot(4,2,6)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumPYR./(binned_struct_US.sumPYR+binned_struct_LS.sumPYR), [], colors_IP_light(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PYRs locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);


subplot(4,2,7)
barwerror(binned_struct_LS.depth,binned_struct_LS.sumINT, [], colors_IP_light(1,:))
hold on
h= barwerror(binned_struct_US.depth,binned_struct_US.sumINT+binned_struct_LS.sumINT);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('INT locked/all spontaneous')
ylabel('count');
xlabel('depth')

subplot(4,2,8)
h= barwerror(binned_struct_US.depth,binned_struct_LS.sumINT./(binned_struct_US.sumINT+binned_struct_LS.sumINT), [], colors_IP_light(1,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('INT locked/all spontaneous')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

figname = '/shirly1/amir1/colors/combinedSps_bars';
save_print(figname);

% induced

fig1 = figure;

subplot(4,2,1)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumbpi, [], colors_NPB(3,:))
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumbpi+binned_struct_LI.sumbpi);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all induced')
ylabel('count');
xlabel('depth')

subplot(4,2,2)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumbpi./(binned_struct_UI.sumbpi+binned_struct_LI.sumbpi), [], colors_NPB(3,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('BPI locked/all induced')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,3)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumpos, [], colors_NPB(2,:))
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumpos+binned_struct_LI.sumpos);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all induced')
ylabel('count');
xlabel('depth')

subplot(4,2,4)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumpos./(binned_struct_UI.sumpos+binned_struct_LI.sumpos), [], colors_NPB(2,:));
set(h,'EdgeColor','k')
set( gca, 'tickdir', 'out', 'box', 'off' )
title('PUNITs locked/all induced')
ylabel('fraction');
xlabel('depth')
ylim ([0 1]);

subplot(4,2,5)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumPYR, [], colors_IP_light(2,:));
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumPYR+binned_struct_LI.sumPYR);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('PYRs locked/all induced');
ylabel('count');
xlabel('depth');

subplot(4,2,6)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumPYR./(binned_struct_UI.sumPYR+binned_struct_LI.sumPYR), [], colors_IP_light(2,:));
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('PYRs locked/all induced');
ylabel('fraction');
xlabel('depth');
ylim ([0 1]);

subplot(4,2,7)
barwerror(binned_struct_LI.depth,binned_struct_LI.sumINT, [], colors_IP_light(1,:));
hold on
h= barwerror(binned_struct_UI.depth,binned_struct_UI.sumINT+binned_struct_LI.sumINT);
set(h,'FaceColor','none');
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('INT locked/all induced');
xlabel('depth');
ylabel('count');

subplot(4,2,8)
h= barwerror(binned_struct_US.depth,binned_struct_LI.sumINT./(binned_struct_UI.sumINT+binned_struct_LI.sumINT), [], colors_IP_light(1,:));
set(h,'EdgeColor','k');
set( gca, 'tickdir', 'out', 'box', 'off' );
title('INT locked/all induced');
xlabel('depth');
ylabel('fraction');
ylim ([0 1]);

figname = '/shirly1/amir1/colors/combinedInd_bars';
save_print(figname);

return