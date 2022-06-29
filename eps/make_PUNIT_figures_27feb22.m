% calls: make_TPR_TPP_figures_bars, make_TPR_TPP_histograms2, make_TPR_TPP_figures_histograms

% cases:
% (1) rate of fields with TPP out of all fields
% (2) field size
% (3) Gain (all fields)
% (4) Phase precession slope (sig)
% (5) Phase precession fit (sig)
% (6) TPP repeats per field (sig)
% (7) Units with one field or more
% (8) Units with exectly one field
% (9) rate of units with TPP out of all units with fields
% (10) histogram of number of fields per cell
% (11) spatial information per cell*direction [bits/s]
% (12) all cells with phase lock (unit phase lock)
% (13) phs lock angle (unit phase lock)
% (14) Units with exectly one field (per direction)

% revisions
% 11-oct-20 (1) adapted for mac as well
%           (2) added CM test for phase locking
% 07-jan-21 (1) added optional argument sst and dilution of fD and uD
% 27-Feb-22 (1) updated to work with the latest dataset


function make_PUNIT_figures( figs, dataname, path, varargin )



%--------------------------------------------------------------------%
% constants
%--------------------------------------------------------------------%
alpha           = 0.05;
pre_slopes      = [-0.1, -0.01];
color_gro       = [1 2 3];
minf_spk        = 30;
max_field_size  = 100; %[cm]

%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%

nfigs        = 1;

if nargin < 1 || isempty(figs)
    figs     = 1:nfigs;
end

if nargin < 2 || isempty(dataname)
    dataname = 'TPR_23-Aug-2021_Z6_extremum.mat';
end


if nargin < 3 || isempty(path)
    if ismac 
        path        = '/Users/eranstark/Documents/da/punits/pf/';
    elseif isunix   % shirly (nanna) Linux
        path     = '/media/shirly/Data/_Shirly/Lab/eps/';
    elseif ispc     % shirly (nanna) Windows
        path     = 'D:\_Shirly\Lab\matlab\eps\';
    end
end

[ sst ...
    ]                       = ParseArgPairs(...
    { 'sst' ...
    }...
    , { [] ...
    }...
    , varargin{ : } );

filename     = sprintf('%s/%s', path, dataname);
load( filename, 'L' );
uD           = L.uD;
fD           = L.fD;
tD           = L.tD;

%--------------------------------------------------------------------%
% dilute sessions and units that appear in L and do not appear in sst
%--------------------------------------------------------------------%
if ~isempty( sst )
    
    fprintf( 1, '%s: diluting ...', mfilename )
    
    filebase1               = sst.filebase;
    shankclu1               = sst.shankclu;

    filebase2               = uD.session_tag;
    shankclu2               = uD.shankclu;
    
    filebase3               = fD.session_tag;
    shankclu3               = fD.shankclu;
    
    % dilute the uD structure
    u1                      = unique( filebase1 );
    u2                      = unique( filebase2 );
    u                       = unique( [ u1; u2 ] );
    [ ~, f1 ]               = ismember( filebase1, u );
    [ ~, f2 ]               = ismember( filebase2, u );
    s1                      = [ f1 shankclu1 ];
    s2                      = [ f2 shankclu2( :, 1 : 2 ) ];
    ns2                     = size( s2, 1 );
    isvalid                 = false( ns2, 1 );
    for i                   = 1 : ns2
        isvalid( i )        = ismember( s2( i, : ), s1, 'rows' );
    end
    uD                      = struct_select( uD, isvalid );
    
    % dilute the fD structure
    u1                      = unique( filebase1 );
    u2                      = unique( filebase3 );
    u                       = unique( [ u1; u2 ] );
    [ ~, f1 ]               = ismember( filebase1, u );
    [ ~, f2 ]               = ismember( filebase3, u );
    s1                      = [ f1 shankclu1 ];
    s2                      = [ f2 shankclu3( :, 1 : 2 ) ];
    ns2                     = size( s2, 1 );
    isvalid                 = false( ns2, 1 );
    for i                   = 1 : ns2
        isvalid( i )        = ismember( s2( i, : ), s1, 'rows' );
    end
    fD                      = struct_select( fD, isvalid );
    
    fprintf( 1, 'done\n' )
    
end

%--------------------------------------------------------------------%
% get group indices
%--------------------------------------------------------------------%

min_spk_cond    = fD.nspk >= minf_spk;
maxField_cond   = fD.stat(:,7) <= max_field_size;
stable          = fD.field_stab(:,2) <= alpha;
min_conds       = min_spk_cond & maxField_cond & stable;


% all pyr / int
Efshankclu                  = [fD.session_idx, fD.shankclu];
Eushankclu            = [uD.session_idx, uD.shankclu];
[~, i_un]             = unique(Eushankclu, 'rows');
uDatype               = logical( uD.shankclu(:,3) );
uDatype_uni           = uDatype(i_un);

% create a t tag for whether each unit has a field
Etshankclu             = [tD.session_idx, tD.shankclu];
nt                     = size(tD.shankclu,1);
tField                 = false(nt,1);
tnFields               = zeros(nt,1);
for i = 1 : nt
    tidx               = Etshankclu(i,1:3);
    idx2               = ( tidx(1) == Efshankclu(:,1) ) & ( tidx(2) == Efshankclu(:,2) ) & ( tidx(3) == Efshankclu(:,3) );
    nspk_field         = min_conds(idx2);
    if any(idx2) && any(nspk_field)
        tField(i,1)    = true;
    end
    tnFields(i,1)      = sum(idx2 & min_conds);
end

% Precession, recession and lock in fields
fTPP_cond             = ( fD.pre(:,1) >= pre_slopes(1) ) & ( fD.pre(:,1) <= pre_slopes(2) );
fTPP_logi             = ( fD.pre(:,3) <= alpha ) & fTPP_cond;

% punits
punits_fields                  = fD.extremum > 0;
fD.shankclu(punits_fields,3)   = 2;
pyr_fields                  = fD.shankclu(:,3) == 1 & min_conds;
int_fields                  = fD.shankclu(:,3) == 0 & min_conds;
pun_fields                  = fD.shankclu(:,3) == 2 & min_conds;
Efshankclu                  = [fD.session_idx, fD.shankclu];

active               = tD.stability(:,2) <= alpha & tD.active2;
punits_tunits               = tD.extremum > 0;
tD.shankclu(punits_tunits,3)= 2;
t_pyr                       = tD.shankclu(:,3) == 1 & active;
t_int                       = tD.shankclu(:,3) == 0 & active;
t_pun                       = tD.shankclu(:,3) == 2 & active;



%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%
for     fig_num = figs
    
    switch fig_num
        
        % rate of fields with TPP out of all fields
        case 1
            xstr                           = {'PYR', 'INT', 'Punits'};
            id                             = [ pyr_fields,  int_fields, pun_fields];
            tstr                           = 'Fields with theta phase precession';
            make_TPR_TPP_figures_bars( fTPP_logi, id, xstr, tstr, alpha, 'color_gro',color_gro  )
            
            %axis tight

        % field size
        case 2   
            tstr                    = 'field size';
            id                      = [ pyr_fields,  int_fields, pun_fields];
            prm                     = fD.stat(:,7);
            xtext                   = 'Field size [cm]';
            ytext                   = 'number of fields';
            grps                    = {'PYR','INT','PUNIT'};
            
            make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext,'logscale',0,'logscale2',0)
            xlim([15 100]);
            
        % Gain (all fields)
        case 3
            tstr                    = 'Gain (all fields)';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields,  int_fields, pun_fields];
            prm                     = fD.stat(:,5);
            xtext                   = 'within field FR gain';
            ytext                   = 'ratio of fields';
            ticks                   = [ 0.02 0.2 2 20 200];
            xlimits                 = [0.2 500];
            
            prm( prm > 500 )        = 500;
            
            make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext,'logscale',1,'logscale2',0, ...
                'ticks',ticks,'xlimits',xlimits)

            
            
        % Phase precession slope (sig)
        case 4
            fig_handel( 4 )         = figure;
            tstr                    = 'Phase precession slope (sig)';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
            prm                     = fD.pre(:,1);
            ytext                   = 'Number of fields';
            xtext                   = 'Slope';
             make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext,'logscale',0,'logscale2',0)

            
        % Phase precession fit (sig)
        case 5
            fig_handel( 5 )         = figure;
            tstr                    = 'Phase precession fit (sig)';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
            prm                     = fD.pre(:,2);
            ytext                   = 'Number of fields';
            xtext                   = 'Fit (R)';
            xlimits                 = [0.06 0.6];
             make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext, ...
                'logscale',0,'logscale2',0, 'xlimits', xlimits)
            
        % TPP repeats per field (sig)
        case 6
            fig_handel( 6 )         = figure;
            tstr                    = 'TPP repeats per field (sig)';
            prm                     = abs (fD.pre(:,7) ) ;
            xtext                   = '#cycles';
            ytext                   = '#fields';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
             make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext, ...
                'logscale',0,'logscale2',0, 'xlimits', [])
            
        % Units with one field or more
        case 7
            id                      = [pyr_fields, int_fields, pun_fields];
            ngroup                  = [t_pyr, t_int, t_pun];
            tstr                    = 'Units with one field or more';
            grps                    = {'PYR','INT','PUNIT'};
            make_TPR_TPP_figures_bars( sum(id), sum(ngroup), grps, tstr, alpha,'color_gro',color_gro )
            %axis tight
            
        % Units with exectly one field    
        case 8
            fig_handel( 8 )         = figure;          
            % to get cells with one exectly
            Efshankclu              = [fD.session_idx, fD.shankclu];
            uEfshankclu             = unique(Efshankclu, 'rows');
            nfields                 = size(uEfshankclu,1);
            one_field               = false( size(Efshankclu,1),1 );
            
            % get data
            for i = 1 : nfields
                u                   = uEfshankclu(i,1:3);
                rev_id              = ( ( Efshankclu(:,1) == u(1) ) & ( Efshankclu(:,2) == u(2) ) & ( Efshankclu(:,3) == u(3) ) );
                if sum(rev_id) == 1
                    one_field(i)    = true(1,1);
                end
            end
            
            % number of cells with fields
            id                      = [pyr_fields & one_field, int_fields & one_field, pun_fields & one_field];
            ngroup                  = [sum(t_pyr), sum(t_int), sum(t_pun)];
            tstr                    = 'Units with exectly one field';
            ylimits                 = [0 0.3];
            make_TPR_TPP_figures_bars( id, ngroup, {'PYR','INT', 'PUNIT'}, tstr, alpha,'color_gro',color_gro,'ylimits',ylimits )               
            %axis tight
            
        % rate of units with TPP out of all units with fields
        case 9
            
            % all pyr with TPP
            pyr_idx                 = logical( pyr_fields & fTPP_logi );
            pyr_cell_id             = [fD.session_idx(pyr_idx), fD.shankclu(pyr_idx,:)];
            [a]                     = unique(pyr_cell_id,'rows');
            pyr_w_TPP               = size(a,1);
            
            % all int with TPP
            int_idx                 = logical( int_fields & fTPP_logi );
            int_cell_id             = [fD.session_idx(int_idx), fD.shankclu(int_idx,:)];
            [a]                     = unique(int_cell_id,'rows');
            int_w_TPP               = size(a,1);
            
                        % all pun with TPP
            int_idx                 = logical( pun_fields & fTPP_logi );
            int_cell_id             = [fD.session_idx(int_idx), fD.shankclu(int_idx,:)];
            [a]                     = unique(int_cell_id,'rows');
            pun_w_TPP               = size(a,1);
            
            %plot
            fig_handel( 9 )         = figure;
            tstr                    = 'Units with TPP out of all units with fields';
            grps                    = {'PYR','INT','PUNIT'};
            make_TPR_TPP_figures_bars( [pyr_w_TPP,int_w_TPP, pun_w_TPP] , [sum(pyr_fields), sum(int_fields), sum(pun_fields)], grps, tstr, alpha, 'color_gro', color_gro  )
            %axis tight
            
        % histogram of number of fields per cell
        case 10
            %plot
            id                      = [ t_pyr, t_int, t_pun];
            prm                     = tnFields;
            xtext                   = 'Number of fields';
            ytext                   = 'Number of units';
            tstr                    = 'Number of fields per unit';
            grps                    = {'PYR', 'INT', 'PUN'};

            make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext, ...
                'logscale',0,'logscale2',0, 'com', 1)
            
%         spatial information per cell
        case 11
            tstr                    = 'spatial information [bits/s]';
            prm                     = tD.info(:,1);
            id                      = [ t_pyr, t_int, t_pun];
            xtext                   = 'SI [bits/s]';
            ytext                   = 'number of units * directions';
            grps                 = {'PYR', 'INT', 'PUN'};
            xlimits                 = [0 20];
                                    ticks           = [0.015 0.15 1.5 15];

            make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
                'xtext',xtext, 'color_gro', color_gro,'ytext', ytext, ...
                'logscale',1,'logscale2',0,'xlimits',xlimits,'ticks',ticks)
            
        % all cells with phase lock (unit phase lock)  
        case 12
            fig_handel( 12 )       = figure;
            id                     = [ t_pyr, t_int, t_pun];
            vals                   = tD.prefPhs(:,4);
            logi_prm1              = vals <= alpha ;
            tstr                   = 'Phase lock ';
            xstr                   = {'PYR','INT','PUN'};
            make_TPR_TPP_figures_bars( logi_prm1 , id, xstr, tstr, alpha, 'color_gro', color_gro )
            
        % phs lock angle (unit phase lock)  
        case 13
            fig_handel( 13 )        = figure;
            tstr                    = 'Phase lock angle (sig)';
            id                     = [ t_pyr, t_int, t_pun];
            param_type              = 'circ';
            prm                   = tD.prefPhs(:,1);
            make_TPR_TPP_figures_histograms( prm, id, tstr, param_type, 'color_gro', color_gro )
            
            
            % test for equal medians using circular non-parametric "ANOVA" (CM test)
            vidx                    = ~isnan( prm );
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

            
            
        % Units with exectly one field (per direction)
        case 14
            fig_handel( 14 )        = figure;
            
            % to get cells with one exectly
            Efshankclu           = [fD.session_idx, fD.shankclu, fD.dir];
            uEfshankclu          = unique(Efshankclu, 'rows');
            nfields              = size(uEfshankclu,1);
            one_field            = false( size(Efshankclu,1),1 );
            
            % get data
            for i = 1 : nfields
                u               = uEfshankclu(i,1:4);
                rev_id          = ( ( Efshankclu(:,1) == u(1) ) & ( Efshankclu(:,2) == u(2) ) & ( Efshankclu(:,3) == u(3) ) & ( Efshankclu(:,4) == u(4) ));
                if sum(rev_id) == 1
                    one_field(i) = true(1,1);
                end
            end
            
            % number of cells with fields
            id                   = [pyr_fields & one_field, int_fields & one_field, pun_fields & one_field];
            ngroup              = [sum(t_pyr), sum(t_int), sum(t_pun)];
            tstr                = 'Units with exectly one field (per direction)';
            ylimits             = [0 0.3];
            make_TPR_TPP_figures_bars( id, ngroup, {'PYR','INT', 'PUNIT'}, tstr, alpha,'color_gro',color_gro, 'ylimits', ylimits)
            
            
    end
end


return;
%EOF
