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

function fig_handel = make_PUNIT_figures( figs, dataname, path, varargin )



%--------------------------------------------------------------------%
% constants
%--------------------------------------------------------------------%
alpha           = 0.05;
pre_slopes      = [-0.005, -0.1];
re_slopes       = [0.04, 0.25];
color_gro       = [1 2 3];

%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%

nfigs        = 1;

if nargin < 1 || isempty(figs)
    figs     = 1:nfigs;
end

if nargin < 2 || isempty(dataname)
    dataname = 'TPR_TPP_data_17-Aug-2020_extremum.mat';
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
%sD           = L.sD;
fD           = L.fD;

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

% all pyr / int
Eushankclu            = [uD.session_idx, uD.shankclu];
[~, i_un]             = unique(Eushankclu, 'rows');
uDatype               = logical( uD.shankclu(:,3) );
uDatype_uni           = uDatype(i_un);



% Precession, recession and lock in fields
fTPP_cond             = ( abs ( fD.pre(:,1)  ) >= abs( pre_slopes(1) ) ) & ( abs ( fD.pre(:,1)  ) <= abs( pre_slopes(2) ) );
fTPP_logi             = ( fD.pre(:,3) <= alpha ) & fTPP_cond;
fTPR_cond             = ( abs ( fD.re(:,1)  ) >= re_slopes(1) ) & ( abs ( fD.re(:,1)  ) <= re_slopes(2) );

% punits
punits_fields                  = fD.extremum > 0;
fD.shankclu(punits_fields,3)   = 2;
pyr_fields                  = fD.shankclu(:,3) == 1;
int_fields                  = fD.shankclu(:,3) == 0;
pun_fields                  = fD.shankclu(:,3) == 2;

punits_fields               = uD.extremum > 0;
uD.shankclu(punits_fields,3)= 2;
unit_pyr                    = uD.shankclu(:,3) == 1;
unit_int                    = uD.shankclu(:,3) == 0;
unit_pun                    = uD.shankclu(:,3) == 2;
Efshankclu                  = [fD.session_idx, fD.shankclu];


%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%
for     fig_num = figs
    
    switch fig_num
        
        % rate of fields with TPP out of all fields
        case 1
            fig_handel( 1 )                = figure;
            xstr                           = {'PYR', 'INT', 'Punits'};
            id                             = [ pyr_fields,  int_fields, pun_fields];
            tstr                           = 'Fields with theta phase precession';
            make_TPR_TPP_figures_bars( fTPP_logi, id, xstr, tstr, alpha, 'color_gro',color_gro  )
            %axis tight

        % field size
        case 2   
            fig_handel( 2 )         = figure;
            tstr                    = 'field size';
            param_type              = 'linear';
            binsize                 = 5;
            id                      = [ pyr_fields,  int_fields, pun_fields];
            prm                     = fD.stat(:,7);
            xtext                   = 'Field size [cm]';
            ytext                   = 'number of fields';
            grps                    = {'PYR','INT','PUNIT'};
            make_TPR_TPP_histograms2( prm, id, tstr, param_type,'binsize',binsize,'xtext',xtext,'ytext', ytext, 'color_gro',color_gro, 'grps', grps)
            
        % Gain (all fields)
        case 3
            fig_handel(3)           = figure;
            tstr                    = 'Gain (all fields)';
            grps                    = {'PYR','INT','PUNIT'};
            param_type              = 'log';
            id                      = [ pyr_fields,  int_fields, pun_fields];
            prm                     = fD.stat(:,5);
            xtext                   = 'within field FR gain';
            ytext                   = 'ratio of fields';
            binsize                 = 0.2;
            ylimits                 = [ 0 0.35];
            make_TPR_TPP_histograms2( prm, id, tstr, param_type,'binsize',binsize,'xtext',xtext,'ytext', ytext, 'color_gro',color_gro, 'grps', grps, 'ratio',1, 'ylimits',ylimits)
            
        % Phase precession slope (sig)
        case 4
            fig_handel( 4 )         = figure;
            tstr                    = 'Phase precession slope (sig)';
            param_type              = 'linear';
            grps                    = {'PYR','INT','PUNIT'};
            binSize                 = 0.005;
            edges                   = pre_slopes(2) : binSize : pre_slopes(1);
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
            prm                     = fD.pre(:,1);
            ytext                   = 'Number of fields';
            xtext                   = 'Slope';
            make_TPR_TPP_histograms2( prm, id, tstr, param_type,'edges',edges,'xtext',xtext,'ytext', ytext, 'color_gro',color_gro, 'grps', grps)
            
        % Phase precession fit (sig)
        case 5
            fig_handel( 5 )         = figure;
            tstr                    = 'Phase precession fit (sig)';
            param_type              = 'linear';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
            prm                     = fD.pre(:,2);
            ytext                   = 'Number of fields';
            xtext                   = 'Fit (R)';
            binsize                 = 0.03;
            make_TPR_TPP_histograms2( prm, id, tstr, param_type,'binsize',binsize,'xtext',xtext,'ytext', ytext, 'color_gro',color_gro, 'grps', grps)
            
        % TPP repeats per field (sig)
        case 6
            fig_handel( 6 )         = figure;
            tstr                    = 'TPP repeats per field (sig)';
            param_type              = 'lin';
            binsize                 = 0.5;
            prm                     = abs (fD.pre(:,7) ) ;
            xtext                   = '#cycles';
            ytext                   = '#fields';
            grps                    = {'PYR','INT','PUNIT'};
            id                      = [ pyr_fields & fTPP_logi,  int_fields & fTPP_logi, pun_fields & fTPP_logi];
            make_TPR_TPP_histograms2( prm, id, tstr, param_type,'binsize',binsize,'xtext',xtext,'ytext', ytext, 'color_gro',color_gro, 'grps', grps)
            
        % Units with one field or more
        case 7
            % get data
            [~, i_un]               = unique(Efshankclu, 'rows');
            pyr_units_fields        = sum(fD.shankclu(i_un,3) == 1);
            int_units_fields        = sum(fD.shankclu(i_un,3) == 0);
            pun_units_fields        = sum(fD.shankclu(i_un,3) == 2);
            [~, i_un]               = unique(Eushankclu, 'rows');
            pyr                     = sum(uD.shankclu(i_un,3) == 1);
            int                     = sum(uD.shankclu(i_un,3) == 0);
            pun                     = sum(uD.shankclu(i_un,3) == 2);
           
            %plot
            fig_handel( 7)          = figure;
            id                      = [pyr_units_fields, int_units_fields, pun_units_fields];
            ngroup                  = [pyr, int, pun];
            tstr                    = 'Units with one field or more';
            grps                    = {'PYR','INT','PUNIT'};
            make_TPR_TPP_figures_bars( id, ngroup, grps, tstr, alpha,'color_gro',color_gro )
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
            ngroup                  = [sum(unit_pyr), sum(unit_int), sum(unit_pun)];
            tstr                    = 'Units with exectly one field';
            make_TPR_TPP_figures_bars( id, ngroup, {'PYR','INT', 'PUNIT'}, tstr, alpha,'color_gro',color_gro )               
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
            
            %get data
            [my_unit, iu]           = unique(Eushankclu,'rows');
            unit_nfields            = NaN(size(my_unit,1),1);
            for i = 1 : size(my_unit, 1)
                unit_i              = my_unit(i,1:3);
                rev_id              = ( ( Efshankclu(:,1) == unit_i(1) ) & ( Efshankclu(:,2) == unit_i(2) ) & ( Efshankclu(:,3) == unit_i(3) ) );
                unit_nfields(i)     = sum(rev_id);
            end
            
            %plot
            fig_handel(10)          = figure;
            id                      = [ unit_pyr(iu), unit_int(iu), unit_pun(iu)];
            prm                     = unit_nfields;
            xtext                   = 'Number of fields';
            ytext                   = 'Number of units';
            tstr                    = 'Number of fields per unit';
            param_type              = 'linear';
            binc                    = 0:3;
            grps                 = {'PYR', 'INT', 'PUN'};
            make_TPR_TPP_histograms2( prm, id, tstr, param_type, 'binc', binc, 'color_gro',color_gro,'xtext',xtext,'ytext',ytext,'mass_cent',1,'grps',grps)

            
        % spatial information per cell*direction [bits/s]
        case 11
            fig_handel(11)          = figure;
            tstr                    = 'spatial information [bits/s]';
            prm                     = uD.info(:,1);
            id                      = [ unit_pyr, unit_int, unit_pun];
            param_type              = 'log';
            xtext                   = 'SI [bits/s]';
            ytext                   = 'number of units * directions';
            binsize                 = 0.2;
            edges                   = -11.1 : binsize : 5.1;
            grps                 = {'PYR', 'INT', 'PUN'};
            make_TPR_TPP_histograms2( prm, id, tstr, param_type, 'color_gro',color_gro,'xtext',xtext,'ytext',ytext,'edges',edges,'grps',grps);           
            
        % all cells with phase lock (unit phase lock)  
        case 12
            fig_handel( 12 )       = figure;
            id                     = [ unit_pyr, unit_int, unit_pun];
            vals                   = uD.lock(:,4);
            logi_prm1              = vals <= alpha ;
            tstr                   = 'Phase lock (units*directions)';
            xstr                   = {'PYR','INT','PUNIT'};
            make_TPR_TPP_figures_bars( logi_prm1 , id, xstr, tstr, alpha, 'color_gro', color_gro )
            
        % phs lock angle (unit phase lock)  
        case 13
            fig_handel( 13 )        = figure;
            tstr                    = 'Phase lock angle (sig)';
            id                      = [ unit_pyr, unit_int, unit_pun];
            param_type              = 'circ';
            prm                     = uD.lock(:,1);
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
            ngroup              = [sum(unit_pyr), sum(unit_int), sum(unit_pun)];
            tstr                = 'Units with exectly one field (per direction)';
            make_TPR_TPP_figures_bars( id, ngroup, {'PYR','INT', 'PUNIT'}, tstr, alpha,'color_gro',color_gro )
            
            
    end
end

fig_handel                          = fig_handel(figs);

return;
%EOF
