params                  = { 'types', {'SINE','RAMP','PSINE', 'TRIANG'} 'durs', [ 0.005 0.03 ], 'vals', [ 1e-6 1e-3 ], 'franges', [], 'dfranges', [], 'times', [] };
[ ~, tims47, durs47, vals47 ] = get_triggers( filebase, 47, 1, 'eq', 0, params );
[ ~, tims42, durs42, vals42 ] = get_triggers( filebase, 42, 1, 'eq', 0, params );
[ ~, tims43, durs43, vals43 ] = get_triggers( filebase, 43, 1, 'eq', 0, params );
[ ~, tims45, durs45, vals45 ] = get_triggers( filebase, 45, 1, 'eq', 0, params );
[ ~, tims44, durs44, vals44 ] = get_triggers( filebase, 44, 1, 'eq', 0, params );
[ ~, tims46, durs46, vals46 ] = get_triggers( filebase, 46, 1, 'eq', 0, params );
[ ~, tims32, durs32, vals32 ] = get_triggers( filebase, 32, 1, 'eq', 0, params );
[ ~, tims40, durs40, vals40 ] = get_triggers( filebase, 40, 1, 'eq', 0, params );
[ ~, tims34, durs34, vals34 ] = get_triggers( filebase, 34, 1, 'eq', 0, params );

save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals47', 'vals47')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals45', 'vals45')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals43', 'vals43')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals42', 'vals42')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs47', 'durs47')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs45', 'durs45')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs43', 'durs43')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs42', 'durs42')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims42', 'tims42')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims43', 'tims43')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims45', 'tims45')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims47', 'tims47')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims44', 'tims44')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs44', 'durs44')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals44', 'vals44')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims46', 'tims46')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs46', 'durs46')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals46', 'vals46')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims32', 'tims32')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs32', 'durs32')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals32', 'vals32')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/tims40', 'tims40')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/durs40', 'durs40')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/vals40', 'vals40')

load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/tims42.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/tims43.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/tims44.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/tims45.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/tims46.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/vals42.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/vals44.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/vals43.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/vals45.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/vals46.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/durs42.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/durs43.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/durs44.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/durs45.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_11/durs46.mat')



% convert durations to samples, and create range matrices
mat42                           = tims42 + [ zeros( length( durs42 ), 1 ) round( durs42 * 20000 ) - 1 ];
mat43                           = tims43 + [ zeros( length( durs43 ), 1 ) round( durs43 * 20000 ) - 1 ];
mat45                           = tims45 + [ zeros( length( durs45 ), 1 ) round( durs45 * 20000 ) - 1 ];
mat47                           = tims47 + [ zeros( length( durs47 ), 1 ) round( durs47 * 20000 ) - 1 ];
mat44                           = tims44 + [ zeros( length( durs44 ), 1 ) round( durs44 * 20000 ) - 1 ];
mat46                           = tims46 + [ zeros( length( durs46 ), 1 ) round( durs46 * 20000 ) - 1 ];
mat40                           = tims40 + [ zeros( length( durs40 ), 1 ) round( durs40 * 20000 ) - 1 ];
mat32                           = tims32 + [ zeros( length( durs32 ), 1 ) round( durs32 * 20000 ) - 1 ];
mat34                           = tims34 + [ zeros( length( durs34 ), 1 ) round( durs34 * 20000 ) - 1 ];

% combine into cell array at set the rule
mats                            = {mat32, mat42, mat47, mat44, mat34, mat45};
seqmat                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;
% call the routine
[ seqs, smat ]                  = seqdetect( mats, seqmat, 1 );

save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/seqs_2', 'seqs')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_14/smat_2', 'smat')

evtfname = [ filebase '.evt.s02' ]; rc = make_evt_file( evtfname, seqs, 20000 );

make_cs_chip_figures_seqplot (filebase, seqs, smat)

resize = '-bestfit'; pstr = '-dpdf';
resize = ''; pstr = '-dpng';
figname = '/media/shirly/C22865A128659567/mice/mC400/figs/seqs/mC400_14_seq2';
fig = gcf;

figure( fig );
figi = gcf;
figi.Renderer = 'painters';
pause( 0.2 )
print( fig, pstr, figname, resize )

PSTH_stimDurs                   = PSTH_stimDurs/2;

periods = (smat (:,1:2))/20000;
seq2 = seqs(1320:end,:);
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_12/seq1/seq1', 'seq1')
periods = seq1/20000;
PSTH_stimVal = [ 1e-5 6e-4 ];
shanknums                       = [];                   % each shank stimulated by its own local source
    [ stimchans, ~, ~, source ]     = get_stimchans( par );
    idx                             = ismember( source, PSTH_sources );
    stimchans                       = stimchans( idx );
        for j                   = 1 : length( PSTH_ilevels )
        ilev                = PSTH_ilevels{ j };
        s                   = load_spikes( filebase, shanknums, ilev );
            stimType        = 'SINE';
%             if j==1
%                 for i = 1:size(s.shankclu,1)
%                     shanknums                       = s.shankclu(i,1:2);                   % each shank stimulated by its own local source
                    multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums, 'stimTypes', 0, 'valRange',...
                    PSTH_stimVal, 'channels', stimchans, 'uflag', PSTH_uflag, 'multi', PSTH_cmp, 'periods',periods);
                    
%                 end
%             else
%                 shanknums                       = [];
%                 multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums, 'stimTypes', 0, 'valRange',...
%                     PSTH_stimVal, 'channels', stimchans, 'uflag', PSTH_uflag, 'multi', PSTH_cmp, 'periods',periods);
%             end
        end
        
 s = celltypeClassification( filebase,'slevel', slevel,'graphics',1, 'sourceType', 'LED', 'Overwrite',1, 'wavRange',...
     [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');

figname = '/media/shirly/C22865A128659567/mice/mC400/dat/mC400_12/seq1/4_15_4';
figure (9);
figi = gcf;
figi.Renderer='painters';
pause (0.2)
print ('-dpng',figname);