filebase                      = '/media/shirly/C22865A128659567/mice/mC400/dat/mC400_27/mC400_27';
ilev                            = 'B';
s                               = load_spikes( filebase, [], ilev );

% detect seqs
    params                        = { 'types', {'SINE','RAMP','PSINE', 'TRIANG','PULSE'}, 'durs', [ 0.005 0.03 ], 'vals', [ 1e-6 1e-3 ], 'franges', [], 'dfranges', [], 'times', [] };

    [ ~, tims47, durs47, vals47 ] = get_triggers( filebase, 47, 1, 'eq', 0, params );
    [ ~, tims42, durs42, vals42 ] = get_triggers( filebase, 42, 1, 'eq', 0, params );
    [ ~, tims43, durs43, vals43 ] = get_triggers( filebase, 43, 1, 'eq', 0, params );
    [ ~, tims45, durs45, vals45 ] = get_triggers( filebase, 45, 1, 'eq', 0, params );
    [ ~, tims44, durs44, vals44 ] = get_triggers( filebase, 44, 1, 'eq', 0, params );
    [ ~, tims46, durs46, vals46 ] = get_triggers( filebase, 46, 1, 'eq', 0, params );
    [ ~, tims32, durs32, vals32 ] = get_triggers( filebase, 32, 1, 'eq', 0, params );
    [ ~, tims40, durs40, vals40 ] = get_triggers( filebase, 40, 1, 'eq', 0, params );
    [ ~, tims34, durs34, vals34 ] = get_triggers( filebase, 34, 1, 'eq', 0, params );
    [ ~, tims31, durs31, vals31 ] = get_triggers( filebase, 31, 1, 'eq', 0, params );
    [ ~, tims33, durs33, vals33 ] = get_triggers( filebase, 33, 1, 'eq', 0, params );
    [ ~, tims39, durs39, vals39 ] = get_triggers( filebase, 39, 1, 'eq', 0, params );
    [ ~, tims41, durs41, vals41 ] = get_triggers( filebase, 41, 1, 'eq', 0, params );
    [ ~, tims35, durs35, vals35 ] = get_triggers( filebase, 35, 1, 'eq', 0, params );


    mat42                         = tims42 + [ zeros( length( durs42 ), 1 ) round( durs42 * 20000 ) - 1 ];
    mat43                         = tims43 + [ zeros( length( durs43 ), 1 ) round( durs43 * 20000 ) - 1 ];
    mat45                         = tims45 + [ zeros( length( durs45 ), 1 ) round( durs45 * 20000 ) - 1 ];
    mat47                         = tims47 + [ zeros( length( durs47 ), 1 ) round( durs47 * 20000 ) - 1 ];
    mat44                         = tims44 + [ zeros( length( durs44 ), 1 ) round( durs44 * 20000 ) - 1 ];
    mat46                         = tims46 + [ zeros( length( durs46 ), 1 ) round( durs46 * 20000 ) - 1 ];
    mat40                         = tims40 + [ zeros( length( durs40 ), 1 ) round( durs40 * 20000 ) - 1 ];
    mat32                         = tims32 + [ zeros( length( durs32 ), 1 ) round( durs32 * 20000 ) - 1 ];
    mat34                         = tims34 + [ zeros( length( durs34 ), 1 ) round( durs34 * 20000 ) - 1 ];
    mat31                         = tims31 + [ zeros( length( durs31 ), 1 ) round( durs31 * 20000 ) - 1 ];
    mat33                         = tims33 + [ zeros( length( durs33 ), 1 ) round( durs33 * 20000 ) - 1 ];
    mat39                         = tims39 + [ zeros( length( durs39 ), 1 ) round( durs39 * 20000 ) - 1 ];
    mat41                         = tims41 + [ zeros( length( durs41 ), 1 ) round( durs41 * 20000 ) - 1 ];
    mat35                         = tims35 + [ zeros( length( durs35 ), 1 ) round( durs35 * 20000 ) - 1 ];

% generate seqs and smat
mats                            = {mat32, mat42, mat33, mat43};
seqmat                          = [14 19 20 28; 14 19 20 28; 14 19 20 28; 14 19 NaN NaN] * 20;
[ seqs, smat ]                  = seqdetect( mats, seqmat, 1 );

seqs0=seqs;
seqs = seqs0(1:900,:);
smat0=smat;
smat=smat0(1:900*length(mats),:,:,:);

save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_27/seqs_2', 'seqs')
save ('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_27/smat_2', 'smat')

evtfname = [ filebase '.evt.p01' ]; rc = make_evt_file( evtfname, tims, 20000 );

% compute peth, gain, and significance
[ r1, seqs1 ] = seqpeth( s, seqs, smat );
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_r1', 'r1')

load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/seqs_2.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/smat_2.mat')
[ r2, seqs2 ] = seqpeth( s, seqs, smat );
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_r2', 'r2')
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_seqs2', 'seqs2')

load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/seqs_3.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/smat_3.mat')
[ r3, seqs3 ] = seqpeth( s, seqs, smat );
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_r3', 'r3')
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_seqs3', 'r3')

load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/seqs_4.mat')
load('/media/shirly/C22865A128659567/mice/mC400/dat/mC400_25/smat_4.mat')
[ r4, seqs4 ] = seqpeth( s, seqs, smat );
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_r4', 'r4')

[ r5, seqs5 ] = seqpeth( s, seqs, smat );
save ('/media/shirly/C22865A128659567/mice/mC400/mat/seq/mC400_26_r5', 'r5')

% plot peths

[ shankclu, fpeth, sgain, pAct ]= seqsort( r1, r2, 'cmp', 'intersect', 'graphics', 1 )
title ('mC400__26__seq5')
% save figure
resize = '-bestfit'; pstr = '-dpdf';
resize = ''; pstr = '-dpng';
figname = '/media/shirly/C22865A128659567/mice/mC400/figs/seqs/mC400_26_part1_seq45';

figname = '/media/shirly/C22865A128659567/mice/mC400/figs/rast/mC400_26_seq4';
fig = gcf;
resize = '-fillpage';
figure( fig );
figi = gcf;
figi.Renderer = 'painters';
pause( 0.2 )
print( fig, pstr, figname, resize )


filename= 'mC400_21';
filebase = filebase_lookup(filename);
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, [], [], [], 'spearman' );


schans1                     = [ 32 46 42 45];
shanknums = [ 1 7 ; 4 4;3 2; 4 2];
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, 'stim', schans1, 'spearman' );
figname = '/media/shirly/Data/_Shirly/Lab/uLED/fig8/correlation_for_review/mC400_25_seq4';
ccname = '/media/shirly/Data/_Shirly/Lab/uLED/fig8/correlation_for_review/mC400_25_seq4_cc';
save (ccname, 'cc')