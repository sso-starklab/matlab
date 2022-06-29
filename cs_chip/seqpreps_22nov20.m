% filebase        = filebaseLookup( 'mC400', -23 );
% seqnums         = [ 2 3 ];

% seqpreps( filebase, seqnums )

% revisions
% 16-nov-20 added template matching including shuffling
%           modified for two seperate sequences

function [ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, seqnums )

[ ~, filename ]                 = fileparts( filebase );

% load( [ fileparts( filebase ) '/r2_new.mat' ] );        % assume that content is r2
% load( [ fileparts( filebase ) '/r3_new.mat' ] );        % assume that content is r3
% L1              = r2;
% L2              = r3;

% pipeline:
% (1) decide on filebase, decide on sequences (in terms of values, channels, order, and duration)

% (2) use get_triggers to load information for the relevant parameters
% here, use all waveforms, 5-30 ms durations, all values
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

% (3) for each sequence: use seqdetect to generate seqs and smat
switch filename
    case 'mC400_23'
        
        % sequence #2 for this day (from google sheet of mC400_s)
        mats1                	= {mat43, mat32, mat41, mat46, mat31, mat44, mat33, mat42, mat45};
        seqmat1              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
        
        % sequence #3 for this day (from google sheet of mC400_s)
        mats2               	= {mat43, mat32, mat41, mat44, mat31, mat46, mat33, mat42, mat45};
        seqmat2              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
        
    case 'mC400_26'
        % sequence #4 for this day (from google sheet of mC400_s)
        mats1                	= {mat34, mat44, mat47, mat42, mat45, mat32, mat46, mat33};
        seqmat1                 = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25;14 22 20 25;14 22 20 25; 14 22 NaN NaN] * 20;
        
        % sequence #5 for this day (from google sheet of mC400_s)
        mats2                 	= {mat32, mat42, mat46, mat44, mat45, mat34, mat47, mat33, mat43};
        seqmat2                 = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25;14 22 20 25;14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;
end

[ seqs1, smat1 ]              	= seqdetect( mats1, seqmat1, 1 );
[ seqs2, smat2 ]              	= seqdetect( mats2, seqmat2, 1 );

% (4) for each sequence: use seqpeth to compute peth, gain, and significance for each unit
ilev                            = 'B';
s                               = load_spikes( filebase, shanknums, ilev );
[ r1, seqs1 ]                   = seqpeth( s, seqs1, smat1 );
[ r2, seqs2 ]                   = seqpeth( s, seqs2, smat2 );

% (5) decide on a pair of sequences, use seqsort to plot peths of the
% (e.g. intersected) units according to their responses during the first sequence
[ shankclu0, fpeth, fbins ]     = seqsort( r1, r2, 'cmp', 'intersect', 'graphics', 0 );

%------------------------------------------------------------------------
% (6) plot rasters for a few trials from each sequence
% decide on: 
% -periods (subset of seqs)
% -unit identities (subset of r.shankclu)
shankclu                        = shankclu0;
clu                             = s.clu;
res                             = s.res;
map                             = s.map;

figs( 1 )                       = figure;

subplot( 2, 1, 1 )
periods1                        = seqs1;
seqrast( clu, res, periods1( 1 : 50 : end, : ), 'map', map, 'shankclu', shankclu, 'stimes', r1.stimes, 'spike_length', '.' );
xlims( 1, : )                   = xlim;

subplot( 2, 1, 2 )
periods2                         = seqs2;
seqrast( clu, res, periods2( 1 : 50 : end, : ), 'map', map, 'shankclu', shankclu, 'stimes', r2.stimes, 'spike_length', '.' );
xlims( 2, : )                   = xlim;

xlims                           = [ min( xlims( : , 1 ) ) max( xlims( : , 2 ) )  ];
for i                           = 1 : 2
    subplot( 2, 1, i )
    xlim( xlims )
end


%------------------------------------------------------------------------
% (7) quantify fit of each and every trial with the template
% rerun seqpeth and seqsort s.t. template will be shorter than data
[ r1, seqs1 ]                   = seqpeth( s, seqs1, smat1, 'nT_calc', [ 0 1 ] );
[ r2, seqs2 ]                   = seqpeth( s, seqs2, smat2, 'nT_calc', [ 0 1 ] );
[ shankclu0, fpeth, fbins ]     = seqsort( r1, r2, 'cmp', 'intersect', 'nT_calc', [ 0 1 ], 'graphics', 0 );
[ mat1, bins1, rh1, ah1 ]       = seqrast( clu, res, periods1, 'map', map, 'shankclu', shankclu, 'graphics', 0 );
[ mat2, bins2, rh2, ah2 ]       = seqrast( clu, res, periods2, 'map', map, 'shankclu', shankclu, 'graphics', 0 );
matsi(1).mats = mat1;
matsi(2).mats = mat2;
binsi(1).bins = bins1;
binsi(2).bins = bins2;
figs( 3 )                       = figure;

for l = 1:2
% template:
fpeth1                          = fpeth{ l }; % bin_number x unit, with firing rate in each bin
fbins1                          = fbins{ l }( 1, : )'; % time at each bin [s]
binsizeT                        = diff( fbins1( 1 : 2 ) );
fpeth1s                         = scale( fpeth1' )';
template                        = fpeth1s( : );
nbins0                          = size( fpeth1, 2 );

% data:
binsizeD                        = diff( binsi(l).bins( 1 : 2 ) );
BS                              = binsizeT / binsizeD;
if BS ~= round( BS )
    if abs( BS - round( BS ) ) < sqrt( eps )
        BS                      = round( BS );
    else
        error( 'check - BS must be an integer' )
    end
end

% bin the data to match the template binsize
nunits                          = size( s.shankclu, 1 );
ntrials                         = size( matsi(l).mats, 2 ) / nunits;
pads                            = ( BS - 1 ) / 2;                  % before binning, pad with zeros to prevent offset
z                               = zeros( pads, ntrials * nunits );
data                            = bincols( [ z; matsi(l).mats; z ], BS );
datat                           = bincols( [ binsi(l).bins( 1 ) * ones( pads, 1 ); binsi(l).bins; binsi(l).bins( end ) * ones( pads, 1 ) ], BS, [], 'mean' );

% filter the data with the Gaussian kernel used for the fpeth
sdGauss                         = 1;
g                               = makegaussfir( sdGauss, 1 );
dataf                           = firfilt( full( data ), g );

% clip the data to match the template duration
slct                            = [ find( datat >= fbins1( 1 ), 1, 'first' ) find( datat >= fbins1( end ), 1, 'first' ) ];
nbins1                          = diff( slct ) + 1;
if nbins0 ~= nbins1 
    error( 'duration mismatch' )
end
sidx                            = slct( 1 ) : slct( 2 );
dataf                           = dataf( sidx, : );
datat                           = datat( sidx, : );
%[ fbins1 - datat ]

% generate a shuffled version of the data
%mat1s = stmix( mat1, [], 1 );            % shuffle each column
nums                            = ceil( rand( length( find( matsi(l).mats ) ), 1 ) * numel( matsi(l).mats ) ); % shuffle the entire matrix
mat1s                           = zeros( size( matsi(l).mats ) );
mat1s( nums )                   = 1; % almost OK (large numbers)
datas                           = bincols( [ z; mat1s; z ], BS ); % bin
datasf                          = firfilt( full( datas ), g ); % filter
%datasf = fliplr( datasf( sidx, : ) );% clip
datasf                          = datasf( sidx, : );% clip

cc                              = NaN( ntrials, 1 );
cc_shuf                         = NaN( ntrials, 1 );
% shirly: think about how to write this without for (hints: use reshape, permute, etc.)
for i                           = 1 : ntrials
    x                           = dataf( :, i : ntrials : nunits * ntrials )'; 
    xs                          = scale( x' )';
    x1                          = xs( : );
    cc( i )                     = calc_pearson( x1, template );
%    cc( i )                     = calc_spearman( x1, template );
    
    x_shuf                      = datasf( :, i : ntrials : nunits * ntrials )';
    x_shufs                    	= scale( x_shuf' )';
    xs1                         = x_shufs( : );
    cc_shuf( i )                = calc_pearson( xs1, template );
%    cc_shuf( i )                = calc_spearman( xs1, template );
end

%------------------------------------------------------------------------
% (7.1) summarize the template matching analysis

% stats 
pval0                           = signrank( cc, 0 ); 
pval0_shuf                      = signrank( cc_shuf, 0 );
pvals                           = utest( cc, cc_shuf ); 

% histograms
binsize                         = 0.02;
edges                           = ( -1 - binsize / 2 ) : binsize : ( 1 + binsize / 2 );
bins                            = ( edges( 1 : end - 1 ) + edges( 2 : end ) )' / 2;
nbins                           = length( bins );
h                               = histc( cc, edges );
h( nbins )                      = sum( h( nbins + [ 0 1 ] ) ); 
h( nbins + 1 )                  = [];
h1                              = h;
h                               = histc( cc_shuf, edges );
h( nbins )                      = sum( h( nbins + [ 0 1 ] ) ); 
h( nbins + 1 )                  = [];
h1_shuf                         = h;

% prep for graphics
hranges                         = [ minmax( find( h1 ~= 0 ) ) minmax( find( h1_shuf ~= 0 ) ) ];
xlims                           = [ -1 1 ] * ( max( abs( bins( minmax( hranges ) ) ) ) + binsize * 2 );
color1                          = [ 1 0 0 ];
color2                          = [ 0 0 0.7 ];

% graphics
figs( 2 )                       = figure;
sh1                             = stairs( bins, h1 );
set( sh1, 'Color', color1 );
hold on
sh2                             = stairs( bins, h1_shuf );
set( sh2, 'Color', color2 );
xlim( xlims )
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'Correlation coefficient' )
ylabel( 'Number of trials' )
alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
alines( median( cc ), 'x', 'color', color1, 'linestyle', '--' );
alines( median( cc_shuf ), 'x', 'color', color2, 'linestyle', '--' );
title( sprintf( 'cc: %0.3g; cc shuf: %0.3g; compare: %0.3g', pval0, pval0_shuf, pvals ) );

%------------------------------------------------------------------------
% (7.2) repeat step 6 with selected trials
% now, select e.g. 10 extremum trials from each side, and plot those
nselected                       = 4;
[ ~, sidx ]                     = sort( cc );
tidx1                           = sidx( 1 : nselected );
tidx2                           = sidx( ( ntrials - nselected + 1 ) : ntrials );
tidx                            = sort( [ tidx1; tidx2 ] );

if l==1
 figs( 3 )    
subplot( 2, 1, l )
periods1                        = seqs1;
seqrast( clu, res, periods1( tidx, : ), 'map', map, 'shankclu', shankclu, 'stimes', r1.stimes );
xlims3( 1, : )                   = xlim;
axis square
else
    figs( 3 ) 
subplot( 2, 1, l )
periods2                         = seqs2;
seqrast( clu, res, periods2( tidx, : ), 'map', map, 'shankclu', shankclu, 'stimes', r2.stimes );
xlims3( 2, : )                   = xlim;
axis square

end
end
xlims                           = [ min( xlims3( : , 1 ) ) max( xlims3( : , 2 ) )  ];
for i                           = 1 : 2
    subplot( 2, 1, i )
    xlim( xlims )
end


return

% EOF


shanknums = [ 1 4; 1 8; 3 10; 3 12; 4 7; 4 9; 4 14 ];
filebase = filebaseLookup( 'mC400', -26 )
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums );

