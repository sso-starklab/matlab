filebase        = filebaseLookup( 'mC400', -23 );
seqnums         = [ 2 3 ];

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

% sequence #2 for this day (from google sheet of mC400_s)
mats2                           = {mat43, mat32, mat41, mat46, mat31, mat44, mat33, mat42, mat45};
seqmat2                         = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
[ seqs2, smat2 ]              	= seqdetect( mats2, seqmat2, 1 );

% sequence #3 for this day (from google sheet of mC400_s)
mats3                        	= {mat43, mat32, mat41, mat44, mat31, mat46, mat33, mat42, mat45};
seqmat3                      	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
[ seqs3, smat3 ]              	= seqdetect( mats3, seqmat3, 1 );

% (4) for each sequence: use seqpeth to compute peth, gain, and significance for each unit
ilev                            = 'B';
s                               = load_spikes( filebase, [], ilev );
[ r2, seqs2 ]                   = seqpeth( s, seqs2, smat2 );
[ r3, seqs3 ]                   = seqpeth( s, seqs3, smat3 );

% (5) decide on a pair of sequences, use seqsort to plot peths of the
% (e.g. intersected) units according to their responses during the first sequence
[ shankclu0, fpeth, sgain, pAct ] = seqsort( r4_part2, r5_part2, 'cmp', 'intersect', 'graphics', 1 );

% (6) plot rasters for a few trials
% decide on: 
% -periods (subset of seqs)
% -unit identities (subset of r.shankclu)
periods                         = seqs;
shankclu                        = shankclu0;
clu                             = s.clu;
res                             = s.res;
map                             = s.map;
binsizeSEC                      = 0.001; % [s]
padSEC                          = [ -0.05 0.05 ]; % [s]
spkFs                           = 20000;


periods0 = periods( 10:35:end, : );


figure, 
clf, 
[ mat, bins, rh, ah ] = seqrast_04nov20( clu, res, periods, 'map', map, 'shankclu', shankclu, 'stimes', r4.stimes );
figure, 
clf, 
[ ~, ~, rh, ah ] = seqrast( clu, res, periods, 'map', map, 'shankclu', shankclu, 'stimes', r4.stimes, 'emptyTrials', 0 );

% now, extract the times of the spikes of these units at these times
% periods, shankclu, clu, res
% binsizeSEC
% padSEC
% 
% % (1) get the rasters as sparse martices for all units in clu
% binsize                         = ceil( binsizeSEC * spkFs );                       % [samples]
% h                               = [ 0 max( diff( periods, [], 2 ) + 1 ) - 1 ];      % [samples]
% h1                              = [ floor( h( 1 ) / binsize ) ceil( h( 2 ) / binsize ) ]; % [bins]
% h2                              = [ floor( padSEC( 1 ) / binsizeSEC ) ceil( padSEC( 2 ) / binsizeSEC ) ]; % [bins]
% halfwin                         = h1 + h2;
% [ rI, bins ]                    = get_rasters( clu, res, [], periods( :, 1 )...
%     , 'binsize', binsize, 'halfwin', halfwin );
% 
% % (2) now, keep only the desired units
% uidx                            = ismember( map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' );
% if isempty( uidx )
%     error( '' )
% end
% r                               = rI( uidx );
% 
% % (3) organize in the requested order
% %map1                            = map( uidx, : )
% [ ~, bb ]                       = sortrows( shankclu( :, 1 : 2 ) );
% [ ~, ridx ]                     = sort( bb );
% %map1( ridx, : )
% r                               = r( ridx );
% 
% % (4) combine into one array
% ntrials                         = size( periods, 1 );
% nbins                           = size( r{ 1 }, 1 );
% nunits                          = length( r );
% if ~isequal( size( r{ 1 }, 2 ), ntrials )
%     error( 'mismatch' )
% end
% if ~isequal( size( r{ 1 }, 1 ), nbins )
%     error( 'mismatch' )
% end
% mat                             = sparse( nbins, nunits * ntrials );
% for i                           = 1 : nunits
%     cidx                        = ntrials * ( i - 1 ) + 1 : ntrials * i; 
%     mat( :, cidx )              = r{ i };
% end
% 
% %         [ r1, bins ]                    = get_rasters( clu, res, [], periods( :, 1 )...
% %     , 'binsize', 1, 'halfwin', h );
