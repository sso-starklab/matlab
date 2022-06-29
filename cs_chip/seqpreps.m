% filebase        = filebaseLookup( 'mC400', -23 );
% seqnums         = [ 2 3 ];

% seqpreps( filebase, seqnums )

% revisions
% 16-nov-20 added template matching including shuffling
% 30-nov-20 added template matching for the stim pattern
% 10-jan-20 added filename cases

function [ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, templateType, schans1, corrtype )

nargs                           = nargin;
if nargs < 3 || isempty( templateType )
    templateType                = 'peth';
end
if nargs < 4 || isempty( schans1 )
    schans1                     = []; % one-based, must be identical to number of units
end
if isempty( schans1 ) && isequal( templateType, 'stim' )
    error( 'missing schans1' )
end
if nargs < 5 || isempty( corrtype )
    corrtype                    = 'pearson';
end

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
    case 'mC400_11'
        % sequence #1 for this day (from google sheet of mC400_s)
%         mats1                            = {mat42, mat45, mat43, mat46};
%         seqmat1                          = [11 17 8 14; 11 17 8 14; 11 17 8 14; 11 17 NaN NaN] * 20;
        % special conditions for seq1: 
%             seqs0=seqs1;
%             seqs1 = seqs0(1:1008,:);
%             smat0=smat1;
%             smat1=smat0(1:1008*length(mats1),:,:,:);

        % sequence #2 for this day (from google sheet of mC400_s)
        mats1                            = {mat42, mat45, mat43, mat46, mat44, mat47};
        seqmat1                          = [11 17 30 35; 11 17 30 35; 11 17 30 35; 11 17 30 35;11 17 30 35;11 17 NaN NaN] * 20;
        
        % duplication of sequence #2 to make the code work
        mats2                            = {mat42, mat45, mat43, mat46, mat44, mat47};
        seqmat2                          = [11 17 30 35; 11 17 30 35; 11 17 30 35; 11 17 30 35;11 17 30 35;11 17 NaN NaN] * 20;
    case 'mC400_13'
        % sequence #3 for this day (from google sheet of mC400_s)
%         mats2                            = {mat42, mat47, mat44, mat45, mat43};
%         seqmat2                          = [11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 NaN NaN] * 20;
%         mats1                            = { mat44, mat45, mat43};
%         seqmat1                          = [11 15 25 30; 11 15 25 30; 11 15 NaN NaN] * 20;
        % special conditions for seq3: 
%             seqs0=seqs1;
%             seqs1 = seqs0(1445:end,:);
%             smat0=smat1;
%             smat1=smat0((1444*length(mats1))+1:2444*length(mats1),:,:,:);
        % sequence #4 for this day (from google sheet of mC400_s)
%         mats1                            = {mat42, mat45, mat43, mat47, mat44};
%         seqmat1                          = [11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 NaN NaN] * 20;
        mats1                            = {mat42, mat45, mat43};
        seqmat1                          = [11 15 25 30; 11 15 25 30; 11 15 NaN NaN] * 20;        
        % duplication of sequence #4 to make the code work
        mats2                            = {mat42, mat45, mat43, mat47, mat44};
        seqmat2                          = [11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 25 30; 11 15 NaN NaN] * 20;
        
    case 'mC400_17'
%         % sequence #1 for this day (from google sheet of mC400_s)
%         mats1                            = {mat44, mat43, mat45, mat47, mat42, mat46};
%         seqmat1                          = [11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20 ;11 15 NaN NaN] * 20;
%         mats1                            = {mat43, mat45,mat47, mat42};
%         seqmat1                          = [11 15 15 20; 11 15 15 20;11 15 15 20; 11 15 NaN NaN] * 20;
        mats2                            = {mat44, mat43, mat45};
        seqmat2                          = [11 15 15 20; 11 15 15 20; 11 15 NaN NaN] * 20;
        
         % sequence #2 for this day (from google sheet of mC400_s)
%         mats1                            = {mat42, mat45, mat43, mat47, mat44, mat46};
%         seqmat1                          = [11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20 ;11 15 NaN NaN] * 20;
                % special conditions for seq2: 
%                 seqs0=seqs1;
%                 seqs1 = seqs0(1:1000,:);
%                 smat0=smat1;
%                 smat1=smat0(1:1000*length(mats1),:,:,:);

        % sequence #3 for this day (from google sheet of mC400_s)
%         mats1                            = {mat32, mat42, mat45, mat43, mat47, mat44, mat46};
%         seqmat1                          = [11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20 ;11 15 NaN NaN] * 20;
%         mats1                            = {mat32, mat42, mat45};
%         seqmat1                          = [11 15 15 20; 11 15 15 20;11 15 NaN NaN] * 20;
%         mats2                            = {mat32, mat42, mat45};
%         seqmat2                          = [11 15 15 20; 11 15 15 20;11 15 NaN NaN] * 20;

        % sequence #4 for this day (from google sheet of mC400_s)
%         mats1                            = {mat32, mat42, mat35, mat47, mat44, mat46, mat43, mat45};
%         seqmat1                          = [11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20 ;11 15 NaN NaN] * 20;

%         
%         % sequence #5 for this day (from google sheet of mC400_s)
        mats1                            = {mat42, mat32, mat47, mat35, mat44, mat46, mat43, mat45};
        seqmat1                          = [11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20; 11 15 15 20 ;11 15 NaN NaN] * 20;

    case 'mC400_19'
%         % duplication of sequence #2 to make the code work
%         mats1                            = {mat32, mat33, mat47, mat44};
%         seqmat1                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
%         mats2                            = {mat32, mat33, mat47, mat44};
%         seqmat2                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;        
    case 'mC400_21'
        % sequence #1 for this day (from google sheet of mC400_s)
        mats2                            = {mat32, mat42, mat47, mat44, mat34, mat45};
        seqmat2                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;
        
        % sequence #2 for this day (from google sheet of mC400_s)
        mats1                            = {mat32, mat42, mat47, mat34, mat44, mat45};
        seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;

        % sequence #3 for this day (from google sheet of mC400_s)
%         mats1                            = {mat46, mat34, mat42, mat47, mat44, mat32, mat45};
%         seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;
% 
%         % sequence #4 for this day (from google sheet of mC400_s)
%         mats1                            = {mat46, mat34, mat42, mat47, mat32, mat44, mat45};
%         seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;
        
    case 'mC400_23'
        % sequence #1 for this day (from google sheet of mC400_s)
        mats1                            = {mat43, mat32, mat41, mat46, mat31, mat44, mat33, mat45, mat42};
        seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
%         mats1                            = {mat41, mat46, mat31, mat44, mat33};
%         seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 NaN NaN] * 20;        
        % sequence #2 for this day (from google sheet of mC400_s)
%         mats1                	= {mat43, mat32, mat41, mat46, mat31, mat44, mat33, mat42, mat45};
%         seqmat1              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
%         mats1                	= {mat41, mat46, mat31, mat44, mat33, mat42, mat45};
%         seqmat1              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
        
        % sequence #3 for this day (from google sheet of mC400_s)
        mats2               	= {mat43, mat32, mat41, mat44, mat31, mat46, mat33, mat42, mat45};
        seqmat2              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
%         mats1               	= {mat41, mat44, mat31, mat46, mat33, mat42, mat45};
%         seqmat1              	= [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;        
%         % sequence #4 for this day (from google sheet of mC400_s)
        mats1                            = {mat41, mat46, mat31, mat44, mat33, mat42, mat45, mat32, mat43};
        seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
        
%         % sequence #5 for this day (from google sheet of mC400_s)
%         mats1                            = {mat41, mat46, mat31, mat44, mat33, mat43, mat32, mat45, mat42};
%         seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
%         mats1                            = {mat41, mat46, mat31, mat44, mat33, mat43};
%         seqmat1                          = [14 20 22 30; 14 20 22 30; 14 20 22 30;14 20 22 30;14 20 22 30;14 20 NaN NaN] * 20;
    
    case 'mC400_25'
        % sequence #2 for this day (from google sheet of mC400_s)
        mats2                            = {mat42, mat47, mat32, mat45, mat34, mat46, mat33, mat44, mat43};
        seqmat2                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
%         mats1                            = {mat42, mat47, mat32, mat45};
%         seqmat1                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;        
        % sequence #3 for this day (from google sheet of mC400_s)
%         mats1                            = {mat42, mat47, mat32, mat45, mat34, mat46, mat44, mat33, mat43};
%         seqmat1                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
%         mats1                            = { mat32, mat45, mat34, mat46,mat44};
%         seqmat1                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
%         % sequence #4 for this day (from google sheet of mC400_s)
        mats1                            = {mat34, mat47, mat32, mat46, mat42, mat45, mat44, mat33, mat43};
        seqmat1                          = [12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
        mats1                            = {mat32, mat46, mat42, mat45};
        seqmat1                          = [ 12 19 20 28; 12 19 20 28; 12 19 20 28; 12 19 NaN NaN] * 20;
        
    case 'mC400_26'
        % sequence #1 for this day (from google sheet of mC400_s)
%         mats1                            = {mat34, mat42, mat47, mat32, mat43, mat46};
%         seqmat1                          = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;
        
        % sequence #2 for this day (from google sheet of mC400_s)
%         mats1                            = {mat34, mat43, mat46, mat32, mat42, mat47};
%         seqmat1                          = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;
%         mats1                            = {mat43, mat46, mat32, mat42, mat47};
%         seqmat1                          = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;
        
        % sequence #3 for this day (from google sheet of mC400_s)
%         mats1                            = {mat34, mat42, mat47, mat44, mat45, mat32};
%         seqmat1                          = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;       
%         mats1                            = {mat42, mat47, mat44, mat45, mat32};
%         seqmat1                          = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;           
%         % sequence #4 for this day (from google sheet of mC400_s)
        mats2                	= {mat34, mat44, mat47, mat42, mat45, mat32, mat46, mat33};
        seqmat2                 = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25;14 22 20 25;14 22 20 25; 14 22 NaN NaN] * 20;
        
        % sequence #5 for this day (from google sheet of mC400_s)
        mats1                 	= {mat32, mat42, mat46, mat44, mat45, mat34, mat47, mat33, mat43};
        seqmat1                 = [14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25; 14 22 20 25;14 22 20 25;14 22 20 25; 14 22 20 25; 14 22 NaN NaN] * 20;

        
        
end

[ seqs1, smat1 ]              	= seqdetect( mats1, seqmat1, 1 );
[ seqs2, smat2 ]              	= seqdetect( mats2, seqmat2, 1 );

% (4) for each sequence: use seqpeth to compute peth, gain, and significance for each unit
ilev                            = 'B';
s                               = load_spikes( filebase, shanknums, ilev );
nunits                          = size( s.shankclu, 1 );
[ r1, seqs1, smat1 ]                   = seqpeth( s, seqs1, smat1 );
[ r2, seqs2, smat2 ]                   = seqpeth( s, seqs2, smat2 );

% (5) decide on a pair of sequences, use seqsort to plot peths of the
% (e.g. intersected) units according to their responses during the first sequence
[ shankclu0, fpeth, fbins ]     = seqsort( r1, r2, 'cmp', 'intersect', 'graphics', 1 );

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

binsizeMS                       = 5;            % [ms] used seqpeth
binsizeSEC                      = 0.001;        % [s] used in seqrast
% seqpeth uses bins of binsizeMS; number of bins depends on mean stimulus
% duration (we use diff( nT_calc ) durations)
% seqrast used bins of binsizeSEC; number of bins depends on max stimulus
% duration, and is padded at each side (e.g. by 50 ms)
[ r1, seqs1 ]                   = seqpeth( s, seqs1, smat1, 'nT_calc', [ 0 1 ], 'binsizeMS', binsizeMS );
[ r2, seqs2 ]                   = seqpeth( s, seqs2, smat2, 'nT_calc', [ 0 1 ], 'binsizeMS', binsizeMS );
[ shankclu0, fpeth, fbins ]     = seqsort( r1, r2, 'cmp', 'intersect', 'nT_calc', [ 0 1 ], 'graphics', 0 );
[ mat1, bins1, ~, ~, bs1, hw1 ] = seqrast( clu, res, periods1, 'map', map, 'shankclu', shankclu, 'graphics', 0, 'binsizeSEC', binsizeSEC );
[ mat2, bins2, ~, ~, bs2, hw2 ] = seqrast( clu, res, periods2, 'map', map, 'shankclu', shankclu, 'graphics', 0, 'binsizeSEC', binsizeSEC );

%------------------------------------------------------------------------
% general template properties (time):
fbins1                          = fbins{ 1 }( 1, : )';                      % time at each bin [s]
templateEdges                   = fbins1( [ 1 end ] );
nbins0                          = length( fbins1 );
binsizeT                        = diff( fbins1( 1 : 2 ) );                  % should be identical to binsizeMS/1000

%------------------------------------------------------------------------
% peth template (firing rates):
fpeth1                          = fpeth{ 1 };                               % bin_number x unit, with firing rate in each bin
templateP                      	= scale( fpeth1' )';
% templateP                       = fpeth1s( : );
%nbins0                          = size( fpeth1, 2 );

%------------------------------------------------------------------------
% stim (applied light) template:

% specific to a given sequence:
schans                          = schans1;
periods                         = seqs1; % @ spkFs

% general for the session
par                             = LoadXml( filebase );
spkFs                           = par.SampleRate;
nchans                          = par.nChannels;

nschans                         = length( schans );
nperiods                        = size( periods, 1 );

% now, use the actual bin sizes in the rasters:
winSamples                      = ( hw1 * bs1 + bs1 / 2 * [ -1 1 ] );
periodsW                        = periods( :, 1 ) * ones( 1, 2 ) + ones( nperiods, 1 ) * winSamples;
mdurs                           = diff( winSamples ) + 1;
% winsize                         = diff( datat( 1 : 2 ) ) * spkFs;
% offsetPre                       = round( datat( 1 ) * spkFs / winsize ) * winsize;
% offsetPost                      = round( datat( end ) * spkFs / winsize ) * winsize - 1;
% periods                         = [ periods( :, 1 ) - offsetPre periods( :, 1 ) + offsetPost ];
% durs                            = diff( periods, [], 2 ) + 1;
% mdurs                           = round( mean( durs ) );
% % periods                         = [ periods( :, 1 ) periods( :, 1 ) + mdurs - 1 ];

% make a pointer to the data
a                               = memmapfile( [ filebase '.dat' ], 'Format', 'int16' );

% extract the data for the ith period
stimKeep                        = NaN( mdurs, nschans, nperiods );
for i                           = 1 : 10 : nperiods
    st                          = periodsW( i, 1 );
    et                          = periodsW( i, 2 );
    idx                         = ( ( st - 1 ) * nchans + 1 ) : et * nchans;
    stim                        = a.Data( idx );
    nsamps                      = length( idx ) / nchans;
    stim                        = reshape( stim, [ nchans nsamps ] );
    stim                        = single( stim( schans, : ) );
    stimKeep( :, :, i )         = stim';
end
% average, scale, and dounsample
mstim                           = nanmean( stimKeep, 3 );
sstim                           = scale( mstim );

% bin into binsizeD sized bins
BS                              = bs1;%winsize;%binsizeT * spkFs;
if BS ~= round( BS )
    if abs( BS - round( BS ) ) < sqrt( eps )
        BS                      = round( BS );
    else
        error( 'check - BS must be an integer' )
    end
end
bstim                           = bincols( sstim, BS );
%tstim                           = ceil( bincols( ( winSamples( 1 ) : winSamples( 2 ) )', bs1, [], 'mean' ) );
tstim                           = ( ( winSamples( 1 ) + bs1/2 ) : bs1 : ( winSamples( 2 ) - bs1/2 ) )' / spkFs; % [s]

% (bstim, tstim) are equivalent to (bins1,mat1). verify:
if ~isequal( bins1, tstim )
    error( 'mismatch between stim template and rasters' )
end



% bin into binsizeT sized bins


% clip the template to match the desired template duration


%------------------------------------------------------------------------
% data (mat1, bins1):
binsizeD                        = diff( bins1( 1 : 2 ) );                   % should be identical to binsizeSEC
BS                              = binsizeT / binsizeD;                      % should be odd
if BS ~= round( BS )
    if abs( BS - round( BS ) ) < sqrt( eps )
        BS                      = round( BS );
    else
        error( 'check - BS must be an integer' )
    end
end

% prepare for binning
nunits                          = size( s.shankclu, 1 );
ntrials                         = size( mat1, 2 ) / nunits;
pads                            = ( BS - 1 ) / 2;                           % before binning, pad with zeros to prevent offset

% bin the data to match the desired template binsize
z                               = zeros( pads, ntrials * nunits );
data                            = bincols( [ z; mat1; z ], BS );
datat                           = bincols( [ bins1( 1 ) * ones( pads, 1 ); bins1; bins1( end ) * ones( pads, 1 ) ], BS, [], 'mean' );

% bin the stim template to match the desired template binsize
z1                              = zeros( pads, nschans );
bstimS                          = bincols( [ z1; bstim; z1 ], BS );
% datat applies also to bstimS

% prepare for clipping
sidx                            = inrange( datat, templateEdges );
nbins1                          = sum( sidx );
if nbins0 ~= nbins1 
%     error( 'duration mismatch' )
        fprintf( 'duration mismatch' )
end

% clip the stim template and the time vector to match the desired template duration
templateS                     	= bstimS( sidx, : )';
datat                           = datat( sidx, : );

%------------------------------------------------------------------------
% filter the data with the Gaussian kernel used for the fpeth
sdGauss                         = 1;
g                               = makegaussfir( sdGauss, 1 );
dataf                           = firfilt( full( data ), g );

% clip the data 
dataf                           = dataf( sidx, : );
%[ fbins1 - datat ]

% decide which template to use
switch templateType
    case 'peth'
        template                = templateP;
    case 'stim'
        if ~isequal( nschans, nunits )
            error( 'number of channels selected for template matching must be identical to number of units' )
        end
        template                = templateS;
end
template                        = template( : );

% generate a shuffled version of the data
%mat1s = stmix( mat1, [], 1 );            % shuffle each column
nums                            = ceil( rand( length( find( mat1 ) ), 1 ) * numel( mat1 ) ); % shuffle the entire matrix
mat1s                           = zeros( size( mat1 ) );
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

    switch corrtype
        case 'pearson'
            cc( i )             = calc_pearson( x1, template );
        case 'spearman'
            cc( i )             = calc_spearman( x1, template );
    end
    
    x_shuf                      = datasf( :, i : ntrials : nunits * ntrials )';
    x_shufs                    	= scale( x_shuf' )';
    xs1                         = x_shufs( : );
    switch corrtype
        case 'pearson'
            cc_shuf( i )     	= calc_pearson( xs1, template );
        case 'spearman'
            cc_shuf( i )        = calc_spearman( xs1, template );
    end
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
nselected                       = 10;
[ ~, sidx ]                     = sort( cc );
tidx1                           = sidx( 1 : nselected );
tidx2                           = sidx( ( ntrials - nselected + 1 ) : ntrials );
tidx                            = sort( [ tidx1; tidx2 ] );

figs( 3 )                       = figure;

subplot( 2, 1, 1 )
periods1                        = seqs1;
seqrast( clu, res, periods1( tidx, : ), 'map', map, 'shankclu', shankclu, 'stimes', r1.stimes, 'spike_length', '.' );
xlims( 1, : )                   = xlim;

subplot( 2, 1, 2 )
periods2                         = seqs2;
seqrast( clu, res, periods2( tidx, : ), 'map', map, 'shankclu', shankclu, 'stimes', r2.stimes, 'spike_length', '.' );
xlims( 2, : )                   = xlim;

xlims                           = [ min( xlims( : , 1 ) ) max( xlims( : , 2 ) )  ];
for i                           = 1 : 2
    subplot( 2, 1, i )
    xlim( xlims )
end


return

% EOF


shanknums = [ 1 4; 1 8; 3 10; 3 12; 4 7; 4 9; 4 14 ];
schans1                     = [ 34 44 47 42 45 32 46 33 ]
filebase = filebaseLookup( 'mC400', -26 )
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, 'peth' );
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, 'stim', schans1( 1 : 7 ) );
[ fpeth, mat1, mat2, bins1, bins2, r1, r2, figs ] = seqpreps( filebase, shanknums, 'stim', schans1( 1 : 7 ), 'spearman' );

