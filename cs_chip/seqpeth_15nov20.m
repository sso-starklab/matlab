% seqpeth                   compute PETH during sequences
%
% call                      r = seqpeth( s, seqs, smat )
%   
% gets                      s                       see load_spikes
%                           seqs, smat              see seqdetect
% 
% optional arguments
% 
%                           padBuffer   {0.01}      [s]     time not to use for baseline before every seqs
%                           ilev        {'B'}               argument 3 to load_spikes 
%                           nT_base     {4}                 number of stim duration before/after trigger for baseRate
%                           nT_calc     {[-1 2 ]}           number of stim duration for fpeth
%
% returns                   r           structure with fields (one row per unit)
%                                           shankclu
%                                           nchans
%                                           nseqs
%                                           totdur
%                                           pAct
%                                           gain
%                                           fpeth
%                                           fbins
%
% calls                     LoadXml (blab)
%                           ParseArgPairs (general)
%                           intersectranges, inrange, isoverlap, sortranges (sets)
%                           load_spikes, multipeth_make, poissonTest (spikes)
%
% see also                  seqdetect, seqsort

% 15-oct-20 SSo & ES

% revisions
% 18-oct-20 actually written from script
% 03-nov-20 added seqs (after pruning) to output
% 15-nov-20 support specific bin sizes

% pipeline:
% (1) decide on filebase, decide on sequences (in terms of values, channels, order, and duration)
% (2) use get_triggers to load information for the relevant parameters
% (3) for each sequence: use seqdetect to generate seqs and smat
% (4) for each sequence: use seqpeth to compute peth, gain, and significance for each unit
% (5) decide on a pair of sequences, use seqsort to plot peths of the
% (e.g. intersected) units according to their responses during the first sequence

function [ r, seqs ] = seqpeth( s, seqs, smat, varargin )

% argument handling
nargs                           = nargin;
if nargs < 3 || isempty( s ) || isempty( seqs ) || isempty( smat )
    return
end
[ padBuffer, nT_base, nT_calc, spkFs, binsizeMS ]            = ParseArgPairs(...
    { 'padBuffer', 'nT_base', 'nT_calc', 'spkFs', 'binsizeMS' }...
    , { 0.01, 4, [ -1 2 ], 20000, 5 }...
    , varargin{ : } );

% check input for proper structure and content
if ~ismatrix( seqs )
    error( 'seqs must be a matrix' )
end
if ~ismatrix( smat )
    error( 'smat must be a matrix' )
end
[ nseqs, c1 ]                   = size( seqs );
[ nsmat, c2 ]                   = size( smat );
if c1 ~= 2 || c2 ~= 4
    error( 'seqs must be a 2-column matrix and smat must be a 4-column matrix' )
end
if ~isequal( sortranges( seqs ), seqs )
    error( 'seqs must be a matrix of ranges' )
end
if ~isequal( sortranges( smat( :, 1 : 2 ) ), smat( :, 1 : 2 ) )
    error( 'first two columns of smat must form a matrix of ranges' )
end
uchans                          = unique( smat( :, 3 ) );
nchans                          = length( uchans );
if nseqs * nchans ~= nsmat
    error( 'smat must have %d times the number of rows that seqs has', nchans )
end
imat                            = intersectranges( seqs, smat( :, 1 : 2 ) ); 
if ~isequal( imat, smat( :, 1 : 2 ) )
    error( 'some extra ranges in smat, or some missing ranges in seqs' )
end
[ ~, loc ]                      = isoverlap( smat( :, 1 : 2 ), seqs );
if ~isequal( unique( uhist( loc ) ), nchans )
    error( 'some ranges in seqs do not encompass %d channels', nchans )
end

% remove sequences with too short pre-sequence intervals
spre                            = smat( 1 : nchans : nsmat, 4 ) / spkFs;
ridx1                           = spre <= padBuffer;
ridx2                           = ismember( loc, find( ridx1 ) );
seqs( ridx1, : )                = [];
smat( ridx2, : )                = [];
[ ~, loc ]                      = isoverlap( smat( :, 1 : 2 ), seqs );
nseqs                           = size( seqs, 1 );
nsmat                           = size( smat, 1 );
spre                            = smat( 1 : nchans : nsmat, 4 ) / spkFs;

% verify that nT before speriods is free from other stims and adapt nT to maximize total base duration
sdur                            = diff( seqs / spkFs, [], 2 );
msdur                           = min( sdur );                              % if get_triggers was tight, and sequences were exact, min is very close to mean 
baseDur1                        = sum( spre >= ( msdur * nT_base ) ) * ( nT_base - padBuffer );     % keep only sequences with long pre
baseDur2                        = length( sdur ) * ( floor( min( spre / msdur ) ) - padBuffer );       % reduce nT but keep all sequences
if baseDur1 > baseDur2
    nTB                         = nT_base;
    ridx1                       = spre < ( msdur * nT_base );
    ridx2                       = ismember( loc, find( ridx1 ) );
    seqs( ridx1, : )            = [];
    smat( ridx2, : )            = [];
    nseqs                       = size( seqs, 1 );
else
    nTB                         = floor( min( spre / msdur ) );
end
if nTB == 0
    error( 'chosen nT is %d, check input', nTB )
end    

% use the remaining seqs/smat for all subsequent computations
speriods                     	= seqs / spkFs;

% 1. compute PETH for the entire seq and extract base rate
[ bpeth, bbins, ~, ~, bdurs, shankclu ]   = multipeth_make( '', 'data', s, 'stimTypes', 0, 'periods', speriods, 'nT', nTB, 'toplot', 0, 'sdGauss', 0 );
brange                          = -[ mean( bdurs ) * nTB padBuffer ];
bidx                            = inrange( bbins, brange );
baseRate                       	= mean( bpeth( bidx, : ), 1 ); % base firing rates

% 2. compute PETH for each invdividual event within the sequence and evaluate gain and significance
elements                        = unique( smat( :, 3 ) );
nelements                       = length( elements );
nunits                          = size( bpeth, 2 );
pAct                            = NaN( nunits, nelements );
gain                            = NaN( nunits, nelements );
for i                           = 1 : nelements
    element                     = elements( i );
    periods                     = smat( smat( :, 3 ) == element, 1 : 2 ) / spkFs;
    [ ipeth, ibins, ~, ~, idurs ] = multipeth_make( '', 'data', s, 'stimTypes', 0, 'periods', periods, 'toplot', 0, 'sdGauss', 0 );
    irange                      = [ 0 mean( idurs ) ];
    bidx                        = inrange( ibins, irange );
    inTime                      = sum( idurs );
    binSize                     = diff( ibins( 1 : 2 ) );
    nTrigs                      = size( periods, 1 );
    inCount                     = sum( ipeth( bidx, : ) * nTrigs * binSize );
    pInc                        = poissonTest( baseRate, inCount, inTime * ones( 1, nunits ) );
    pAct( :, i )                = pInc;
    gain( :, i )                = ( inCount / inTime ) ./ baseRate;
end

% 3. compute PETH for the entire seq (as in step 1) but for a narrower range and with smoothing
%[ fpeth, fbins ]                = multipeth_make( '', 'data', s, 'stimTypes', 0, 'periods', speriods, 'toplot', 0, 'nT', nT_calc );
binsize                         = round( binsizeMS / 1000 * spkFs );        % [samples]
mdur                            = mean( diff( round( speriods * spkFs ), [], 2 ) + 1 ); % [samples]
ntmax                           = max( abs( nT_calc ) );
halfwin                         = ntmax * ceil( mdur / binsize );
[ fpeth, fbins ]                = multipeth_make( '', 'data', s, 'stimTypes', 0, 'periods', speriods, 'toplot', 0, 'binsize', binsize, 'halfwin', halfwin );
drange                          = halfwin * binsizeMS / ntmax / 1000 * nT_calc; % [s]
kidx                            = inrange( fbins, drange );
fbins                           = fbins( kidx );
fpeth                           = fpeth( kidx, : );
fbins                           = round( fbins / ( binsizeMS / 1000 ) ) * binsizeMS / 1000;

% 5. find stimulus patterns [s]
nblocks                         = size( smat, 1 ) / nelements;
t0                              = repmat( smat( 1 : nelements : ( nblocks * nelements ), 1 ), [ 1 nelements ] )';
d                               = smat( :, 1 : 2 ) - t0( : ) * ones( 1, 2 );
m1                              = ones( nelements, 1 ) * ( 1 : nelements : ( nblocks * nelements ) );
m2                              = ( 0 : 1 : ( nelements - 1 ) )' * ones( 1, nblocks );
m                               = m1 + m2;
d1                              = d( :, 1 );
d2                              = d( :, 2 );
stimes                          = [ mean( d1( m ), 2 ) mean( d2( m ), 2 ) ] / spkFs;

% 5. assign output
uv                              = ones( nunits, 1 );
r.shankclu                      = shankclu;
r.nchans                        = uv * nchans;
r.nseqs                         = uv * nseqs;
r.totdur                        = uv * sum( bdurs );
r.stimes                        = repmat( permute( stimes, [ 3 1 2 ] ), [ nunits 1 1 ] );

r.pAct                          = pAct;
r.gain                          = gain;
r.fpeth                         = fpeth';
r.fbins                         = uv * fbins;

return

% EOF

filebase                        = filebaseLookup( 'mC400', -15 );
ilev                            = 'B';
s                               = load_spikes( filebase, [], ilev );

% 1 is bad, use 2/3/4
load('/Volumes/slab1/mC400/dat/mC400_15/seqs_2.mat')
load('/Volumes/slab1/mC400/dat/mC400_15/smat_2.mat')
r2 = seqpeth( s, seqs, smat );

load('/Volumes/slab1/mC400/dat/mC400_15/seqs_3.mat')
load('/Volumes/slab1/mC400/dat/mC400_15/smat_3.mat')
r3 = seqpeth( s, seqs, smat );

load('/Volumes/slab1/mC400/dat/mC400_15/seqs_4.mat')
load('/Volumes/slab1/mC400/dat/mC400_15/smat_4.mat')
r4 = seqpeth( s, seqs, smat );
