% seqdetect             detect sequence of events according to a rule
%
% call                  [ seqs, smat ] = seqdetect( mats, seqmat )
%
% gets                  mats        cell array of matrices of ranges.
%                                       one matrix per sequence element
%
%                       seqmat      matrix describing the rule. 
%                                       each row corresponds to one sequence
%                                       element. 
%                                   syntax of each row: [ mindur maxdur minlag maxlag ]
%                                   lags are from onset of element i to onset of element i+1
%
% returns               seqs        [ nseqs by 2 ] matrix of all ranges
%                                   each row is a range of one sequence
%
%                       smat        [ nseqs x nelements x 4 ] matrix
%                                   each row contains: 
%                                       -the range of one event
%                                       -event number
%                                       -the interval from the previous event
%
% example:              seqmat = [ 5 15 25 30; 10 20 NaN NaN ];
% 
%                       this corresponds to a 5-15 sample event in mats{ 1 }, followed by a 10-20
%                       sample event in mats{ 2 }. The mats{ 2 } event must lag the mats{ 1 }
%                       event by 25-30 samples.
%
% calls                 sortranges              (sets)
%                       tempmatch               (ssp)

% 04-oct-20 ES

function [ seqs, smat ] = seqdetect( mats, seqmat, dt )


%------------------------------------------------------------------
% check arguments
seqs                            = [];
smat                            = [];

%------------------------------------------------------------------
% check arguments
nargs                           = nargin;
if nargs < 2 || isempty( mats ) || isempty( seqmat )
    error( 'missing arguments' )
end
if nargs < 3 || isempty( dt )
    dt                          = 0;
end
if ~isa( mats, 'cell' )
    error( 'input type mismatch: mats must be cell array' )
end
mats                            = mats( : );
if ~isa( seqmat, 'numeric' ) || ~ismatrix( seqmat )
    error( 'seqmat must be a matric of numbers' )
end
ne                              = length( mats );
if size( seqmat, 1 ) ~= ne
    error( 'input size mismatch: seqmat must have same number of rows as the number of matrices in mats' )
end
for i                           = 1 : ne
    mat                         = mats{ i };
    if ~isa( mat, 'numeric' ) || ~ismatrix( mat ) || size( mat, 2 ) ~= 2 
        fprintf( 1, 'element %d of mats must be a matrix of ranges', i )
        return
    end
    if ~isequal( sortranges( mat ), mat )
        fprintf( 1, 'element %d of mats must be a matrix of non-overlapping ranges', i )
        return
    end
end

%------------------------------------------------------------------
% (1) prune elements of each matrix by durations
nev                             = zeros( ne, 1 );
for i                           = 1 : ne
    mat                         = mats{ i };
    durs                        = diff( mat, [], 2 ) + dt;
    drange                      = seqmat( i, 1 : 2 );
    didx                        = durs >= drange( 1 ) & durs <= drange( 2 );
    mats{ i }                   = mat( didx, : );
    nev( i )                    = sum( didx );
end

%------------------------------------------------------------------
% (2) combine all events
ei                              = cumsum( nev );
tev                             = ei( ne );
si                              = [ 1; ei( 1 : ne - 1 ) + 1 ];
emat                            = zeros( tev, 4 );
for i                           = 1 : ne
    ridx                        = si( i ) : ei( i );
    emat( ridx, 1 : 3 )         = [ mats{ i } ones( nev( i ), 1 ) * i ];
end
emat                            = sortrows( emat, 1 );
difs                            = [ 0; emat( 2 : tev, 1 ) - emat( 1 : tev - 1, 2 ) ];
emat( :, 4 )                    = difs;

%------------------------------------------------------------------
% (3) detect sequences, regardless of lags
tmpl                            = ( 1 : ne )';
c                               = tempmatch( tmpl, emat( :, 3 ) );
c                               = abs( c );
tmpj                            = tmpl( [ 1 : ne - 1 ne - 1 ] );
z1                              = ( tmpl - mean( tmpl ) ) / std( tmpl );
z2                              = ( tmpj - mean( tmpj ) ) / std( tmpj );
th                              = sum( z1 .* z2 ) / ( ne - 1 );
th                              = ( th + 1 ) / 2;
tidx                            = find( c >= th );
ns                              = length( tidx );
m1                              = tidx * ones( 1, ne );
m2                              = ones( ns, 1 ) * ( 0 : ne - 1 );
m                               = ( m1 + m2 )';
smat0                           = emat( m( : ), : );

%------------------------------------------------------------------
% (4) keep the sequences with the proper lags
dlags                           = [ 0 inf; seqmat( 1 : ne - 1, 3 : 4 ) ];
lmat                            = repmat( dlags, [ ns 1 ] );
tf                              = smat0( :, 4 ) >= lmat( :, 1 ) & smat0( :, 4 ) <= lmat( :, 2 );
tf                              = reshape( tf, [ ns ne ] );
tf                              = sum( tf, 2 ) == ne;
sidx                            = 1 + ( find( tf ) - 1 ) * ne;
nseqs                           = length( sidx );
m1                              = sidx * ones( 1, ne );
m2                              = ones( nseqs, 1 ) * ( 0 : ne - 1 );
m                               = ( m1 + m2 )';
smat                            = smat0( m( : ), : );
seqs                            = [ smat0( sidx, 1 ) smat0( sidx + ne - 1, 2 ) ];

return

% EOF

%------------------------------------------------------------------
% 04-oct-20 example run: 

% load the data (durations and onset times)
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs42.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs43.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs45.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs47.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims42.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims43.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims45.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims47.mat')

% convert durations to samples, and create range matrices
mat42                           = tims42 + [ zeros( length( durs42 ), 1 ) round( durs42 * 20000 ) - 1 ];
mat43                           = tims43 + [ zeros( length( durs43 ), 1 ) round( durs43 * 20000 ) - 1 ];
mat45                           = tims45 + [ zeros( length( durs45 ), 1 ) round( durs45 * 20000 ) - 1 ];
mat47                           = tims47 + [ zeros( length( durs47 ), 1 ) round( durs47 * 20000 ) - 1 ];

% combine into cell array at set the rule
mats                            = { mat42, mat45, mat43, mat47 };
%seqmat                          = [ 9 12 29 31; 9 12 29 31; 9 12 29 31; 9 12 NaN NaN ] * 20;
% 08-oct-20 correction
seqmat                          = [ 7 15 25 35; 7 15 25 35; 7 15 25 35; 7 15  25 35 ] * 20;

% call the routine
[ seqs, smat ]                  = seqdetect( mats, seqmat, 1 );

%------------------------------------------------------------------
% 08-oct-20 - still missed events, do not appear in mats{2} - 
% may be missed detection / classification

% 4 missed:
%    122727937   122728112           4         625
%    122898801   122898976           1       41761
 

%------------------------------------------------------------------
% 08-oct-20: PETH plotting
alfa                            = 0.05;
spkFs                           = 20000;

ilev                            = 'B';
filebase                        = filebaseLookup( 'mC400', -12 );
s                               = load_spikes( filebase, [], ilev );

% compute PSTH for entire sequence - to determine baseline rate for each unit
speriods                     	= seqs / spkFs;
nT                              = 4; % number of stim duration before/after trigger
[ bpeth, bbins, btrigs, btims, bdurs, shankclu ] = multipeth_make( filebase, 'data', s, 'stimTypes', 0, 'periods', speriods, 'nT', nT, 'toplot', 0 );
brange                          = mean( bdurs ) * [ -4 -0.5 ];
bidx                            = inrange( bbins, brange );
baseRate                       	= mean( bpeth( bidx, : ), 1 ); % base firing rates

% compute PSTH for each invdividual event within the sequence and evaluate gain and significance
elements                        = unique( smat( :, 3 ) );
nelements                       = length( elements );
nunits                          = size( bpeth, 2 );
pAct                            = NaN( nunits, nelements );
gain                            = NaN( nunits, nelements );
for i                           = 1 : nelements
    element                     = elements( i );
    periods                     = smat( smat( :, 3 ) == element, 1 : 2 ) / spkFs;
    [ ipeth, ibins, itrigs, itims, idurs, ishankclu ] = multipeth_make( filebase, 'data', s, 'stimTypes', 0, 'periods', periods, 'toplot', 0, 'sdGauss', 0 );
    irange                      = [ 0 mean( idurs ) ];
    bidx                        = inrange( ibins, irange );
    inTime                      = sum( idurs );
    binSize                     = diff( ibins( 1 : 2 ) );
    nTrigs                      = size( periods, 1 );
    inCount                     = sum( ipeth( bidx, : ) * nTrigs * binSize );
    [ pInc, pDec, surprise ]    = poissonTest( baseRate, inCount, inTime * ones( 1, nunits ) );
    pAct( :, i )                = pInc;
    gain( :, i )                = ( inCount / inTime ) ./ baseRate;
end


% % plot PETH only for the units that are sign. activated by a single light source
% pAct( sum( pAct <= alfa, 2 ) == 1, : ) <= alfa
% shanknums = shankclu( sum( pAct <= alfa, 2 ) == 1, : )
% multipeth_make( filebase, 'data', s, 'stimTypes', 0, 'shanknums', shanknums, 'periods', speriods );

% plot PETH for the units sorted by the time of maximal gain
sgain                           = gain;
[ maxval, maxidx ]              = max( sgain, [], 2 );
[ ~, sidx ]                     = sort( maxidx ); 
sgain                           = sgain( sidx, : ); 

% plot PETH for a subset of the units (only those sig. modulated by at
% least 1 light source), sorted by the time of maximal gain
[ fpeth, fbins, ftrigs, ftims, fdurs ] = multipeth_make( filebase, 'data', s, 'stimTypes', 0, 'periods', speriods, 'toplot', 0, 'nT', [ -1 2 ] );

uidx                            = sum( pAct <= alfa / 4, 2 ) >= 1;
shankcluG                       = shankclu( uidx, : );
sgainG                          = gain( uidx, : );
pActG                           = pAct( uidx, : );
fpethG                          = fpeth( :, uidx );
[ maxval, maxidx ]              = max( sgainG, [], 2 );
[ ~, sidx ]                     = sort( maxidx ); 
sgainG                          = sgainG( sidx, : ); 
pActG                           = pActG( sidx, : );
fpethG                          = fpethG( :, sidx );

% find stimulus pattern (s)
nblocks                         = size( smat, 1 ) / nelements;
t0                              = repmat( smat( 1 : nelements : ( nblocks * nelements ), 1 ), [ 1 nelements ] )';
d                               = smat( :, 1 : 2 ) - t0( : ) * ones( 1, 2 );
m1                              = ones( nelements, 1 ) * ( 1 : nelements : ( nblocks * nelements ) );
m2                              = ( 0 : 1 : ( nelements - 1 ) )' * ones( 1, nblocks );
m                               = m1 + m2;
d1                              = d( :, 1 );
d2                              = d( :, 2 );
stimes                          = [ mean( d1( m ), 2 ) mean( d2( m ), 2 ) ] / spkFs;


% plot everything
figure, colormap( mygray )
subplot( 2, 2, 1 ), imagesc( 1 : nelements, 1 : sum( uidx ), scale( sgain' )' ), axis xy
set( gca, 'box', 'off', 'tickdir', 'out' )

% subplot( 2, 2, 2 ), imagesc( scale( sgainG' )' ), axis xy
% set( gca, 'box', 'off', 'tickdir', 'out' )
% 
% subplot( 2, 2, 3 ), imagesc( bbins, 1 : sum( uidx ), scale( bpethG )' ), axis xy
% set( gca, 'box', 'off', 'tickdir', 'out' )
% alines( [ 0 mean( bdurs ) ], 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
% 
subplot( 2, 2, 2 )
imagesc( fbins, 1 : sum( uidx ), scale( fpethG )' ), axis xy
set( gca, 'box', 'off', 'tickdir', 'out' )
alines( [ 0 mean( fdurs ) ], 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
xlim( mean( fdurs ) * [ -1 2 ] )
for i = 1 : nelements
    xe                              = stimes( i, : );
    ye                              = ylim;
    ph                              = patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), [ 0 0 0.7 ] );
    set( ph, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3 )
end

for ct                          = 0 : 1
    subplot( 2, 2, 3 + ct )
    cidx                        = shankcluG( :, 3 ) == ct;
    unums                       = shankcluG( cidx, 1 : 2 );
    ustr                        = cell( 1, sum( cidx ) );
    for i                       = 1 : sum( cidx )
        ustr{ i }               = sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
    end
    imagesc( fbins, 1 : sum( cidx ), scale( fpethG( :, cidx ) )' )
    axis xy
    set( gca, 'ytick', 1 : sum( cidx ), 'YTickLabel', ustr )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    alines( [ 0 mean( fdurs ) ], 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    xlim( mean( fdurs ) * [ -1 2 ] )
    for i = 1 : nelements
        xe                              = stimes( i, : );
        ye                              = ylim;
        ph                              = patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), [ 0 0 0.7 ] );
        set( ph, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3 )
    end
end

