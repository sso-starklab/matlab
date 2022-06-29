% make_cs_chip_figures_seqplot
%
% calls                     load_spikes, multipeth_make, inrange,
%                           poissonTest, mygray, alines

% 08-oct-20 SSo & ES

function make_cs_chip_figures_seqplot( filebase, seqs, smat )

% 08-oct-20: PETH plotting
alfa                            = 0.05;
spkFs                           = 20000;
ilev                            = 'B';

%filebase                        = filebaseLookup( 'mC400', -12 );
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

return

% EOF

