% spikes_stats_supplement   add geometric dispersion and signed extremum to sst structure
%
% call                      sst = spikes_stats_supplement_bpi( sst )
%
% gets and returns          sst structure
%
% calls                     calc_spatial_waveform_features

% 30-mar-21 SSo


function sst = spikes_stats_supplement_bpi( sst )

% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( sst )
    return
end

% allocate space
n                               = size( sst.shankclu, 1 );
nsites_all                      = NaN( n, 1 );
for i                           = 1 : n
    nsites_all( i )             = size( sst.mean{ i }, 1 );
end
nsites                          = max( nsites_all );
z0                              = NaN( n, nsites );
vA                              = z0;
vT                              = z0;
bpi                             = NaN( n, 1 );
tV                              = z0;
vB                              = z0;
pA                              = z0;
nA                              = z0;
pT                              = z0;
nT                              = z0;

% go over units and compute stats
for i               = 1 : n
    w              = sst.mean{ i };
    [ vA0, vT0, bpi0, tV0, vB0, pA0, nA0, pT0, nT0 ] = calc_spatial_waveform_features( w, [], [], 0 );

    % update the structure
    cidx                        = 1 : length( vA0 );
    vA( i, cidx )               = vA0';
    vT( i, cidx )               = vT0';
    bpi( i )                    = bpi0;
    tV( i, cidx )              	= tV0';
    vB( i, cidx )              	= vB0';
    pA( i, cidx )             	= pA0';
    nA( i, cidx )            	= nA0';
    pT( i, cidx )           	= pT0';
    nT( i, cidx )             	= nT0';
    

end

% now expand the structure  with the new fields
sst.vA                          = vA;
sst.vT                          = vT;
sst.bpi                          = bpi;
sst.tV                          = tV;
sst.vB                          = vB;
sst.pA                          = pA;
sst.nA                          = nA;
sst.pT                          = pT;
sst.nT                          = nT;


return

% EOF


% use example #1 (look at a specific dataset, do not save anything):

% load sst for B and above units
[ shankclu, sst ]   = determine_units( filebase, [], 'B' );
[ imat, i1, i2 ]    = intersect( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' );
idx1                = false( size( sst.shankclu, 1 ), 1 );
idx1( i1 )          = 1;
sst                 = struct_select( sst, idx1 );

% add the new fields
sst2                = spikes_stats_supplement( sst );

% calculate some statistics...

% use example #2 (update the sst on the disk):
sstfname            = [ filebase '.sst' ];
load( sstfname, '-mat' )
sst                 = spikes_stats_supplement_bpi( sst );
save( sstfname, 'sst' )

% 30-mar-21: example run with addition of weighted-mean bpi
sst1                = spikes_stats_supplement_bpi( sst );
bpi                 = calc_com( sst1.vB', abs( sst1.vA )' )';
%bpi( all( isnan( sst1.vB ), 2 ) ) = NaN; % some are zero because there are no waveforms with positive peak before negative trough
rmv = all( isnan( sst1.vB ) | isnan( sst1.vA ), 2 ); % some are zero becuase there are no waveforms with postive peak before negative trough and extremum above TH
bpi( rmv ) = NaN;

figure
subplot( 2, 2, 1 )
plot( sst1.bpi, bpi, '.b' )
xlim( [ 0 1 ] ), ylim( [ 0 1 ] ), line( [ 0 1 ], [ 0 1 ], 'linestyle', '--', 'color', [ 0 0 0 ] );
xlabel( 'BPI (max)' ), ylabel( 'BPI (weighted' ), set( gca, 'tickdir', 'out', 'box', 'off' )

hold on, 
plot( 0, bpi( isnan( sst1.bpi ) & ~isnan( bpi ) ), '.b' )
plot( sst1.bpi( ~isnan( sst1.bpi ) & isnan( bpi ) ), zeros( sum( ~isnan( sst1.bpi ) & isnan( bpi ) ), 1 ), '.b' )

subplot( 2, 2, 2 )
hist( sst1.bpi( ~isnan( sst1.bpi ) ), 50 )


% sst1.max( :, ~isnan( sst1.bpi  ) )
% 
bidx = find( ~isnan( sst1.bpi  ) );

figure
for i = 67 : length( bidx )
    plot( sst1.max( :, bidx ( i ) ), 'b' )
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx ), sst1.filebase{ bidx( i ) }, sst1.shankclu( bidx( i ) , 1 ) ...
        , sst1.shankclu( bidx( i ) , 2 ), sst1.bpi( bidx( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    pause, 
end






