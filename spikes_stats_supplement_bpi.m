% spikes_stats_supplement   add geometric dispersion and signed extremum to sst structure
%
% call                      sst = spikes_stats_supplement_bpi( sst )
%
% gets and returns          sst structure
%
% calls                     calc_spatial_waveform_features

% 30-mar-21 SSo

% revisions
% 28-jun-21 (1) changed tV to be unit level
%           (2) added upol (unit level metric)
% 13-sep-21 added csd (unit level metric)
% 15-feb-22  call calc_spatial_waveform with sd per unit, instead of unitary TH
% 19-jun-22 return abvTH from calc_spatial_waveform

function sst = spikes_stats_supplement_bpi( sst, TH, spkFs )

% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( sst )
    return
end
if nargs < 2 || isempty( TH )
    TH                          = 40;                                       % uV
end
if nargs < 3 || isempty( spkFs )
    spkFs                       = 20000;                                    % [Hz]
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
tV                              = NaN( n, 1 );
upol                            = NaN( n, 1 );
vB                              = z0;
pA                              = z0;
nA                              = z0;
pT                              = z0;
nT                              = z0;
csd                             = cell( n, 1 );
abvTHi                          = false( n, 1 );

% go over units and compute stats
for i                           = 1 : n
    
    w                           = sst.mean{ i };
    [ vA0, vT0, bpi0, tV0, vB0, pA0, nA0, pT0, nT0, upol0 ,tp3, rA3, abvTH] = calc_spatial_waveform_features_28jun22( w, TH{i,1}, spkFs, 0, 'THflag', 0);
    [ csd0 ]                                                 = csd_mat (w);
    % update the structure
    cidx                        = 1 : length( vA0 );
    vA( i, cidx )               = vA0';
    vT( i, cidx )               = vT0';
    bpi( i )                    = bpi0;
    tV( i )                     = tV0';
    upol( i )                   = upol0';
    vB( i, cidx )              	= vB0';
    pA( i, cidx )             	= pA0';
    nA( i, cidx )            	= nA0';
    pT( i, cidx )           	= pT0';
    nT( i, cidx )             	= nT0';
    csd {i}                     = num2cell(csd0);
    tpbpi(i)                    = tp3;
    rAbpi(i)                    = rA3;
    abvTHi (i)                  = abvTH;
end

% now expand the structure  with the new fields
sst.pA                          = pA;
sst.nA                          = nA;
sst.pT                          = pT;
sst.nT                          = nT;
sst.vA                          = vA;
sst.vT                          = vT;
sst.vB                          = vB;
sst.bpi                         = bpi;
sst.tV                          = tV;
sst.upol                        = upol;
sst.csd                         = csd;
sst.tpbpi                       = tpbpi;
sst.rAbpi                       = rAbpi;
sst.abvTH                       = abvTHi;

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
[ res, sst, ff ] = shirly_eps_analysis( [], 'onlygather', 1 );  % added 06-apr-21
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

% bidx = find( ~isnan( sst1.bpi  )& sst1.bpi >0.6);

figure
for i = 1 : length( bidx )
    subplot( 1, 2, 1 )
    plot( sst1.max( :, bidx ( i ) ), 'b' )
    alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( [ 8.5 25.5 ], 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 20 );
    
    subplot( 1, 2, 2 )
    plot( deslope( sst1.max( :, bidx ( i ) ), 1 ), 'b' )
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx ), sst1.filebase{ bidx( i ) }, sst1.shankclu( bidx( i ) , 1 ) ...
        , sst1.shankclu( bidx( i ) , 2 ), sst1.bpi( bidx( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( [ 8.5 25.5 ], 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 20 );
    pause, 
end


idx =~ismember(bidx, bidx5)
rm = find(idx)

figure
for i = 1 : length(rm)
    plot( sst1.max( :, bidx(rm ( i )) ), 'b' )
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx(rm) ), sst1.filebase{ bidx(rm( i )) }, sst1.shankclu( bidx(rm( i )) , 1 ) ...
        , sst1.shankclu( bidx(rm( i )) , 2 ), sst1.bpi( bidx(rm( i )  )) );
    title( replacetok( str, '\_', '_' ) )
    alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( [ 8.5 25.5 ], 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );
    pause, 
end




% work on example bidx element 336
num                             = bidx( 336 );
nunits                          = length( sst.filebase );
vec                             = false( nunits, 1 );
vec( num )                      = 1;
ssti                            = struct_select( sst, vec );

%------------------------------------------------------------------------
% 23-jun-21 upgrade the BPI algorithm s.t.:
%
% for each channel:
% (1) detect all local extrema
% (2) keep only valid extrema (absolute value above threshold, e.g. 50 uV)
% (3) if local maximum before local minimum, keep the max/min as channel-specific BPI
% (4) determine which local extremum is larger (for the channel), and keep
% that (and maybe the time sample) 
%
% at the unit level:
% (1) determine whether all extrema are positive (Punit), negative (Nunit), or both (Dunit)
% (2) determine the BPI according to the channel with the maximal p2p

[ res, sst, ff ] = shirly_eps_analysis( [], 'onlygather', 1 );  % added 06-apr-21
%sst1                = spikes_stats_supplement_bpi( sst );

uidx = ismember( sst.filebase, 'mC41_53' ) & ismember( sst.shankclu, [ 2 29 ], 'rows' );
sst1 = struct_select( sst, uidx );

sst2                = spikes_stats_supplement_bpi( sst1 );


% 28-jun-21 checks after upgrading completed:
[ res, sst, ff ] = shirly_eps_analysis( [], 'onlygather', 1 );  % added 06-apr-21

uhist( sst.upol( ~isnan( sst.upol ) ) ) % nunits, dunits, punits
sum( sst.upol( ~isnan( sst.upol ) ) == 0 ) % dunits

figure
plot( sort( sst.bpi( ~isnan( sst.bpi ) & ~isinf( sst.bpi ) ) ) )

figure
plot( sort( sst.bpi( ~isnan( sst.bpi ) )))