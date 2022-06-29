% classify_waveform_plot 
%
% CALL              [ fig shat ] = classify_waveform_plot( s, tstr, jitflag, softlim, reclass )
% 
% CALLS             classify_waveform, minmax, replacetok

% See also wnAnalysisCelltypeClassificationPlot

% 12-sep-13 ES

% revisions
% 30-aug-18 (1) handling of ph to allow compatibility with R2018a
%           (2) hpfUsed temporary workaround
% 17-aug-19 cleaned up
% 23-mar-21 edge conditions handled
% 20-apr-21 legend call update for R2020

function [ fig, shat ] = classify_waveform_plot( s, tstr, jitflag, softlim, reclass, circtagged, hpfUsed )

% constants
colors                          = [ 0 0 0.7; 1 0 0 ];
colorsLight                     = [ 1 0 1 ; 0 0.7 0 ];
XLIM                            = [ 0.1 0.9 ];
YLIM                            = [ 0.6 1.4 ];
pTH                             = 0.05;

% initialize output
fig                             = [];

% arguments
nargs                           = nargin;
if nargs < 2 || isempty( tstr )
    tstr                        = '';
end
if nargs < 3 || isempty( jitflag )
    jitflag                     = 1;
end
if nargs < 4 || isempty( softlim )
    softlim                     = 0;
end
if nargs < 5 || isempty( reclass )
    reclass                     = 1;
end
if nargs < 6 || isempty( circtagged )
    circtagged                  = 1;
end
if nargs < 7 || isempty( hpfUsed )
    hpfUsed                     = 0;
end

if isempty( s )
    return
end
nclu                            = size( s.shankclu, 1 );
if nclu == 0
    return
end

%-------------------------------------------------------------------------
% preparations
%-------------------------------------------------------------------------
% check fields
if all( ismember( { 'd', 't', 'inh', 'act', 'exc' }, fields( s ) ) )
elseif all( ismember( { 'tp', 'fmax', 'filebase' }, fields( s ) ) )
    d                           = s.tp;
    t                           = 1 ./ s.fmax * 1000;
    nul                         = false( size( d ) );
    if ismember( { 'opttag' }, fields( s ) )
        act                     = s.opttag;
    else
        act                     = nul;
    end
    if isempty( tstr )
        tstr                    = s.filebase;
    end
    if ismember( 'pyr', fields( s ) )
        pyr                     = s.pyr;
    else
        pyr                     = [];
    end
    s                           = struct( 'd', d, 't', t, 'act', act, 'inh', nul, 'exc', nul );
    if ~isempty( pyr )
        s.pyr                   = pyr;
    end
else
    error( 'input type mismatch' )
end

% jitter points to improve visualization
if jitflag
    dd                          = sort( diff( sort( s.d ) ) );
    dd                          = dd( dd > 0 );
    tt                          = sort( diff( sort( s.t ) ) );
    tt                          = tt( tt > 0 );
    if isempty( dd )
        jx                      = 0;
    else
        dd                      = dd( ceil(  length( dd ) * 0.05 ) );
        jx                      = dd / 2 * randn( nclu, 1 );
    end
    if isempty( tt )
        jy                      = 0;
    else
        tt                      = tt( ceil(  length( tt ) * 0.05 ) );
        jy                      = tt * randn( nclu, 1 );
    end
else
    jx                          = 0;
    jy                          = 0;
end

% reclassify if required
shat                            = s;
if reclass
    [ shat.pidx, pval, fsep ]   = classify_waveform( [ s.d s.t ], 0, hpfUsed );
    shat.pidx                   = double( shat.pidx );
end
if ~reclass
    if ismember( { 'pyr' }, fields( s ) )
        shat.pidx               = s.pyr;
    end
    pval                        = zeros( size( shat.pidx ) );
    fsep                        = [];
end
shat.pval                       = pval;
shat.d                          = shat.d + jx;
shat.t                          = shat.t + jy;

% tags
c0                              = shat.inh;
c2                              = shat.exc;
cact( :, 1 )                    = shat.act == 1;
cact( :, 2 )                    = shat.act == 2;
if sum( cact( :, 2 ) ) > 0
    colorsLight                 = flipud( colorsLight );
end
gt                              = c0 | c2;                                  % inh/exc
cls                             = zeros( size( shat.pidx ) );               % by waveform, if any classification
cls( shat.pidx == 1 )           = 2;                                        % excitatory
cls( shat.pidx == 0 | shat.pidx == 2 ) = 1;                                 % inhibitory
cls( pval > pTH )               = 0;                                        % not classified
cls( isnan( shat.pidx ) )       = NaN;
x                               = [ shat.d shat.t ];

if softlim
    xlims                       = minmax( [ x( :, 1 ); XLIM( : ) ] );
    ylims                       = minmax( [ x( :, 2 ); YLIM( : ) ] );
else
    xlims                       = XLIM;
    ylims                       = YLIM;
end

%-------------------------------------------------------------------------
% plot
%-------------------------------------------------------------------------
if ~isempty( fsep )
    fh                          = fimplicit( fsep, [ xlims ylims ] );
    set( fh, 'color', [ 0 0 0 ] );
end
hold on

% narrow (not mono)
if sum( cls == 1 & ~gt  )
    ph( 5 )                     = plot( x( cls == 1 & ~gt, 1 ), x( cls == 1 & ~gt, 2 ), '.' ); 
    set( ph( 5 ), 'color', colors( 1, : ) );
    set( ph( 5 ), 'marker', 'o', 'markersize', 4 )
else
    ph( 5 )                     = 0;
end

% wide (not mono)
if sum( cls == 2 & ~gt  )
    ph( 6 )                     = plot( x( cls == 2 & ~gt, 1 ), x( cls == 2 & ~gt, 2 ), '.' );
    set( ph( 6 ), 'color', colors( 2, : ) );
    set( ph( 6 ), 'marker', 'o', 'markersize', 4 )
else
    ph( 6 )                     = 0;
end

% not classified
if sum( cls == 0 & ~gt )
    ph( 4 )                     = plot( x( cls == 0 & ~gt, 1 ), x( cls == 0 & ~gt, 2 ), '.k' );
else
    ph( 4 )                     = 0;
end

% excitatory
if sum( c2 )
    ph( 2 )                     = plot( x( c2, 1 ), x( c2, 2 ), '.' ); 
    set( ph( 2 ), 'color', colors( 2, : ) );
else
    ph( 2 )                     = 0;
end

% inhibitory
if sum( c0 )
    ph( 3 )                     = plot( x( c0, 1 ), x( c0, 2 ), '.' ); 
    set( ph( 3 ), 'color', colors( 1, : ) );
else
    ph( 3 )                     = 0;
end

% light tagged
if sum( cact( :, 1 ) )
    ph( 1 )                     = plot( x( cact( :, 1 ), 1 ), x( cact( :, 1 ), 2 ), '.' ); 
    set( ph( 1 ), 'color', colorsLight( 1, : ) );
    if circtagged
        set( ph( 1 ), 'marker', 'o', 'markersize', 8 )
    end
else
    ph( 1 )                     = 0;
end
if sum( cact( :, 2 ) )
    ph( 7 )                     = plot( x( cact( :, 2 ), 1 ), x( cact( :, 2 ), 2 ), '.' ); 
    set( ph( 7 ), 'color', colorsLight( 2, : ) );
    if circtagged
        set( ph( 7 ), 'marker', 'o', 'markersize', 8 )
    end
else
    ph( 7 )                     = 0;
end

% now make legend:
ndots                           = sum( ph( 1 : 7 ) ~= 0 );
if ndots == 7
    legend( ph( [ 1 7 2 : 3 6 5 4 ]), { 'Light activated1', 'Light activated2', 'Excitatory', 'Inhibitory'...
        , 'Wide', 'Narrow', sprintf( 'Unclassified (%d)', sum( cls == 0 & ~gt ) ) }, 'Location', 'NorthEast' )
elseif ndots == 6 && ph( 4 ) == 0
    legend( ph( [ 1 7 2 : 3 6 5 ]), { 'Light activated', 'Excitatory', 'Inhibitory'...
        , 'Wide', 'Narrow' }, 'Location', 'NorthEast' )
elseif ndots == 4  && ph( 4 ) == 0 && ph( 1 ) == 0
    legend( ph( [ 2 3 6 5 ]), { 'Excitatory', 'Inhibitory'...
        , 'Wide', 'Narrow' }, 'Location', 'NorthEast' )
elseif ndots == 3  && ph( 4 ) == 0 && ph( 1 ) == 0 && ph( 3 ) == 0
    legend( ph( [ 2 6 5 ]), { 'Excitatory'...
        , 'Wide', 'Narrow' }, 'Location', 'NorthEast' )
elseif ndots == 3  && ph( 4 ) == 0 && ph( 2 ) == 0 && ph( 3 ) == 0
    legend( ph( [ 1 6 5 ]), { 'Light tagged'...
        , 'Wide', 'Narrow' }, 'Location', 'NorthEast' )
elseif ndots == 4  && ph( 4 ) == 0 && ph( 5 ) == 0 && ph( 6 ) == 0
    legend( ph( [ 1 7 2 3 ]), { 'Light activated1', 'Light activated2'...
        , 'Excitatory', 'Inhibitory' }, 'Location', 'NorthEast' )
end

xlim( xlims )
ylim( ylims ),
xlabel( 'Trough-to-peak [ms]' ),
ylabel( 'Duration [ms]' ),
set( gca, 'tickdir', 'out', 'box', 'off' )
axis square

th                              = text( 0.2, 1.2, 'INT' );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', colors( 1, : ) );
th                              = text( 0.8, 1.3, 'PYR' );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', colors( 2, : ) );

title( replacetok( tstr, '\_', '_' ) );

return

% EOF

