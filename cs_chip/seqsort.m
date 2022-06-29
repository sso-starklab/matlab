% seqsort
%
% call              [ shankclu, fpeth, fbins, sgain, pAct ] = seqsort( r1, r2 )
%
% cmp         'intersect', 'unite', 'xor'
%               intersect       r1 & r2
%               unite           r1 | r2
%               xor             r1 & ~r2
%
% calls                     load_spikes, multipeth_make, inrange,
%                           poissonTest, mygray, alines
%
% see also                  make_cs_chip_figures_seqplot, seqdetect

% 15-oct-20 SSo & ES

% revisions
% 18-oct-20 actually written from script
% 21-oct-20 expanded ustr also for all units
% 15-nov-20 added fbins to output arguments

% pipeline:
% (1) decide on filebase, decide on sequences (in terms of values, channels, order, and duration)
% (2) use get_triggers to load information for the relevant parameters
% (3) for each sequence: use seqdetect to generate seqs and smat
% (4) for each sequence: use seqpeth to compute peth, gain, and significance for each unit
% (5) decide on a pair of sequences, use seqsort to plot peths of the
% (e.g. intersected) units according to their responses during the first sequence

function [ shankclu, fpeth, fbins, sgain, pAct ] = seqsort( r1, r2, varargin )

% argument handling
nargs                           = nargin;
if nargs < 2 || isempty( r1 ) || isempty( r2 )
    return
end
[ cmp, alfa ...
    , graphics, nT_calc ...
    , fAlpha, eAlpha, pColor ]  = ParseArgPairs(...
    { 'cmp', 'alfa' ...
    , 'graphics', 'nT_calc' ...
    , 'fAlpha', 'eAlpha', 'pColor' }...
    , { 'intersect', 0.05 ...
    , 1, [ -1 2 ] ...
    , 0.3, 0.3, [ 0 0 0.7 ] } ...
    , varargin{ : } );
if ~ismember( cmp, { 'intersect', 'unite', 'xor' } )
    error( 'unsupported logic for cmp' )
end
if ~isa( r1, 'struct' ) || ~isa( r2, 'struct' ) || ~isequal( fields( r1 ), fields( r2 ) )
    error( 'r1 and r2 must be structures with identical fields' )
end
if ~isfield( r1, 'shankclu' ) || ~isfield( r2, 'shankclu' ) || ~isequal( r1.shankclu, r2.shankclu )
    error( 'r1 and r2 must both have identical shankclu fields' )
end

% get some general parameters
nunits                          = size( r1.shankclu, 1 );

% 1. arrange structures in a cell array
nr                              = 2;
r                               = cell( nr, 1 );
for ri                          = 1 : nr
    cmd                         = sprintf( 'r{ %d } = r%d;', ri, ri );
    eval( cmd )
end

% 2. determine which units were active in each structure
uidxs                           = false( nunits, nr );
for ri                          = 1 : nr
    nBonf                      	= r{ ri }.nchans( 1 );
    uidxs( :, ri )             	= sum( r{ ri }.pAct <= alfa / nBonf, 2 ) >= 1;
end

% 3. apply the comparison to the uidx
switch cmp
    case 'intersect'
        uidx                    = all( uidxs, 2 );
    case 'unite'
        uidx                    = any( uidxs, 2 );
    case 'xor'
        uidx                    = uidxs( :, 1 ) & ~all( uidxs( :, 2 : end ), 2 );
end

% keep only the relevant units
for ri                          = 1 : nr
    r{ ri }                  	= struct_select( r{ ri }, uidx );
end

% 4. sort units by the time of maximal gain of the first structure (r1)
[ ~, maxidx ]                   = max( r{ 1 }.gain, [], 2 );
[ ~, sidx ]                     = sort( maxidx ); 
shankclu                        = r{ 1 }.shankclu( sidx, : );
sgain                           = cell( nr, 1 );
pAct                            = cell( nr, 1 );
fpeth                           = cell( nr, 1 );
fbins                           = cell( nr, 1 );
for ri                          = 1 : nr
    sgain{ ri }                 = r{ ri }.gain( sidx, : );
    pAct{ ri }                  = r{ ri }.pAct( sidx, : );
    fpeth{ ri }                 = r{ ri }.fpeth( sidx, : );
    fbins{ ri }                 = r{ ri }.fbins( sidx, : );
end

if ~graphics
    return
end

% plot everything
figure, colormap( mygray )
for ri                          = 1 : nr
    
    rdur                        = mean( r{ ri }.totdur ./ r{ ri }.nseqs );
    stimes                      = permute( r{ ri }.stimes( 1, :, : ), [ 2 3 1 ] );
    nchans                      = r{ ri }.nchans( 1 );
    tims                        = fbins{ ri }( 1, : );
    peth                        = fpeth{ ri };
    
    subplot( 3, 2, ri )
    imagesc( tims, 1 : sum( uidx ), scale( peth' )' ), axis xy
    set( gca, 'box', 'off', 'tickdir', 'out' )
    alines( [ 0 rdur ], 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    xlim( rdur * nT_calc )
    for i                       = 1 : nchans
        xe                    	= stimes( i, : );
        ye                    	= ylim;
        ph                      = patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), pColor );
        set( ph, 'FaceAlpha', fAlpha, 'EdgeAlpha', eAlpha )
    end
    cidx                    = true (size(shankclu,1),1);
    unums                	= shankclu( cidx, 1 : 2 );
    ustr                  	= cell( 1, sum( cidx ) );
    for i                	= 1 : sum( cidx )
        ustr{ i }         	= sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
    end
    set( gca, 'ytick', 1 : sum( cidx ), 'YTickLabel', ustr )

    for ct                   	= 0 : 1
        subplot( 3, 2, nr * ( ct + 1 ) + ri )
        cidx                 	= shankclu( :, 3 ) == ct;
        unums                	= shankclu( cidx, 1 : 2 );
        ustr                  	= cell( 1, sum( cidx ) );
        for i                	= 1 : sum( cidx )
            ustr{ i }         	= sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
        end
        imagesc( tims, 1 : sum( cidx ), scale( peth( cidx, : )' )' )
        axis xy
        set( gca, 'ytick', 1 : sum( cidx ), 'YTickLabel', ustr )
        set( gca, 'box', 'off', 'tickdir', 'out' )
        alines( [ 0 rdur ], 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
        xlim( rdur * nT_calc )
        for i                   = 1 : nchans
            xe                  = stimes( i, : );
            ye                  = ylim;
            ph                	= patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), pColor );
            set( ph, 'FaceAlpha', fAlpha, 'EdgeAlpha', eAlpha )
        end
    end % ct

end % ri

return

% EOF

filebase                        = filebaseLookup( 'mC400', -15 );
ilev                            = 'B';
s                               = load_spikes( filebase, [], ilev );

load('/Volumes/slab1/mC400/dat/mC400_15/seqs_2.mat')
load('/Volumes/slab1/mC400/dat/mC400_15/smat_2.mat')
r2 = seqpeth( s, seqs, smat );

load('/Volumes/slab1/mC400/dat/mC400_15/seqs_3.mat')
load('/Volumes/slab1/mC400/dat/mC400_15/smat_3.mat')
r3 = seqpeth( s, seqs, smat );

seqsort( r2, r3, 'cmp', 'intersect' )



