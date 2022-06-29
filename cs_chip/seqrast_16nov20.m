% seqrast
%
% CALL                      [ mat, bins, rh, ah ] = seqrast( clu, res, periods, varargin )
%
% GETS                      clu
%                           res 
%                           periods
%
% OPTIONAL ARGUMENTS        shankclu            specific clusters to take into analysis
%                           map                 
%                           binsizeSEC
%                           padSEC 
%                           spkFs 
%                           stimes 
%                           graphics 
%                           fAlpha
%                           eAlpha
%                           pColor
%                           emptyTrials         {1}, take into account trials with no spikes
%
% RETURNS                   mat              matrix of rasters of the units 
%                                            ( nbins, nunits * ntrials )
%                           bins             
%                           rh               raster handle
%                           ah               seperates the clus
%   
% DOES                      plots rasters specific periods of time
%
% CALLS                     get_rasters, plot_raster (spikes)
%               
% see also                  make_cs_chip_figures_seqplot, seqdetect, seqpeth, seqsort


% 04-nov-20 ES + SSo
%
% revisions:
% 10-nov-20                 added emptyTrials to varargin
% 15-nov-20                 added spike_length to allow plotting dots

function [ mat, tim, rh, ah ] = seqrast( clu, res, periods, varargin )

% initialize output
rh                                      = [];
ah                                      = [];

% default values
binsizeSEC_DEFAULT                      = 0.001;            % [s]
padSEC_DEFAULT                          = [ -0.05 0.05 ];   % [s]
spkFs_DEFAULT                           = 20000;            % [Hz]

% argument handling
nargs                           = nargin;
if nargs < 3 || isempty( clu ) || isempty( res ) || isempty( periods )
    return
end
[ shankclu, map ...
    , binsizeSEC, padSEC, spkFs ...
    , stimes, graphics ...
    , fAlpha, eAlpha, pColor, spike_length, emptyTrials ]  = ParseArgPairs(...
    { 'shankclu', 'map' ...
    , 'binsizeSEC', 'padSEC', 'spkFs' ...
    , 'stimes', 'graphics' ...
    , 'fAlpha', 'eAlpha', 'pColor', 'spike_length', 'emptyTrials' }...
    , { [], [] ...
    , binsizeSEC_DEFAULT, padSEC_DEFAULT, spkFs_DEFAULT ...
    , [], 1 ...
    , 0.1, 0.1, [ 0 0 0.7 ], [], 1 } ...
    , varargin{ : } );

% (1) get the rasters as sparse martices for all units in clu
binsize                         = ceil( binsizeSEC * spkFs );                       % [samples]
h                               = [ 0 max( diff( periods, [], 2 ) + 1 ) - 1 ];      % [samples]
h1                              = [ floor( h( 1 ) / binsize ) ceil( h( 2 ) / binsize ) ]; % [bins]
h2                              = [ floor( padSEC( 1 ) / binsizeSEC ) ceil( padSEC( 2 ) / binsizeSEC ) ]; % [bins]
halfwin                         = h1 + h2;
[ r, bins ]                     = get_rasters( clu, res, [], periods( :, 1 )...
    , 'binsize', binsize, 'halfwin', halfwin );
tim                             = bins / spkFs;

    
    
% (2) now, keep only the desired units and organize in the requested order
if ~isempty( shankclu )
    uidx                            = ismember( map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' );
    if isempty( uidx )
        error( '' )
    end
    r                               = r( uidx );
    
    [ ~, bb ]                       = sortrows( shankclu( :, 1 : 2 ) );
    [ ~, ridx ]                     = sort( bb );
    r                               = r( ridx );
end



% (3) combine into one array
ntrials                         = size( periods, 1 );
nbins                           = size( r{ 1 }, 1 );
nunits                          = length( r );
if ~isequal( size( r{ 1 }, 2 ), ntrials )
    error( 'mismatch' )
end
if ~isequal( size( r{ 1 }, 1 ), nbins )
    error( 'mismatch' )
end
mat                             = sparse( nbins, nunits * ntrials );
for i                           = 1 : nunits
    cidx                        = ntrials * ( i - 1 ) + 1 : ntrials * i;
    mat( :, cidx )              = r{ i };
end

if ~emptyTrials
    j=0;
    f=0;
    l=1;
    trail_kept = [];
    mat_full_final = [];
    mat_full = [];
    mat1 = [];
    matf = [];
    mat0 = [];
    
    mat_full = full(mat);
    for i = 1:nunits
        mati(i).units = mat_full(:,l:l+ntrials-1);
        l = l+ntrials;
    end
    for k =1:nunits
    mat0 = mati(k);
    j=0;
    for i = 1:ntrials
        idxmat = mat0.units(:,i) == 0;
        idxzero = sum(idxmat) ==length(idxmat); %if 0 than there was a spike
        if ~idxzero
            j= j+1;
            trail_kept(i,k) = 1;
%             mat1(:,j) = mat0.units(:,i);
        end
    end
    matf(k).units = mat0.units;
    end
    trail_full = sum(trail_kept,2) == nunits;
    if sum(trail_full)<15
        trail_full = sum(trail_kept,2) >= nunits-1;
    end
    % now we choose 20 trials in random out of the trials for each unit
    num_of_trials = sum(trail_full); % to be changed by Simcky
    l=1;
    for k = 1:nunits
%         aa = size(matf(k).units,2);
%         vecaa = 1:aa;
%         iwant = vecaa(randperm(aa,num_of_trials));
%         iget = sort(iwant);
%         mat_full_final(:,l:l+num_of_trials-1) = matf(k).units(:,iget);
%         l= l+num_of_trials;
            mat_full_final(:,l:l+num_of_trials-1) = matf(k).units(:,trail_full);
            l= l+num_of_trials;
    end
end


% (4) plot
if ~graphics
    return
end
newplot
% plot the raster
%rh                              = plot_raster( mat, tim );
if ~emptyTrials
    rh                              = plot_raster( mat_full_final, tim, [], spike_length );
ah                              = alines( num_of_trials * ( 1 : ( nunits - 1 ) ) + 0.5, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' )
else
 rh                              = plot_raster( mat, tim, [], spike_length );   
ah                              = alines( ntrials * ( 1 : ( nunits - 1 ) ) + 0.5, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
end
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'Time [s]' )
ylabel( 'Trial number x unit' )

% add patches
if ~isempty( stimes )
    stimes                      = permute( stimes( 1, :, : ), [ 2 3 1 ] );
    nchans                      = size( stimes, 1 );
    for i                       = 1 : nchans
        xe                    	= stimes( i, : );
        ye                    	= ylim;
        ph                      = patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), pColor );
        set( ph, 'FaceAlpha', fAlpha, 'EdgeAlpha', eAlpha )
    end
end


if ~emptyTrials
if ~isempty( shankclu )
    ye                          = num_of_trials / 2 : num_of_trials : num_of_trials * ( nunits - 0.5 );
    unums                       = shankclu( :, 1 : 2 );
    ustr                        = cell( 1, nunits );
    for i                       = 1 : nunits
        ustr{ i }               = sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
    end
    set( gca, 'ytick', ye, 'YTickLabel', ustr )
end
else
% add ustr
if ~isempty( shankclu )
    ye                          = ntrials / 2 : ntrials : ntrials * ( nunits - 0.5 );
    unums                       = shankclu( :, 1 : 2 );
    ustr                        = cell( 1, nunits );
    for i                       = 1 : nunits
        ustr{ i }               = sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
    end
    set( gca, 'ytick', ye, 'YTickLabel', ustr )
end
end
return

% EOF
