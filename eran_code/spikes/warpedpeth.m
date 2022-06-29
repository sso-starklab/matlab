% warpedpeth            by duration and center
% 
% call                  [ stats fig ] = warpedpeth( Clu, Res, trigs, periods )
%
% gets                  Clu, Res, trigsm periods
%
% optional arguments (given as name/value pairs):
%                       nbinsCenter   {4}       (will be 9)
%                       nbinsPre      {12}      (will also be 12 post)
%                       Fs            {20000}, [Hz]
%                       sts           {[]}      baseline periods (@Fs)
%                       vals          {[]}      excluded baseline periods (@Fs)
%                       graphics      {0}
%
% returns               stats
%                       fig
%
% calls                 ParseArgPairs, sortcols             (general)
%                       imagescbar, alines, myjet           (graph) 
%                       geteventsinranges, setdiffranges    (sets)
%                       get_rasters                         (spikes)
%                       firfilt                             (ssp)
%                       calc_index                          (stats)

% 23-dec-13 ES

% revisions
% 22-oct-14 modified to match nclu in case some units do not spike in any period
% 31-oct-14 modified to compute surprise
% 03-nov-14 modified to return dummy structure if no intra-period spikes
% 23-mar-21 cleaned up

%-----------------------------------------------------------------%
% compute ripple-warped PSTH:
% the raw data:
% rast: nbins x nevents (nbins)
% periods: nevents x 2 (onset/offset)
%   ripsNew.edges - spec.periods{ 1 }
% ripple time for each event (ripsNew.trigs)
% the sampling frequency for all here is the same (eegFs)

% arguments:
% trigs, edges, sts, vals : all @ spkFs!!
% spkFs
% s: shankclu, clu, res, map

function [ stats, fig ] = warpedpeth( Clu, Res, trigs, periods, varargin )

% initialize
stats                           = [];
fig                             = [];

% arguments
nargs                           = nargin;
if nargs < 4 || isempty( Clu ) || isempty( Res ) || isempty( trigs ) || isempty( periods )
    return
end
if ~isequal( size( Clu ), size( Res ) )
    return
end
if ~isequal( size( trigs, 1 ), size( periods, 1 ) )
    return
end
[ nbinsCenter, nbinsPre, Fs, map, sts, vals, graphics ...
    ]                           = ParseArgPairs(...
    { 'nbinsCenter', 'nbinsPre', 'Fs', 'map', 'sts', 'vals', 'graphics',...
    }...
    , { 4, 12, 20000, [], [], [], 0 ...
    }...
    , varargin{ : } );
if ~isempty( vals ) && size( vals, 2 ) ~= 2
    return
end
if ~isempty( sts ) && size( sts, 2 ) ~= 2
    return
end

% get parameters
if isempty( map )
    map                         = unique( Clu );
end
ridx                            = ~ismember( Clu, map );
Clu( ridx )                     = [];
Res( ridx )                     = [];
nclu                            = length( map );
halfwin                         = nbinsCenter + nbinsPre;
nbins                           = 2 * halfwin + 1;
durs                            = diff( periods, 1, 2 ) + 1;
durPre                          = trigs - periods( :, 1 );                  % @ spkFs
durPost                         = periods( :, 2 ) - trigs;                  % @ eegFs -> spkFs
nevents                         = length( durs );
% focus on relevant data
mat                             = [ -durPre durPost ] * ( nbinsPre + nbinsCenter ) / nbinsCenter + trigs * [ 1 1 ];
[ out, idx ]                    = geteventsinranges( Res, mat, 0 );
res                             = out;
clu                             = Clu( idx );
if isempty( clu )
    binsHat                     = -( nbinsPre + nbinsCenter ) / nbinsCenter : 1 / nbinsCenter : ( nbinsPre + nbinsCenter ) / nbinsCenter; % [wbins]
    wbins                       = binsHat * ( mean( diff( periods, 1, 2 ) + 1 ) / Fs );
    nans                        = NaN * ones( length( map ), 1 );
    stats.map                   = map;
    stats.wgain                 = NaN * ones( length( map ), length( wbins ) );
    stats.wbins                 = wbins;
    stats.phi                   = nans;
    stats.plo                   = nans;
    stats.mi                    = nans;
    stats.gain                  = nans;
    stats.surp                  = nans;
    return
end

% compute
rr                              = cell( nclu, 1 );
uClu                            = unique( Clu );
bbins                           = zeros( 2 * halfwin + 1, nevents );
binsizes                        = bbins;
for i                           = 1 : nevents
    binsize0                    = round( durPre( i ) / nbinsCenter );       % centeral bin will be as pre
    binsize1                    = round( durPost( i ) / nbinsCenter );
    [ rr0, bbins0 ]             = get_rasters( clu, res, 1, trigs( i ), 'binsize', binsize0, 'halfwin', halfwin, 'clunums', uClu );
    [ rr1, bbins1 ]             = get_rasters( clu, res, 1, trigs( i ), 'binsize', binsize1, 'halfwin', halfwin, 'clunums', uClu );
    % rearrange the rasters
    for j                       = 1 : nclu                                  % @ spks/bin
        jin                     = find( uClu == map( j ) );
        if isempty( jin )
            continue
        end
        rr{ j }( :, i )         = [ rr0{ jin }( 1 : halfwin + 1 ); rr1{ jin }( ( halfwin + 2 ) : ( 2 * halfwin + 1 ) ) ];
    end
    bbins( :, i )               = [ bbins0( 1 : halfwin + 1 ); bbins1( ( halfwin + 2 ) : ( 2 * halfwin + 1 ) ) ];
    binsizes( :, i )            = [ binsize0 * ones( halfwin + 1, 1 ); binsize1 * ones( halfwin, 1 ) ]; % samples/bin
end
wbins                           = mean( bbins, 2 ) / Fs;                    % time [s]
binsHat                         = -( nbinsPre + nbinsCenter ) / nbinsCenter : 1 / nbinsCenter : ( nbinsPre + nbinsCenter ) / nbinsCenter; % [wbins]

% accumulate histograms
rsum                            = NaN * ones( nbins, nclu );
for j = 1 : nclu
    if isempty( rr{ j } )
        continue
    end
    rsum( :, j )                = sum( rr{ j }, 2 );                        % spks/bin
end
wpeth                           = bsxfun( @rdivide, rsum, sum( binsizes, 2 ) / Fs ); % convert to spks/s

% also get the baseline rate (e.g. during SWS, no ripples, no light):
swsnorip                        = setdiffranges( sts, intersectranges( periods, vals ) ); % @spkFs
[ ~, idx ]                      = geteventsinranges( Res, swsnorip, 0 );
cluhat                          = Clu( idx );
maxClu                          = max( map );
cidx                            = ismember( 1 : maxClu, map );
count                           = accumarray( cluhat, 1, [ maxClu 1 ] );
count                           = count( cidx );
Tbase                           = sum( diff( swsnorip, 1, 2 ) + 1 ) / Fs;
baserates                       = count / Tbase;

% convert to gain histograms:
wgain                           = bsxfun( @rdivide, wpeth, baserates' );

% stats for the PSTH:
tout                            = Tbase;
tin                             = sum( sum( binsizes( ( nbinsPre + 1 ) : ( nbins - nbinsPre ), : ) ) ) / Fs;
nin                             = sum( rsum( ( nbinsPre + 1 ) : ( nbins - nbinsPre ), : ) )';
nout                            = count( : );
lambda                          = nout / tout * tin;                        % expected count (total over entire tin)
pAct                            = 1 - poisscdf( nin - 1, lambda );          % prob to get %% NIN OR MORE %%
pSup                          	= poisscdf( nin, lambda );
rin                             = nin ./ tin;
rout                            = nout ./ tout;
mi                              = calc_index( rin, rout );
mgain                           = mean( wgain( ( nbinsPre + 1 ) : ( nbins - nbinsPre ), : ) )';

% summarize
stats.map                       = map;
stats.wgain                     = wgain';
stats.wbins                     = wbins';
stats.phi                       = pAct;
stats.plo                       = pSup;
stats.mi                        = mi;
stats.gain                      = mgain;
stats.surp                      = log10( ( pSup + eps ) ./ ( pAct + eps ) ); 


% plot
if ~graphics
    return
end

fig                             = figure;
mat                             = wgain;
mat                             = firfilt( mat, [ 0.25 0.5 0.25 ] ); % 3-bin smoother
mathat                          = sortcols( mat );
[ ~, ah ]                       = imagescbar( wbins, 1 : size( mathat, 2 ), mathat );
subplot( ah( 1 ) ),
ylabel( 'Cell #' )
xlabel( 'Time [s]' )
onoff                       = wbins( abs( binsHat ) == 1 ).';
alines( [ onoff + diff( wbins( 1 : 2 ) ) * [ -1 1 ] / 2 0 ], 'x', 'color', [ 1 1 1 ], 'linestyle', '--' );
subplot( ah( 2 ) )
xlabel( 'Gain' )
alines( 1, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
colormap( myjet )

return

% EOF

