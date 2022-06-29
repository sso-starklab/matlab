% spikes2spikes_asg         add ASG and STC to s2s structure
%
% call                      s2s = spikes2spikes_asg( filebase )
%
% gets                      filebase
% 
% returns                   s2s structure
%
% optional arguments (given as name/value pairs):
%
%               CCH computation parameters:
%
%                           isDeadTime          {1}         assume deadtime in same-shank detection
%                           deadTimeMS          {0.4}       [ms], ignored if isDeadTime is 0
%                           BinSizeMS           {1}         [ms], CCH bin size 
%                           halfWidthMS         {50}        [ms], CCH half-width 
%                           jitWindowMS         {5}         [ms], filtering half-width
%                           convType            {'gauss'}   argument to cch_conv
%                           roiMS               {[NaN 5]}   [ms], for significance testing 
%                           alfa                {0.001}     significance testing (global bands)
%                           supportEdges        {[-1 3]}    jitWindowMS multiples i.e. keep
%                                                               h(t) only for the bins -5 : 15 
%
%               Flow control:
%                           verbose             {1}
%                           suffix              {'s2s'}
%                           Overwrite           {-2}        1: compute and overwrite 
%                                                           0: only compute
%                                                           -1: load/compute but do not write
%                                                           -2: load/compute and write
%               Data selection parameters:
%
%                           shanknums           {[]}
%                           periods             {-1}
%
% calls                     ParseArgPairs, load_spikes, LoadXml
%                           LoadStims, uniteranges, 
%                           get_hfo_times, sortranges, resampleranges, inranges
%                           calc_asg
%
% see also                  spikes2spikes.m

% 03-mar-21 ES

% to do:
% (1) check with a few sessions
% (2) speed up using local computations

function s2s = spikes2spikes_asg( filebase, varargin )

% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ isDeadTime ...
    , deadTimeMS, BinSizeMS, halfWidthMS, jitWindowMS, roiMS ...
    , convType, alfa ...
    , shanknums, ilevel, periods ...
    , verbose, suffix, Overwrite ...
    , supportEdges ]            = ParseArgPairs( ...
    { 'isDeadTime' ...
    , 'deadTimeMS', 'BinSizeMS', 'halfWidthMS', 'jitWindowMS', 'roiMS' ...
    , 'convType', 'alfa' ...
    , 'shanknums', 'ilevel', 'periods' ...
    , 'verbose', 'suffix', 'Overwrite' ...
    , 'supportEdges' } ...
    , { 1 ...
    , 0.4, 1, 50, 5, [ NaN 5 ] ...
    , 'gauss', 0.001 ...
    , [], 'E', -1 ...
    , 1, '', -2 ...
    , [ -1 3 ] }...
    , varargin{ : } );

% i/o
if isempty( suffix )
    filename                 	= [ filebase '.s2s' ];
else
    filename                	= [ filebase '.' suffix ];
end
if exist( filename, 'file' ) && Overwrite < 0 && nargout > 0
    if verbose
        fprintf( 1, 'Loading %s\n', filename )
    end
    load( filename, '-mat', 's2s' );
    return
end

%------------------------------------------------------------------------
% (1) load all the spikes
spk                             = load_spikes( filebase, shanknums, ilevel );  % use
par                             = LoadXml( filebase );
SpikesFs                        = par.SampleRate;

% keep only the relevant spikes
if periods( 1 ) < 0
    
    % remove spikes during stimulation
    [ Vals, Trigs ]             = LoadStims( filebase );
    if periods( 1 ) == -1 || periods( 1 ) == -3                             % any trigger
        tidx                    = true( size( Trigs ) );
    elseif periods( 2 ) == -2                                               % specific subset
        tidx                    = ismember( Trigs, abs( periods ) );
    else
        tidx                    = [];
    end
    if isempty( tidx )
        uvals1                  = [];
    else
        uvals1                  = uniteranges( Vals( tidx, 1 : 2 ) );       % combine all the segments
    end
    if periods( 1 ) == -2 || periods( 1 ) == -3                             % remove HFO times
        hperiods                = get_hfo_times( filebase );
        hperiods                = sortranges( hperiods );
        eegFs                   = par.lfpSampleRate;
        uvals2                  = resampleranges( hperiods, SpikesFs, eegFs );
    else
        uvals2                  = [];
    end
    uvals                       = uniteranges( uvals1, uvals2 );
    if ~isempty( uvals )
        ridx                    = inranges( spk.res, uvals );               % remove any spike that is in any segment
        spk.res( ridx )         = [];
        spk.clu( ridx )         = [];
    end
    Bsec                        = sum( diff( uvals, [], 2 ) + 1 ) / SpikesFs;
    Rsec                        = spk.res( end ) / SpikesFs;
    Tsec                        = Rsec - Bsec;

elseif ~isempty( periods ) && size( periods, 2 ) == 2
    
    % keep only spikes during specified periods
    kidx                        = inranges( spk.res, periods );
    spk.res                     = spk.res( kidx );
    spk.clu                     = spk.clu( kidx );
    Tsec                        = sum( diff( periods, [], 2 ) + 1 ) / SpikesFs;
    
end

%------------------------------------------------------------------------
% (2) prepare fields for s2s structure

% compute CCH time lag vector
nBins                           = halfWidthMS / BinSizeMS;               	% number of bins on each side of the CCH
t                               = ( -nBins : nBins )' * BinSizeMS;      	% [ms]
idx                             = ( supportEdges( 1 ) * jitWindowMS ) : ( supportEdges( 2 ) * jitWindowMS );
sidx                            = find( ismember( t, idx ) );
nsupport                        = length( sidx );

% allocate space
shankclu                        = spk.shankclu( :, 1 : 2 );
n                               = size( shankclu, 1 );
m                               = length( t );
nspks                           = NaN( n, 1 );
cch                             = NaN( m, n, n );
pred                            = NaN( m, n, n );
pvalsUpper                      = NaN( m, n, n );
pvalsLower                      = NaN( m, n, n );
cc_act                          = false( n, n );
cc_sil                          = false( n, n );
g1mat                           = NaN( n, n );
g2mat                           = NaN( n, n );
g1tidx                          = NaN( nsupport, n, n );
g1gcch                          = NaN( nsupport, n, n );
g2tidx                          = NaN( nsupport, n, n );
g2gcch                          = NaN( nsupport, n, n );

%------------------------------------------------------------------------
% (3) go over pairs and compute CCH, pred, pvals, STCs, and ASGs
for n1                          = 1 : n
    
    % select the spike train of n1
    shankclu1                   = shankclu( n1, : );
    clunum1                     = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
    st1                         = spk.res( spk.clu == clunum1 );
    if verbose
        fprintf( 1, '%d ', n1 )
    end
    if isempty( st1 )
        continue
    end
    
    for n2                      = 1 : n
        
        if verbose
            fprintf( 1, '.' )
        end
        if n1 == n2
            continue
        end
        
        % select the spike train of n2
        shankclu2               = shankclu( n2, : );
        clunum2                 = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
        st2                     = spk.res( spk.clu == clunum2 );
        if isempty( st2 )
            continue
        end

        % determine dead time
        if shankclu1( :, 1 ) == shankclu2( :, 1 ) && isDeadTime
            dtMS                = deadTimeMS;
        else
            dtMS                = 0;
        end
        
        % compute CCH and ASG
        [ g1, g2, act, sil, s ] = calc_asg( st1, st2 ...
            , 'BinSizeMS', BinSizeMS, 'halfWidthMS', halfWidthMS ...
            , 'SpikesFs', SpikesFs, 'jitWindowMS', jitWindowMS ...
            , 'roiMS', roiMS, 'deadTimeMS', dtMS ...
            , 'convType', convType, 'alpha0', alfa ... 
            , 'graphics', 0 );        
        
        cc_act( n1, n2 )        = act;
        cc_sil( n1, n2 )        = sil;
        
        % collect CCH, pred, and pvals
        nspks( n1, : )          = s.nspks1;
        nspks( n2, : )          = s.nspks2;
        pvalsUpper( :, n1, n2 ) = s.pvals;
        pvalsLower( :, n1, n2 ) = 1 - s.pvals;
        if n2 == ( n1 + 1 )
            cch( :, n1, n1 )    = s.ach1;
            cch( :, n2, n2 )    = s.ach2;
        end
        cch( :, n1, n2 )        = s.cch;
        pred( :, n1, n2 )       = s.pred;
        
        % collect ASG
        g1mat( n1, n2 )         = g1;
        g2mat( n1, n2 )         = g2;
        
        % determine the support of the ASGs
        idx1t                   = ismember( sidx, s.g1base );
        idx2t                   = ismember( sidx, s.g2base );
        [ ~, idx1s ]            = intersect( s.g1base, sidx );
        [ ~, idx2s ]            = intersect( s.g2base, sidx );
        
        % assign the STC for ASGe and ASGi
        g1tidx( idx1t, n1, n2 ) = s.g1base( idx1s );
        g1gcch( idx1t, n1, n2 ) = s.gcch( s.g1base( idx1s ) );
        g2tidx( idx2t, n1, n2 ) = s.g2base( idx2s );
        g2gcch( idx2t, n1, n2 ) = s.gcch( s.g2base( idx2s ) );
        
    end         % n2
    
    if verbose
        fprintf( 1, '\n' )
    end
    
end             % n1

%------------------------------------------------------------------------
% (4) compute global bands and detect high/low bins

cch1                            = reshape( cch, [ m n * n ] );
pred1                           = reshape( pred, [ m n * n ] );
roiMS( isnan( roiMS ) )         = 0;
t_ROI                           = t >= roiMS(1) & t <= roiMS(2);
nBonf                           = sum( t_ROI );

% global limits, with a conservative Bonferroni correction, on the tested range:
gbUpper                         = poissinv( 1 - alfa / nBonf, max( pred1( t_ROI, : ), [], 1 ) );
gbLower                         = poissinv( alfa / nBonf, min( pred1( t_ROI, : ), [], 1 ) );

% detect the bins with high/low global counts
hiBins                          = false( size( cch1 ) );
loBins                          = false( size( cch1 ) );
hiBins( t_ROI, : )              = bsxfun( @gt, cch1( t_ROI, : ), gbUpper ) & ( cch1( t_ROI, : ) > 0 );
loBins( t_ROI, : )              = bsxfun( @lt, cch1( t_ROI, : ), gbLower ) & ( ones( nBonf, 1 ) * gbLower ) > 0;

% reshape back
pvalsUpper                      = reshape( pvalsUpper, [ m n n ] );
pvalsLower                      = reshape( pvalsLower, [ m n n ] );
hiBins                          = reshape( hiBins, [ m n n ] );
loBins                          = reshape( loBins, [ m n n ] );
gbUpper                         = reshape( gbUpper, [ n n ] );
gbLower                         = reshape( gbLower, [ n n ] );

%------------------------------------------------------------------------
% (5) assign output

% 'standard' fields (conforms to old spikes2spikes output)
s2s.filebase                    = filebase;
s2s.shankclu                    = spk.shankclu( :, 1 : 2 );
s2s.nspks0                      = nspks;
s2s.nspks                       = nspks;
s2s.Tsec                        = Tsec;
s2s.ccg                         = cch;
s2s.t                           = t;
s2s.State                       = [];
s2s.cutoff                      = 0;
s2s.window                      = jitWindowMS;
s2s.convtype                    = convType;
s2s.alpha                       = alfa;
s2s.t_ROI                       = t_ROI;
s2s.periods                     = periods;
s2s.pred                        = pred;
s2s.pvalsUpper                  = pvalsUpper;
s2s.pvalsLower                  = pvalsLower;
s2s.hiBins                      = hiBins;
s2s.loBins                      = loBins;
s2s.gbUpper                     = gbUpper;
s2s.gbLower                     = gbLower;

% STC and ASG fields (new)
s2s.g1mat                       = g1mat;
s2s.g2mat                       = g2mat;
s2s.g1tidx                      = g1tidx;
s2s.g1gcch                      = g1gcch;
s2s.g2tidx                      = g2tidx;
s2s.g2gcch                      = g2gcch;
s2s.act                         = cc_act;
s2s.sil                         = cc_sil;

% save to disk
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( filename, 'file' ) )
    if verbose
        fprintf( 1, 'Saving %s\n', filename )
    end
    save( filename, 's2s', '-v6' );
end

return

% EOF
