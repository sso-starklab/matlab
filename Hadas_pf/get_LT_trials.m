% get_LT_trials             partition linear track behavior into trials
%
% call:                     [ trials tmov ] = get_LT_trials( filebase )
% 
% gets - optional arguments as parameter/value pairs:
%
%                           Overwrite       {-2}                    Overwrite only if does not exist
%                           graphics        {1}
%                           savetype        {'png'}
%                           durRange        { 1 120 ]               for temporal pruning
%                           alphaLevel      {0.05}                  for position pruning
%                           pt_bits         { [ 9 10 11 12 ] }      PT mapping
%                           sol_bits        { [ 14 15 ] }           Solenoid mapping
%                           extractMode     { 'byPTs' }             alternatives include: bySOLs, bySPD, byPOS 
%                           useMov          { 1 }    
%
%                           LP              { 1 }                   [Hz]; cutoff for FIR used in bySPD/byPOS
%                           spdTH           { [ 1 20 ] };           [cm/s]; thresholds for parseSchmitt used in bySPD/byPOS
%                           trackEnds       { [ 30 280 ] }          [cm], used in bySPD/byPOS see below
%                           putCOMx         { [ 65 250 ] }          [cm], used in bySPD/byPOS see below
%                           runEnds         { [ 1 inf ] }           [1/mov.Fs], used in bySPD/byPOS see below
%                           maxExpansion    { 60 };                 [s], used in bySPD/byPOS (should be half the upper value of durRange); see below
%
% returns:                  trials, a 5-column matrix
%                           tmov, a 9-field structure
%
%   trials: 
%       [ start_time_whlFs end_time_whlFs start_side start_time_spkFs end_time_spkFs ]
%       start_side: right=1; left=2
%
%   tmov:
%       tnum                serial number of the trial associated with the sample
%       tidx                time index (in whlFs) in the file
%       ttyp                trial type (same as start_side)
%       pos, spd, dir, ang  mov parameters (xy position, spd and dir, orientation)
%       whl                 name of whl file
%       Fs                  whlFs
%
% calls:                    ParseArgPairs, mindist          (general)
%                           LoadXml                         (blab)
%                           get_stimchans, parseDigital     (formats)
%                           eeg2whl                         (movement)
%                           circ_mean                       (circs)
%                           makefir, firfilt, parseSchmitt, local_max (ssp)
%
% does:
% (1) determines the trials
%       -loads the digital events
%       -prunes events:
%               (1) simultaneous low-to-high transitions of two different PTs
%               (2) sequential low-to-high transitions of the same PT/Sol
%       -defines trials as sol1 then sol2 or vice versa
%       -prunes trials: 
%               (1) too short/long
%               (2) with reward too far from the platform (based on x-pos)
%     organizes in a matrix (trials)
% (2) obtains the movement parameters
%   for movement along the X-axis, load only 
%       x, speed, movement direction, and head orientation
%   for each trial determine the movement features
%   concatenate and arrange in a structure tmov
%
% NOTEs:
% (1) the mapping of the digital channels in room 540 is, as of Aug 19, 2018:
% 1     ROI1            MC
% 2     ROI2            MC
% 3     ROI3            MC
% 4     ROI4            MC
% 5     Sync            MC
% 6     Mains           MC
% 7     DigOut 6/7      DSP
% 8     DigOut 7/7      DSP
% 9     PT1 (right)     BC
% 10    PT2             BC
% 11    PT3             BC
% 12    PT4 (left)      BC
% 13    Manual          BC
% 14    S1 (right)      BC
% 15    S2 (left)       BC
% 16    NIU             BC  (constantly high)
%
% (2)  bySPD or byPOS are not rigorous algorithms, and are only
% intended to fine tune the initial guess of the user. Thus, before running
% this routine with any of these extractMode, the user must visualize the
% data using e.g. 
% 
% >> figure, subplot( 2, 1, 1 ), hist( mov.pos( :, 1 ), 50 ), subplot( 2, 1, 2 ), plot( mov.pos( :, 1 ), 'b' )
% 
% and use these data to define initial guesses in the form of:
%
%       -putCOMx        spatial guesses of center of platforms
%       -trackEnds      spatial limits of xpos to consider
%       -runEnds        temporal limits of LT session [mov.Fs]

% 02-sep-18 ES

% revisions
% 11-sep-18 genericalized
% 30-oct-18 (1) digChan added as a possible external input
%           (2) useMov flag added (and made default)
%           (3) initialize output
%           (4) byPTs method added (and made default)
% 18-nov-18 (1) logic for digChan amended
% 08-sep-19 (1) added extractMode, defaults to byPTs, alt. include bySOL and byMOV
%           (2) changed internal checks of byPTs
% 09-sep-19 (1) added bySPD and determination of borders
% 10-sep-19 (1) finalized bySPD and automatic borders refinement
%           (2) fixed graphics
%           (3) added runEnds in units of mov.Fs
%           (4) written byPOS
% 12-sep-19 (1) identified problem with sampling
% 03-oct-19 (1) modified sampling to remain at movFs throughout until the
%               end (to prevent overlap between end of one trial n and
%               beginning of trial n+1)
%           note: this results in two things: 
%           (i) all events are shifted forward in time by no more than 1/sf samples
%           (ii) the events in trials( :, 4 : 5 ) are no longer identical to those in output of parseDigital
% 27-dec-19 (1) byPOS/bySPD: modified the borders (instead of shifting by SD/2, shifted by 2*SD)
%           (2) byPOS/bySPD: modified the trial definition: instead of time
%           from border crossing-to-border crossing, defined a trial from
%           directional border crossing to directional border crossing
%           (3) byPTs: modified the trial definition according to the above
% 21-Jan-20 added check that mov contains data
% 03-Jun-20 HS - fixed legends problems in the first figure
% 13-jul-20 (1) byPOS: kept only xpos in range
%           (2) byPOS, bySPD: modified trial durations by expanding trials
%               towards COM. this is required since the detection is biased
%               towards the center of the track to minimize false
%               positives, but has the cost of shortening the effective
%               length of the track
%               the expansion can be done in space (to COM), in time (e.g.
%               1 second), or composite (e.g. to COM, but limit expansion 1
%               second). This is a parameter maxExpansion, defaulted to inf
%           (3) modified durRange to be [ 0.5 120 ] (some animals are really fast)
% 19-jul-20 (1) limited edge expansion to onset of next trial / end of
%               previous trial to prevent spill-over

% to do:
%       * If in PTs mode, try to expand periods up to the position of
%       limits (currently, using PTs mode means trials are cropped at the
%       PTs, which are not the end of the linear track)

function [ trials, tmov, tparams, mat ] = get_LT_trials( filebase, varargin )

%--------------------------------------------------------------------%
% initialize output
%--------------------------------------------------------------------%
trials                      = zeros( 0, 5 );
tparams                     = zeros( 0, 5 );
tmov.tnum                   = [];
tmov.tidx                   = [];
tmov.ttyp                   = [];
tmov.pos                    = [];
tmov.spd                    = [];
tmov.dir                    = [];
tmov.ang                    = [];
tmov.whl                    = [];
tmov.Fs                     = [];

%--------------------------------------------------------------------%
% arguments
%--------------------------------------------------------------------%
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ Overwrite, graphics, savetype, durRange...
    , alphaLevel, pt_bits, sol_bits, digChan ...
    , useMov, extractMode ...
    , LP, spdTH, maxExpansion ...
    , trackEnds, putCOMx, runEnds ] = ParseArgPairs(...
    { 'Overwrite', 'graphics', 'savetype', 'durRange' ...
    , 'alphaLevel', 'pt_bits', 'sol_bits', 'digChan' ...
    , 'useMov', 'extractMode' ...
    , 'LP', 'spdTH', 'maxExpansion' ...
    , 'trackEnds', 'putCOMx', 'runEnds' }...
    , { -2, 1, 'png', [ 0.5 120 ] ...
    , 0.05, [ 9 10 11 12 ], [ 14 15 ], [] ...
    , 1, 'byPTs' ...
    , 1, [ 1 20 ], 60 ...
    , [ 30 280 ], [ 65 250 ], [ 1 inf ] 
    } ...
    , varargin{ : } );

if length( pt_bits ) ~= 4 || any( pt_bits ~= round( pt_bits ) ) || any( pt_bits < 0 )
    error( 'must have 4 PT bits' )
end
if length( sol_bits ) ~= 2 || any( sol_bits ~= round( sol_bits ) ) || any( sol_bits < 0 )
    error( 'must have 2 SOL bits' )
end
if length( durRange ) ~= 2 || any( durRange < 0 )
    error( 'durRange must have 2 non-negative numbers [s]' )
end
if ~ismember( extractMode, { 'byPTs', 'bySOLs', 'bySPD', 'byPOS' } )
    error( 'unsupported extractMode' )
end
trackEnds                   = clipmat( trackEnds, [ 0 inf ] );

%--------------------------------------------------------------------%
% parameters
%--------------------------------------------------------------------%
pt1bit                      = pt_bits( 1 );
pt2bit                      = pt_bits( 2 );
pt3bit                      = pt_bits( 3 );
pt4bit                      = pt_bits( 4 );
sol1bit                     = sol_bits( 1 );
sol2bit                     = sol_bits( 2 );
minDur                      = durRange( 1 );
maxDur                      = durRange( 2 );

par                         = LoadXml( filebase );
spkFs                       = par.SampleRate;
whlFs                       = par.lfpSampleRate / 32;
sf                          = whlFs / spkFs;

% load digital events and position data
if ~isempty( digChan ) && ( digChan > par.nChannels || round( digChan ) ~= digChan || digChan < 1 )
    digChan                 = [];
end
if isempty( digChan )
    digChan                 = get_stimchans( par, [], 'digitalin' );
end
movfname                    = [ filebase '.mov.mat' ];
digfname                    = [ filebase '.dig.' num3str( digChan ) ];
[ ~, filename, extname ]    = fileparts( filebase ); 
filename                    = [ filename extname ];
if exist( digfname, 'file' )
    load( digfname, '-mat', 'mat' )
else
    mat                     = parseDigital( filebase, 'chan', digChan );
end
if useMov
    if exist( movfname, 'file' )
        load( movfname, '-mat', 'mov' )
    else
        mov                     = eeg2whl( filebase );
    end
    if all( diff(mov.pos)==0 )
        error('Mov file is empty');
    end
end

% output file handling
trlfname                    = sprintf( '%s.trl', filebase );
cpathi                      = strfind( filebase, 'dat' );
cpath                       = filebase( 1 : ( cpathi - 1 ) );
fpath                       = sprintf( '%sfigs/', cpath );
if ~exist( fpath, 'dir' )
    eval( sprintf( '!mkdir %s', fpath ) )
end
figfname                    = sprintf( '%s%s.trl', fpath, filename );

% load the trl file
if Overwrite < 0 && exist( trlfname, 'file' )
    fprintf( '%s: Loading *.trl.mat file %s\n', mfname, trlfname )
    load( trlfname, '-mat', 'trials', 'tmov' )
    return
end

%--------------------------------------------------------------------%
% Step 1. Divide into trials.
%--------------------------------------------------------------------%
% A trial is ideally defined as the time between crossing PT1 and PT4,
%   or at the reverse direction between PT4 and PT1. 
% When PT measurements are imprecise, we can operationally define
%   a trial as the time between rewards, i.e. Sol1 and Sol2

if ismember( extractMode, { 'bySPD', 'byPOS' } )
    
    SHIRLY_RUFF = 0;
    if SHIRLY_RUFF
           xpos                = mov.pos( :, 1 );
           xpos(xpos <200) = 0;
    else
    % get the movement parameters
    xpos                = mov.pos( :, 1 );
    spd                 = mov.spd;
    end
    % focus on the relevant times
    ir                  = intersectranges( [ 1 length( xpos ) ], runEnds );
    zidx                = true( size( xpos ) );
    zidx( ir( 1 ) : ir( 2 ) ) = 0;
    xpos( zidx )        = 0;
    spd( zidx )         = 0;
    
    % determine track borders from movement data
    switch extractMode
        case 'bySPD'
            
            % filter the speed
            h                   = makefir( [ 0 LP ], mov.Fs, [], 'low' ); % FIR filter
            h                   = h / sum( h );
            spdf                = firfilt( spd, h ); % this can generate negative speed
            
            % detect segments of high filtered speed
            segs                = parseSchmitt( spdf, spdTH( 2 ), spdTH( 1 ) );
            nsegs               = size( segs, 1 );
            
            % for each such stroke, determine start and end points in time (by local minima)
            eidx                = local_max( spdf, 'min' );
            ne                  = length( eidx );
            segsNew             = zeros( size( segs ) );
            j                   = 1;
            for i               = 1 : nsegs
                si              = segs( i, 1 );
                ei              = segs( i, 2 );
                while eidx( j ) < si
                    j = j + 1;
                    if j > ne
                        break
                    end
                end
                segsNew( i, 1 ) = eidx( j - 1 );
                while eidx( j ) < ei
                    j           = j + 1;
                    if j > ne
                        break
                    end
                end
                segsNew( i, 2 ) = eidx( j );
            end
            
            % unite overlapping segments of high speed
            segs                = uniteranges( segsNew );
            
            % accumulate the x-positions for the start and end points of each stroke
            xbord               = xpos( segs );
            xb                  = xbord( : );
            
        case 'byPOS'
            
            xb                  = xpos( ir( 1 ) : ir( 2 ) );
            
    end
    
    % this part is the same for bySPD and byPOS
    
    % keep only data in the desired range
    idx                 = inrange( xb, trackEnds );
    xb                  = xb( idx );
    
    % partition these numbers into 3 groups: two borders and the rest (in the middle)
    ngroups             = 3;
    idx                 = kmeans( xb, ngroups );
    
    % fit Gaussian models to these data
    gm                  = zeros( ngroups, 3 );
    for i = 1 : ngroups
        xi              = xb( idx == i );
        gm( i, : )      = [ mean( xi ) std( xi ) i ];
    end
    
    % choose the two closest to the initial guess
    gidx                = mindist( gm( :, 1 ), putCOMx );
    com                 = gm( gidx, : );
    
    % to minimize false negatives, expand the borders towards the center of the track
    bords               = zeros( 2, 1 );
    [ ~, midx ]         = min( com( :, 1 ) );
    bords( 1 )          = com( midx, 1 ) + com( midx, 2 ) * 2;
    [ ~, midx ]         = max( com( :, 1 ) );
    bords( 2 )          = com( midx, 1 ) - com( midx, 2 ) * 2;
    
    % generate a new mat by crossings of the track borders
    fprintf( 'Defined COM %s: %0.3g and %0.3g cm\n', extractMode, com( 1 ), com( 2 ) );
    fprintf( 'Defined borders %s: %0.3g and %0.3g cm\n', extractMode, bords( 1 ), bords( 2 ) );
    %extractMode = 'byPTs';
    
    % detect crossings of each border
    matN                    = [];
    for i                   = 1 : 2
        idx                 = xpos > bords( i );
        switch i
            case 1
                solbit      = pt1bit;                                           % right side (larger x value)
            case 2
                solbit      = pt4bit;                                           % left side (smaller x value)
        end
        didx                = diff( idx );
        ev1                 = find( didx == 1 ) + 1;                            % low-to-high
        ev2                 = find( didx == -1 ) + 1;                           % high-to-low
        v1                  = ones( length( ev1 ), 1 );
        v2                  = ones( length( ev2 ), 1 );
        mat1                = [ ev1 solbit * v1 1 * v1 ];
        mat2                = [ ev2 solbit * v2 -1 * v2 ];
        matN                = [ matN; mat1; mat2 ];
    end
    
    % sort matN and assign it a new name (mat)
    mat                     = sortrows( matN, 1 );
    
else
    
    % downsample the digital events to movFs
    mat( :,1 )              = floor( mat( :, 1 ) * sf ) + 1;
    
end


% process digital events
switch extractMode
    
    case 'byPTs'
        
        % get the event times and arrange
        pt1                     = mat( mat( :, 2 ) == pt1bit, : );
        pt4                     = mat( mat( :, 2 ) == pt4bit, : );
        mat                     = sortrows( [ pt1; pt4 ], 1 );
        
    case 'bySOLs'
        
        % get the event times and arrange
        sol1                    = mat( mat( :, 2 ) == sol1bit, : );          % right side
        sol2                    = mat( mat( :, 2 ) == sol2bit, : );          % left side
        sol1on                  = sol1( sol1( :, 3 ) == 1, : );
        sol2on                  = sol2( sol2( :, 3 ) == 1, : );
        sols                    = sortrows( [ sol1on; sol2on ], 1 );

        % make some sanity checks (do not prune since these are actual rewards and not PT crossings)
        strangeSol1             = sum( abs( diff( sol1( :, 3 ) ) ) ~= 2 );
        strangeSol2             = sum( abs( diff( sol2( :, 3 ) ) ) ~= 2 );
        nsame                   = sum( abs( diff( sols( :, 2 ) ) ) ~= 1 );
        if strangeSol1 || strangeSol2
            fprintf( 1, 'Strange1: %d; Strange2: %d\n', strangeSol1, strangeSol2 )
        end
        if nsame
            fprintf( 1, '%d consecutive same-solenoid activations\n', nsame )
        end
        
        % differences
        r2l                     = sol2bit - sol1bit;
        
end

% combine into trials
switch extractMode
    
    case { 'bySOLs' }
        
        idx                     = find( diff( sols( :, 2 ) ) == r2l );        % right-to-left
        r0                      = sols( idx, : ) + 1;
        r1                      = sols( idx + 1, : );
        nt                      = size( r0, 1 );
        trials1                 = [ r0( :, 1 ) r1( :, 1 ) ones( nt, 1 ) * 1 ];
        idx                     = find( diff( sols( :, 2 ) ) == -r2l );       % left-to-right
        r0                      = sols( idx, : ) + 1;
        r1                      = sols( idx + 1, : );
        nt                      = size( r0, 1 );
        trials2                 = [ r0( :, 1 ) r1( :, 1 ) ones( nt, 1 ) * 2 ];
        trials                  = sortrows( [ trials1; trials2 ], 1 );      % by definition interleaved, R2L/L2R
        
    case { 'byPTs', 'byPOS', 'bySPD' }
        
        % identify consecutive pairs of events:
        % byPTs:
        % [pt4bit -1]->[pt1bit 1], are high->low, thus R->L: name these "1"
        % [pt1bit -1]->[pt4bit 1], are low->high, thus L->R: name these "2"
        
        % byPOS/bySPD:
        % [pt4bit -1]->[pt1bit -1], are high->low, thus R->L: name these "1"
        % [pt1bit 1]->[pt4bit 1], are low->high, thus L->R: name these "2"

        if isequal( extractMode, 'byPTs' )
            dval                = 2;
            tnums               = [ 2 1 ];
        else
            dval                = 0;
            tnums               = [ 1 2 ];
        end
        row                     = find( ismember( mat( 2 : end, 2 : 3 ) - mat( 1 : end - 1, 2 : 3 ), [ pt1bit - pt4bit dval ], 'rows' ) );
        trials1                 = [ mat( row, 1 ) mat( row + 1, 1 ) ones( length( row ), 1 ) * tnums( 1 ) ];
        row                     = find( ismember( mat( 2 : end, 2 : 3 ) - mat( 1 : end - 1, 2 : 3 ), [ pt4bit - pt1bit dval ], 'rows' ) );
        trials2                 = [ mat( row, 1 ) mat( row + 1, 1 ) ones( length( row ), 1 ) * tnums( 2 ) ];
        trials                  = sortrows( [ trials1; trials2 ], 1 );      % by definition interleaved, R2L/L2R
        
        pt1                     = mat( mat( :, 2 ) == pt1bit, : );
        pt4                     = mat( mat( :, 2 ) == pt4bit, : );
        pt1on                   = pt1( pt1( :, 3 ) == 1, : );
        pt4on                   = pt4( pt4( :, 3 ) == 1, : );
        pts                     = sortrows( [ pt1on; pt4on ], 1 );
        sols                    = pts;
        sol1                    = pt1;
        sol2                    = pt4;
end

% expand trials towards edges of track
switch extractMode

    case { 'byPOS', 'bySPD' }
        % for each crossing, detect the time of first prior crossing
        % of relevant COM
        t0                      = trials( :, 1 : 2 );
        x0                      = xpos( t0 );
        difs                    = [ abs( x0( : ) - com( 1, 1 ) )  abs( x0( : ) - com( 2, 1 ) ) ];
        [ ~, minidx ]           = min( difs, [], 2 );
        minidx                  = reshape( minidx, size( x0 ) );
        n                       = size( x0, 1 );
        t1                      = NaN( size( t0 ) );
        mmxpos                  = ir;
        % go over trials
        % for each trial, determine whether running from 1->2 or 2->1
        % if from 1->2, then for 1, go backwards in time; and for 2, go forward
        % if from 2->1, then for 1, go forward in time; for 2 go backwards
        % when starting from 1, look for a smaller value
        % when starting from 2, look for a larger value
        for j                   = 1 : n
            
            if isequal( minidx( j, : ), [ 1 2 ] )
                % onset
                tf              = 0;
                k               = t0( j, 1 );
                if j > 1
                    tprev       = t0( j - 1, 1 );
                else
                    tprev       = 1;
                end
                while ~tf && k > mmxpos( 1 ) && k > tprev
                    k           = k - 1;                            % go backwards 
                    tf          = xpos( k ) < com( 1, 1 );          % look for smaller
                end
                if k == tprev
                    [ ~, mi ]   = min( xpos( k : t0( j, 1 ) ) );
                    t1( j, 1 )  = k + mi - 1;                       % add (instead of removing!!) 1 to prevent overlapping trials
                else
                    t1( j, 1 )  = k;
                end
                
                % offset
                tf              = 0;
                k               = t0( j, 2 );
                if j < n
                    tnext       = t0( j + 1, 1 );
                else
                    tnext       = inf;
                end
                while ~tf && k < mmxpos( 2 ) && k < tnext
                    k           = k + 1;                            % go forward
                    tf          = xpos( k ) > com( 2, 1 );          % look for larger
                end
                if k == tnext
                    [ ~, mi ]   = max( xpos( t0( j, 2 ) : k ) );
                    t1( j, 2 )  = t0( j, 2 ) + mi - 1;
                else
                    t1( j, 2 )  = k;
                end
                
            else
                
                % onset
                tf              = 0;
                k               = t0( j, 1 );
                if j > 1
                    tprev       = t0( j - 1, 1 );
                else
                    tprev       = 1;
                end
                while ~tf && k > mmxpos( 1 ) && k > tprev
                    k           = k - 1;                            % go backwards
                    tf          = xpos( k ) > com( 2, 1 );          % look for larger
                end
                t1( j, 1 )      = k;
                if k == tprev
                    [ ~, mi ]   = max( xpos( k : t0( j, 1 ) ) );
                    t1( j, 1 )  = k + mi - 1;                       % add (instead of removing!!) 1 to prevent overlapping trials
                else
                    t1( j, 1 )  = k;
                end
                
                % offset
                tf              = 0;
                k               = t0( j, 2 );
                if j < n
                    tnext       = t0( j + 1, 1 );
                else
                    tnext       = inf;
                end
                while ~tf && k < mmxpos( 2 ) && k < tnext
                    k           = k + 1;                            % go forward
                    tf          = xpos( k ) < com( 1, 1 );          % look for smaller
                end
                if k == tnext
                    [ ~, mi ]   = min( xpos( t0( j, 2 ) : k ) );
                    t1( j, 2 )  = t0( j, 2 ) + mi - 1;
                else
                    t1( j, 2 )  = k;
                end
                
            end % direction

        end % trials
        
        % limit expansion by time
        dt                      = t1 - t0;
        kidx                    = abs( dt ) > ceil( maxExpansion * mov.Fs ); 
        t1( kidx )              = t0( kidx );

        % now replace the trial edge time with the expanded edges
        trials( :, 1 : 2 )      = t1;
        
    otherwise
        
end
        
% upsample the digital events back to spkFs
trials( :, 4 )              = ( trials( :, 1 ) - 1 ) / sf + 1;
trials( :, 5 )              = trials( :, 2 )  / sf;

% prune ultra-short/long trials
tdur                        = diff( trials( :, 4 : 5 ), [], 2 ) / spkFs; % [s]
rmv1                        = tdur < minDur;
rmv2                        = tdur >= maxDur;
rmv                         = rmv1 | rmv2;
if sum( rmv ) > 0
    fprintf( 1, 'Removed %d trials: %d shorter than %0.3g s, %d longer than %0.3g s\n' ...
        , sum( rmv ), sum( rmv1 ), minDur, sum( rmv2 ), maxDur )
    trials( rmv, : ) = [];
end
nt                          = size( trials, 1 );

% prune by position during onset/offset (fit a MOG)
if useMov
    xpos                    = mov.pos( :, 1 );
    xtrials                 = [ xpos( trials( :, 1 : 2 ) ) ( 1 : nt )' ];
    x                       = [ xtrials( :, 1 ); xtrials( :, 2 ) ];
    idx                     = kmeans( x, 2 );
    x1                      = x( idx == 1 );
    x2                      = x( idx == 2 );
    gm1                     = [ mean( x1 ) std( x1 ) ];
    gm2                     = [ mean( x2 ) std( x2 ) ];
    p1                      = 2 * normcdf( -abs( ( x - gm1( 1 ) ) / gm1( 2 ) ), 0, 1 );
    p2                      = 2 * normcdf( -abs( ( x - gm2( 1 ) ) / gm2( 2 ) ), 0, 1 );
    valid                   = p1 > alphaLevel | p2 > alphaLevel;
    valid                   = all( reshape( valid, [ nt 2 ] ), 2 );
    trials( ~valid, : )     = [];
    if sum( ~valid )
        fprintf( 1, 'Removed %d/%d trials: positions far from %0.3g or %0.3g\n' ...
            , sum( ~valid ),nt, gm1( 1 ), gm2( 1 ) )
    end
end

% report and summarize
fprintf( 1, '\n\n%s: Pruning summary; extractMode: %s\n', mfilename, extractMode )
fprintf( 1, 'initial counts (low-to-high + high-to-low): Sol1: %d; Sol2: %d\n' ...
    , size( sol1, 1 ), size( sol2, 1 ) )
fprintf( 1, '%d alternating solenoid trials; %d after duration pruning\n' ...
    , size( sols, 1 ), nt )
fprintf( 1, '%d after position pruning: 1-to-2: %d; 2-to-1: %d trials.\n'...
    , size( trials, 1 ), sum( trials( :, 3 ) == 1 ), sum( trials( :, 3 ) == 2 ) )

if ~useMov
    % never plot since requires position data
    % never save since the format will not match
    return
end


if graphics( 1 )
    
    xmat                    = [ xtrials( :, 1 ) xtrials( :, 3 )  xtrials( :, 2 )  xtrials( :, 3 ) ];
    xplot                   = xmat';  
    xplot                   = xplot( : );
    figs( 1 )               = figure;
    
    %----------------------------------------------------

    % consider: make sure that the plots on the r.h.s. are of the raw data 
    
    matP                    = mat;
    matP( :, 1 )            = ( matP( :, 1 ) - 1 ) / sf;
    pt1                     = matP( matP( :, 2 ) == pt1bit, : );
    pt2                     = matP( matP( :, 2 ) == pt2bit, : );
    pt3                     = matP( matP( :, 2 ) == pt3bit, : );
    pt4                     = matP( matP( :, 2 ) == pt4bit, : );
    
    pt1( :, 1 )             = round( pt1( :, 1 ) * sf );
    pt2( :, 1 )             = round( pt2( :, 1 ) * sf );
    pt3( :, 1 )             = round( pt3( :, 1 ) * sf );
    pt4( :, 1 )             = round( pt4( :, 1 ) * sf );
    
    % position at onset/offset of sol1:
    psol1on                 = mov.pos( sol1( sol1( :, 3 ) == 1, 1 ), : );
    psol1off                = mov.pos( sol1( sol1( :, 3 ) == -1, 1 ), : );
    psol2on                 = mov.pos( sol2( sol2( :, 3 ) == 1, 1 ), : );
    psol2off                = mov.pos( sol2( sol2( :, 3 ) == -1, 1 ), : );
    
    ppt1on                  = mov.pos( pt1( pt1( :, 3 ) == 1, 1 ), : );
    ppt2on                  = mov.pos( pt2( pt2( :, 3 ) == 1, 1 ), : );
    ppt3on                  = mov.pos( pt3( pt3( :, 3 ) == 1, 1 ), : );
    ppt4on                  = mov.pos( pt4( pt4( :, 3 ) == 1, 1 ), : );
    
    subplot( 2, 3, 1 )
    hold on
    if ~isempty( psol1on )
        ph(1) = plot( psol1on( :, 1 ), psol1on( :, 2 ), '.' );
        set( ph(1), 'color', [ 1   0 1 ] )
    end
    if ~isempty( psol1off )
        ph(2) = plot( psol1off( :, 1 ), psol1off( :, 2 ), '.' );
        set( ph(2), 'color', [ 1 0.5 1 ] )
    end
    if ~isempty( psol2on )
        ph(3) = plot( psol2on( :, 1 ), psol2on( :, 2 ), '.' );
        set( ph(3), 'color', [ 0   1 1 ] )
    end
    if ~isempty( psol2off )
        ph(4) = plot( psol2off( :, 1 ), psol2off( :, 2 ), '.' );
        set( ph(4), 'color', [ 0.5 1 1 ] )
    end
    lims                    = minmax( mov.pos ) + [ -10 10 ];
    set( gca, 'xlim', lims, 'ylim', lims )
    axis square
    legend( ph, { sprintf( 'Sol1 on (%0.3g)', mean( psol1on( :, 1 ) ) ) ...
        , sprintf( 'Sol1 off (%0.3g)', mean( psol1off( :, 1 ) ) ) ...
        , sprintf( 'Sol2 on (%0.3g)', mean( psol2on( :, 1 ) ) ) ...
        , sprintf( 'Sol2 off (%0.3g)', mean( psol2off( :, 1 ) ) ) } ...
        , 'Location', 'NorthWest' )
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'X position [cm]' )
    ylabel( 'Y position [cm]' )
    title( extractMode )
    
    
    subplot( 2, 3, 4 )
    hold on
    PT_mid_flag = 0;
    if ~isempty( ppt1on )
        ph(1) = plot( ppt1on( :, 1 ), ppt1on( :, 2 ), '.' );
        set( ph(1), 'color', [ 1 0 0 ] )
    end
    if ~isempty( ppt2on )
        ph(2) = plot( ppt2on( :, 1 ), ppt2on( :, 2 ), '.' );
        set( ph(2), 'color', [ 0 0.7 0 ] )
        PT_mid_flag = 1;
    end
    if ~isempty( ppt3on )
        ph(3) = plot( ppt3on( :, 1 ), ppt3on( :, 2 ), '.' );
        set( ph(3), 'color', [ 0 0 1 ] )
    end
    if ~isempty( ppt4on )
        ph(4) = plot( ppt4on( :, 1 ), ppt4on( :, 2 ), '.' );
        set( ph(4), 'color', [ 1 1 0 ] )
    end
    if PT_mid_flag
    legend( ph, { sprintf( 'PT1 on (%0.3g)', mean( ppt1on( :, 1 ) ) ) ...
        , sprintf( 'PT2 on (%0.3g)', mean( ppt2on( :, 1 ) ) ) ...
        , sprintf( 'PT3 on (%0.3g)', mean( ppt3on( :, 1 ) ) ) ...
        , sprintf( 'PT4 on (%0.3g)', mean( ppt4on( :, 1 ) ) ) } ...
        , 'Location', 'NorthWest' )
    else
        legend( [ph(1), ph(4)], { sprintf( 'PT1 on (%0.3g)', mean( ppt1on( :, 1 ) ) ) ...
        , sprintf( 'PT4 on (%0.3g)', mean( ppt4on( :, 1 ) ) ) } ...
        , 'Location', 'NorthWest' )
    end
    set( gca, 'xlim', lims, 'ylim', lims )
    axis square
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'X position [cm]' )
    ylabel( 'Y position [cm]' )
    
    %----------------------------------------------------
    
    subplot( 1, 3, 2 )
    plot( xplot( 1 : 2 : end ), xplot( 2 : 2 : end ), '.-' )
    axis tight
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'X position [cm]' )
    ylabel( 'Trials number' )
    title( sprintf( '%s: before pruning (%d trials)', replacetok( filename, '\_', '_' ), nt ) )

    xlims                   = xlim;
    ylims                   = ylim;
    vx                      = linspace( xlims( 1 ), xlims( 2 ), 200 );
    gm                      = gm1;
    vy1                     = 1 / ( 2 * pi * gm( 2 ) ) * exp( -( ( vx - gm( 1 ) ).^2 / 2 / gm( 2 ).^2 ) );
    gm                      = gm2;
    vy2                     = 1 / ( 2 * pi * gm( 2 ) ) * exp( -( ( vx - gm( 1 ) ).^2 / 2 / gm( 2 ).^2 ) );
    line( vx, vy1 / max( vy1 ) * diff( ylims ) + ylims( 1 ), 'color', [ 1 1 1 ] * 0.2 )
    line( vx, vy2 / max( vy2 ) * diff( ylims ) + ylims( 1 ), 'color', [ 1 1 1 ] * 0.2 )
    
    %----------------------------------------------------

    % prune by position during onset/offset
    nt                      = size( trials, 1 );
    xtrials                 = [ xpos( trials( :, 1 : 2 ) ) ( 1 : nt )' ];
    xmat                    = [ xtrials( :, 1 ) xtrials( :, 3 )  xtrials( :, 2 )  xtrials( :, 3 ) ];
    xplot                   = xmat';
    xplot                   = xplot( : );
    
    subplot( 1, 3, 3 )
    plot( xplot( 1 : 2 : end ), xplot( 2 : 2 : end ), '.-' )
    axis tight
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'X position [cm]' )
    ylabel( 'Trials number' )
    title( sprintf( '%s: after position pruning (%d trials)', replacetok( filename, '\_', '_' ), nt ) )
   
end

%--------------------------------------------------------------------%
% Step 2. Get movement data for each trial.
%--------------------------------------------------------------------%
% (2) then load the movement data
%   for movement along the X-axis, load only 
%       x, speed, movement direction (and head orientation, if available)
%   for each trial determine the movement features
%   concatenate and arrange in a 6-column matrix: 
%       [ trial index, time, x, speed, direction, orientation ]
%   here, time is in units of the whl file

% get time indices
nt                          = size( trials, 1 );
ns                          = diff( trials( :, 1 : 2 ), [], 2 ) + 1;
tidx                        = single( zeros( sum( ns ), 1 ) );
tnum                        = tidx;
ttyp                        = tidx;
p                           = [ 0; cumsum( ns ) ]; 
for i                       = 1 : nt
    idx                     = trials( i, 1 ) : trials( i, 2 );
    pidx                    = ( p( i ) + 1 ) : ( p( i + 1 ) );
    tidx( pidx )            = idx( : );
    tnum( pidx )            = i;
    ttyp( pidx )            = trials( i, 3 );
end
tmov.tnum                   = tnum;
tmov.tidx                   = tidx;
tmov.ttyp                   = ttyp;
tmov.pos                    = mov.pos( tidx, : );
tmov.spd                    = mov.spd( tidx, : );
tmov.dir                    = mov.dir( tidx, : );
tmov.ang                    = mov.ang( tidx, : );
tmov.whl                    = mov.whl;
tmov.Fs                     = mov.Fs;

% compute some statistics for each trial:
xrange                      = zeros( nt, 2 );
tdur                        = zeros( nt, 1 );
mspd                        = zeros( nt, 1 );
mdir                        = zeros( nt, 1 );
mang                        = zeros( nt, 1 );
for i                       = 1 : nt
    idx                     = tmov.tnum == i;
    xrange( i, : )          = minmax( tmov.pos( idx, 1 ) );
    tdur( i )               = sum( idx ) / tmov.Fs;
    mspd( i )               = mean( tmov.spd( idx, : ) );
    mdir( i )               = circ_mean( tmov.dir( idx, : ) );
    mang( i )               = circ_mean( tmov.ang( idx, : ) );
end
tparams                     = [ xrange tdur mspd mdir mang ];

if graphics
    
    pbins                   = -pi / 20 : pi / 20 : ( 2 * pi + pi / 20 );
    bins                    = ( pbins( 1 : end - 1 ) + pbins( 2 : end ) ) / 2;
    sbins                   = 2.5 : 5 : 152.5;
    tbins                   = minDur : 2 : maxDur;
    
    figs( 2 ) = figure;
    subplot( 4, 2, 1 ), hist( mdir( trials( :, 3 ) == 1 ), bins ), xlim ( [ 0 2 * pi ] ), title( 'Right to left' ), xlabel( 'Direction [rad]' )
    subplot( 4, 2, 2 ), hist( mdir( trials( :, 3 ) == 2 ), bins ), xlim ( [ 0 2 * pi ] ), title( 'Left to right' )
    subplot( 4, 2, 3 ), hist( mang( trials( :, 3 ) == 1 ), bins ), xlim ( [ 0 2 * pi ] ), xlabel( 'Orienation [rad]' )
    subplot( 4, 2, 4 ), hist( mang( trials( :, 3 ) == 2 ), bins ), xlim ( [ 0 2 * pi ] ), 
    subplot( 4, 2, 5 ), hist( mspd( trials( :, 3 ) == 1 ), sbins ), xlabel( 'Speed [cm/s]' )
    subplot( 4, 2, 6 ), hist( mspd( trials( :, 3 ) == 2 ), sbins ), 
    subplot( 4, 2, 7 ), hist( tdur( trials( :, 3 ) == 1 ), tbins ), xlabel( 'Duration [s]' )
    subplot( 4, 2, 8 ), hist( tdur( trials( :, 3 ) == 2 ), tbins ), 
        xlabel( sprintf( '%s (%d trials)', replacetok( filename, '\_', '_' ), nt ) )

    for i = 1 : 8
        subplot( 4, 2, i )
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    
end

%--------------------------------------------------------------------%
% save the results and figure
%--------------------------------------------------------------------%
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( trlfname, 'file' ) )
    fprintf( '%s: Saving MATLAB %s\n', mfname, trlfname )
    save( trlfname, 'trials', 'tmov' )
    if graphics && ~isempty( savetype )
        fprintf( '%s: Saving %s\n', mfname, figfname )
        for i = 1 : length( figs )
            print( figs( i ), '-dpng', [ figfname '.part' num2str( i ) '.' savetype ] )
        end
    end
end

return

% EOF

% evaluation and visualization for dataset with all options
% was evaluated on mP23_37, with 
% timeLimits = [ 5.76e4 2.434e5 ]; % [mov.Fs]

[ TRIALS{ 1 }, TMOV{ 1 } ]          = get_LT_trials( filebase, 'Overwrite', 0, 'graphics', 1, 'extractMode', 'byPTs' );
[ TRIALS{ 2 }, TMOV{ 2 } ]        = get_LT_trials( filebase, 'Overwrite', 0, 'graphics', 1, 'extractMode', 'bySOLs' );
[ TRIALS{ 3 }, TMOV{ 3 } ]        = get_LT_trials( filebase, 'Overwrite', 0, 'graphics', 1, 'extractMode', 'bySPD', 'runEnds', timeLimits );

colors = [ 1 0 0; 0 0.7 0; 0 0 1 ];
styles = { '-', '--', '--' }

figure, hold on,

for j = [ 1 3 ]
    trials = TRIALS{ j };
    tmov = TMOV{ j };
    utrials = unique( tmov.tnum );
    for i = 1 : length( utrials )
        tnum = utrials( i );
        idx = tmov.tnum == tnum;
        ph = plot( tmov.tidx( idx, 1 ), tmov.pos( idx, 1 ) );
        set( ph, 'color', colors( j, : ), 'linestyle', styles{ j } );
    end
end


%legend( ph, { 'byPTs', 'bySOLs' }, 'Location', 'NorthWest' )

