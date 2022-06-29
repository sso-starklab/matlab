% CALC_FIELD        compute a place field, phase map etc (supports 1D or 2D positions and periods)
%
% CALL              [ F, X ] = CALC_FIELD( ST, POS, BSIZE, FS )
%                   [ ...,, C, T, V, SIDX, EDGES ] = CALC_FIELD( ..., PERIODS, MINSPD, VALS, FUNC, SPD, SD )
%
% GETS              st          spike times (sec), relative to pos( 1, : )
%                   pos         position (1D or 2D; cm), first position at time 1/Fs
%                   bsize       bin size (cm) (or bin edges, number of bins.. see BINDATA)
%                   Fs          of positions (Hz)
%                   periods     optional (sec); mx2 matrix, start and end time of each period
%                   minspd      {10}; minimum speed (cm/sec)
%                   vals        [] additional values (e.g. phases) sampled at the same instants as the spikes
%                   func        {'circ_mean'} function to operate on the additional data
%                   spd         {[]:derived from position signal}; optional speed (cm/sec) 
%                   SD          SD (cm); default - 1 spatial bin
%
% RETURNS           F           smoothed rate map. by convention - NaN for places never visited
%                   X           bin centers
%                   C           spike count (per bin)
%                   T           time spent in each bin (sec)
%                   V           additional map (e.g. phase map); NaN for places w/o spikes
%                   SIDX        indices of spikes considered for the analysis (filtered by periods & speed)
%                   EDGES       bin edges
%
% CALLS             DERIVE, CAR2POL, GAUSSKERNEL, FIRFILT, BINDATA

% 14-dec-10 ES

% revisions
% 13-may-11 much faster algorithm for counting spikes/position bin
% 24-may-11 speed threshold
% 26-may-11 vals argument added
% 07-aug-11 round up st0 to prevent index 0
% 21-mar-21 SDx and SDy forced to double for compatibility with imfilter

% algorithm:
% (1) keep only the spikes (and positions) of the relevant times (periods)
% (2) translate spike times to the same sampling frequency as the position sampling
% (3) bin the positions 
% (4) compute time/position bin
% (5) count number of spikes/position bin
% (6) smooth each spatially (spike count, time)
% (7) divide
%
% in case of an additional argument VALS, the values of that argument
% during each spike are averaged (or operated upon) using the function func 
%
% discard: 
% (1) spikes when the animal's speed is <10 cm/sec
% (2) trials w/ < 2 theta cycles, < 3 spikes, and mean speed >10 cm/sec
%
% place fields:
% (1) peak rate > 2 Hz
% (2) size - rate > 10% of peak

function [ f, x, c, t, v, uidx, edges ] = calc_field( st, pos, bsize, Fs, periods, minSpd, vals, func, spd, SD )

% arguments
nargs = nargin;
if nargs < 4 || isempty( st ) || isempty( pos ) || isempty( bsize ) || isempty( Fs )
    error( 'missing arguments' )
end
dims = size( pos, 2 );
if dims > 2 
    error( 'input size mismatch: pos' )
end
if nargs < 5 || isempty( periods ), periods = []; end
if ~isempty( periods )
    if size( periods, 2 ) ~= 2
        error( 'input size mismatch: periods' )
    end
end
if nargs < 6 || isempty( minSpd ), minSpd = 0; end % cm/sec
if nargs < 7 || isempty( vals ), vals = []; end
if ~isempty( vals ) && ~isequal( size( vals ), size( st ) ), error( 'input size mismatch: vals/st' ), end
if nargs < 8 || isempty( func ), func = 'circ_mean'; end % assume that the vals is circular (e.g. phase)
if ismember( func, { 'mean', 'nanmean' } )
    vtype = 'linear'; 
else
    vtype = 'circular';
end
if nargs < 9 || isempty( spd )
    if dims == 1
        spd = abs( derive( pos, 1 / Fs ) );
    else
        dx = derive( pos, 1 / Fs ); 
        spd = car2pol( dx( :, 1 ), dx( :, 2 ) );
    end
end
if size( spd, 1 ) ~= size( pos, 1 )
    error( 'input size mismatch: spd/pos' )
end
if nargs < 10 || isempty( SD ), SD = []; end
if ~isempty( SD )
    if length( SD ) ~= 1
        error( 'inpute size mismatch: SD' )
    end
end

% keep only the data in the requested periods
if ~isempty( periods )
    pidx = ~true( size( st ) );
    for i = 1 : size( periods, 1 ),
        pidx = pidx | st >= periods( i, 1 ) & st <= periods( i, 2 );
    end
    if ~isempty( vals )
        vals = vals( pidx, : );
    end
    pidx = find( pidx );
    periods2 = round( periods * Fs );
    offset = cumsum( [ 0; diff( periods, [], 2 ) + 1 / Fs ] );
    idx = [];
    stkeep = [];
    for i = 1 : size( periods, 1 ),
        idx = [ idx periods2( i, 1 ) : periods2( i, 2 ) ];
        stkeep = [ stkeep; st( st >= periods( i, 1 ) & st <= periods( i, 2 ) ) - periods( i, 1 ) + offset( i ) ];
    end
    st = stkeep;
    pos = pos( idx( : ), : );
    spd = spd( idx( : ), : );
else
    pidx = [ 1 : size( pos, 1 ) ]';
end
nsamples = size( pos, 1 );

% convert inputs (st, vals) to vector of same length as the position signal:
% (1) spikes falling in [0 1/Fs] -> bin1; [1/Fs+eps 2/Fs]->bin2;...
% (2) remove spikes when the speed is low
uidx = [];
%if length( st ) ~= nsamples
    % expand st
    sidx = sparse( zeros( nsamples, 1 ) );
    st0 = ceil( st * Fs );
    st0( st0 == 0 ) = 1; % round up
    ust = unique( st0 );
    if minSpd > 0
        kidx = spd >= minSpd;
        for i = 1 : length( ust )
            if kidx( ust( i ) )
                sidx( ust( i ) ) = sum( st0 == ust( i ) );
                uidx = [ uidx; pidx( st0 == ust( i ) ) ];
            end
        end
    else
        for i = 1 : length( ust )
            sidx( ust( i ) ) = sum( st0 == ust( i ) );
            uidx = [ uidx; pidx( st0 == ust( i ) ) ];
        end
    end
    st = sidx;
    % expand vals over the selected spikes
    if ~isempty( vals )
        sdata = sparse( zeros( nsamples, 1 ) );
        if minSpd > 0
            kidx = spd >= minSpd;
            for i = 1 : length( ust )
                if kidx( ust( i ) )
                    sdata( ust( i ) ) = feval( func, vals( st0 == ust( i ) ) );
                end
            end
        else
            for i = 1 : length( ust )
                sdata( ust( i ) ) = feval( func, vals( st0 == ust( i ) ) );
            end
        end
        vals = sdata;
    end
% elseif minSpd > 0
%     ridx = spd < minSpd;
%     st( ridx ) = 0;
%     uidx( ridx ) = [];
%     if ~isempty( vals )
%         vals( ridx ) = 0;
%     end
% end

% bin positions and compute time/bin
[ bdata x bidx edges ] = bindata( pos, bsize );
t = bdata / Fs;

% count spikes in each position bin
if dims == 1
    nbins = length( x );
    c = zeros( nbins, 1 );
    for i = 1 : nbins
        c( i ) = full( sum( st( bidx == i ) ) );
    end
    if ~isempty( vals )
        d = zeros( nbins, 1 );
        for i = 1 : nbins
            d( i ) = feval( func, full( vals( st > 0 & bidx == i ) ) );
        end
    end
else
    nbinsx = length( x{ 1 } );
    nbinsy = length( x{ 2 } );
    c = zeros( nbinsx, nbinsy );
    for i = 1 : nbinsx
        for j = 1 : nbinsy
            c( i, j ) = full( sum( st( bidx( :, 1 ) == i & bidx( :, 2 ) == j ) ) );
        end
    end
    if ~isempty( vals )
        d = zeros( nbinsx, nbinsy );
        for i = 1 : nbinsx
            for j = 1 : nbinsy
                d( i, j ) = feval( func, full( vals( st > 0 & bidx( :, 1 ) == i & bidx( :, 2 ) == j ) ) );
            end
        end
    end
end
if ~isempty( vals )
     % fill with the grand mean (filter2 cannot handle NaNs)
     % suboptimal; should fill according to nearest neighbor
    nans = isnan( d );
    if strcmp( vtype, 'linear' )
        d( nans ) = mean( d( ~nans ) );
    elseif strcmp( vtype, 'circular' )
        d( nans ) = circ_mean( d( ~nans ) );
    end
end

% smooth and divide
if dims == 1
    %w = gausskernel( SD, 0, 6 * SD + 1, 1 );
    if ~isempty( SD )
        SDx = SD / mean( diff( x ) );
    else
        SDx = 1;
    end
    nx = ceil( 6 * SDx ) + mod( ceil( 6 * SDx ) + 1, 2 );
    w = gausskernel( SDx, 0, nx, 1 );
    w = w / sum( w );
    if isnan( w ), w = 1; end
    cf = firfilt( c, w );
    tf = firfilt( t, w );
    if ~isempty( vals )
        if strcmp( vtype, 'linear' )
            v = firfilt( d, w );
        elseif strcmp( vtype, 'circular' )
            xx = firfilt( cos( d ), w );
            yy = firfilt( sin( d ), w );
            v = atan2( yy, xx );
            v = mod( v, 2 * pi );
        end
    end
else
    %w = gausskernel( SD, SD, 6 * SD + 1, 6 * SD + 1 );
    if ~isempty( SD )
        SDx = double( SD / mean( diff( x{ 1 } ) ) );
        SDy = double( SD / mean( diff( x{ 2 } ) ) );
    else
        SDx = 1;
        SDy = 1;
    end
    nx = ceil( 6 * SDx ) + mod( ceil( 6 * SDx ) + 1, 2 );
    ny = ceil( 6 * SDy ) + mod( ceil( 6 * SDy ) + 1, 2 );
    w = gausskernel( SDx, SDy, nx, ny );
    w = w / sum( sum( w ) );
    if isnan( w ), w = 1; end
    cf = imfilter( c, w, 'symmetric' );
    tf = imfilter( t, w, 'symmetric' );
    if ~isempty( vals )
        if strcmp( vtype, 'linear' )
            v = imfilter( d, w, 'symmetric' );
        elseif strcmp( vtype, 'circular' )
            xx = imfilter( cos( d ), w, 'symmetric' );
            yy = imfilter( sin( d ), w, 'symmetric' );
            v = atan2( yy, xx );
            v = mod( v, 2 * pi );
        end
    end
end
f = cf ./ tf;
% prevent spilling of sparse data
f( t == 0 ) = NaN; 
if ~isempty( vals )
    if SD == 0
        v = d;
        v( nans ) = NaN;
    else
        if dims == 1
            % make sure there is a continuous stretch of data
            fnans = [ 1 : ( find( ~nans , 1 ) - 1 ) ( find( ~nans , 1, 'last' ) + 1 ) : nbins ];
            v( fnans ) = NaN;
        else
            v( nans ) = NaN;
        end
    end
else
    v = NaN;
end

return

% EOF

[ ff x1 c2 t1 ] = calc_field( s1 / Fs, mov.pos( :, 1 ), edges, mov.Fs, trials + ones( size( trials, 1 ), 1 ) * [ -1 1 ] );
[ ff x1 c2 t1 ] = calc_field( s1 / Fs, mov.pos( :, 1 ), edges, mov.Fs, trials );
[ ff x1 c1 t1 ] = calc_field( s1 / Fs, mov.pos( :, 1 ), edges, mov.Fs );
