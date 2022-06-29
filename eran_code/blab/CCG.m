function [out, t, Pairs, Bias ] = CCG(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
% constructs multiple cross and Auto correlogram
% usage: [ccg, t, pairs] = CCG(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
%
% T gives the time of the events, (events need not be sorted by TIME)
% G says which one is in which group
% BinSize gives the size of a bin in input units (i.e. not scaled by SampleRate)
% HalfBins gives the number of bins on each side of 0 - so the total is 1+2*HalfBins
% SampleRate is for x-axis scaling only.  It defaults to 20000
% GSubset says which groups to plot the CCGS of (defaults to all but group 1)
% Normalization indicates the type of y-axis normalization to be used.
% 'count' indicates that the y axis should show the raw spike count in each bin.
% 'hz' will normalize to give the conditional intensity of cell 2 given that cell 1 fired a spike (default)
% 'hz2' will give the joint intensity, measured in hz^2.
% 'scale' will scale by both firing rates so the asymptotic value is 1.  This gives you
%   the ratio of the coincidence rate to that expected for uncorrelated spike trains
%
% optional input Epochs allows you to compute from only spikes in certain epochs
% and it will bias-correct so you don't see any triangular shape.
% Warning: if gaps between epochs are shorter than total CCG length, this will mess
% up the edges.
%
% The output array will be 3d with the first dim being time lag and the second two
% specifying the 2 cells in question (in the order of GSubset)
% If there is no output specified, it will plot the CCGs
%
% This file calls a C program so your CCG is computed fast.
% to use it, you need to compile mex file from CCGHeart.c
% run : mex -v CCGHeart.c
% different architectures will produce different extensions of mex file,
% also different versions of matlab link mex file to different libraries
% they are mostly taken into account in the code of CCG.m, but if your
% version or archicture is different from those - modify CCGFun string to
% match the name of the mex file your compiler generated.
% optional output t gives time axis for the bins in ms
% optional output argument pairs gives a nx2 array with the indices of the spikes
% in each train that fall in the CCG.

% 15-sep-01 Ken Harris

% revisions
% 2009      Anton Sirota: Epochs, Bias 
% 2018      Eran Stark: prevent memory explosion by indexing the sparse

%----------------------------------------------------------------------%
% preps
%----------------------------------------------------------------------%

nargs               = nargin;
if nargs < 5
    SampleRate      = 20000;
end
if nargs < 6
    GSubset         = unique( G );
end
if nargs < 7
    Normalization   = 'hz';
end
if nargs < 8
    Epochs          = [];
end

if length( G ) == 1
    G               = ones( length( T ), 1 );
    GSubset         = 1;
    nGroups         = 1;
else
    nGroups         = length( GSubset );
end

% Prepare Res and Clu arrays.
G                   = G( : );
T                   = T( : );

if isempty( Epochs )
    Included = find( ismember( G, GSubset ) & isfinite( T ) );
    Epochs = [];
else
    Included = find( ismember( G, GSubset ) & isfinite( T ) & WithinRanges( T, Epochs ) );
    % check gaps between epochs are not too short
    GapLen          = Epochs( 2 : end, 1 ) - Epochs( 1 : ( size( Epochs, 1 ) - 1 ), 2 );
    TooShort        = find( GapLen < BinSize * ( HalfBins + .5 ) );
    if ~isempty( TooShort )
        fprintf( 'WARNING: Epochs ' );
        fprintf( '%d ', TooShort );
        fprintf( 'are followed by too-short gaps.\n' );
    end
end

Res                 = T( Included );

% if no spikes, return nothing
if length( Res ) <= 1
    nBins           = 1 + 2 * HalfBins;
    out             = zeros( nBins, nGroups, nGroups );
    t               = 1000 * ( -HalfBins : HalfBins ) * BinSize / SampleRate;
    Pairs           = [];
    return
end

%----------------------------------------------------------------------%
% compute
%----------------------------------------------------------------------%

% Indexing array (index before full to prevent memory explosion)
G2Clu               = sparse( GSubset, 1, 1 : nGroups );
Clu                 = full( G2Clu( G( Included ) ) );

% sort by time
[ Res, ind ]        = sort( Res );
Clu = Clu(ind);

% Now call the C program
CCGFun              = 'CCGHeart';
if nargout >= 3
    [ Counts, RawPairs ] = feval( CCGFun, Res, uint32( Clu ), BinSize, uint32( HalfBins ) );
    rsRawPairs      = reshape( RawPairs, [ 2 length( RawPairs ) / 2 ] )';
else
    Counts          = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
end

% shape the results
nBins               = 1 + 2 * HalfBins;
% if there are no spikes in the top cluster, CCGEngine will produce a output the wrong size
nPresent            = max( Clu );
Counts              = double( reshape( Counts, [ nBins nPresent nPresent ] ) );
if nPresent < nGroups
    % extent array size with zeros
    Counts( nBins, nGroups, nGroups ) = 0;
end

if nargout >= 3
    Pairs           = Included( ind( double( rsRawPairs ) + 1 ) );
end

%----------------------------------------------------------------------%
% scale
%----------------------------------------------------------------------%

% OK so we now have the bin counts.  Now we need to rescale it.
% remove bias due to edge effects - this should be vectorized
if isempty( Epochs )
    Bias            = ones( nBins, 1 );
else
    nTerm           = [ HalfBins : -1 : 1 , 0.25 , 1 : HalfBins ];
    Bias            = zeros(nBins,1);
    TotLen          = 0;
    for e           = 1 : size( Epochs, 1 )
        EpochLen    = Epochs( e, 2 ) - Epochs( e, 1 );
        EpochBias   = clip( EpochLen - nTerm * BinSize, 0, inf ) * BinSize;
        Bias        = Bias + EpochBias';
        TotLen      = TotLen + EpochLen;
    end
    Bias = Bias / TotLen / BinSize;
end

if isempty( Epochs)
    Trange          = max( Res ) - min( Res ); % total time
else
    Trange          = sum( diff ( Epochs, [], 2 ) );
end
t                   = 1000 * ( -HalfBins : HalfBins ) * BinSize / SampleRate;

% count the number of spikes in each group:
nSpikesPerGroup     = zeros( 1, nGroups );
for g               = 1 : nGroups
    nSpikesPerGroup( g ) = sum( Clu == g );
end

% normalize each group
for g1              = 1 : nGroups
    for g2          = g1 : nGroups
        switch Normalization
            case 'hz'
                Factor      = SampleRate / ( BinSize * nSpikesPerGroup( g1 ) );
                AxisUnit    = '[Hz]';
            case 'hz2'
                Factor      = SampleRate * SampleRate / ( Trange * BinSize );
                AxisUnit    = '[Hz^2]';
            case 'count'
                Factor      = 1;
                AxisUnit    = '[Spikes]';
            case 'scale'
                Factor      = Trange / (BinSize * nSpikesPerGroup( g1 ) * nSpikesPerGroup( g2 ) );
                AxisUnit    = '[Scaled]';
            otherwise
                warning( [ 'Unknown Normalization method ', Normalization ] );
        end
        ccg( :, g1, g2 )    = flipud( Counts( :, g1, g2 ) ) * Factor ./ Bias;
        ccg( :, g2, g1 )    = Counts( :, g1, g2 ) * Factor ./ Bias;
        
        % now plot, if there is no output argument
        if nargout == 0
            FigureIndex     = g1 + nGroups * ( nGroups - g2 );
            subplot( nGroups, nGroups, FigureIndex );
            bar( t, ccg( :, g1, g2 ) );
            if g1 == g2
                FiringRate  = SampleRate * nSpikesPerGroup( g1 ) / Trange;
                xlabel( 'ms' );
                ylabel( [ 'ACG ', AxisUnit ] )
                Ttitle      = sprintf( '%d (~%5.2fHz)', GSubset( g1 ), FiringRate );
                title( Ttitle );
            else
                ylabel( ['CCG ', AxisUnit ] )
                Ttitle      = sprintf( '%d vs %d', GSubset( g1 ), GSubset( g2 ) );
                title( Ttitle );
            end
            if g1 == 1
                ylabel( sprintf( '%d', GSubset( g2 ) ) );
            end
            if g2 == nGroups
                Ttitle      = sprintf( '%d (~%5.2fHz)', GSubset( g1 ), FiringRate );
                title( Ttitle );
                ylabel( AxisUnit )
            end
            axis tight
        end         % nargout
        
    end             % g2
end                 % g1

% only make an output argument if its needed (to prevent command-line spew)
if nargout > 0
    out             = ccg;
end

return

% EOF
