% LOCAL_MAX     detect all local extrema in a matrix
%
% call          [ idx, vals, etype ] = local_max( x, mode, dt )
%
% gets          x           data matrix
%               mode        {'max'} string; options: 'ext', 'min', 'max'
%               dt          {0} number [samples]; closer extrema are pruned
%
% returns       idx         indices of local extrema
%               vals        values
%               etype       1: local maximum; -1: local minimum
%
% calls         nothing
%
% does
% (1) computes the direction of the derivative
% (2) detects points of direction change (local extrema), classifies as
%       maxima/minima
% (3) optionally, prunes adjacent maxima to be no more than dt samples
%       apart. after pruning, satisfies (for 'max'/'min'):
% 
%           sum( diff( idx ) < dt ) == 0
%
%       for 'ext', satisfies:
% 
%           sum( diff( idx( etype == 1 ) ) < dt ) == 0
%           sum( diff( idx( etype == -1 ) ) < dt ) == 0
%           sum( diff( etype ) == 0 ) == 0
%
% note: if X is a matrix, this routine returns inds as a 2-column matrix
% (index of extremum, column index) and vals is a matched vector

% 15-jan-13 ES

% revisions
% 30-mar-21 (1) added 3rd output argument etype
%           (2) added 3rd input argument dt
% 31-mar-21 (1) exhaustive pruning

function [ idx, vals, etype ] = local_max( x, mode, dt )

idx                             = [];
vals                            = [];
etype                           = [];

% arguments
nout                            = nargout;
nin                             = nargin;
if nin < 1 || isempty( x )
    return
end
sx                              = size( x );
if length( sx ) == 2 && sum( sx == 1 ) >= 1 && sx( 2 ) > 1
    x                           = x';
end
if length( sx ) > 2
    error( 'input size mismatch: x must be a vector or a matrix' )
end
m                               = sx( 1 );
n                               = sx( 2 );

if nin < 2 || isempty( mode )
    mode                        = 'max';
end
mode                            = lower( mode( 1 : 3 ) );
if ~ismember( mode, { 'ext', 'min', 'max' } )
    error( 'input mismatch: mode must be ''ext''/''min''/''max''' )
end

if nin < 3 || isempty( dt )
    dt                          = 0;
end
if numel( dt ) ~= 1 || dt < 0 || imag( dt ) ~= 0
    error( 'input mismatch: dt must be a non-negative real scalar' )
end

% compute second 'derivative'
d2                              = diff( sign( diff( x ) ) );

% identify extrema
switch mode
    case 'ext'
        [ rowMin, colMin ]      = find( d2 > 1 );
        [ rowMax, colMax ]      = find( d2 < -1 );
        etype                   = [ -1 * ones( length( rowMin ), 1 ); ones( length( rowMax ), 1 ) ];
        mat                     = [ [ [ rowMin colMin ]; [ rowMax colMax ] ] etype ];
        smat                    = sortrows( mat, [ 2 1 ] );
        row                     = smat( :, 1 );
        col                     = smat( :, 2 );
        etype                   = smat( :, 3 );
    case 'min'
        [ row, col ]            = find( d2 > 1 );
        etype                   = -1 * ones( length( row ), 1 );
    otherwise
        [ row, col ]            = find( d2 < -1 );
        etype                   = ones( length( row ), 1 );
end
row                             = row + 1;

% prune nearby extrema
if dt > 0
    
    if isequal( mode, 'ext' )
        
        % prune the max
        maxidx                  = find( etype == 1 );
        rowMax                  = row( maxidx );
        colMax                  = col( maxidx );
        d                       = diff( rowMax );
        cls                     = find( d < dt & d > 0 );                   % two nearby maxima of the same column
        n0                      = inf;
        while ~isempty( cls ) && n0 > length( cls )
            n0                  = length( cls );
            idx1                = sub2ind( [ m n ], rowMax( cls ), colMax( cls ) );
            idx2                = sub2ind( [ m n ], rowMax( cls + 1 ), colMax( cls + 1 ) );
            mat               	= [ x( idx1 ) x( idx2 ) ] ;
            %mat               	= [ x( rowMax( cls ) ) x( rowMax( cls + 1 ) ) ] ;
            [ ~, rmv ]       	= min( mat, [], 2 );                        % remove the smaller maximum
            rmv                 = maxidx( unique( cls + rmv - 1  ) );
            row( rmv )        	= [];
            col( rmv )        	= [];
            etype( rmv )      	= [];
            maxidx           	= find( etype == 1 );
            rowMax           	= row( maxidx );
            colMax           	= col( maxidx );
            d                   = diff( rowMax );
            cls                	= find( d < dt & d > 0 );
        end
        
        % prune the min
        minidx                  = find( etype == -1 );
        rowMin                  = row( minidx );
        colMin                  = col( minidx );
        d                       = diff( rowMin );
        cls                     = find( d < dt & d > 0 );
        n0                      = inf;
        while ~isempty( cls ) && n0 > length( cls )
            n0                  = length( cls );
            idx1                = sub2ind( [ m n ], rowMin( cls ), colMin( cls ) );
            idx2                = sub2ind( [ m n ], rowMin( cls + 1 ), colMin( cls + 1 ) );
            mat               	= [ x( idx1 ) x( idx2 ) ] ;
            %mat              	= [ x( rowMin( cls ) ) x( rowMin( cls + 1 ) ) ] ;
            [ ~, rmv ]          = max( mat, [], 2 );
            rmv                 = minidx( unique( cls + rmv - 1 ) );
            row( rmv )          = [];
            col( rmv )          = [];
            etype( rmv )        = [];
            minidx              = find( etype == -1 );
            rowMin              = row( minidx );
            colMin              = col( minidx );
            d                   = diff( rowMin );
            cls                 = find( d < dt & d > 0 );
        end
        
        % prune until no two consecutive max/min
        cls                     = find( diff( etype ) == 0 );               % two consecutive maxima/minima
        n0                      = inf;
        while ~isempty( cls ) && n0 > length( cls )
            n0                  = length( cls );
            idx1                = sub2ind( [ m n ], row( cls ), col( cls ) );
            idx2                = sub2ind( [ m n ], row( cls + 1 ), col( cls + 1 ) );
            mat               	= [ x( idx1 ) x( idx2 ) ] ;
            %mat                 = [ x( row( cls ) ) x( row( cls + 1 ) ) ];  % values
            mat1                = mat( etype( cls ) == 1, : );              % max
            [ ~, rmv1 ]         = min( mat1, [], 2 );
            rmv1                = unique( cls( etype( cls ) == 1 ) + rmv1 - 1 );
            mat2                = mat( etype( cls ) == -1, : );           	% min
            [ ~, rmv2 ]         = max( mat2, [], 2 );
            rmv2                = unique( cls( etype( cls ) == -1 ) + rmv2 - 1 );
            rmv                 = sort( [ rmv1; rmv2 ] );
            row( rmv )       	= [];
            col( rmv )        	= [];
            etype( rmv )        = [];
            cls                 = find( diff( etype ) == 0 );
        end
        
    else
        
        % prune max or min
        d                       = diff( row );
        cls                   	= find( d < dt & d > 0 );
        n0                      = inf;
        while ~isempty( cls ) && n0 > length( cls )
            n0                  = length( cls );
            idx1                = sub2ind( [ m n ], row( cls ), col( cls ) );
            idx2                = sub2ind( [ m n ], row( cls + 1 ), col( cls + 1 ) );
            mat               	= [ x( idx1 ) x( idx2 ) ] ;
            %mat               	= [ x( row( cls ) ) x( row( cls + 1 ) ) ] ;
            switch mode
                case 'max'
                    [ ~, rmv ] 	= min( mat, [], 2 );
                    rmv      	= unique( cls + rmv - 1 );
                    
                case 'min'
                    [ ~, rmv ] 	= max( mat, [], 2 );
                    rmv     	= unique( cls + rmv - 1 );
            end
            row( rmv )      	= [];
            col( rmv )      	= [];
            etype( rmv )     	= [];
            d                   = diff( row );
            cls                 = find( d < dt & d > 0 );
        end

    end
    
end

% organize output
if n == 1
    idx                         = row;
else
    idx                         = [ row col ];
end
if nout > 1
    vals                        = x( row + ( col - 1 ) * m );
end

return

% EOF