% medfilt               median filtering with edge mirroring and optional hollowing
%
% call                  y = medfilt( x, W, hf )
%
% gets                  x       matrix, columns are filtered
%                       W       scalar (filter order)
%                       hf      {0}, flag 
%
% returns               y       median filtered version of x
%
% does
% (1) forces W to be odd
% (2) mirrors edges of x
% (3) median filters
%       if hf is false, same as medfilt1
%       if hf is true, then central sample of W is ignored
%
% notes
% (1) algorithm includes only one sort:
%   -sort the first W items in the list
%   -then remove the 1st item, add the (W+1)'th item, and update the sorted array
%   -repeat this until reaching the end of the list
%   -at each step, compute the median:
%       either by simply taking the middle element of the sorted array
%       or, if hollowing, remove the middle element of the unsorted array
%       and then average the two middle elements of the sorted array
% (2) this is all in MATLAB using C-style writing
%   run time is about two times slower than medfilt1
%   to speed this up, can implement as a mex file
%
% see also              medfilt1

% 20-nov-21 ES

function y = medfilt( x, W, hf )

nargs                           = nargin;
if nargs < 1 || isempty( x )
    y                           = [];
    return
end
if ~ismatrix( x )
    error( 'x must be a matrix' )
end
[ m, n ]                        = size( x );

if nargs < 2 || isempty( W )
    y                           = [];
    return
end
if nargs < 3 || isempty( hf )
    hf                          = 0;
end
if hf > 0.5
    hf                          = 1;
else
    hf                          = 0;
end

% (1) force W to be odd
W                               = floor( W / 2 ) * 2 + 1;
if W == 1
    y                           = x;
    return
elseif W < 1
    error( 'W must be non-negative integer' )
end
if m < W
    error( 'W cannot be larger than number of x''s rows' )
end

% (2) mirror edges
hw                              = floor( W / 2 );
x1                              = x( hw : -1 : 1, : );
x2                              = x( m : -1 : m - hw + 1, : );
x                               = [ x1; x; x2 ];

% (3) allocate space
y                               = zeros( m, n );
new_sidx                        = zeros( W - 1, 1 );
new_sx                          = zeros( W - 1, 1 );
hsx                             = zeros( W - 1, 1 );

% (4) compute the running median
for c                           = 1 : n
    % sort every column independently of the others

    for r                       = 1 : m
        % rows in the same columns are linked
        
        if r == 1
            
            % sort the first W samples
            xidx                = x( 1 : W, c );
            [ sx, sidx ]        = sort( xidx );
            
        else
            
            % remove the first item from the previous segment and shift indices
            j                   = 0;
            for i               = 1 : W
                if sidx( i ) ~= 1
                    j           = j + 1;
                    new_sidx( j ) = sidx( i ) - 1;
                    new_sx( j ) = sx( i );
                end
            end

            % find the rank of the new item relative to the sorted array
            nx                  = x( r + W - 1, c );
            j                   = 1;
            if nx > new_sx( 1 )
                while j < W
                    if nx <= new_sx( j )
                        break
                    else
                        j       = j + 1;
                    end
                end
            end
            
            % actually add the new item to the sorted array
            k                   = 1;
            for i               = 1 : W
                if i == j
                    sx( i )     = nx;
                    sidx( i )   = W;
                else
                    sx( i )     = new_sx( k );
                    sidx( i )   = new_sidx( k );
                    k           = k + 1;
                end
            end
            
        end
        
        % compute the median
        if hf
            % remove the middle item of the unsorted array from the sorted array
            j                   = 0;
            for i               = 1 : W
                if sidx( i ) ~= ( hw + 1 )
                    j           = j + 1;
                    hsx( j )    = sx( i );
                end
            end
            % then compute the median (mean of two if W is odd)
            mx                  = ( hsx( hw ) + hsx( hw + 1 ) ) / 2;
        else
            mx                  = sx( hw + 1 );
        end
        y( r, c )               = mx;
        
    end
    
end

return

% EOF
