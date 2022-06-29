% BINDATA           bin 1 or 2D data
%
% CALL              [ BDATA, BINS, BIDX, EDGES ] = BINDATA( DATA, BSIZE )
%
% GETS              data        data (1D or 2D)
%                   bsize       bin size (scalar), vector of bin edges (for 1D data), or cell
%                                   array of 2 vectors of bin edges (for 2D data)
%                                   to specify NUMBER of bins instead, use negative number/s
%
% RETURNS           bdata       the binned data (vector/matrix) of nxbins x nybins bins
%                   bins        the bin centers (1D: vector; 2D: cell array of two vectors)
%                   bidx        map from data->bdata (indicates bin number)
%                   edges       bin edges
%
% CALLS             nothing

% 24-may-11 ES

% revisions
% 07-jun-11 added buffer of 100*eps to upper edges
% 20-jan-13 corrected for 2D, 1 element case
% 18-aug-19 cleaned up
% 14-sep-19 added default values for bsize and argument checks

function [ bdata, bins, bidx, edges ] = bindata( data, bsize )

bdata                       = [];
bins                        = [];
bidx                        = [];
edges                       = [];

nargs                       = nargin;
if nargs < 1 || isempty( data )
    return
end
if nargs < 2 || isempty( bsize )
    bsize                   = 10 * ones( 1, size( data, 2 ) );
end

[ nn, dims ]                = size( data );
if nn < dims && dims ~= 2
    data                    = data';
    [ nn, dims ]            = size( data );
end
if dims > 2
    error( 'dimensionality not supported' )
end

if dims == 1

    % generate the edges array
    if ~isa( bsize, 'numeric' ) || length( bsize( : ) ) ~= numel( bsize )
        error( 'input size mismatch' )
    elseif numel( bsize ) == 1      % scalar
        if bsize < 0 % number of bins
            nbins           = -bsize;
            bsize           = ( max( data ) - min( data ) ) / nbins;
            edges           = min( data ) : bsize : max( data );
        else % bin size
            if bsize > range( data )
                edges       = [ min( data ) max( data ) ];
            else
                mid         = mean( [ min( data ) max( data ) ] );
                edges       = unique( [ fliplr( mid : -bsize : min( data ) - bsize / 2 ) ...
                    mid : bsize : max( data ) + bsize / 2 ] );
            end
        end
    else                            % vector
        edges               = bsize;
    end
    edges                   = edges( : );
    edges( end )            = edges( end ) + 100 * eps;
    
    % bin data
    [ bdata, bidx ]         = histc( data, edges ); % bidx==0 for data values outside edges
    bdata( end )            = [];
    % bin centers
    bins                    = filter( ones( 2, 1 ) / 2 , 1, edges ); 
    bins( 1 )                = [];
    
elseif dims == 2

    % generate the edges arrays
    if isa( bsize, 'cell' ) && length( bsize ) ~= 2 || isa( bsize, 'numeric' ) && length( bsize( : ) ) ~= numel( bsize )
        error( 'input size mismatch' )
    elseif isa( bsize, 'cell' ) % cell array w/ two vectors - use one per dimension
        edgesx              = bsize{ 1 };
        edgesy              = bsize{ 2 };
    elseif numel( bsize ) <= 2 % scalar/s (positive: bin size; negative: -number of bins)
        xdata               = data( :, 1 );
        ydata               = data( :, 2 );
        if numel( bsize ) == 1
            bsizex          = bsize;
            bsizey          = bsize;
        else
            bsizex          = bsize( 1 );
            bsizey          = bsize( 2 );
        end
        if bsizex < 0
            nbinsx          = -bsizex;
            bsizex          = ( max( xdata ) - min( xdata ) ) / nbinsx;
            edgesx          = min( xdata ) : bsizex : max( xdata );
        else
            if bsizex > range( xdata )
                edgesx      = [ min( xdata ) max( xdata ) ];
            else
                mid         = mean( [ min( xdata ) max( xdata ) ] );
                edgesx      = unique( [ fliplr( mid : -bsizex : min( xdata ) - bsizex / 2 ) ...
                    mid : bsizex : max( xdata ) + bsizex / 2 ] );
            end
        end
        if bsizey < 0
            nbinsy          = -bsizey;
            bsizey          = ( max( ydata ) - min( ydata ) ) / nbinsy;
            edgesy          = min( ydata ) : bsizey : max( ydata );
        else
            if bsizey > range( ydata )
                edgesy      = [ min( ydata ) max( ydata ) ];
            else
                mid         = mean( [ min( ydata ) max( ydata ) ] );
                edgesy      = unique( [ fliplr( mid : -bsizey : min( ydata ) - bsizey / 2 ) ...
                    mid : bsizey : max( ydata ) + bsizey / 2 ] );
            end
        end
    else % 1 vector - use the same for both dimensions
        edgesx              = bsize;
        edgesy              = bsize;
    end
    edgesx                  = unique( edgesx( : ) );
    edgesy                  = unique( edgesy( : ) );
    edgesx( end )           = edgesx( end ) + 100 * eps;
    edgesy( end )           = edgesy( end ) + 100 * eps;
    nbinsx                  = length( edgesx ) - 1;
    nbinsy                  = length( edgesy ) - 1;

    % bins centers
    binsx                   = filter( ones( 2, 1 ) / 2 , 1, edgesx ); 
    binsx( 1 )              = [];
    binsy                   = filter( ones( 2, 1 ) / 2 , 1, edgesy ); 
    binsy( 1 )              = [];
    bins                    = { binsx binsy };
    edges                   = { edgesx edgesy };

    % remove irrelevant data
    kidx                    = data( :, 1 ) >= edgesx( 1 ) & data( :, 1 ) < edgesx( nbinsx + 1 ) & data( :, 2 ) >= edgesy( 1 ) & data( :, 2 ) < edgesy( nbinsy + 1 );
    data                    = data( kidx, : );

    % fast binning (can be made faster if use sort to hist too)
    bdata                   = zeros( nbinsx, nbinsy );
    idx                     = zeros( size( data ) );
    [ sdata, sidx ]         = sortrows( data, 1 );
    if size( data, 1 ) == 1
        cn                  = [ 0; cumsum( histc( data( :, 1 ), edgesx ) )' ];
    else
        cn                  = [ 0; cumsum( histc( data( :, 1 ), edgesx ) ) ];
    end
    for i                   = 1 : nbinsx
        yidx                = cn( i ) + 1 : cn( i + 1 );
        [ y2, y2idx ]       = histc( sdata( yidx, 2 ), edgesy );
        bdata( i, : )       = y2( 1 : nbinsy );
        idx( sidx( yidx ), 1 ) = i;
        idx( sidx( yidx ), 2 ) = y2idx;
    end
    
    % expand indices to original dataseet
    bidx = zeros( nn, 2 );
    bidx( kidx, : ) = idx;

end

return

% EOF