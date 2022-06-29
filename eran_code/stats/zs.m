% ZS                z-score for matrix columns.
%
% call              [ Z M S N ] = ZS( X, TVEC, FLAG )
%
% gets              X       observations in rows, variables in columns
%                   TVEC    optional column vector, same for all columns in
%                           X. if supplied, each column is standardized in
%                           groups, according to the pointer vector TVEC.
%                   FLAG    {0}; 1 normalizes by N and thus allows
%                           computation of correlations
%
% returns           Z       centered and scaled observations
%                   M       means of each column (if TVEC, then a matrix of
%                               length( unique( TVEC ) ) by size( X, 2 )
%                   S       standard deviations of each column (")
%                   N       number of observations (per group)
%
% NOTES             (1) when no NaNs and only Z output, a fast C routine is used.
%                   (2) there is a machine error (means are not precisely zero)
%
% calls             ZSC.

% 23-may-04 ES

% revisions
% 03-dec-04 handle NaN input correctly
% 03-oct-05 1. TVEC supported
%           2. zero variance -> zero Z-score
% 09-sep-06 TVEC loop corrected
% 22-sep-06 N output added
% 16-mar-07 C code used for rapid computations; for large input the C
%           routine is ~3 times faster; with safe-guarding against NaNs etc
%           it is ~twice as fast and easily handles huge data sets
% 19-mar-07 FLAG added, to allow correct computation of correlations
% 16-jul-12 double check removed
% 17-aug-19 cleaned up

function [ z, m, s, n ] = zs( x, tvec, flag )

nargs                       = nargin;
nout                        = nargout;
if nargs < 1 || isempty( x ) || ndims( x ) ~= 2
    error( 'matrix data' )
end
if nargs < 2 || isempty( tvec )
    tvec                    = [];
end
if nargs < 3
    flag                    = 0;
end
[ n, ncols ]                = size( x );
z                           = zeros( n, ncols );
nans                        = isnan( x );

if isempty( tvec )
    
    % standardize each column as a whole
    nanf                    = sum( sum( nans ) );
    if nout <= 1 && ~nanf
        % call a fast routine
        z                   = zsc( double( x ), flag );
    else
        % otherwise, do everthing local
        d                   = ones( n, 1 );
        if nanf
            s               = nanstd( x );
            m               = nanmean( x );
            y               = x - d * m;
        else
            m               = sum( x, 1 ) / n;
            y               = x - d * m;
            s               = sqrt( sum( y .* y, 1 ) / ( n - ~flag ) );
        end
        clear x
        zidx                = s == 0;
        den                 = d * s;
        z( :, ~zidx )       = y( :, ~zidx ) ./ den( :, ~zidx );
    end
    
else
    
    % standardize each column by parts, according to pointer vector tvec
    tvec                    = tvec( : ).';
    if n ~= length( tvec )
        error( 'input size mismatch' )
    end
    uidx                    = unique( tvec );
    
    % go over subsets
    for i                   = uidx( : ).'
        
        idx                 = tvec == i;
        n( i, : )           = sum( idx );
        nanf                = sum( sum( nans( idx ) ) );
        if nout <= 1 && ~nanf
            z( idx, : )     = zsc( x( idx, : ), flag );
        else
            di              = ones( n( i ), 1 );
            xi              = x( idx, : );
            if nanf
                si          = nanstd( xi );
                mi          = nanmean( xi );
                yi          = xi - di * mi;
            else
                mi          = sum( xi, 1 ) / n( i );
                yi          = xi - di * mi;
                si          = sqrt( sum( yi .* yi, 1 ) / ( n( i ) - ~flag ) );
            end
            clear xi
            zidx            = si == 0;
            den             = di * si;
            z( idx, ~zidx ) = yi( :, ~zidx ) ./ den( :, ~zidx );
            m( i, : )       = mi;
            s( i, : )       = si;
        end
        
    end
    
end

% correct machine error
%z = z + randn( n, ncols ) * eps;
%c = z' * z / ( n - 1 );    % gives Pearson's correlation coefficients

return

% EOF
