% calc_stc_asg      from gain CCH
% 
% call              calc_stc_asg( gcch, t_ROI, dt, asgMode, roiHardness )
%
% gets              gcch            [spks/s], gain CCH
%                   t_ROI           logical vector, same size as cch rows
%                   dt              [s], CCH bin size
%                   asgMode         {1}         0   detects global maxima in ROI (always exists, may be non-causal)
%                                               1   detects local maxima in ROI (may not exist, but preserves causality)
%                   roiHardness     {[1 0]}     hard/soft edges (limits/does not limit ASG base by ROI)
% 
% returns           g1              ASG of positive extremum in the ROI
%                   g2              ASG of negative extremum in the ROI
%                   g1base          support for positive extremum [samples]
%                   g2base          support for negative extremum [samples]
% 
% calls             nothing

% 07-mar-21 ES

function [ g1, g2, g1base, g2base ] = calc_stc_asg( gcch, t_ROI, dt, asgMode, roiHardness )

nargs                           = nargin;
if nargs < 4 || isempty( asgMode )
    asgMode                     = 1;
end
if nargs < 5 || isempty( roiHardness )
    roiHardness                 = [ 1 0 ];
end

% prepare
ft_ROI                          = find( t_ROI );
n                               = size( gcch, 1 );
gcch( isnan( gcch ) )           = 0;

% check input size

% support of positive extremum
switch asgMode
    case 0
        x                       = gcch( t_ROI );
        [ ~, maxidx ]           = max( x );
    case 1
        pidx                    = ft_ROI( [ 1 end ] ) + [ -1 1 ]';
        xidx                    = [ pidx( 1 ); ft_ROI; pidx( 2 ) ];
        nROI                    = length( ft_ROI ) + 1;
        x                       = gcch( xidx );
        sx                      = find( diff( x ) == 0 );                   % same values of x
        ux                      = x;
        ix                      = ( 1 : length( x ) )';
        ix( sx )                = [];
        ux( sx )                = [];
        tmpidx                  = find( diff( sign( diff( ux ) ) ) < -1 ) + 1;
        maxidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) > ux( end - 1 )
            maxidx              = sx( end ) - 1;
        end
end
if length( maxidx ) > 1
    [ ~, subidx ]               = max( x( maxidx ) );
    maxidx                      = maxidx( subidx );
end
if isempty( maxidx )
    sidx                        = [];
else
    pidx                    	= ft_ROI( maxidx );
    si                        	= 1;                                    	% first preceding negative value
    for i                     	= pidx : -1 : 1
        if gcch( i ) < 0
            si                 	= i + 1;
            break
        end
    end
    ei                         	= n;                                    	% first proceeding negative value
    for i                      	= pidx : n
        if gcch( i ) < 0
            ei                 	= i - 1;
            break
        end
    end
    sidx                       	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
g1base                          = sidx;
% compute the integral
a1                              = nanmean( gcch( sidx ) );
b1                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g1                              = a1 * b1;

% support of negative extremum
switch asgMode
    case 0
        [ ~, minidx ]           = min( x );
    case 1
        tmpidx                  = find( diff( sign( diff( ux ) ) ) > 1 ) + 1;
        minidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) < ux( end - 1 )
            minidx              = sx( end ) - 1;
        end
end
if length( minidx ) > 1
    [ ~, subidx ]               = min( x( minidx ) );
    minidx                      = minidx( subidx );
end
if isempty( minidx )
    sidx                        = [];
else
    pidx                     	= ft_ROI( minidx );
    si                         	= 1;                                       	% first preceding positive value
    for i                      	= pidx : -1 : 1
        if gcch( i ) > 0
            si                	= i + 1;
            break
        end
    end
    ei                       	= n;                                    	% first proceeding positive value
    for i                      	= pidx : n
        if gcch( i ) > 0
            ei               	= i - 1;
            break
        end
    end
    sidx                     	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
g2base                          = sidx;
% compute the integral
a2                              = nanmean( gcch( sidx ) );
b2                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g2                              = a2 * b2;

return

% EOF