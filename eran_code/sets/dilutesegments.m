% dilutesegments        of short duration/separation
%
% CALL                  mat = dilutesegments( mat, mindur, minisi )
% 
% GETS                  mat             ranges in rows
%                       mindur          {0}, samples
%                       minisis         {0}, samples
%
% DOES
%                       merges segments with separation <minisi
%                       removes (combined) segments with duration <mindur
% 
% CALLS                 sortranges
%
% See also 
%                       getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, uniteranges 

% 16-jul-12 ES

% revisions
% 31-mar-13 empty mat returned as empty
% 17-aug-19 cleaned up

function [ mat, loc ] = dilutesegments( mat, mindur, minisi )

nargs           = nargin;
if nargs < 1 || isempty( mat )
    return
end
if nargs < 2 || isempty( mindur )
    mindur      = 0;
end
if nargs < 3 || isempty( minisi )
    minisi      = 0;
end

mat             = sortranges( mat );
mindur          = mindur( 1 );
minisi          = minisi( 1 );
loc             = ( 1 : size( mat, 1 ) )';

% gap must be higher than minisi - otherwise combine/remove
if ~isempty( mat )
    if minisi < 0
        isiMode             = 'remove';
        minisi              = abs( minisi );
    else
        isiMode             = 'combine';
    end
    switch isiMode
        case 'combine'
            sgap            = find( ( mat( 2 : end, 1 ) - mat( 1 : end - 1, 2 ) ) < minisi );
            for gi          = 1 : length( sgap )
                mat( sgap( gi ) + 1, 1 ) = mat( sgap( gi ), 1 );
            end
            mat( sgap, : )  = [];
            loc( sgap, : )  = [];
        case 'remove'
            gaps            = mat( :, 1 ) - [ -inf; mat( 1 : end - 1,2 ) ]; 
            sidx            = find( gaps < minisi ); 
            mat( sidx, : )  = [];
            loc( sidx, : )  = [];
    end
end

% combined duration must be above minDuration - otherwise remove
ridx                        = diff( mat, [], 2 ) < mindur;
mat( ridx, : )              = [];
loc( ridx, : )              = [];

return

% EOF