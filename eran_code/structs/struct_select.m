% struct_select         select entries from fields of a structure
%
% CALL                  [ out rc ] = struct_select( in, sidx, dim )
%
% GETS                  in        input structure
%                       sidx      logical vector, n x 1
%                       dim       [] optional (prevent variability in case of two identical
%                                   dimensions)
% 
% RETURNS               out       output structure
%                       rc        return code map:
%                                     -1  the field is a cell array, with no dimensions of size n
%                                             (no nested checking though)
%                                     0   no dimension of the field is of size n; nothing done
%                                             (e.g. headers)
%                                     1   the field had exactly 1 dimension of size n; diluted
%                                     m   the field had m dimension of size n; nothing done
%                                 NaN - if mismatching arguments
%
% NOTE                  that this will work only with fields with 1-4D arrays
%
% CALLS                 nothing
%
% See also              struct_cat, struct_compare, struct_extract, struct_merge, struct_reorder, struct_select
%
% example: 
% st.header = 'somename'; 
% st.field1 = rand( 10, 1 ); 
% st.field2 = rand( 10, 2 ); 
% st.field3 = cell( 2, 10, 3 ); 
% st.field4 = cell( 1 ); 
% st.field4{ 1 } = rand( 10, 1);
% idx = 2 : 2 : 10; 
% slct = false( 10, 1 ); 
% slct( idx ) = 1;
% [ out rc ] = struct_select( st, slct )
% 
% in this case the header and field4 are not modified (rc 0 and -1 respectively), whereas all the rest are diluted.

% 22-mar-13 ES

% revisions
% 17-nov-13 added support for specified dimension (quick modification, should organize..)
% 17-aug-19 cleaned up

function [ s, rc, fieldnames ] = struct_select( in, sidx, dim )

% initialize output
s                           = [];
rc                          = NaN;

% arguments
nargs                       = nargin;
if nargs < 2 || ~isa( in, 'struct' ) || ~isa( sidx, 'logical' )
    return
end
if nargs < 3 || isempty( dim )
    dim                     = [];
end

% prepare
s                           = in;
fieldnames                  = fields( s );
sidx                        = sidx( : ); 
n                           = length( sidx );
m                           = length( fieldnames );
rc                          = zeros( m, 1, 'single' );
for i                       = 1 : m
    idx                     = size( s.( fieldnames{ i } ) ) == n;
    
    if isempty( dim )               % no dimension specified
        
        if sum( idx ) == 0          % nothing of size n: do nothing
            if isa( s.( fieldnames{ i } ), 'cell' ) % no nested checking for now
                rc( i )     = -1;
            else
                rc( i )     = 0;
            end
            
        elseif sum( idx ) == 1      % one field of size n: dilute it
            rc( i )         = 1;
            fidx            = find( idx );
            switch length( idx )
                case 2
                    if fidx == 1
                        s.( fieldnames{ i } )( ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, ~sidx ) = [];
                    end
                case 3
                    if fidx == 1
                        s.( fieldnames{ i } )( ~sidx, :, : ) = [];
                    elseif fidx == 2
                        s.( fieldnames{ i } )( :, ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, :, ~sidx ) = [];
                    end
                case 4
                    if fidx == 1
                        s.( fieldnames{ i } )( ~sidx, :, :, : ) = [];
                    elseif fidx == 2
                        s.( fieldnames{ i } )( :, ~sidx, :, : ) = [];
                    elseif fidx == 3
                        s.( fieldnames{ i } )( :, :, ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, :, :, ~sidx ) = [];
                    end
                otherwise
                    rc = -length( idx );
            end
            
        else                        % more than 1 field of size n: do nothiing
            rc( i ) = sum( idx );
        end
        
    else
        
        if length( idx ) >= dim && idx( dim ) == 1
            switch length( idx )
                case 2
                    if dim == 1
                        s.( fieldnames{ i } )( ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, ~sidx ) = [];
                    end
                case 3
                    if dim == 1
                        s.( fieldnames{ i } )( ~sidx, :, : ) = [];
                    elseif dim == 2
                        s.( fieldnames{ i } )( :, ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, :, ~sidx ) = [];
                    end
                case 4
                    if dim == 1
                        s.( fieldnames{ i } )( ~sidx, :, :, : ) = [];
                    elseif dim == 2
                        s.( fieldnames{ i } )( :, ~sidx, :, : ) = [];
                    elseif dim == 3
                        s.( fieldnames{ i } )( :, :, ~sidx, : ) = [];
                    else
                        s.( fieldnames{ i } )( :, :, :, ~sidx ) = [];
                    end
                otherwise
                    rc = -length( idx );
            end
        else
            rc( i ) = -1;
        end
    end
end

return

% EOF
