% struct_cat         concat entries in similar structures
%
% [ out rc ] = struct_cat( in1, in2 )
%
% in1, in2  input structures
% sidx      logical vector, n x 1
%
% out       output structure
% rc        return code map:
%               -1  the field is a cell array, with no dimensions of matching sizes
%               0   no dimension of the field is of matching size; nothing done
%                       (e.g. headers)
%               1   the field had exactly 1 dimension of mismatching size; merged
%               m   the field had m dimensions of mismatching size; nothing done
%           NaN - if mismatching structures
% 
% calls     nothing
%
% see also      struct_cat, struct_compare, struct_extract, struct_merge, struct_reorder, struct_select
%
% note
% idetical-size fields will be skipped (no way to know which dimension to
% concatenate along). Thus, struct_cat( s, s ) is the same as struct_cat(
% s, [] ). To bypass this, use struct_cat( s, s, 1 ) to concatenate rows
% etc

% 18-apr-13 ES

% revisions
% 28-apr-13 allow same fields but different order
% 12-may-13 recursion for vectors
% 28-oct-13 default concatenation dimension added

function [ out rc ] = struct_cat( in1, in2, dim )
rc = NaN;
out = [];
if nargin < 2 || ~isa( in1, 'struct' ) || ~isa( in2, 'struct' )
    % not structures
    if isa( in1, 'struct' )
        if length( in1 ) > 1
            in2 = in1( 1 );
            in1 = struct_cat( in1( 2 : end ) );
        else
            out = in1;
            return
        end
    end
end
if nargin < 3 || isempty( dim )
    dim = [];
end
if length( in1 ) ~= 1 || length( in2 ) ~= 1
    % not single-element strucutres
    return
end
out = in1;
fieldnames = fields( in1 );
if ~isequal( sort( fieldnames ), sort( fields( in2 ) ) )
    % dissimilar structures
    return
end
m = length( fieldnames );
rc = zeros( m, 1, 'single' );

for i = 1 : m
    s1 = size( in1.( fieldnames{ i } ) );
    if all( s1 == 0 )
        rc( i ) = 1;
        out.( fieldnames{ i } ) = in2.( fieldnames{ i } );
        continue
    end
    s2 = size( in2.( fieldnames{ i } ) );
    ns1 = length( s1 );
    ns2 = length( s2 );
    if ns1 < ns2
        s1 = [ s1 ones( 1, ns2 - ns1 ) ];
    elseif ns1 > ns2
        s2 = [ s2 ones( 1, ns1 - ns2 ) ];
    end
    
    idx = s1 == s2;
    if sum( idx ) == 0 % nothing of identical size: do nothing
        if isa( in1.( fieldnames{ i } ), 'cell' ) % no nested checking for now
            rc( i ) = -1;
        else
            rc( i ) = 0;
        end
    elseif sum( ~idx ) == 1; % one field of dissimilar size: concat along it
        rc( i ) = 1;
        out.( fieldnames{ i } ) = cat( find( ~idx ), in1.( fieldnames{ i } ), in2.( fieldnames{ i } ) );
    elseif sum( idx ) > 1 && ~isempty( dim ) && idx( dim ) == 1 % default dimension
        rx( i ) = 1;
        out.( fieldnames{ i } ) = cat( dim, in1.( fieldnames{ i } ), in2.( fieldnames{ i } ) );
    else % more than 1 field of dissimilar : do nothiing (keep the value from in1)
        rc( i ) = sum( idx );
    end
end
return

% EOF
