% struct_merge          merge dissimilar structures
%
% CALL                  [ out rc ] = struct_merge( in1, in2 )
%
% GETS                  in1, in2  input structures
%
% RETURNS 
%                       out       output structure
%                       rc        return code:
%                                     1       successful merge
%                                     0       similar structures
%                                     NaN     error( not structures, empty arguments...)
% 
% CALLS                 nothing
%
% See also              struct_cat, struct_compare, struct_extract, struct_merge, struct_reorder, struct_select

% 18-apr-13 ES

% revisions
% 02-mar-15 in case of shared fields, in1 takes precedence
% 10-jun-20 cleaned up

function [ out, rc ] = struct_merge( in1, in2 )

rc                      = NaN;
out                     = [];
if nargin < 2 || ~isa( in1, 'struct' ) || ~isa( in2, 'struct' )
    return
end

fieldnames1                     = fields( in1 );
fieldnames2                     = fields( in2 );
sharedfields                    = intersect( fieldnames1, fieldnames2 );
if ~isempty( sharedfields )
    in2                         = rmfield( in2, sharedfields );
    if isempty( in2 )
        return
    else
        fieldnames2             = setdiff( fieldnames2, sharedfields ); 
    end
end

out                             = in1;
m                               = length( fieldnames2 );

for i                           = 1 : m
    out.( fieldnames2{ i } )    = in2.( fieldnames2{ i } );
end
rc                              = 1;

return

% EOF
