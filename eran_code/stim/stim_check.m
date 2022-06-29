% stim_check        check an existing strucutre for format
%
% tf = stim_check( s )
% 
% see also          stim_make, stim_get

% 18-jan-13


function tf = stim_check( s )

tf = false( 1 );
if isa( s, 'struct' ) && isequal( fieldnames( s ), fieldnames( stim_make ) )
    tf = 1;
end

% addendum:
if ~isa( [ s.generator ], 'cell' )
    tf = 0;
end

return