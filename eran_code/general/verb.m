% VERB      string to screen
%
% verb( txt, flag )

% 17-jul-12 ES

function verb( txt, flag )
if isa( txt, 'char' )
    if flag == 1
        fprintf( 1, '%s\n', txt )
    elseif flag == -1
        fprintf( 1, '%s', txt )
    end
end
return