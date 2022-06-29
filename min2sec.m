% min2sec           minutes to seconds
%
% convert minutes to seconds. works on each element of an array
%
% sec = min2sec( mins )
%
% see also sec2min

% 08-apr-12 ES

function sec = min2sec( mins )
sec = floor( mins ) * 60 + mod( mins, 1 ) * 100;
return