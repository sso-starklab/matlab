% xmlfindchild      find a tagged child
%
% idx = xmlfindchild( s, tag )
%
% s         nested strucutre
% tag       string

% this a callback from xmlupdate and xmlgetchans

% 02-jan-12 ES


function idx = xmlfindchild( s, tag )
ok = 0;
idx = [];
for i = 1 : length( s.child )
    if strcmp( s.child( i ).tag, tag )
        ok = 1;
        idx = [ idx i ];
    end
end
if ~ok
    idx = NaN;
end
return