
filebase = 'G:\mice\mF84\dat\mF84_03\mF84_03';
s2sfname = [filebase, '.s2s'];
load (s2sfname, 's2s', '-mat');

cch_tag = s2s.ccg - s2s.pred;
nunits = length (s2s.nspks);

for i=1:nunits
    cch_tag (:, i, :) = cch_tag (:, i, :) / s2s.  (i);
end

stg = cch_tag*1000;

idx = s2s.t>=-20 & s2s.t<=20;
figure, bar (s2s.t(idx), stg(idx, 3, 9));
