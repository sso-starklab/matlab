% BURSTF            burst filter spike train
%
% GETS              X           spike times (vector)
%                   CUTOFF      cutoff number of samples for a burst
%
% RETURNS           X           spike time vector with all ISIs > CUTOFF
%                   IDX         of removed spikes

% 19-may-10 ES

% 28-feb-12 corrected idx output

function [ x, idx ] = burstf( x, cutoff )

isis = diff( x ); 
idx = find( isis <= cutoff ) + 1; 
x( idx ) = [];

return