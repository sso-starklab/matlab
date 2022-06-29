% region code:	function output 
% Neocortex     1
% before CA1	2
% CA1           3
% after CA1     4
% Unknown       0
% TBD           6

% in the txt file:
% same as function output, with the exception of linear probes: 
%                                      in the linear probe sessions, the txt file contains the ripple channel, 
%                                      and the function then caculates for each shankclu whether it is
%                                      above, below or on the pyramidal layer

function [regions] = add_region(sessions, locations ,regions_table_path)
    
    regions = zeros(length(sessions),1);
    var_names = {'Animal', 'Session', 'Region'};
    tab = readtable(regions_table_path);
    tab.Properties.VariableNames = var_names;
    
    locations = num2cell(locations);
    regions = cellfun(@(x,y) reg(x, y, tab), sessions, locations, 'UniformOutput', false);
    regions = cell2mat(regions);
    
end

function region = reg(session, location, tab)

    if strcmp(session(1:end-3), 'mK01') || strcmp(session(1:end-3), 'mO251') || strcmp(session(1:end-3), 'mV99')
        pyr = tab.Region(strcmp(tab.Session, session));
        if pyr == 0                         % Neocortex
            region = 1;
        elseif pyr == 33                    % After CA1
            region = 4;
        else
            if location > (pyr + 10)        % Neocortex
                region = 1;
            elseif location < (pyr - 10)    % After CA1
                region = 4;
            else
                region = 3;                 % CA1
            end
        end
    elseif sum(strcmp(tab.Session, session)) == 0
        region = 6;                         % not affiliated yet 
    else
        region = tab.Region(strcmp(tab.Session, session));
    end    
end

% example to run the code:
% sst.region = add_region(sst.filebase, sst.geo_com ,'region_table.txt');