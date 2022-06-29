%function [xml, rxml] = LoadXml(FileBase)
%loads the xml file using xmltools (have to have it in the path)
% rxml returns it's original layout - very messy structure but contains all
% the xml file contents.
% xml - is the ouput structure which is backwards compatible to LoadPar
% output, so you can use it instead ..also loads some usefull stuff -
% Anatomoical groups with Skips , Spike electrode groups
% more can be added later (e.g. parameters of the process scripts)
% this script is written for xml version 1.1 .. older version doesn't work.
% additions are welcome
%
% 02-jan-13 small modification, also check for
%           type/voltageRange/source/target/location/reward fields and output those
%           type/voltage range exist for all channels
%           location/reward are mapped to source/target


function [xml, rxml] = LoadXml( FileBase )
xml = struct;

%if xml was in the filebase by chance
FileBase0 = FileBase;
xmli = strfind(FileBase,'.xml');
if ~isempty(xmli)
    FileBase = FileBase(1:xmli-1);
end
if strcmp( [ FileBase '.xml' ], FileBase0 )
    FileBase = FileBase0;
    rxml = xmltools( FileBase );
elseif exist( [ FileBase '.prm.xml' ], 'file' )
    rxml = xmltools([FileBase '.prm.xml']);
else
    rxml = xmltools([FileBase '.xml']);
end

rxml = rxml.child(2);

% from this level all children are the different parameters fields
xmli = strfind(FileBase,'.prm');
if ~isempty(xmli)
    FileBase = FileBase(1:xmli-1);
end

xml.FileName = FileBase;

for i=1:length(rxml.child)

    switch rxml.child(i).tag
        
        case 'generalInfo'
            xml.Date = rxml.child(i).child(1).value; % date of xml file creation?

        case 'acquisitionSystem'
            xml.nBits = str2num(rxml.child(i).child(1).value); % number of bits of the file
            xml.nChannels = str2num(rxml.child(i).child(2).value);
            xml.SampleRate = str2num(rxml.child(i).child(3).value);
            xml.SampleTime = 1e6/xml.SampleRate; %to make backwards compatible
            xml.VoltageRange = str2num(rxml.child(i).child(4).value);
            xml.Amplification = str2num(rxml.child(i).child(5).value);
            xml.Offset =  str2num(rxml.child(i).child(6).value);
            
        case 'fieldPotentials'
            xml.lfpSampleRate = str2num(rxml.child(i).child.value);
            
        case 'anatomicalDescription'
            tmp = rxml.child(i).child.child;
            for grpI =1:length(tmp)
                nI = length(tmp(grpI).child);
                xml.AnatGrps(grpI).Channels = zeros( 1, nI );
                xml.AnatGrps(grpI).Skip = zeros( 1, nI );
                xml.AnatGrps(grpI).Type = cell( 1, nI );
                xml.AnatGrps(grpI).VoltageRange = zeros( 1, nI );
                xml.AnatGrps(grpI).MaxVoltage = zeros( 1, nI );
                xml.AnatGrps(grpI).Source = cell( 1, nI );
                xml.AnatGrps(grpI).Target = zeros( 1, nI );
                xml.AnatGrps(grpI).Wavelength = zeros( 1, nI );
                xml.AnatGrps(grpI).Power = zeros( 1, nI );
                xml.AnatGrps(grpI).Location = cell( 1, nI );
                xml.AnatGrps(grpI).Reward = cell( 1, nI );
                for chI=1:length(tmp(grpI).child)
                    xml.AnatGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(chI).value);
                    xml.AnatGrps(grpI).Skip(chI) = str2num(tmp(grpI).child(chI).attribs.value);
                    s = tmp( grpI ).child( chI ).child;
                    for cI = 1 : length( s )
                        if strcmpi(  s( cI ).tag , 'type' )
                            xml.AnatGrps(grpI).Type{chI} =  s( cI ).value;
                        end
                        if strcmpi( s( cI ).tag, 'voltageRange' )
                            xml.AnatGrps(grpI).VoltageRange(chI) =  str2num( s( cI ).value );
                        end
                        if strcmpi( s( cI ).tag, 'voltage' )
                            xml.AnatGrps(grpI).MaxVoltage(chI) =  str2num( s( cI ).value );
                        end
                        if strcmpi( s( cI ).tag, 'source' )
                            xml.AnatGrps(grpI).Source{chI} =  s( cI ).value;
                        end
                        if strcmpi( s( cI ).tag, 'target' )
                            xml.AnatGrps(grpI).Target(chI) =  str2num( s( cI ).value );
                        end
                        if strcmpi(  s( cI ).tag , 'wavelength' )
                            xml.AnatGrps(grpI).Wavelength(chI) =  str2num( s( cI ).value );
                        end
                        if strcmpi(  s( cI ).tag , 'power' )
                            xml.AnatGrps(grpI).Power(chI) =  str2num( s( cI ).value );
                        end
                        if strcmpi( s( cI ).tag, 'location' )
                            xml.AnatGrps(grpI).Location{chI} =  s( cI ).value;
                        end
                        if strcmpi( s( cI ).tag, 'reward' )
                            xml.AnatGrps(grpI).Reward{chI} =  s( cI ).value;
                        end
                    end
                end
            end
            
        case 'spikeDetection'
            if ~isempty(rxml.child(i).child)
                tmp =rxml.child(i).child.child;
                for grpI =1:length(tmp)
                    for chI=1:length(tmp(grpI).child(1).child)
                        xml.SpkGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(1).child(chI).value);
                    end
                    if length(tmp(grpI).child)>1
                        xml.SpkGrps(grpI).nSamples = str2num(tmp(grpI).child(2).value);
                        xml.SpkGrps(grpI).PeakSample = str2num(tmp(grpI).child(3).value);
                    end
                    if length(tmp(grpI).child)>3
                        xml.SpkGrps(grpI).nFeatures = str2num(tmp(grpI).child(4).value);
                    end
                    %backwards compatibility
                    xml.nElecGps = length(tmp);
                    xml.ElecGp{grpI} = xml.SpkGrps(grpI).Channels;
                end
            else
                xml.nElecGps = 0;
            end


        case 'programs'
            tmp = rxml.child(i).child;
            for i=1:length(tmp)
                if strcmp(tmp(i).child(1).value,'process_mhipass')
                    for j=1:length(tmp(i).child(2).child )
                        if strcmp(tmp(i).child(2).child(j).child(1).value,'frequency')
                            xml.HiPassFreq = str2num(tmp(i).child(2).child(j).child(2).value);
                            break
                        end
                    end
                end
            end
    end


end


% general recursive parsing will have to wait.
