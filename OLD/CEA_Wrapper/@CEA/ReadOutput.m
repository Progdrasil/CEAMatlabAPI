function [ data ] = ReadOutput(obj, fileStr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %Debug variables
%     clear all;
%     close all;
%     clc;
%     currentpath = which('runWrapCEA.m');
%     [pathstr,name,ext] = fileparts(currentpath);
%     fileStr = strcat(pathstr,'/wrapper.dat');
    
    IOout = fopen(fileStr,'r');
    tline = fgets(IOout) ;
    test = false;
    perfo = false;
    while ischar(tline)
            
        %Finds the line containing the O/F number
        if any(strfind(tline,'O/F'),1)
            partstr = textscan(tline,'%s');
            data.OF = str2num(partstr{end}{2}) ;
            test = true;
        elseif test
            parts = textscan(tline,'%s');
%             celldisp(parts)
            if any(strfind(tline,'MOLE'),1)
                test = false;
            else
                if  ~isempty(parts{1})
                    if any(strfind(tline,'PERFORMANCE'),1)
                            perfo = true;
                    end
                    if length(parts{1})>4
                        if perfo
                            sub = 2;
                        else
                            sub = 3;
                        end
                            
                        N = length(parts{1}) -sub;
                        parts{1}{1} = regexprep(parts{1}{1},'[,/()*]','');
                        j = 1;
                        for i = N:length(parts{1})
                            data.(parts{1}{1})(j) = str2num(parts{1}{i});
                            j = j+1;
                        end
                    end
                end
            end
        end
        tline = fgets(IOout);
    end

    fclose(IOout);
    obj.data = data;
end

