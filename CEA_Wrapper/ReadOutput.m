function [ data ] = ReadOutput( fileStr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    IOout = fopen(fileStr,'r');
    tline = fgets(IOout) ;
    disp(tline)
    while ischar(tline)
            
        %Finds the line containing the O/F number
        if any(strfind(tline,'NOHARM(I)'),1)
            disp(tline)
            found = tline ;
            ii = 0 ;
            % Separates the line in strings
            tempcell = strread(found, '%s', 'delimiter','=') ;
            out.noharm = tempcell{end}(1) ;
        end

%         if any(strfind(tline,'NO.          OMEGA          RAD./SEC.      HERTZ'),1)
%             tline = fgets(fid);
%             tline = fgets(fid);
%             C_data = textscan(fileID,,'%d %f %f %f');
%             C_data(1);
%         end

        tline = fgets(IOout);
        disp(tline)
    end

    data = fclose(IOout);
    
end

