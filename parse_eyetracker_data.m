%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse EyeTracker Data file 
%
% Program name: parse_eyetracker_data.m
%
% This program parses eyetracker data file, removes first 26 text lines
% and replaces sync pulse column with 0,1,2 and outputs a file which is
% used by program "eyetrack_process.m". 
%
% Entry Requirements:
%    Input file:  modify variable "fileStr" Example: fileStr='cg021810';
%    Output file: "fileStr" string appended with "_mod" string and mat file generated.
%           Example: output file: cg021810_mod.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%%%%%%%%%%%%%%%%%%%%% User input ****************************
%fileStr='cg021810';
%fileStr='mas02142011';
%fileStr='CG040710';
%fileStr='rnd021810';
%fileStr='rnd040710_test';

%fileStr=''; %Enter new file string without extension

N1=26; %Make sure first 26 lines of eyetracker data file have to be removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[filename, pathname, filterindex] = uigetfile('*.txt', 'Pick an eye-tracker file');
inFileStr=strcat(pathname,filename);
outFileStr=strcat(pathname,filename(1:find(filename=='.')-1),'mod.txt');
matFileStr=strcat(pathname,filename(1:find(filename=='.')-1),'_mod.mat'); 

inFile=fopen(inFileStr,'r');
outFile=fopen(outFileStr,'w');


%Read first 4 lines of XML file which do not contain parameters and values
for lineNum=1:26,
    fgets(inFile);
end
				
eofFlag=0;lineNum=1;
while eofFlag==0, 
    
    temp1=fscanf(inFile, '%d', 1);
    if feof(inFile)==0,

    fprintf(outFile,'%d',temp1);
    fprintf(outFile,'\t');

    for count=1:4, %TotalTime DeltaTime	X_Gaze	Y_Gaze
    temp2=fscanf(inFile, '%f', 1);
    fprintf(outFile,'%f',temp2);
    fprintf(outFile,'\t');
    end
    
    temp1=fscanf(inFile, '%d', 1); %Region
    fprintf(outFile,'%d',temp1);
    fprintf(outFile,'\t');

    for count=1:2, %PupilWidth	PupilHeight
    temp2=fscanf(inFile, '%f', 1);
    fprintf(outFile,'%f',temp2);
    fprintf(outFile,'\t');
    end
    
    temp1=fscanf(inFile, '%d', 1); %Quality	
    fprintf(outFile,'%d',temp1);
    fprintf(outFile,'\t');
    
    
    temp2=fscanf(inFile, '%f', 1); %Fixation
    fprintf(outFile,'%f',temp2);
    fprintf(outFile,'\t');

    temp1=fscanf(inFile, '%d', 1); %	Count	
    fprintf(outFile,'%d',temp1);
    fprintf(outFile,'\t');
    
    str1=fgets(inFile); %Marker

    if strncmp(str1(2),'S',1)==1,
      fprintf(outFile,'1');   %replace "S" with "1"
    else
        if strncmp(str1(2),'s',1)==1,
          fprintf(outFile,'2');  %replace "s" with "2"
        else
            fprintf(outFile,'0'); %replace everything else with "0"
        end
    end
    
    fprintf(outFile,'\n');
    
    lineNum=lineNum+1;
    if mod(lineNum,10000)==0,
        lineNum
    end
    else
       eofFlag=1; 
    end
end

fclose(outFile);


array1=dlmread(outFileStr);
size(array1)

eyetrack_data=array1;
save(matFileStr,'eyetrack_data');

fclose all;


