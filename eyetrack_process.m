%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process EyeTracker Data file 
%
% Program name: eyetrack_process.m
%
% This program processes eyetracker data file as follows.
%   1. Identifies the start and end of each scan.
%   2. Removes blink portions based on pupil height criteria 
%     provided percentage of blink segments < a certain threshold.
%       in some data, pupil height and aspect ratio are highly unreliable.
%   3. Runs a 2-pole IIR High pass Filter at cut-off frequency 0.025Hz to remove
%      slow drift due to head movements.
%   4. Runs a low pass filter with 7 Hz  cutoff frequency: FIR filter with order 42.
%       to cut off high frequency noise.
%   5. Performs Saccade detection. [200ms for saccade onset; 20-200ms
%   saccade duration assumed.]
%   6. Prints X,Y data to file along with blink/saccade information [
%    _scanInfo_ files]
%
% Entry Requirements:
%    Run "parse_eyetracker_data.m" on desired data file first and generate output
%    files with extension "_mod.mat" . Then run this program after specifying
%     name of the file in string "fileStr0" and setting other user options.
%     Then use "parse_eyetracker_outdata.m" to get an example of how to
%     scan output files with "_scanInfo" string generated by this program.
%
%    Run the program as follows: eyetrack_process(fileNum)
%      fileNum=0   => enter new input file name in variable "newFileStr"
%      fileNum=1   =>  input file name='cg_eyetrack_data_021810_mod.mat'
%      fileNum=2   =>  input file name='cg_eyetrack_data_021810_saccade_4.mat'
%                   saccades simulated
%      fileNum=3   =>  input file name='mas02142011_mod.mat'
%                       Pupil height and aspect ration are highly noisy
%                       hence blink detection unreliable
%      fileNum=4   =>  input file name='CG040710_mod.mat'
%                       aspect ratio noisy in some segments hence aspect
%                       ration not used for blink detection
%      fileNum=5   =>  input file name='rnd021810_mod.mat'                
%      fileNum=6   =>  input file name='rnd040710_mod.mat'
%           lots of adjacent blinks in some segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eyetrack_process(fileNum)


close all; %closes all figures

%%%%%%%%%%%%%%%%%%%%%% User input required ****************************

params.sampleRate = 60; %eyetracker sample rate in Hz

%Sync pulse information
syncPulsePeriodMin=5/6; %5/6 seconds
syncPeriodMin=syncPulsePeriodMin*params.sampleRate; 
syncPulsePeriod=2; % 2 seconds for fileNum=2
syncPeriodMax=params.sampleRate*syncPulsePeriod*4; %User input
numPulses=378; %160; %number of sync pulses per scan

xyWindow=1;  %User specified fixation windowl within 1 degree of visual arc
BlinkOverride=0; % set this to 1 to perform blink removal always irrepsective of noise/jitter in pupilheight
newFileStr=''; %enter the new file name without extension. example 'cg_eyetrack_data_021810'

%%%%%%%%%%%%%%%%%%%%%% End of User input required****************************

%%%%%%%%%%%%%%%%%%%%%% User input optional ****************************

plotFigure=1; %plot X,Y every scan optionally
plotFigure_2=1; %plot blink related plots and misc
printFile=0; %set this to 1 to print data to file

includeLPF=1; %include low pass filter
includeHPF=1; %include high pass filter

aspectRatioThr=1.25; % ratio of pupil width/pupil height
aspectRatioThr2=0.8;
pupilHeightThr=0.75;
aspectRatioMemory=10;
blinkMemory=60; % to capture several adjacent blinks
blinkThr=0.5;  %0.1-0.2 for normal segments; 0.9 for abnormal data
nPtsRemove = 5*2; %number of data points before and after each eyeblink to remove

timeIndex=2;
xDataIndex=4;
yDataIndex=5;
pupilWidthIndex=7; %data column number for pupil width
pupilHeightIndex=8;%data column number for pupil height

sampleDiffMean=4;
sampleDiffMax=4+2;
numGuardSamples=60; %used by sync pulse detection

%%%%%%%%%%%%%%%%%%%%%% End of User input optional****************************



%check Operating system; doesn't seem necessary
%if isunix==1,
%  filePathStr='dat/';
%else
%    filePathStr='dat\';
%end

filePathStr='dat/'; %this works for both unix and windows

%%New File
if fileNum==0,
fileStr0=strcat(filePathStr,newFileStr);
end

%%%%%%% Test with existig files
if fileNum==1,
fileStr0=strcat(filePathStr,'cg_eyetrack_data_021810_mod');
numPulses=378; 
end

if fileNum==2,
fileStr0=strcat(filePathStr,'cg_eyetrack_data_021810_saccade_4'); 
%same as cg_eyetrack_data_021810 but with saccades artificially simulated
numPulses=378; 
end

if fileNum==3,
fileStr0=strcat(filePathStr,'mas02142011_mod'); %.mat';
numPulses=160; 
end

if fileNum==4,
fileStr0=strcat(filePathStr,'CG040710_mod'); % due to aspectRatio jitters
numPulses=378;
end
if fileNum==5,
fileStr0=strcat(filePathStr,'rnd021810_mod');
numPulses=378;
end
if fileNum==6,
fileStr0=strcat(filePathStr,'rnd040710_mod'); %scan 3,9 blinkMemory and genuine saccades clash
numPulses=378;
end
%genuine saccadic movements in scan beginnings affecting long term average?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Processing data')

fileStr=strcat(fileStr0,'.mat'); %'dat\cg_eyetrack_data_021810.mat';
startNum=0;

load(fileStr); %loads array eyetrack_data which had been stored in file

[numRows numCols]=size(eyetrack_data);
eyetrack_data_orig=eyetrack_data;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Identifies the start and end of each scan.
% Sync pulse detection and periodicity and substitutes missing pulses

startSeg = find(eyetrack_data(:, 12)> 0, 1, 'first');
[syncUpArray] = find(eyetrack_data(:, 12)==1);
[syncDownArray] = find(eyetrack_data(:, 12)==2);


N1=min([length(syncUpArray) length(syncDownArray) ]);
sampleDiffArray=syncDownArray(1:N1)-syncUpArray(1:N1);

%substitute missing pulses

N2=max([length(syncUpArray) length(syncDownArray) ]);
count1=1; upCount=1; downCount=1;
for count=1:length(syncDownArray),
    if abs(syncDownArray(downCount)-syncUpArray(upCount))<=sampleDiffMax,
        syncUpDownArray(count1)=syncDownArray(downCount);
        count1=count1+1;upCount=upCount+1;downCount=downCount+1;
    else
        sampleDiff=syncDownArray(downCount)-syncUpArray(upCount);
        if sampleDiff<0, %up pulse missing
            syncUpDownArray(count1)=syncDownArray(downCount);
                sampleIndex1=syncUpDownArray(count1); %sampleIndex1  
                eyetrack_data_mod(sampleIndex1, 12)=1;
            count1=count1+1;             downCount=downCount+1;
    
        else %down pulse missing
            syncUpDownArray(count1)=syncUpArray(upCount)+sampleDiffMean;
                sampleIndex1=syncUpDownArray(count1); %sampleIndex1  
                eyetrack_data_mod(sampleIndex1, 12)=2;
            upCount=upCount+1;
            count1=count1+1;
        end
    end

                   
end


syncPeriodArray=syncUpDownArray(2:end)-syncUpDownArray(1:end-1);
disp('Average Sync Pulse Period in seconds')
mean(syncPeriodArray(1:numPulses/2))/params.sampleRate

if plotFigure_2==1,
figure(11)
hold off
plot(syncPeriodArray(1:numPulses/2)/params.sampleRate)
title('Sync Pulse Period for first scan');
xlabel(' sync pulse number');
end
    disp('Enter 1 if  you wish to change the sync pulse period.');
    disp('Enter 0 if  you do not wish to change the sync pulse period.');
    disp('Default is 1-2 seconds ');
    inChar=input('');
    if inChar==1,
        disp('Enter new sync pulse period');
        syncPeriodMax=input('');
        syncPeriodMax=syncPeriodMax*params.sampleRate*4;
    end
    

[scanEndIndexArray] = find(syncPeriodArray>syncPeriodMax);

%Identify start and end of each scan

syncPeriodArrayMod=[];start1=1;
scanStartEndArray(1)=syncUpDownArray(1);

for count=1:length(scanEndIndexArray), %syncUpDownArray),
    
    syncPeriodArrayMod=[syncPeriodArrayMod syncPeriodArray(start1:scanEndIndexArray(count)-1)];
    start1=scanEndIndexArray(count)+1;
    
    sampleIndex1=syncUpDownArray(scanEndIndexArray(count))+syncPeriodMax;
    endSampleXCurrentScan=eyetrack_data_mod(sampleIndex1, xDataIndex);
    sampleIndex2=syncUpDownArray(scanEndIndexArray(count)+1)-numGuardSamples;
    startSampleXNextScan=eyetrack_data_mod(sampleIndex2, xDataIndex);
    slope=(startSampleXNextScan-endSampleXCurrentScan)/(sampleIndex2-sampleIndex1);

    scanStartEndArray(count*2)=sampleIndex1;
    scanStartEndArray(count*2+1)=sampleIndex2;
       
    endSampleYCurrentScan=eyetrack_data_mod(sampleIndex1, yDataIndex);
    startSampleYNextScan=eyetrack_data_mod(sampleIndex2, yDataIndex);
    slope=(startSampleYNextScan-endSampleYCurrentScan)/(sampleIndex2-sampleIndex1);
    
end

scanStartEndArray((length(scanEndIndexArray)+1)*2)=syncUpDownArray(end)+syncPeriodMax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check if blinking removal needs to be done
% remove blinks only if % of blinks is below a certain threshold
%   2. Removes blink portions based on pupil height criteria 
%     provided percentage of blink segments < a certain threshold.
%       in some data, pupil height and aspect ratio are highly unreliable.

aspectRatio=eyetrack_data(:,pupilWidthIndex)./eyetrack_data(:,pupilHeightIndex);

%concatenate relevant data from scans
pupilHeightArray=[];
aspectRatioArray=[];
for scanNum=1: length(scanStartEndArray)/2,
   pupilHeightArray=[pupilHeightArray; eyetrack_data(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2), pupilHeightIndex)];
   aspectRatioArray=[aspectRatioArray; aspectRatio(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))];
end

meanPupilHeightArray=mean(pupilHeightArray);

array1=find(pupilHeightArray < meanPupilHeightArray*pupilHeightThr);
array2=find(aspectRatioArray > aspectRatioThr); %Blink

array3=find(aspectRatioArray < aspectRatioThr); %valid fixation
% resolve aspect ratio jitters by including memory
array4=array3(2:end)-array3(1:end-1);
count1=0;count2=0;
for count=1:length(array4),
   if array4(count)==1,
      count1=count1+1;   
   else
       if count1>aspectRatioMemory,
           array5(count2+(1:count1+1))=array3(count-count1:count);
           count2=count2+count1+1;
       end
        count1=0;

   end
end
meanPupilHeightArray=(mean(pupilHeightArray(array5))+max(pupilHeightArray(array5)))/2;
array6=find(pupilHeightArray < meanPupilHeightArray*pupilHeightThr); %0.5;
maxVal=max([length(array6)/length(pupilHeightArray) (length(aspectRatioArray)-length(array5))/length(aspectRatioArray)]);

disp('Fraction of Data classified as blinks')
maxVal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%CODE FOR EYEBLINK REMOVAL%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [blinkArray]=find(aspectRatio>aspectRatioThr);
    [blinkArray2]=find(aspectRatio<aspectRatioThr2);
    numSamples=length(eyetrack_data(:,pupilHeightIndex))-length(blinkArray);

    meanPupilHeight=mean(eyetrack_data(:,pupilHeightIndex));
    meanPupilHeightThr=meanPupilHeight*pupilHeightThr;
    blinkArray0=find(eyetrack_data(:,pupilHeightIndex) < meanPupilHeightThr);

    %removeArray=sort([ blinkArray; blinkArray2 blinkArray0]); 
    %Aspect ratio  not used because aspect ratio is jittery and unreliable in some
    %segments in CG040710.txt
    removeArray=sort([ blinkArray0]);
    removePts=unique(removeArray);

    array1=removePts(2:end)-removePts(1:end-1);
    array2=find(array1>blinkMemory);
    
    if maxVal >  blinkThr && BlinkOverride==0,
            disp('No blinking removal due to noisy pupil height/aspect ratio');
    else
        disp('Blinking removal being done');
        eyetrack_data(removePts(1)-nPtsRemove:removePts(array2(1))+nPtsRemove,xDataIndex) = eyetrack_data(removePts(1)-nPtsRemove-1,xDataIndex);
        eyetrack_data(removePts(1)-nPtsRemove:removePts(array2(1))+nPtsRemove,yDataIndex) = eyetrack_data(removePts(1)-nPtsRemove-1,yDataIndex);
    end
    blinkStartEndArray(1)=removePts(1)-nPtsRemove;
    blinkStartEndArray(2)=removePts(array2(1))+nPtsRemove;
            
    for count=1:length(array2)-1,

        if maxVal >  blinkThr && BlinkOverride==0,
        else
            eyetrack_data(removePts(array2(count)+1)-nPtsRemove:removePts(array2(count+1))+nPtsRemove,xDataIndex) = eyetrack_data(removePts(array2(count)+1)-nPtsRemove-1,xDataIndex);
            eyetrack_data(removePts(array2(count)+1)-nPtsRemove:removePts(array2(count+1))+nPtsRemove,yDataIndex) = eyetrack_data(removePts(array2(count)+1)-nPtsRemove-1,yDataIndex);
        end       
        blinkStartEndArray(count*2+1)=removePts(array2(count)+1)-nPtsRemove;
        blinkStartEndArray(count*2+2)=removePts(array2(count+1))+nPtsRemove;
            
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if plotFigure_2==1,

    %dataIndex=xDataIndex;
    dataIndex=yDataIndex;
    
    array11=eyetrack_data(:,dataIndex);

    array12=zeros(1,length(eyetrack_data));
    for count=1:length(blinkStartEndArray)/2,
        array12(blinkStartEndArray(count*2-1):blinkStartEndArray(count*2))=4;
    end

    N1=syncUpDownArray(1); %39100
    N0=N1; N3=1000;
    
    figure(10)
    hold off
    for N0=N1:1000:N1, %length(eyetrack_data),
        hold off
        plot(eyetrack_data_orig(N0+(1:N3), dataIndex),'r-')
        hold on
        plot(array11(N0+(1:N3)),'g-')
        plot(array12(N0+(1:N3)),'k-')

        plot(eyetrack_data_orig(N0+(1:N3), 8),'b-');
        plot(aspectRatio(N0+(1:N3)),'m-');

        title('Red: original; Green: Blink Removed; Black: Blink pulse; Blue: PupilHeight; Magenta:AspectRatio');
        xlabel('First 1000 samples in Scan 1');

        disp('Hit Enter to continue');

        pause
    end

end %end if plotFigure_2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   3. Runs a 2-pole High pass Filter at cut-off frequency 0.025Hz to remove
%      slow drift due to head movements.
if includeHPF==1,
%High Pass Filter at 0.025Hz [0.01-0.1 Hz]

%Need Signal Processing toolbox
%fs=60;
%N=2; %10; %5; similar to 40 for FIR %20; %50; %30;
%Wn=[0.025]*2/fs; 
%ftype = 'low';
%[B,A] = butter(N,Wn,ftype);

    B=[0.171030589091181   0.342061178182362   0.171030589091181]*1.0e-005;
    A =[1.000000000000000  -1.996297601769122   0.996304442992686];

    delay=0;
    sig=eyetrack_data(:,xDataIndex)';
    temp = [sig zeros(1,delay)];
    temp=filter(B,A,temp);
    filtered_xValues=temp(delay+1:delay+length(sig));
    x_hpf=sig-filtered_xValues;

    sig=eyetrack_data(:,yDataIndex)';
    temp = [sig zeros(1,delay)];
    temp=filter(B,A,temp);
    filtered_yValues=temp(delay+1:delay+length(sig));
    y_hpf=sig-filtered_yValues;

else
    x_hpf=eyetrack_data(:,xDataIndex);
    y_hpf=eyetrack_data(:,yDataIndex); 
end


%   4. Runs a low pass filter with 7 Hz  cutoff frequency: FIR filter with order 42.
%       to cut off high frequency noise.
if includeLPF==1,
%Low Pass Filter at 7 Hz [6-8 Hz] using FIR

%LPF
%Need Signal Processing toolbox
%fs=60;
%[n,fo,mo,w] = firpmord([6 8], [1 0], [.1 .005], fs);
%if mod(n,2)==0, n=n+1;end
%delay=(n-1)/2;
%B2=firls(n,fo,mo,w); %firls is better


    A2=1;
    delay=21;
    B2 = [0.003017468846531   0.005237721033048   0.005739582407630 ...
       0.003254020008419  -0.002108963973915  -0.008372877677971 ...
       -0.012229747971581  -0.010590536102917  -0.002513599310528 ...
      0.009532344978106   0.020002280318416   0.022459377270665 ...
       0.012889590534292  -0.007355469052629  -0.030614242712175 ...
       -0.044885534448784  -0.038326091507919  -0.004581379533159 ...
      0.053391771347861   0.122941757733578   0.185338011012221 ...
       0.222195582649439   0.222195582649439   0.185338011012221 ...
       0.122941757733578   0.053391771347861  -0.004581379533159...
       -0.038326091507919  -0.044885534448784  -0.030614242712175...
      -0.007355469052629   0.012889590534292   0.022459377270665 ...
      0.020002280318416   0.009532344978106  -0.002513599310528 -0.010590536102917  -0.012229747971581  ...
       -0.008372877677971  -0.002108963973915   0.003254020008419   ...
       0.005739582407630   0.005237721033048 0.003017468846531];
   
    sig=x_hpf;
    temp = [sig zeros(1,delay)];
    temp=filter(B2,A2,temp);
    x_lpf=temp(delay+1:delay+length(sig));

    sig=y_hpf;
    temp = [sig zeros(1,delay)];
    temp=filter(B2,A2,temp);
    y_lpf=temp(delay+1:delay+length(sig));

else
   x_lpf=x_hpf;
   y_lpf=y_hpf;
end

if plotFigure_2==1,
    %Plot high pass and low pass filtered data
    figure(12)
    hold off
    plot(eyetrack_data(18000:100000,4),'r-')
    hold on; plot(x_hpf(18000:100000),'b-');
    plot(x_lpf(18000:100000),'g-');
    title('Red: Original with blinks removed; Blue: High Pass Filtered; Green: Low pass Filtered');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Saccade detection

%Short Term  Low pass filter
numSamplesShortTerm=6; % 6 samples at Ts=1/60 sec = 100ms
alpha=1/numSamplesShortTerm; 

B3=alpha;
A3=[1 -(1-alpha)];
x_lpf_short_term=filter(B3,A3,x_lpf);
y_lpf_short_term=filter(B3,A3,y_lpf);

shortTermThr=2;
saccadeMemory=10;

cols=max(scanStartEndArray(2:2:end)-scanStartEndArray(1:2:end-1));
rows=length(scanStartEndArray)/2;
xData=zeros(rows,cols);
yData=zeros(rows,cols);
fixWinArray=zeros(length(scanStartEndArray)/2,100);
saccadeCount=1;

for scanNum=1:length(scanStartEndArray)/2,
    array_x=abs(x_lpf_short_term(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(x_lpf));
    array_y=abs(y_lpf_short_term(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(y_lpf));

    array2_x=find(array_x>shortTermThr);
    array2_y=find(array_y>shortTermThr);
    array2=unique(sort([array2_x array2_y]));
    
    if length(array2)>0,
        array3=array2(2:end)-array2(1:end-1);
        array4=find(array3>saccadeMemory);
        saccadeStart(scanNum,1)=scanStartEndArray(scanNum*2-1)+array2(1);
        saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+array2(1);
        saccadeCount=saccadeCount+1;
        count1=1;
        
        for count=1:length(array4),
            saccadeEnd(scanNum,count)=scanStartEndArray(scanNum*2-1)+array2(array4(count));
            saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+array2(array4(count));
            saccadeCount=saccadeCount+1;
            if array4(count)+1 <= length(array2),
                saccadeStart(scanNum,count+1)=scanStartEndArray(scanNum*2-1)+array2(array4(count)+1);
                saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+array2(array4(count)+1);
                saccadeCount=saccadeCount+1;
            end
            count1=count1+1;
        end

        saccadeEnd(scanNum,count1)=scanStartEndArray(scanNum*2-1)+array2(end);
        saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+array2(end);
        saccadeCount=saccadeCount+1;

     end %end if array2
 

    xyNum(scanNum)=length([scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2)]);
    xData(scanNum,1:xyNum(scanNum))=x_lpf(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(x_lpf);
    yData(scanNum,1:xyNum(scanNum))=y_lpf(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(y_lpf);
    
    %Statistics
    meanX(scanNum)=mean(xData(scanNum,1:xyNum(scanNum))); 
    meanY(scanNum)=mean(yData(scanNum,1:xyNum(scanNum))); 
    stdX(scanNum)=std(xData(scanNum,1:xyNum(scanNum))); 
    stdY(scanNum)=std(yData(scanNum,1:xyNum(scanNum))); 
    
    array10=find(abs(xData(scanNum,1:xyNum(scanNum)))>xyWindow);
    array11=find(abs(yData(scanNum,1:xyNum(scanNum)))>xyWindow);
    array12=unique([array10 array11]);
    if length(array12)<=100,
        fixWinArray(scanNum,1:length(array12))=array12;
    else
        fixWinArray(scanNum,1:100)=array12(1:100);
    end
    fracPtsOutsideWindow(scanNum)=length(array12)*100/xyNum(scanNum);
    
    if plotFigure==1,
    %Plot figures for each scan    
    figure(21);
    hold off
    plot(xData(scanNum,1:xyNum(scanNum)),yData(scanNum,1:xyNum(scanNum)),'x');
    title('Scatterplot: X, Y')
    figure(22);
    hold off
    subplot(4,1,1);
    hist(xData(scanNum,1:xyNum(scanNum)),500);
    title('X Histogram')
    subplot(4,1,2);
    hist(yData(scanNum,1:xyNum(scanNum)),500);
    title('Y Histogram')
    subplot(4,1,3);
    plot(xData(scanNum,1:xyNum(scanNum)));
    title('X')
    subplot(4,1,4);
    plot(yData(scanNum,1:xyNum(scanNum)));
    title('Y')

    if fileNum==2 && scanNum==1,
        figure(3)
        hold off
        plot(eyetrack_data_orig(19532+(-100:100),4));
        title('Simulated saccade')
    end

    disp('Mean of X data')
    meanX(scanNum)
    disp('Mean of Y data')
    meanY(scanNum)
    disp('Standard Deviation of X data')
    stdX(scanNum)
    disp('Standard Deviation of Y data')
    stdY(scanNum)
    disp('Fraction of data outside fixation window')
    fracPtsOutsideWindow(scanNum)
    scanNum
    disp('Hit Enter to continue'); 
    pause
    end %end if plotFigure

end %scanNum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print to file

if printFile==1,
    disp('Printing to File. This will take a few minutes');
    outFileStr0=strcat(fileStr0,'_scanInfo_');
    textStr1='Time';textStr2='X';textStr3='Y';textStr4='Blink';textStr5='Saccade';textStr6='Outside fixation window';
    N10=length(blinkStartEndArray);
    blinkArrayIndex=1;saccadeArrayIndex=1;
    [N11, N12]=size(fixWinArray);

    meanX(scanNum)
    meanY(scanNum)
    stdX(scanNum)
    stdY(scanNum)
    fracPtsOutsideWindow(scanNum)


    outFileStr=strcat(outFileStr0, 'numscans', '.txt');
    outFile=fopen(outFileStr, 'w');
    fprintf(outFile, '%d', length(scanStartEndArray)/2);    
    fclose(outFile);

    for scanNum=1:length(scanStartEndArray)/2,
    
    scanNumStr=sprintf('%d',scanNum);
    outFileStr=strcat(outFileStr0, scanNumStr, '.txt');
    outFile=fopen(outFileStr, 'w');

    fprintf(outFile,'%s\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n','NumSamples','Mean(X)','Mean(Y)','Std(X)','Std(Y)','% samples_outside_window');  
    fprintf(outFile,'%d\t\t',xyNum(scanNum));
    fprintf(outFile,'%f\t',meanX(scanNum));
    fprintf(outFile,'%f\t',meanY(scanNum));
    fprintf(outFile,'%f\t',stdX(scanNum));
    fprintf(outFile,'%f\t',stdY(scanNum));
    fprintf(outFile,'%f\n',fracPtsOutsideWindow(scanNum));

    
    fprintf(outFile,'%s\t\t\t%s\t\t\t%s\t\t%s\t%s\t%s\n', textStr1,textStr2,textStr3,textStr4,textStr5,textStr6);
    
    fixWinIndex=1;
    for count=1:xyNum(scanNum),
        dataIndex=scanStartEndArray(scanNum*2-1)+(count-1);
        fprintf(outFile,'%f\t',eyetrack_data_orig(dataIndex,timeIndex));
        fprintf(outFile,'%f\t',xData(scanNum,count));
        fprintf(outFile,'%f\t',yData(scanNum,count));
    
      if N10>0,
        blinkDataIndex=blinkStartEndArray(blinkArrayIndex);
        while blinkDataIndex<dataIndex && blinkArrayIndex<N10,
           blinkArrayIndex=blinkArrayIndex+1; 
           blinkDataIndex=blinkStartEndArray(blinkArrayIndex);
        end

        if  mod(blinkArrayIndex,2)==1, %blinkArrayIndex==1 ||
              fprintf(outFile,'0\t'); %no blink  
        else
           fprintf(outFile,'1\t'); %in the middle of blink 
        end
      else
            fprintf(outFile,'0\t'); %no blink
      end
        fprintf(outFile,'\t');
        
      if saccadeCount>1,
        saccadeDataIndex=saccadeStartEndArray(saccadeArrayIndex);
        while saccadeDataIndex<dataIndex && saccadeArrayIndex<saccadeCount-1,
           saccadeArrayIndex=saccadeArrayIndex+1; 
           saccadeDataIndex=saccadeStartEndArray(saccadeArrayIndex);
        end

        if  mod(saccadeArrayIndex,2)==1, %end of saccade saccadeArrayIndex==1 || 
           fprintf(outFile,'0\t'); %no saccade
        else
            fprintf(outFile,'1\t'); %in the middle of saccade
        end    
      else
            fprintf(outFile,'0\t'); %no saccade
      end
        
        fprintf(outFile,'\t');
        
        if fixWinArray(scanNum,fixWinIndex)>0 && fixWinIndex<100,
            if count==fixWinArray(scanNum,fixWinIndex),
                fprintf(outFile,'1\n');
                fixWinIndex=fixWinIndex+1;
            else
                fprintf(outFile,'0\n');
            end

        else
            fprintf(outFile,'0\n');
        end
        
     end %end count
    
    
    fclose(outFile);
    scanNum
    
    end %end scanNum
end %end if printFile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all;

return

