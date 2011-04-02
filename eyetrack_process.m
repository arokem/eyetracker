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
%    files with extension "_mod.mat" . Then run this program after setting user options.
%
%
%    Run this program as follows:
%       Run "clear all; eyetrack_process(TR_array, numTR_array);"
%
%   TR_array: Array containing TR(Time of repetition) for each scan in seconds.
%   numTR_array: Array containing number of TR pulses for each scan.
%   NOTE: wrong values specified in TR_array will give wrong results.
%
%   Input filename is  specified using "uigetFile" command 
%   when it pops a window asking for input file name.
%
%    Example:  For file: cg021810.txt
%
%   clear all; 
%   TR_array=ones(1,13)*1;
%   numTR_array = ones(1,13)*378;
%   eyetrack_process(TR_array, numTR_array);
%
%   Issues: Array Indexing checks included. More checks may be needed with
%   new data files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eyetrack_process(TR_array, numTR_array)


close all; %closes all figures

%%%%%%%%%%%%%%%%%%%%%% User input required ****************************

params.sampleRate = 60; %eyetracker sample rate in Hz

%This variable is used in filtering 
fs=params.sampleRate;

xyWindow=1;  %User specified fixation window within 1 degree of visual arc

blinkThr=0.5;  %Threshold for percentage of classified Blinks in a given scan; 0.1-0.2 for normal segments; 0.9 for abnormal data
  %If Blink percentage is less than blinkThr, do blink removal per scan. Else dont.
BlinkOverride=0; % set this to 1 to perform blink removal always irrespective of noise/jitter in pupilheight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For really noisy data, if a very high percentage of samples are misclassified as Blinks, 
%or if too many Blinks go unsuppressed,  use this option below.
suppressBlinksOutsideSaccades=0; %set this to 1 to hard-limit X,Y amplitudes above blinkLimitThr, 
                            %outside saccades
blinkLimitThr=0.5;   %related to suppressBlinksOutsideSaccades option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SignalProcessingToolbox=0; % If user has signal processing toolbox, set this to 1 to change filter cutoffs

lpfBlinks=0; %set to 1 to Low pass filter Pupil height and Width. No better performance.


minScanIntervalSamples=params.sampleRate*10; %Minimum number of samples in between 2 scans
%ttlPulse period information; Period: TR-Time of Repetition; Used to catch
%errors in TR array arguments.
TR_min= 0.1; % Min value of TR: 0.1 seconds; 
TR_max=100; % Max value of TR: 10 seconds

%%%%%%%%%%%%%%%%%%%%%% End of User input required****************************



%%%%%%%%%%%%%%%%%%%%%% User input optional ****************************

printFile=0; %set this to 1 to print data to file
saveFigure=0; %set this to 1 to save Figures

pupilHeightThr=0.75; %If pupil Height is less than mean(PupilHeight)*pupilHeightThr, then it is blink.
aspectRatioThr=1.25; % ratio of pupil width/pupil height. If aspect ratio>aspectRatioThr, then it is blink
aspectRatioThr2=0.8; %Lower threshold for aspect ratio
blinkMemory=10; %If adjacent blinks are separated by < blinkMemory, consider it as a single blink.
        %Increase to 60  to capture several adjacent blinks..
        
nPtsRemove = 5*2; %number of data points before and after each eyeblink to remove

shortTermSaccadeThr=2;  
 %If short term average exceeds long term average by 2 degrees, then consider it as a valid saccade
saccadeMemory=10; %If saccades are separated by less than 10 samples, consider it as a single saccade
numSamplesShortTerm=6; % 6 samples at Ts=1/60 sec = 100ms

ttlPulsePeriod=2;  %assume 2 second period TR in automatic derivation

plotFigure=0; %plot X,Y every scan optionally
plotFigure_2=0; %plot blink related plots and misc
plotFigure_TTL=0; %set to 1 to plot TTL period in a scan

deriveTTLpulseInfo=0; % Set this to 1 to automatically derive TTL pulse information.
                      %User need not input TR arrays. Make sure TTL
                      %period(TR) input are reasonable values between TR_min and TR_max


maxNumMissingTTLpulses=4; % Maximum number of missing TTL pulses. Used for automatic derivation only
  %For some data, more than 4 pulses may be missing, in which case, scan
  %start and end will be identified incorrectly. Blink plots have scan
  %start plotted and can indicate if this is the case. Then increase
  %maxNumMissingTTLpulses to a higher number.[Example:arokem_20110110_mod.txt]
  %The tradeoff is: extra samples equal to maxNumMissingTTLpulses 
  % after the end of each scan will be considered as valid samples inside scan, 
  % which may increase fraction of samples outside fixation window.
  % NOTE that user provided number of TTL pulses in each scan cannot be
  % relied upon because of missing TTL pulses.
  
removeSamplesBetweenScans=1; %removes samples between scans
includeLPF=1; %include low pass filter
includeHPF=1; %include high pass filter

sampleDiffMean=4; %number of samples between TTL Up and Down pulse
sampleDiffMax=4+2;
numGuardSamples=60; %used by ttlPulse  detection

timeIndex=2; %data column number for time
xDataIndex=4; %data column number for X data
yDataIndex=5; %data column number for Y data
pupilWidthIndex=7; %data column number for pupil width
pupilHeightIndex=8;%data column number for pupil height
qualityIndex=9; %data column number for quality
ttlPulseIndex=12; %data column number for TTL pulse

%%%%%%%%%%%%%%%%%%%%%% End of User input optional****************************



%%%%%%%%%%%%%%%%%%%%%% Section 1: Check input arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorMsg=nargchk(2,2,nargin);
if isempty(errorMsg)==0,

    disp('Usage: eyetrack_process(TR_array, numTR_array)');
    disp('TR_array: Array containing TR(Time of repetition) for each scan in seconds');
    disp('numTR_array: Array containing number of TR pulses for each scan');
    disp('         ');
    disp('Enter 1 if you want the program to automatically derive TR information');
    inChar=input('');
    if inChar==1,
        deriveTTLpulseInfo=1;
    else
      return;
    end
else
    if length(TR_array)~=length(numTR_array),
       disp('Input same number of elements for both TR_array and numTR_array ');
       return;
    else
        numScans=length(TR_array);
        TR_max_array=find(TR_array> TR_max);
        TR_min_array=find(TR_array < TR_min);
        if isempty(TR_max_array)==0 || isempty(TR_min_array)==0,
            disp('TR array values outside min and max values');
            return;
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%% End of Section 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%% User input file name ****************************
[fileStr, filePathStr, filterindex] = uigetfile('*mat','MultiSelect','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename_end=find(fileStr=='.');
fileStr0=strcat(filePathStr, fileStr(1:filename_end-1));
fileStr=strcat(filePathStr,fileStr);



%%%%%%%%%%%%%%%%%%%%%% Start Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Processing data')


load(fileStr); %loads array eyetrack_data which had been stored in file



[numRows numCols]=size(eyetrack_data);

if numRows==0 || numCols < 12,
   disp('eyetrack_data wrong size.');
   return;
end

eyetrack_data_orig=eyetrack_data;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Section 2: Identifies the start and end of each scan.
    %option 1: user provided 
    % option 2: automatic TTL pulse detection and periodicity and substitutes missing TTL pulses
%          
%
% ttlPulseUpDownArray: stores sample index of TTL Down pulse after missing pulse substitution
% scanStartEndArray: sample index of start and end of each scan
%   scanStartEndArray[1]: start of Scan 1
%   scanStartEndArray[2]: end of Scan 1
%   scanStartEndArray[3]: start of Scan 2
%   scanStartEndArray[4]: end of Scan 2
%
% scanEndIndexArray: sample index of end of each scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[ttlPulseUpArray] = find(eyetrack_data(:, 12)==1); % sample index of TTL UP pulse "S"
[ttlPulseDownArray] = find(eyetrack_data(:, 12)==2); % sample index of TTL Down pulse "s"

%substitute missing TTL pulses

ttlPulseUpDownArray=[];
ttlPulseCount=1; upCount=1; downCount=1;

for count=1:length(ttlPulseDownArray),
  if upCount<=length(ttlPulseUpArray),
    if abs(ttlPulseDownArray(downCount)-ttlPulseUpArray(upCount))<=sampleDiffMax,
        ttlPulseUpDownArray(ttlPulseCount)=ttlPulseDownArray(downCount);
        ttlPulseCount=ttlPulseCount+1;upCount=upCount+1;downCount=downCount+1;
    else %missing TTL up or down pulse
        sampleDiff=ttlPulseDownArray(downCount)-ttlPulseUpArray(upCount);
        if sampleDiff<0, %up pulse missing
            ttlPulseUpDownArray(ttlPulseCount)=ttlPulseDownArray(downCount);
            sampleIndex1=ttlPulseUpDownArray(ttlPulseCount); 
            eyetrack_data(sampleIndex1, ttlPulseIndex)=1;
            ttlPulseCount=ttlPulseCount+1;             downCount=downCount+1;
    
        else %down pulse missing
            ttlPulseUpDownArray(ttlPulseCount)=ttlPulseUpArray(upCount)+sampleDiffMean;
            sampleIndex1=ttlPulseUpDownArray(ttlPulseCount); %sampleIndex1  
            eyetrack_data(sampleIndex1, ttlPulseIndex)=2;
            upCount=upCount+1;
            ttlPulseCount=ttlPulseCount+1;
        end
    end
  end
                   
end



if length(ttlPulseUpDownArray)<2, %No TTL pulses detected!
    disp('No TTL pulses in eyetracker File ')
    return;
end
    
eyetrack_data_length=length(eyetrack_data(:,1));

if  deriveTTLpulseInfo==0, %Use user inputs; 
    %scan end=first TTL pulse+ numTR_array(scanNum)*TR_array(scanNum)*params.sampleRate
    startCurrentScan=ttlPulseUpDownArray(1); upDownCount=1;

   for scanNum=1:numScans,
      scanStartEndArray(scanNum*2-1)=startCurrentScan;
      endIndex=min(round(startCurrentScan+numTR_array(scanNum)*TR_array(scanNum)*params.sampleRate),eyetrack_data_length);
      scanStartEndArray(scanNum*2)= endIndex;
      
      while ttlPulseUpDownArray(upDownCount)<endIndex+ minScanIntervalSamples && upDownCount<length(ttlPulseUpDownArray),
          upDownCount=upDownCount+1;
      end
      startCurrentScan=ttlPulseUpDownArray(upDownCount);
   end


else %automatic derivation of TTL pulse period; needs setting maxNumMissingTTLpulses adjusted for each enw data

    %TR: TTL Pulse period array
    ttlPulsePeriodArray=ttlPulseUpDownArray(2:end)-ttlPulseUpDownArray(1:end-1);


    if plotFigure_TTL==1, %Optionally Plot TTL Pulse period

        disp('Average ttlPulse Pulse Period in seconds')
        mean(ttlPulsePeriodArray(1:floor(ttlPulseCount/2)))/params.sampleRate
    
        figure(11)
        hold off
        plot(ttlPulsePeriodArray(1:ttlPulseCount/2)/params.sampleRate)
        title('TTL Pulse Period for first scan');
        xlabel(' TTL pulse number');

        disp('Enter 1 if  you wish to change the TTL pulse period.');
        disp('Enter 0 if  you do not wish to change the TTL pulse period.');
        disp('Default is 1-2 seconds ');
        inChar=input('');
        if inChar==1,
        disp('Enter new TTL pulse period');
        ttlPulsePeriodMax=input('');
        ttlPulsePeriodMax=ttlPulsePeriodMax*params.sampleRate*maxNumMissingTTLpulses;
        end
    end    


    
    ttlPulsePeriodMax=params.sampleRate*ttlPulsePeriod*maxNumMissingTTLpulses; %10; %User input
    [scanEndIndexArray] = find(ttlPulsePeriodArray>ttlPulsePeriodMax);
    numScans=length(scanEndIndexArray)+1;


       
    %Identify start and end of each scan

    scanStartEndArray(1)=ttlPulseUpDownArray(1);


    for count=1:length(scanEndIndexArray), 
        
        sampleIndex1=min(ttlPulseUpDownArray(scanEndIndexArray(count))+ttlPulsePeriodMax, eyetrack_data_length);
        sampleIndex2=max(1,ttlPulseUpDownArray(scanEndIndexArray(count)+1)-numGuardSamples);
 
        scanStartEndArray(count*2)=sampleIndex1;
        scanStartEndArray(count*2+1)=sampleIndex2;
           
    end

    scanStartEndArray((length(scanEndIndexArray)+1)*2)=min(ttlPulseUpDownArray(end)+ttlPulsePeriodMax,eyetrack_data_length);


end %deriveTTLpulseInfo



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Remove segments in between scans just in case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is technically NOT required for following reasons: 
%  1. 7 Hz Low pass filter is a 43-tap FIR filter which means after 50
%     samples, there is no memory. What happens between the scans for several
%      minutes is irrelevant.
%  2. 0.025 Hz High Pass Filter and Saccade detection also have a time constant (decay time) of less 
%     than 100 samples(memory) and hence what happens between the scans for several minutes is irrelevant.
%  3. Blink removal has no memory and hence what happens between the scans
%      for several minutes is irrelevant.
%  
% However, segments between scans are removed just in case.This does not
% result in any noticeable performance improvement. But it reduces size of
% matlab arrays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if removeSamplesBetweenScans==1,
    
    disp('Removing data in between scans ...')
    extraSamples=5000; 
    dataOffset=extraSamples;
    
    samplesEachScanArray=scanStartEndArray(2:2:end)-scanStartEndArray(1:2:end-1);
    numSamplesInScans=sum(samplesEachScanArray);
    numSamplesInScans=numSamplesInScans+extraSamples;
    
    eyetrack_data=zeros(numSamplesInScans,12);

    scanStartEndArrayCopy=scanStartEndArray;
    sampleIndex=ones(1, extraSamples)*scanStartEndArrayCopy(1);
    %copy 1st sample of 1st scan to first 5000 samples,to stabilize filters
    
    eyetrack_data(1:extraSamples,:)=eyetrack_data_orig(sampleIndex,:);

    for count=1:length(scanStartEndArrayCopy)/2,
        %copy only samples inside scans
        sampleIndex=scanStartEndArrayCopy(count*2-1):scanStartEndArrayCopy(count*2);
        eyetrack_data(extraSamples+(1:length(sampleIndex)),:)=eyetrack_data_orig(sampleIndex,:);
        
        %re-adjust scanStartEndArray
        scanStartEndArray(count*2-1)=extraSamples+1;
        scanStartEndArray(count*2)=extraSamples+length(sampleIndex);
        extraSamples=extraSamples+length(sampleIndex);
    end
    eyetrack_data_orig2=eyetrack_data;
    eyetrack_data_length=length(eyetrack_data(:,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 4 Blink removal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   1. Runs a low pass filter with 14 Hz  cutoff frequency: FIR filter with order 22.
%       to cut off high frequency noise in pupil height and width.
if lpfBlinks==1,
%Low Pass Filter at 14 Hz [12-16 Hz] using FIR

if SignalProcessingToolbox==1,
%Need Signal Processing toolbox

    [n,fo,mo,w] = firpmord([12 16], [1 0], [.1 .005], fs);
    if mod(n,2)==0, n=n+1;end
    delay=(n-1)/2;
    B4=firls(n,fo,mo,w); %firls is better

else
    delay=10;
    B4=[   0.009633	0.010737	-0.010706	-0.024894 ...	
    0.006676	0.045477	0.007453	-0.079602	...
    -0.049732	0.175804	0.415678	0.415678	...
    0.175804	-0.049732	-0.079602	0.007453	...
    0.045477	0.006676	-0.024894	-0.010706	...
    0.010737	0.009633];

end

    A4=1;   

    sig=eyetrack_data(:,pupilHeightIndex);
    temp = [sig; zeros(delay,1)];
    temp=filter(B4,A4,temp);
    eyetrack_data(:,pupilHeightIndex)=temp(delay+1:delay+length(sig));

    sig=eyetrack_data(:,pupilWidthIndex);
    temp = [sig; zeros(delay,1)];
    temp=filter(B4,A4,temp);
    eyetrack_data(:,pupilWidthIndex)=temp(delay+1:delay+length(sig));

    %filtered values which are negative adjusted. else apsect ratio will be
    %negative
    negHeight=find(eyetrack_data(:,pupilHeightIndex)<=0);
    if isempty(negHeight)==0,
        eyetrack_data(negHeight,pupilHeightIndex)=0.001;
    end
    negWidth=find(eyetrack_data(:,pupilWidthIndex)<=0);
    if isempty(negWidth)==0,
        eyetrack_data(negWidth,pupilWidthIndex)=0.001;
    end    

end


%Check if blinking removal needs to be done
% remove blinks only if percentage of blinks is below a certain threshold
%   2. Removes blink portions based on pupil height criteria 
%     provided percentage of blink segments < a certain threshold.
%       in some data, pupil height and aspect ratio are highly unreliable.


aspectRatio=eyetrack_data(:,pupilWidthIndex)./eyetrack_data(:,pupilHeightIndex);

for scanNum=1: length(scanStartEndArray)/2,
   pupilHeightArray=[ eyetrack_data(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2), pupilHeightIndex)];
   aspectRatioArray=[ aspectRatio(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))];

    meanPupilHeightArray=mean(pupilHeightArray);
    array6=find(pupilHeightArray < meanPupilHeightArray*pupilHeightThr); 
    blinkFraction(scanNum)=length(array6)/length(pupilHeightArray);
    
end %end scanNum

disp('Fraction of Data classified as blinks')
blinkFraction 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%CODE FOR EYEBLINK REMOVAL%%%%%%%%%%%%%%%%%%%%%%%%%%%


    meanPupilHeight=mean(eyetrack_data(:,pupilHeightIndex));
    meanPupilHeightThr=meanPupilHeight*pupilHeightThr;
    blinkArray0=find(eyetrack_data(:,pupilHeightIndex) < meanPupilHeightThr);
    %Declare a Blink if pupil height < meanPupilHeight * 0.75
    
    %[blinkArray]=find(aspectRatio>aspectRatioThr);
    %[blinkArray2]=find(aspectRatio<aspectRatioThr2);
    %removeArray=sort([ blinkArray; blinkArray2 blinkArray0]); 
    %Aspect ratio  not used because aspect ratio is jittery and unreliable in some
    %segments in CG040710.txt
    
    blinkStartEndArray=[];
    
  if isempty(blinkArray0)==0,
        
    removeArray=sort([ blinkArray0]);
    removePts=unique(removeArray);

    scanNum=1;

%If Blinks are separated by less than Blin memory, consider it as a single
%Blink and substitute blink portions with last good sample.
    if length(removePts)>2,
        removePts(end+1)=removePts(end)+blinkMemory*2;
        array1=removePts(2:end)-removePts(1:end-1);
        blinkEndArray=find(array1>blinkMemory); %Locate end of blinks
    end
    
    if length(removePts)>2 && length(blinkEndArray)>0 ,
        startSample=max(removePts(1)-nPtsRemove,2);
        endSample=min(removePts(blinkEndArray(1))+nPtsRemove,eyetrack_data_length);
        while  endSample > scanStartEndArray(scanNum*2), %startSample > scanStartEndArray(scanNum*2-1) ||
             scanNum=scanNum+1; 
        end   
        if blinkFraction(scanNum) >  blinkThr && BlinkOverride==0,
            disp('No blinking removal due to noisy pupil height/aspect ratio');
        else
            disp('Blinking removal being done');
            eyetrack_data(startSample:endSample,xDataIndex) = eyetrack_data(startSample-1,xDataIndex);
            eyetrack_data(startSample:endSample,yDataIndex) = eyetrack_data(startSample-1,yDataIndex);
        end
        blinkStartEndArray(1)= startSample;
        blinkStartEndArray(2)= endSample;
    
        for count=1:length(blinkEndArray)-1,
            startSample=max(2,removePts(blinkEndArray(count)+1)-nPtsRemove);
            endSample=min(removePts(blinkEndArray(count+1))+nPtsRemove,eyetrack_data_length);
            while  endSample > scanStartEndArray(scanNum*2), 
                scanNum=scanNum+1; 
            end
            
            if blinkFraction(scanNum) >  blinkThr && BlinkOverride==0,
                %Don't suppress blinks
            else
                eyetrack_data(startSample:endSample,xDataIndex) = eyetrack_data(startSample-1,xDataIndex);
                eyetrack_data(startSample:endSample,yDataIndex) = eyetrack_data(startSample-1,yDataIndex);
            end       
            blinkStartEndArray(count*2+1)=max(2,removePts(blinkEndArray(count)+1)-nPtsRemove);
            blinkStartEndArray(count*2+2)=min(removePts(blinkEndArray(count+1))+nPtsRemove,eyetrack_data_length);
            
        end
      end    
  end  
%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blink Plots for whole data 

    dataIndex=yDataIndex;
    %dataIndex=xDataIndex;
    
    array11=eyetrack_data(:,dataIndex);

    array12=zeros(1,eyetrack_data_length);
  
    for count=1:length(blinkStartEndArray)/2,
        array12(blinkStartEndArray(count*2-1):blinkStartEndArray(count*2))=4; %Mark Blink Portions
    end

    array13=zeros(1,eyetrack_data_length);
    array13(scanStartEndArray(1:2:end))=5; %TTL pulse
    
    array14=find(abs(aspectRatio)>5);
    aspectRatioCopy=aspectRatio;
    if isempty(array14)==0,
        aspectRatioCopy(array14)=5; %clip at 5
    end
    
    N1=scanStartEndArray(1); 
    
    h_fig10=figure(10);

        hold off
        plot(eyetrack_data_orig2(N1:end, dataIndex),'k-')
        hold on
        plot(array11(N1:end),'g-')
        plot(array12(N1:end),'c-')
        plot(array13(N1:end),'r-')

        plot(eyetrack_data_orig2(N1:end, 8),'b-');
        plot(aspectRatioCopy(N1:end),'m-');

        
        title('Black: original; Green: Blink Removed; Cyan: Blink pulse; Blue: PupilHeight; Magenta:AspectRatio;');
        xlabel( 'samples in Scans; Red: TTL pulse start ');

        if saveFigure==1,
            outFileStr=strcat(fileStr0, '_figure10');
            saveas(h_fig10,outFileStr,'jpeg');
        end
        
if plotFigure_2==1, %Optional Blink plots every 1000 samples

    %dataIndex=xDataIndex;
    dataIndex=yDataIndex;
    
    array11=eyetrack_data(:,dataIndex);

    array12=zeros(1,eyetrack-data_length);
  
    for count=1:length(blinkStartEndArray)/2,
        array12(blinkStartEndArray(count*2-1):blinkStartEndArray(count*2))=4;
    end

    array14=find(abs(aspectRatio)>5);
    aspectRatioCopy=aspectRatio;
    if isempty(array14)==0,
        aspectRatioCopy(array14)=5; %clip at 5
    end
    
    N1=scanStartEndArray(1); 
    N0=N1; N3=1000;
    
    figure(40)
    hold off
    for N0=N1:1000:length(eyetrack_data),
        hold off
        plot(eyetrack_data_orig2(N0+(1:N3), dataIndex),'k-')
        hold on
        plot(array11(N0+(1:N3)),'g-')
        plot(array12(N0+(1:N3)),'c-')

        plot(eyetrack_data_orig2(N0+(1:N3), 8),'b-');
        plot(aspectRatioCopy(N0+(1:N3)),'m-');
        
        title('Black: original; Green: Blink Removed; Cyan: Blink pulse; Blue: PupilHeight; Magenta:AspectRatio;');
        xlabel('First 1000 samples in Scan 1');

        disp('Hit Enter to continue');

        pause
    end

end %end if plotFigure_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% x_hpf, y_hpf: High Pas Filtered Data
% x_lpf, y_lpf: Low Pas Filtered Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   3. Runs a 2-pole High pass Filter at cut-off frequency 0.025Hz to remove
%      slow drift due to head movements.
if includeHPF==1,
%High Pass Filter at 0.025Hz [0.01-0.1 Hz]

if SignalProcessingToolbox==1,
%Need Signal Processing toolbox

N=2;  %order of filter 2-pole
Wn=[0.025]*2/fs; 
ftype = 'low'; %lowpass filter
[B,A] = butter(N,Wn,ftype); %Butterworth IIR filter 
else
    B=[0.171030589091181   0.342061178182362   0.171030589091181]*1.0e-005;
    A =[1.000000000000000  -1.996297601769122   0.996304442992686];
end

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
if SignalProcessingToolbox==1,
%Need Signal Processing toolbox

[n,fo,mo,w] = firpmord([6 8], [1 0], [.1 .005], fs); %FIR filter order
if mod(n,2)==0, n=n+1;end
delay=(n-1)/2;
B2=firls(n,fo,mo,w); %firls is better
else

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
end
    A2=1;
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


%Plot high pass and low pass filtered data
    h_fig12=figure(12);
    hold off
    plot(eyetrack_data(dataOffset:end,4),'r-')
    hold on; plot(x_hpf(dataOffset:end),'b-');
    plot(x_lpf(dataOffset:end),'g-');
    title('Red: Original with blinks removed; Blue: High Pass Filtered; Green: Low pass Filtered');
     
    if saveFigure==1,
        outFileStr=strcat(fileStr0, '_figure12');
        saveas(h_fig12,outFileStr,'jpeg');
    end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Saccade detection
%
% saccadeStartEndArray: sample index of start and end of saccades.
%       saccadeStartEndArray[1]: start of saccade 1
%       saccadeStartEndArray[2]: end of saccade 1
%       saccadeStartEndArray[3]: start of saccade 2
%       saccadeStartEndArray[4]: end of saccade 2
%
%   saccadeCount: saccade count in all scans
%   saccadeStart(scanNum,:) : Start of saccades in current scan
%   saccadeEnd(scanNum,:) : End of saccades in current scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Short Term  Low pass filter

alpha=1/numSamplesShortTerm; 

B3=alpha;
A3=[1 -(1-alpha)];
x_lpf_short_term=filter(B3,A3,x_lpf); %1 pole IIR filter short term filtered output
y_lpf_short_term=filter(B3,A3,y_lpf);


cols=max(scanStartEndArray(2:2:end)-scanStartEndArray(1:2:end-1));
rows=length(scanStartEndArray)/2;
xData=zeros(rows,cols);
yData=zeros(rows,cols);
fixWinArray=zeros(length(scanStartEndArray)/2,100);
saccadeCount=1; %Total saccade Count for all scans

for scanNum=1:length(scanStartEndArray)/2,
    array_x=abs(x_lpf_short_term(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(x_lpf));
    array_y=abs(y_lpf_short_term(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(y_lpf));

    array2_x=find(array_x>shortTermSaccadeThr);
    %If short term average exceeds long term average by 2 degrees, then consider it as a valid saccade
    array2_y=find(array_y>shortTermSaccadeThr);
    saccadeIndexCurrentScan=unique(sort([array2_x array2_y])); %Location of saccades


    saccadeCountCurrentScan=1;
    
    if length(saccadeIndexCurrentScan)>1,
        array3=saccadeIndexCurrentScan(2:end)-saccadeIndexCurrentScan(1:end-1);
        saccadeEndIndex=find(array3>saccadeMemory);
        %If saccades are separated by less than 10 samples, consider it as a single saccade
        saccadeStart(scanNum,1)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(1);
        saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(1);
        saccadeCount=saccadeCount+1;
        
        for count=1:length(saccadeEndIndex),
            saccadeEnd(scanNum,count)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(saccadeEndIndex(count));
            saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(saccadeEndIndex(count));
            saccadeCount=saccadeCount+1;
            if saccadeEndIndex(count)+1 <= length(saccadeIndexCurrentScan),
                saccadeStart(scanNum,count+1)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(saccadeEndIndex(count)+1);
                saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(saccadeEndIndex(count)+1);
                saccadeCount=saccadeCount+1;
            end
            saccadeCountCurrentScan=saccadeCountCurrentScan+1;
        end

        saccadeEnd(scanNum,saccadeCountCurrentScan)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(end);
        saccadeStartEndArray(saccadeCount)=scanStartEndArray(scanNum*2-1)+saccadeIndexCurrentScan(end);
        saccadeCount=saccadeCount+1;
        saccadeCountCurrentScan=saccadeCountCurrentScan+1;
    end %end if saccadeIndexCurrentScan
    
    saccadeCountCurrentScan=saccadeCountCurrentScan-1;
      
    %If data is highly noisy like mas02142011.txt, use this option to
    %suppress all amplitudes above blinkLimitThr outside valid saccades.
    %This is hard Limiting
    if suppressBlinksOutsideSaccades==1,

        array_x=abs(x_lpf(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(x_lpf));
        array3_x=find(array_x>blinkLimitThr); %Location of amplitudes above threshold
            
        if saccadeCountCurrentScan>0,
            count3=1;count4=1;
            for count2=1:length(array3_x),
                if count3<saccadeCountCurrentScan &&  array3_x(count2)+scanStartEndArray(scanNum*2-1)>saccadeEnd(scanNum,count3),
                    count3=count3+1; 
                end
                if  array3_x(count2)+scanStartEndArray(scanNum*2-1)>=saccadeStart(scanNum,count3) && array3_x(count2)+scanStartEndArray(scanNum*2-1)<=saccadeEnd(scanNum,count3),
                    %Middle of a saccade; don't suppress
                    
                else
                    array4_x(count4)=array3_x(count2); count4=count4+1;
                end
            end

        else
                array4_x=array3_x; count4=length(array4_x);
                %Location of amplitudes above threshold, other than
                %saccades
        end

        nPtsRemove_2=1;
            
        if length(array4_x)>0,
            removePts=array4_x;
            
            if length(removePts)>2,
                removePts(end+1)=removePts(end)+20;
                array5_x=removePts(2:end)-removePts(1:end-1);
                arrayBlinkLimitThr=find(array5_x> blinkMemory);
            end
           
    
            if length(removePts)>2 && length(arrayBlinkLimitThr)>0 ,
                startSample=max(removePts(1)-nPtsRemove_2,2);
                endSample=removePts(arrayBlinkLimitThr(1))+nPtsRemove_2;
    
                x_lpf(scanStartEndArray(scanNum*2-1)-1+(startSample:endSample)) = x_lpf(scanStartEndArray(scanNum*2-1)-1+startSample-1);
    
                for count=1:length(arrayBlinkLimitThr)-1,
                    startSample=max(2,removePts(arrayBlinkLimitThr(count)+1)-nPtsRemove_2);
                    endSample=min(removePts(arrayBlinkLimitThr(count+1))+nPtsRemove_2,length(x_lpf));
                    x_lpf(scanStartEndArray(scanNum*2-1)-1+(startSample:endSample)) = x_lpf(scanStartEndArray(scanNum*2-1)-1+startSample-1);
                end
            end
        end


        array_y=abs(y_lpf(scanStartEndArray(scanNum*2-1):scanStartEndArray(scanNum*2))-mean(y_lpf));            
        array3_y=find(array_y>blinkLimitThr);
 
        if saccadeCountCurrentScan>0,
            count3=1;count4=1;
            for count2=1:length(array3_y),
                if count3<saccadeCountCurrentScan && array3_y(count2)+scanStartEndArray(scanNum*2-1)>saccadeEnd(scanNum,count3),
                    count3=count3+1; 
                end
                if array3_y(count2)+scanStartEndArray(scanNum*2-1)>=saccadeStart(scanNum,count3) && array3_y(count2)+scanStartEndArray(scanNum*2-1)<=saccadeEnd(scanNum,count3),
                    %Middle of a saccade; don't suppress
                    
                else
                    array4_y(count4)=array3_y(count2); count4=count4+1; %outside saccades
                end
            end       
        else
                array4_y=array3_y;count4=length(array4_y);
                %Location of amplitudes above threshold outside saccades
        end        

        if length(array4_y)>0,
            removePts=array4_y;
            if length(removePts)>2,
                removePts(end+1)=removePts(end)+20;
                array5_y=removePts(2:end)-removePts(1:end-1);
                arrayBlinkLimitThrY=find(array5_y> blinkMemory);
            end
           
    %Substitute data amplitudes above threshold, outside saccades, with last valid sample 
            if length(removePts)>2 && length(arrayBlinkLimitThrY)>0 ,
                startSample=max(removePts(1)-nPtsRemove_2,2);
                endSample=removePts(arrayBlinkLimitThrY(1))+nPtsRemove_2;
                y_lpf(scanStartEndArray(scanNum*2-1)-1+(startSample:endSample)) = y_lpf(scanStartEndArray(scanNum*2-1)-1+startSample-1);
    
                for count=1:length(arrayBlinkLimitThrY)-1,
                    startSample=max(2,removePts(arrayBlinkLimitThrY(count)+1)-nPtsRemove_2);
                    endSample=min(removePts(arrayBlinkLimitThrY(count+1))+nPtsRemove_2,length(y_lpf));
                    y_lpf(scanStartEndArray(scanNum*2-1)-1+(startSample:endSample)) = y_lpf(scanStartEndArray(scanNum*2-1)-1+startSample-1); 
                end   
            end
        end     
        
    end %if suppressBlinksOutsideSaccades
            

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
    if length(array12)>0,
        if length(array12)<=100,
            fixWinArray(scanNum,1:length(array12))=array12;
        else
            fixWinArray(scanNum,1:100)=array12(1:100);
        end
    end
    fracPtsOutsideWindow(scanNum)=length(array12)*100/xyNum(scanNum);
    
    if plotFigure==1,
    %Plot figures for each scan    
    figure(41);
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


%Plot Figure for whole data
   
    array13=zeros(1,eyetrack_data_length);
    array13(scanStartEndArray(1:2:end))=2; %TTL pulse
    
    h_fig21=figure(21);
    hold off
    plot(x_lpf(dataOffset:end),y_lpf(dataOffset:end),'x');
    title('Scatterplot: X, Y')
    h_fig22=figure(22);
    hold off
    subplot(4,1,1);
    hist(x_lpf(dataOffset:end),500);
    title('X Histogram')
    subplot(4,1,2);
    hist(y_lpf(dataOffset:end),500);
    title('Y Histogram')
    subplot(4,1,3);
    hold off
    plot(x_lpf(dataOffset:end));
    hold on; plot(array13(dataOffset:end),'r-')
    title('X')
    subplot(4,1,4);
    plot(y_lpf(dataOffset:end));
    title('Y')

    if saveFigure==1,
        outFileStr=strcat(fileStr0, '_figure21');
        saveas(h_fig21,outFileStr,'jpeg');
        outFileStr=strcat(fileStr0, '_figure22');
        saveas(h_fig22,outFileStr,'jpeg');
    end
    
    disp('Mean of X data')
    mean(meanX)
    disp('Mean of Y data')
    mean(meanY)
    disp('Standard Deviation of X data')
    mean(stdX)
    disp('Standard Deviation of Y data')
    mean(stdY)
    disp('Fraction of data outside fixation window')
    mean(fracPtsOutsideWindow)



    

%%%%%%%%%%%%%%%%%%%%%%%%%%End of Section 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section  7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print to file

if printFile==1,
    disp('Printing to File. This will take a few minutes');
    outFileStr0=strcat(fileStr0,'_scanInfo_');
    textStr1='Time';textStr2='X';textStr3='Y';textStr4='Blink';textStr5='Saccade';textStr6='Outside fixation window';
    N10=length(blinkStartEndArray);
    blinkArrayIndex=1;saccadeArrayIndex=1;
    [N11, N12]=size(fixWinArray);

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
        fprintf(outFile,'%f\t',eyetrack_data_orig2(dataIndex,timeIndex));
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


