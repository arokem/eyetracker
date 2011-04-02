
README FILE FOR PROCESSING EYE TRACKER DATA 
===========================================


1. Run "parse_eyetracker_data.m" program.

 This program parses eyetracker data file, removes first 26 text lines
 and replaces sync pulse column with 0,1,2 and outputs a file which is
 used by program "eyetrack_process.m".

2. Run "clear all; eyetrack_process(TR_array, numTR_array);"

   TR_array: Array containing TR(Time of repetition) for each scan in seconds.
   numTR_array: Array containing number of TR pulses for each scan.

   NOTE: wrong values specified in TR_array will give wrong results.

   Input filename is  specified using "uigetFile" command 
   when it pops a window asking for input file name.

 Example:  For file: cg021810.txt

   clear all; 
   TR_array=ones(1,13)*1;
   numTR_array = ones(1,13)*378;
   eyetrack_process(TR_array, numTR_array);


 User Options:
   A) For really noisy data, if a very high percentage of samples are misclassified 
     as Blinks, or if too many Blinks go unsuppressed, use this option below.

     suppressBlinksOutsideSaccades=1; %set this to 1 to hard-limit X,Y amplitudes above blinkLimitThr, 
			 %outside the region of valid saccades

   B) BlinkOverride=1; % set this to 1 to perform blink removal always irrepsective of noise/jitter in pupilheight

   C) blinkMemory=10; %If adjacent blinks are separated by < blinkMemory, consider it as a single blink.
        %Increase to 60  to capture several adjacent blinks. 

   D)minScanIntervalSamples=params.sampleRate*10; %Minimum number of samples in between 2 scans

   E) maxNumMissingTTLpulses=4; 
    % Maximum number of missing TTL pulses. Used for automatic derivation
    %For some data, more than 4 pulses may be missing, in which case, scan
    %start and end will be identified incorrectly. Blink plots have scan
    %start plotted and can indicate if this is the case. Then increase
    %maxNumMissingTTLpulses to a higher number.[Example:arokem_20110110_mod.txt] 
    %The tradeoff is: extra samples equal to maxNumMissingTTLpulses after the end of each scan will be
    %considered as valid samples inside scan , which may increase fraction of samples outside fixation window.

   % NOTE that user provided number of TTL pulses in each scan cannot be
   % relied upon because of missing TTL pulses.

3. "parse_eyetracker_outdata.m" provides an example of how to 
scan output files printed by "eyetrack_process.m"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Issues:
========
 1.  Array Indexing and array size checks: 
   most issues have been fixed, but there could be residual bugs
  which get exposed by new data files. If they do exist, feel free to fix them.

 2. Noisy Data:

    If the eyetracker is not able to locate the pupil correctly,
    it may print wrong values in pupil height, pupil width, X and Y data.
    Which will result in misclassified blinks and poorer peformance.

    Example plots in plots/ directory:
          anna_blink_data_scan_1 :  good data  cg021810.txt
	  ariel_blink_data_scan_1 :  noisy data  mas02142011

    It may be a good idea to analyze the difference between
    Anna's eyetracker setup and present setup to identify the
   source of this eyetracker noise.

   A) It seems specific steps need to be done to help the eyetracker 
      locate the pupil correctly.Viewpoint eyetracker user guide  chapter 5
      talks about "Locating Pupil and Glint". 

  B) Thomas: former postdoc who helped Anna setup eyetracker might be able to help.
       tlauritzen@gmail.com


NOTE:
=======
these matlab programs can also be run by open-source GNU Octave.
 Download Octave from
http://www.gnu.org/software/octave/download.html

==========================================

