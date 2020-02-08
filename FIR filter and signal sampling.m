## begining - clearing al previous figures, and windows etc...
clear all;
close all;
clf;
pkg load signal;


## loading of an audiofile to project
## s- signal## fp - frequenzy of a file - for us it is 44100 Hz
[s, fp] = auload('bird.wav');


##  audio file is ~1,78s long
## 1,78 * 44100 ~78.5 kHz

t_length = length(s);

## vector of next samples from 0 to ~78500 
t = linspace(0, t_length/fp, length(s) );


## amplification factor
amp = 10.0;
## multiply figure by amp factor
s(1:t_length) = amp * s(1:t_length);


## number of bits of filter - iterative loop
BITS = 5;

##for loop describing good enough resolution of sampling
for n = 1:BITS;
  
  
  a =2^n - 1;
  
  ## rounding a signal quantized value using "a" parameter
  signal_quantized = round( s(1:t_length) * a) / a;
  
  
  ## noise signal
  noise = s(1:t_length) - signal_quantized;
  
  
  ## ploting signals: original, quantized and error of sampling
  ## -b: blue - original signal
  ## -g: green - quantized signal
  ## -r: red - noise signal
  
  plot(t, s(1:t_length), '-b; Original signal;',
        t, signal_quantized, '-g; Quantized signal;',
        t, noise, '-r; Error of sampling - noise;');    
      
  grid on;
  xlabel ('Time [s]');
  ylabel('Amplitude [-]');
  title({"Singular step of quantisation number:", int2str(n)});  
        
      
      
  ## (n, :) - all "n" rows will be read in their entirety
  ## std - standard deviation
  ## logarythmic scale multiplied by 20 - voltage scale??  
  signal_deviation(n,:) = 20 * log10( std( s(1:t_length) ) / std(noise) ); 
 


  ##int2str - converts integers to strings 
  print(["Quantized_signal ", int2str(n), ".png"], "-dpng", "-color");
  
  signal_rounded = round(s * a) / a;
  
  ## audio save function from folder
  ausave(['bird', int2str(n), '.wav'], signal_rounded, fp);
  
  
end
##
        
        
        
################################################################################
################################################################################
##SECOND FRAME WITH PLOTS  


## plot of signal in time multiplied by amp
## figure = new window with plot

 subplot(3,3,1);
 plot(t, s(1:t_length), '-b');
 
 title('Original signal');
 xlabel('Time [s]');
 ylabel('Amplitude [-]');
 grid on;
 hold on;

 
 
 
 subplot(3,3,3);
 plot(signal_deviation, 'r*');
 hold on;
 plot(signal_deviation ,'-b');
 hold on;
 grid on;
 
 title({'SNR signal to noise ratio', 'dependent on bits of filter'});
 xlabel({'Bits of filter -', 'resolution of quantization'});
 ylabel({'value of SNR', 'in logarytmic scale [dB]'});
 
 
 
 n_poly = 1:1:BITS; 

 
 ## method polyfit returns roots of equation with best fit
 ## symbol ' reverts horizontal vector to vertical matrix
 ## 1 means logiczal vector setting usage of roots to true
 polyfit(n_poly', signal_deviation, 1);
 
 
 
 ## overwriting (copying) final audio file to next one
 ## one with o number higher of 1
 file1 = ['bird' , int2str(BITS), '.wav'];
 file2 = ['bird' , int2str(BITS + 1), '.wav'];
 copyfile(file1, file2); 


 ## spectrum analisys of audio file after quantization and sampling
 [s_final, fp] = auload(file2);
 
 ## first half of Fourier spectrum
 N = fp / 2;

   
 ## function determines  to which top to the power of raise a number to be greater than N
 ## you could use logarythm in this place 
    Nf_help = 1;
    help = 2;
    while(help < N)
        help = help* 2;
        Nf_help++;
    endwhile
    ##  
 
 ## power two to a number - FFT is faster that way
 Nf = 2^Nf_help;
 
 ## Nf2~22kHz
 Nf2 = Nf/2 + 1;
 
 
 ## vector of equaly set probes
 f = linspace(0, fp/2, Nf2);
        
      
 ##subplot of only a signal after iterative sampling  
 subplot(3,3,4);
 plot(t, s_final(1:t_length));
 title({'Signal quantized','-without noise'});
 xlabel('Time [s]');
 ylabel('Amplitude [-]');
 hold on;
 grid on;
 
 
 ## fourier transfomate of final signal on frequanzy from 0 to 44.1 kHz
 s_fft = fft(s_final, Nf);
 
 
 ## absolute value of our final signal after fft
 s_fft_abs = abs(s_fft);
 
 
 subplot(3,3,5);
 plot(f, s_fft_abs(1:Nf2) );
 title({"Module of a signal after", "iterative sampling"}); 
 xlabel('Frequenzy [Hz]');
 ylabel({'Module of a signal after', 'iterative sampling [-]'});
 box off; grid on;
 
 
 
 ## phase of signal iterative sampling
 s_fft_angle = angle(s_fft);

 subplot(3,3,6);
 plot(f, s_fft_angle(1:Nf2));
 title({"Phase of signal after", "iterative sampling"});
 xlabel('Frequenzy [Hz]');
 ylabel('Phase angle of a signal [rad]');
 box off; grid on; axis tight;      
        
        
        
 ## noise signal plot
 subplot(3,3,7);
 plot(t,  noise(1:t_length) );
 title('Noise of sampling'); 
 xlabel('Time [s]'); 
 ylabel('Amplitude [-]');
 hold on; grid on;    
        
        
## FFT for noise signal
noise_fft = fft(noise, Nf);

noise_fft_abs = abs(noise_fft);   
        
        
## Module of signal spectrum
subplot(3,3,8);
plot(f, noise_fft_abs(1:Nf2));
title("Module of noise spectrum");
xlabel('Frequenzy [Hz]');
ylabel('Module of noise spectrum [-]');
grid on; hold on;      
        
        
## angle of noise signal
noise_ftt_angle = angle(noise_fft);

subplot(3,3,9);
plot(f, noise_ftt_angle(1:Nf2));
title("Phase angle of noise");
xlabel('Frequenzy [Hz]');
ylabel('Phase angle of signal [rad]');
grid on; hold on;
        
        
        
## saving picture with 1600x900 pixels
print(['Fig. 1 - Plots of signals original sampled and noise.png'], "-dpng", '-S1600,900');   
        
        
        
        
        
#############################################################################
#############################################################################
##THIRD PICTURE WITH PLOTS



## FILTR FIR
freqCut = [2200 2900];

## frequenzy normalized, fp ~44.1kHz
wc = freqCut / (fp/2);

## FIR filter coefficients
firCoeff = fir1(511, wc, 'bandpass');



figure;

subplot(2,3,4);
## filar-type plot
stem(firCoeff);
title("Filter coefficients");
hold on;


## signal audio after sampling - filtered with delay
s_filtered = filter(firCoeff, 1, s_final);


subplot(2,3,1);
plot(t, s_filtered, '-r', 
     t, s_final, '-b');

title("Signal audio - filter with delay: \n red - sig. original, \n blue - sig. filtered]");
xlabel('Time [s]'); 
ylabel('Signal [-]');
grid on; hold on;


## FFT of filtered signal
s_filtered_fft = fft(s_filtered, Nf);


## signal module after filter FIR
s_filtered_fft_abs = abs(s_filtered_fft);

subplot(2,3,2);
plot(f, s_filtered_fft_abs(1:Nf2));
title({"Module of signal spectrum after" , "using bandpass filer"});
xlabel('Frequenzy [Hz]');
ylabel('Module of signal spectrum after reduction [-]');
grid on; hold on;



## Signal phase angle after filter FIR
s_filtered_fft_angle = angle(s_filtered_fft);

subplot(2,3,3);
plot(f, s_filtered_fft_angle(1:Nf2));
title({"Signal phase angle", "after filter FIR"});
xlabel('Frequenzy [Hz]');
ylabel('Signal phase angle [rad]');
grid on; hold on;


ausave(['bird_filtered.wav'], s_filtered, fp); 




## frequenzy response of filter FIR
[H,f] = freqz(firCoeff, 1, 2^Nf_help, fp);

subplot(2,3,5);
plot(f, abs(H));
title({"Frequenzy", "response"});
xlabel('Frequenzy [Hz]')'
ylabel('Frequenzy response [-]');
grid on; hold on;

## frequenzy response phase angle of filter FIR
subplot(2,3,6);
plot(f, angle(H));
title({"Frequenzy response of", "phase angle"});
xlabel('Frequenzy [Hz]');
ylabel('Phase angle of signal [rad]');
grid on; hold on;


## save plots as a picture
print(['Fig. 2 - Plots of filter FIR.png'],"-dpng", '-S1600,900'); 





##############################################################################
##############################################################################
## LAST FRAME WITH PLOTS


figure;
freqz(firCoeff, 1, 2^Nf_help, fp);

print(['Fig. 3 - Plots of frequenzy response of filter FIR.png'],"-dpng", '-S1600,900'); 

     
        
