clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;
res = 1;
Rmax = 200;
fc = 77e9;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
v = -20;
r = 150;
 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells


Bsweep = c/(2*res);
Tchirp = 5.5 * 2 * Rmax/c;
a = Bsweep/Tchirp;

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); %transmitted signal
Rx = zeros(1,length(t)); %received signal
Mix = Tx.*Rx; % beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = r + v*t(i);
    td(i) = 2 * r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + a*(t(i)^2)/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + a*((t(i)-td(i))^2)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);  
end

%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
reshape(Mix, [Nr, Nd]);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) 
signal = fft(Mix, Nr, 2);

 % *%TODO* :
% Take the absolute value of FFT output
signal = abs(signal);


 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
N = length(signal);
signal = signal(1:N/2 + 1);

%calculate range axis
Fs = 1/mean(diff(t));
Bfreq = (0:1:length(signal)-1)*Fs/Nr;
range = (c * Tchirp * Bfreq)/(2*Bsweep);

%normalize signal
signal = signal/length(signal);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
plot(range, signal)
 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 3;
Tc = 3;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 2;
Gc = 2;

% *%TODO* :
% offset the threshold by SNR value in dB
threshold = 5;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
[rows, columns] = size(RDM);
RDM_truncated = zeros(rows, columns);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
for i = 1:rows -2*Tr-2*Gr
    for j= 1:columns-2*Tc-2*Gc
       
        %cell under test indices
        CUT = [i + Tr + Gr, j + Tc + Gc];
        
        outer_row_begin = i;
        outer_row_end = i + 2*Tr + 2*Gr -1;
        
        outer_col_begin = j;
        outer_col_end = j + 2*Tc + 2*Gc -1;
        
        inner_row_begin = i + Tr;
        inner_row_end =  i + Tr + 2*Gr -1;
        inner_col_begin = j + Tc;
        inner_col_end = j + Tc + 2*Gc -1;
        
        count = 0;
        noise_level_sum = 0;
        for k = outer_row_begin:outer_row_end
            for l = outer_col_begin:outer_col_end
                
                if (k >=  inner_row_begin && k <= inner_row_end && l >= inner_col_begin && l <= inner_col_end)
                    count = count + 1;
                    noise_level_sum = db2pow(RDM(k,l)) + noise_level_sum;
                end
            end
        end
        noise_level_avg = noise_level_sum/count;
        noise_level_avg_log = pow2db(noise_level_avg);
        
        cfar_thresh = threshold + noise_level_avg_log;
        
        %supress output based on CFAR threshold
        if RDM(CUT(1), CUT(2)) >= cfar_thresh
            RDM_truncated(CUT(1), CUT(2)) = 1;
        end

    end
end

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM_truncated);
colorbar;


 
 