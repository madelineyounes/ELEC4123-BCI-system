// DP Elective Topic | Analogue Electronics/DSP | Task 2
// EEG Signal Processing

// Assume peripherals are correctly set up

// Load in words

#include "AudioFile.h"

#include "fft_settings.h"
#include "error_handling.hpp"
#include "fft_impl.hpp"
#include "copy_array.hpp"
#include "check_fft.hpp"
#include "fft.hpp"
#include "fft.h"

#include <vector>
#include <math.h>
#include <chrono>
#include <iostream>
#include <ctime>

int make_decision(double* frame, int N, double* noise_mean, double* noise_std);

int main()
{
    // Read in word files 
    // AudioFile setup
    AudioFile<double> alphaFile;
    AudioFile<double> betaFile;
    AudioFile<double> deltaFile;
    AudioFile<double> thetaFile;

    // Load in .wav files as objects
    alphaFile.load("../data/Word wav files/alphaWord.wav");
    betaFile.load("../data/Word wav files/betaWord.wav");
    deltaFile.load("../data/Word wav files/deltaWord.wav");
    thetaFile.load("../data/Word wav files/thetaWord.wav");

    // Obtain length of audio file in seconds
    double lengthInSeconds_alpha = alphaFile.getLengthInSeconds();
    double lengthInSeconds_beta = betaFile.getLengthInSeconds();
    double lengthInSeconds_delta = deltaFile.getLengthInSeconds();
    double lengthInSeconds_theta = thetaFile.getLengthInSeconds();

    // Assume mono channel
    int channel = 0;

    // Assume sampling frequency constant across .wav files
    int fs_word = alphaFile.getSampleRate();

    // Initialize pointers to .wav files
    std::vector<double> alphaWord = alphaFile.samples[channel];
    std::vector<double> betaWord = betaFile.samples[channel];
    std::vector<double> deltaWord = deltaFile.samples[channel];
    std::vector<double> thetaWord = thetaFile.samples[channel];

    // Read in EEG signal
    AudioFile<double> EEGFile;
    EEGFile.load("../data/EEGdata/Synthetic EEG Fs512Hz/Noise/alpha1.wav");
    std::vector<double> EEG = EEGFile.samples[channel];

    int L = EEGFile.getNumSamplesPerChannel();  // Number of samples
    int fs = EEGFile.getSampleRate();           // Sampling frequency
    int N = 1 * fs;                             // Frame size in samples
    int hop = 0.125 * fs;                     // Hop size
    int num_frames = ceil((L - N) / (float)hop);       // Number of frames
    int* decision = new int[num_frames + 1]();  // Allocate memory for decision array

    // Obtain current time in ms
    using std::cout; using std::endl;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::system_clock;

    auto time_now = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    // Check time_now
    // cout << "milliseconds since epoch: " << time_now << endl;
    
    // Noise Parameters
    double noise_mean = 0;
    double noise_std = 0;

    double* EEG_p = EEG.data(); // Obtain pointer to first element in EEG
    for (int i = 1; i <= num_frames; i++) {
        auto tic = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

        // Choose frame
        int j = (i - 1) * hop;
        decision[i] = make_decision(EEG_p + j, N, &noise_mean, &noise_std);
    }
    
    // [Last Step] Send samples of .wav file to speaker
    // Output the chosen .wav file as an array
}



 
int make_decision(double* frame, int N, double* noise_mean, double* noise_std) {
    
    // Noise Learning Rate
    #define MEAN_ETA 0.2
    #define STD_ETA 0.2

    // Windows 
    float WIN_WIDE[29] = {0.8160, 0.8183, 0.8251, 0.8361, 0.8506, 0.8681, 0.8875,
                          0.9080, 0.9285, 0.9479, 0.9654, 0.9799, 0.9909, 0.9977,
                          1.0000, 0.9977, 0.9909, 0.9799, 0.9654, 0.9479, 0.9285,
                          0.9080, 0.8875, 0.8681, 0.8506, 0.8361, 0.8251, 0.8183,
                          0.8160};
    float WIN_NARROW[5] = {0.8160, 0.9080, 1.0000, 0.9080, 0.8160};

    float FFT_TEST_SINE[32] = { 0.0000, 0.2027, 0.4054, 0.6081, 0.8107, 1.0134, 1.2161, 1.4188,
                                1.6215, 1.8242, 2.0268, 2.2295, 2.4322, 2.6349, 2.8376, 3.0403,
                                3.2429, 3.4456, 3.6483, 3.8510, 4.0537, 4.2564, 4.4590, 4.6617,
                                4.8644, 5.0671, 5.2698, 5.4725, 5.6751, 5.8778, 6.0805, 6.2832 };

    double* Frame = new double[N];

    // FFT notes
    // Input array must be a power of 2
    // For 1D FFT, use in the form b = FFT(A,B,n,error);
    // Initialize error, const char * error = NULL, the function returns an error descriptor
    // User needs to define two types called real_type and complex_type
    const char* error = NULL;
    auto b = FFT(FFT_TEST_SINE,Frame,32,error);


};



