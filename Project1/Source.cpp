// DP Elective Topic | Analogue Electronics/DSP | Task 2
// EEG Signal Processing

// Assume peripherals are correctly set up

// Load in words

#include "AudioFile.h"
#include "fft.h"
#include "fft.hpp"
#include "fft_impl.hpp"
#include "fft_settings.h"
#include "copy_array.hpp"
#include "error_handling.hpp"
#include "check_fft.hpp"
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

    double* Frame = new double[N];
    bool b = FFT(); 


};



