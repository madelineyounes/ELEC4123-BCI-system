// DP Elective Topic | Analogue Electronics/DSP | Task 2
// EEG Signal Processing

// Assume peripherals are correctly set up

// Load in words

#include "AudioFile.h"
// https://github.com/adamstark/AudioFile
//#include "csv.h"
//// https://github.com/ben-strasser/fast-cpp-csv-parser
#include "dj_fft.h"
// https://github.com/jdupuy/dj_fft

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <ctime>
#include <math.h>
#include <complex>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <string>

#define HOP_IN_SECS 0.125

using namespace std;

int make_decision(vector<complex<float>> frame, int N, double* noise_mean, double* noise_std);

int main()
{
    // Read in word files
    // AudioFile setup
    AudioFile<double> deltaFile;
    AudioFile<double> thetaFile;
    AudioFile<double> alphaFile;
    AudioFile<double> betaFile;

    // Load in .wav files as objects
    deltaFile.load("../data/Word wav files/deltaWord.wav");
    thetaFile.load("../data/Word wav files/thetaWord.wav");
    alphaFile.load("../data/Word wav files/alphaWord.wav");
    betaFile.load("../data/Word wav files/betaWord.wav");

    // Obtain length of audio file in seconds
    //double lengthInSeconds_alpha = alphaFile.getLengthInSeconds();
    //double lengthInSeconds_beta = betaFile.getLengthInSeconds();
    //double lengthInSeconds_delta = deltaFile.getLengthInSeconds();
    //double lengthInSeconds_theta = thetaFile.getLengthInSeconds();

    // Assume mono channel
    int channel = 0;

    // Assume sampling frequency constant across .wav files
    int fs_word = alphaFile.getSampleRate();

    // Initialize pointers to .wav files
    vector<double> deltaWord = deltaFile.samples[channel];
    vector<double> thetaWord = thetaFile.samples[channel];
    vector<double> alphaWord = alphaFile.samples[channel];
    vector<double> betaWord = betaFile.samples[channel];

    // Pad audio samples to the next integer number of hops
    double samples_per_hop = HOP_IN_SECS * fs_word;
    deltaWord.resize(ceil(deltaFile.getNumSamplesPerChannel() / samples_per_hop) * samples_per_hop, 0.0);
    thetaWord.resize(ceil(thetaFile.getNumSamplesPerChannel() / samples_per_hop) * samples_per_hop, 0.0);
    alphaWord.resize(ceil(alphaFile.getNumSamplesPerChannel() / samples_per_hop) * samples_per_hop, 0.0);
    betaWord.resize(ceil(betaFile.getNumSamplesPerChannel() / samples_per_hop) * samples_per_hop, 0.0);

    //// Read in EEG signal as .wav file
    //AudioFile<double> EEGFile;
    //EEGFile.load("../data/EEGdata/Synthetic EEG Fs512Hz/Noise/delta2.wav");
    //vector<double> EEG = EEGFile.samples[channel];
    //int L = EEGFile.getNumSamplesPerChannel();  // Number of samples
    //int fs = EEGFile.getSampleRate();           // Sampling frequency

    // Read in EEG signal as .csv file
    vector<double> EEG = vector<double>();  // Empty vector
    string line;
    string eeg_name = "beta2";
    ifstream file_in;
    file_in.open("../data/EEGdata/Synthetic EEG Fs512Hz/Noise/" + eeg_name + ".csv");
    while (file_in >> line) {
        EEG.push_back(stod(line));
    }
    int L = EEG.size();
    int fs = 512;       // Synthetic EEG signals are sampled at 512 Hz

    int N = 1 * fs;                                         // Frame size in samples
    int hop = HOP_IN_SECS * fs;                             // Hop size
    int num_frames = ceil((L - N) / (float)hop);            // Number of frames
    vector<int> decision = vector<int>(num_frames + 1, 0);  // Allocate memory for decision array

    // Obtain current time in ms
    using std::cout; using std::endl;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::system_clock;

    auto time_now = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    // Check time_now
    // cout << "milliseconds since epoch: " << time_now << endl;
    int i_finish = 0;
    vector<double>audio_samples = vector<double>(); // Empty vector
    
    // Noise Parameters
    double noise_mean = 0;
    double noise_std = 0;

    for (int i = 1; i <= num_frames; i++) {
        auto tic = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

        // Choose frame
        int j = (i - 1) * hop;
        vector<complex<float>> frame = vector<complex<float>>(EEG.begin() + j, EEG.begin() + j + N);
        decision[i] = make_decision(frame, N, &noise_mean, &noise_std);
        frame.clear(); frame.shrink_to_fit();

        // 'Play' word
        if (i >= i_finish) {
            if (decision[i] == decision[i-1]) {
                // Decisions agree
                if (decision[i] == 1) {
                    audio_samples.insert(audio_samples.end(), deltaWord.begin(), deltaWord.end());
                    i_finish = i + (deltaWord.size() / samples_per_hop);

                } else if (decision[i] == 2) {
                    audio_samples.insert(audio_samples.end(), thetaWord.begin(), thetaWord.end());
                    i_finish = i + (thetaWord.size() / samples_per_hop);

                } else if (decision[i] == 3) {
                    audio_samples.insert(audio_samples.end(), alphaWord.begin(), alphaWord.end());
                    i_finish = i + (alphaWord.size() / samples_per_hop);

                } else if (decision[i] == 4) {
                    audio_samples.insert(audio_samples.end(), betaWord.begin(), betaWord.end());
                    i_finish = i + (betaWord.size() / samples_per_hop);

                } else if (decision[i] == 0) {
                    audio_samples.resize(audio_samples.size() + samples_per_hop, 0.0);
                }
            } else {
                // Decisions disagree
                audio_samples.resize(audio_samples.size() + samples_per_hop, 0.0);
            }
        }
    }

    // Write decision to csv file
    // i = 1 because decision[0] is always 0
    ofstream file_out_decisions;
    file_out_decisions.open("../data/Test Results/" + eeg_name + "_decisions.csv");
    for (int i = 1; i < num_frames + 1; i++) {
        if (i == num_frames) {
            // Last frame
            file_out_decisions << to_string(decision[i]);
        } else {
            file_out_decisions << to_string(decision[i]) + ",";
        }
    }
    file_out_decisions.close();

    // Write audio samples to csv file
    ofstream file_out_audio_samples;
    file_out_audio_samples.open("../data/Test Results/" + eeg_name + "_audio_samples.csv");
    for (int i = 0; i < audio_samples.size(); i++) {
        if (i == audio_samples.size() - 1) {
            file_out_audio_samples << to_string(audio_samples[i]);
        } else {
            file_out_audio_samples << to_string(audio_samples[i]) + ",";
        }
    }
    file_out_audio_samples.close();

    // Free memory
    alphaWord.clear(); alphaWord.shrink_to_fit();
    betaWord.clear(); betaWord.shrink_to_fit();
    deltaWord.clear(); deltaWord.shrink_to_fit();
    thetaWord.clear(); thetaWord.shrink_to_fit();
    EEG.clear(); EEG.shrink_to_fit();
    decision.clear(); decision.shrink_to_fit();
    audio_samples.clear(); audio_samples.shrink_to_fit();
}



 // Note: 'frame' is deleted by the function caller
int make_decision(vector<complex<float>> frame, int N, double* noise_mean, double* noise_std) {
    
    // Noise Learning Rate
    #define MEAN_ETA 0.2
    #define STD_ETA 0.2
    #define NOISE_LENGTH 215
    #define BAND_NARROW 5
    #define BAND_WIDE 29

    // Windows 
    vector<double> WIN_WIDE{0.8160, 0.8183, 0.8251, 0.8361, 0.8506, 0.8681, 0.8875,
                          0.9080, 0.9285, 0.9479, 0.9654, 0.9799, 0.9909, 0.9977,
                          1.0000, 0.9977, 0.9909, 0.9799, 0.9654, 0.9479, 0.9285,
                          0.9080, 0.8875, 0.8681, 0.8506, 0.8361, 0.8251, 0.8183,
                          0.8160};
    vector<double> WIN_NARROW{0.8160, 0.9080, 1.0000, 0.9080, 0.8160};

    // Perform FFT
    vector<complex<float>> Frame_complex = dj::fft1d(frame, dj::fft_dir::DIR_FWD);

    // Calculate power spectrum
    vector<double> Frame = vector<double>(N);
    for (int i = 0; i < N; i++) {
        Frame[i] = pow(Frame_complex[i].real(), 2.0) + pow(Frame_complex[i].imag(), 2.0);
    }

    // Single sided power spectrum
    vector<double> Frame_delta = vector<double>(5);
    vector<double> Frame_theta = vector<double>(5);
    vector<double> Frame_alpha = vector<double>(5);
    vector<double> Frame_beta = vector<double>(29);
    
    // The last element is not included, [first1, last1)
    transform(Frame.begin(), Frame.begin() + 5,
                WIN_NARROW.begin(),
                Frame_delta.begin(),
                multiplies<double>());

    transform(Frame.begin() + 4, Frame.begin() + 9,
                WIN_NARROW.begin(),
                Frame_theta.begin(),
                multiplies<double>());

    transform(Frame.begin() + 8, Frame.begin() + 13,
                WIN_NARROW.begin(),
                Frame_alpha.begin(),
                multiplies<double>());

    transform(Frame.begin() + 12, Frame.begin() + 41,
                WIN_WIDE.begin(),
                Frame_beta.begin(),
                multiplies<double>());

    // Estimate noise
    double noise_mean_curr = accumulate(Frame.begin() + 41, Frame.begin() + 256, 0.0) / NOISE_LENGTH;
    vector<double> noise_squared = vector<double>(NOISE_LENGTH);
    transform(Frame.begin() + 41, Frame.begin() + 256,
        Frame.begin() + 41,
        noise_squared.begin(),
        multiplies<double>()
    );

    double noise_squared_mean_curr = accumulate(noise_squared.begin(), noise_squared.end(), 0.0) / NOISE_LENGTH;
    // std = sqrt(E[X^2] - (E[X])^2)
    double noise_std_curr = sqrt(noise_squared_mean_curr - pow(noise_mean_curr, 2.0));

    if (*noise_mean == 0 && *noise_std == 0) {
        *noise_mean = noise_mean_curr;
        *noise_std = noise_std_curr;
    }
    else {
        *noise_mean = *noise_mean * MEAN_ETA + noise_mean_curr * (1 - MEAN_ETA);
        *noise_std = *noise_std * STD_ETA + noise_std_curr * (1 - STD_ETA);
    }

    double noise_threshold = *noise_mean + 2 * (*noise_std);
    double delta_avg = accumulate(Frame_delta.begin(), Frame_delta.end(), 0.0) / BAND_NARROW - noise_threshold;
    double theta_avg = accumulate(Frame_theta.begin(), Frame_theta.end(), 0.0) / BAND_NARROW - noise_threshold;
    double alpha_avg = accumulate(Frame_alpha.begin(), Frame_alpha.end(), 0.0) / BAND_NARROW - noise_threshold;
    double beta_avg = accumulate(Frame_beta.begin(), Frame_beta.end(), 0.0) / BAND_WIDE - noise_threshold;

    // Free memory
    WIN_WIDE.clear(); WIN_WIDE.shrink_to_fit();
    WIN_NARROW.clear(); WIN_NARROW.shrink_to_fit();
    Frame_complex.clear(); Frame_complex.shrink_to_fit();
    Frame.clear(); Frame.shrink_to_fit();
    Frame_delta.clear(); Frame_delta.shrink_to_fit();
    Frame_delta.clear(); Frame_delta.shrink_to_fit();
    Frame_alpha.clear(); Frame_alpha.shrink_to_fit();
    Frame_beta.clear(); Frame_beta.shrink_to_fit();
    noise_squared.clear(); noise_squared.shrink_to_fit();

    // Make decision on this frame
    if (delta_avg <= 0 && theta_avg <= 0 && alpha_avg <= 0 && beta_avg <= 0) {
        return 0;
    } else {
        vector<double> avgs{delta_avg, theta_avg, alpha_avg, beta_avg};
        return distance(avgs.begin(), max_element(avgs.begin(), avgs.end())) + 1;
    }
};



