#include <stdio.h>
#include <vector>
#include <random>
#include <direct.h>

#include "pcg/pcg_basic.h"

#define DETERMINISTIC() false

static const size_t c_numSamples = 1000;
static const size_t c_numtests = 100;

// 1/goldenRatio
// or goldenRatio-1
static const float c_goldenRatioConjugate = 0.61803398875f;

static const float c_pi = 3.14159265359f;

inline pcg32_random_t GetRNG()
{
    pcg32_random_t rng;
#if DETERMINISTIC()
    pcg32_srandom_r(&rng, 0x1337FEED, 0);
#else
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<uint32_t> dist;
    pcg32_srandom_r(&rng, dist(generator), 0);
#endif
    return rng;
}

inline float RandomFloat01(pcg32_random_t& rng)
{
    return float(pcg32_random_r(&rng)) / 4294967295.0f;
}

std::vector<float> Sequence_WhiteNoise(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    for (size_t index = 0; index < numSamples; ++index)
        ret[index] = RandomFloat01(rng);
    return ret;
}

std::vector<float> Sequence_Regular(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    for (size_t index = 0; index < numSamples; ++index)
        ret[index] = float(index) / float(numSamples);
    return ret;
}

std::vector<float> Sequence_Stratified(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    for (size_t index = 0; index < numSamples; ++index)
        ret[index] = (float(index) + RandomFloat01(rng)) / float(numSamples);
    return ret;
}

std::vector<float> Sequence_GoldenRatio(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    ret[0] = 0.0f;
    for (size_t index = 1; index < numSamples; ++index)
        ret[index] = fmodf(ret[index - 1] + c_goldenRatioConjugate, 1.0f);
    return ret;
}

float Function_Sine(float x)
{
    return sin(c_pi * x);
}

float Function_Step(float x)
{
    return (x < 0.3f) ? 1.0f : 0.0f;
}

float Function_Triangle(float x)
{
    return x;
}

float Lerp(float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

// TODO: golden ratio sequence.
// TODO: spit out CSVs and graph with python or open office.

struct Sequence
{
    const char* name = nullptr;
    std::vector<float>(*fn)(size_t numSamples) = nullptr;
    std::vector<float> sequence;
    std::vector<float> results;
    std::vector<float> error;
};

struct Function
{
    const char* name = nullptr;
    float(*fn)(float x) = nullptr;
    float actualValue = 0.0f;
};

int main(int argc, char** argv)
{
    _mkdir("out");

    Sequence sequences[] =
    {
        {"Uniform", Sequence_WhiteNoise},
        {"Regular", Sequence_Regular},
        {"Stratified", Sequence_Stratified},
        {"GoldenRatio", Sequence_GoldenRatio},
    };

    Function functions[] =
    {
        {"Sine", Function_Sine, 2.0f / c_pi},
        {"Step", Function_Step, 0.3f},
        {"Triangle", Function_Triangle, 0.5f},
    };

    // for each function
    for (size_t functionIndex = 0; functionIndex < _countof(functions); ++functionIndex)
    {
        // TODO: do multiple tests and keep track of mean and stddev.
        //c_numtests

        // for each sequence
        for (size_t sequenceIndex = 0; sequenceIndex < _countof(sequences); ++sequenceIndex)
        {
            sequences[sequenceIndex].results.resize(c_numSamples, 0.0f);
            sequences[sequenceIndex].error.resize(c_numSamples, 0.0f);

            // for each number of samples to test
            for (size_t sampleCount = 0; sampleCount < c_numSamples; ++sampleCount)
            {
                // generate the samples
                sequences[sequenceIndex].sequence = sequences[sequenceIndex].fn(sampleCount + 1);

                // integrate
                float y = 0.0f;
                for (size_t index = 0; index < sampleCount; ++index)
                {
                    float x = sequences[sequenceIndex].sequence[index];
                    y = Lerp(y, functions[functionIndex].fn(x), 1.0f / float(index + 1));
                }

                // store the result
                sequences[sequenceIndex].results[sampleCount] = y;

                // and the error
                sequences[sequenceIndex].error[sampleCount] = std::abs(y - functions[functionIndex].actualValue);
            }
        }

        // Open the CSV file for writing
        char fileName[256];
        sprintf_s(fileName, "out/%s.csv", functions[functionIndex].name);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "wb");

        // write the CSV column headers
        fprintf(file, "\"Samples\"");
        for (size_t sequenceIndex = 0; sequenceIndex < _countof(sequences); ++sequenceIndex)
            fprintf(file, ",\"%s\"", sequences[sequenceIndex].name);
        fprintf(file, "\n");

        // write the results to the csv file
        for (int reportingIndex = 0; reportingIndex < 100; ++reportingIndex)
        {
            float percent = float(reportingIndex) / 99.0f;
            size_t sampleCount = std::min(size_t(percent * float(c_numSamples)), c_numSamples - 1);
            fprintf(file, "\"%i\"", (int)sampleCount+1);
            for (size_t sequenceIndex = 0; sequenceIndex < _countof(sequences); ++sequenceIndex)
                fprintf(file, ",\"%f\"", sequences[sequenceIndex].error[sampleCount]);
            fprintf(file, "\n");
        }

        // all done
        fclose(file);
    }

	return 0;
}

/*
TODO:
- convergence graphs of stratified, white noise, golden ratio 1d. on sine, step, triangle.
  - this to verify golden ratio is good.
  - white and stratified probably need a mean and variance (std dev) on the plot, so should happen many times.

- even if you do 100 million samples, maybe only report 100 results.

Blog:
* motivation: golden ratio is better sampling. having it in sorted order is better for cache, and better for ray marching, to composite things properly.
* Show the right way to do it, show that it works.
* show your failed attempt and the things you found along the way

*/