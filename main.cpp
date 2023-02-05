#include <stdio.h>
#include <vector>
#include <random>
#include <direct.h>

#include "pcg/pcg_basic.h"

#define DETERMINISTIC() true

static const size_t c_numSamples = 100;
static const size_t c_numtests = 10000;

// 1/goldenRatio
// or goldenRatio-1
static const float c_goldenRatioConjugate = 0.61803398875f;

static const float c_pi = 3.14159265359f;

inline pcg32_random_t GetRNG()
{
    static uint64_t s_sequence = 0;
    s_sequence++;

    pcg32_random_t rng;
#if DETERMINISTIC()
    pcg32_srandom_r(&rng, 0x1337FEED, s_sequence);
#else
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<uint32_t> dist;
    pcg32_srandom_r(&rng, dist(generator), s_sequence);
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
    std::vector<float> ret(numSamples);
    for (size_t index = 0; index < numSamples; ++index)
        ret[index] = (float(index) + 0.5f) / float(numSamples);
    return ret;
}

std::vector<float> Sequence_RegularRandomoffset(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    float offset = RandomFloat01(rng);
    std::vector<float> ret(numSamples);
    for (size_t index = 0; index < numSamples; ++index)
        ret[index] = std::fmodf(offset + float(index) / float(numSamples), 1.0f);
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
    std::vector<float> ret(numSamples);
    ret[0] = 0.0f;
    for (size_t index = 1; index < numSamples; ++index)
        ret[index] = fmodf(ret[index - 1] + c_goldenRatioConjugate, 1.0f);
    return ret;
}

std::vector<float> Sequence_GoldenRatioRandomOffset(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    ret[0] = RandomFloat01(rng);
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

struct Sequence
{
    const char* name = nullptr;
    std::vector<float>(*fn)(size_t numSamples) = nullptr;
    std::vector<float> sequence;

    std::vector<float> avgError;
    std::vector<float> avgSquaredError;
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
        //{"Regular", Sequence_Regular},
        {"RegularRandomOffset", Sequence_RegularRandomoffset},
        {"Stratified", Sequence_Stratified},
        //{"GoldenRatio", Sequence_GoldenRatio},
        {"GoldenRatioRandomOffset", Sequence_GoldenRatioRandomOffset},
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
        printf("%s\n", functions[functionIndex].name);

        // for each sequence
        for (size_t sequenceIndex = 0; sequenceIndex < _countof(sequences); ++sequenceIndex)
        {
            sequences[sequenceIndex].avgError.resize(c_numSamples, 0.0f);
            sequences[sequenceIndex].avgSquaredError.resize(c_numSamples, 0.0f);

            // for each test
            int lastPercent = -1;
            for (size_t testIndex = 0; testIndex < c_numtests; ++testIndex)
            {
                int percent = int(100.0f * float(testIndex) / float(c_numtests - 1));
                if (lastPercent != percent)
                {
                    lastPercent = percent;
                    printf("\r    %s: %i%%", sequences[sequenceIndex].name, percent);
                }

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

                    // update the avgError and avgSquaredError;
                    float error = y - functions[functionIndex].actualValue;
                    sequences[sequenceIndex].avgError[sampleCount] = Lerp(sequences[sequenceIndex].avgError[sampleCount], error, 1.0f / float(testIndex + 1));
                    sequences[sequenceIndex].avgSquaredError[sampleCount] = Lerp(sequences[sequenceIndex].avgSquaredError[sampleCount], error * error, 1.0f / float(testIndex + 1));
                }
            }
            printf("\r    %s: 100%%\n", sequences[sequenceIndex].name);
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
                fprintf(file, ",\"%f\"", std::sqrt(sequences[sequenceIndex].avgSquaredError[sampleCount]));
            fprintf(file, "\n");
        }

        // all done
        fclose(file);
    }

	return 0;
}

/*
Blog:
* show the functions too. have google or wolfram graph em.
* NOTE: csv has root mean squared error, per this blog post: https://blog.demofox.org/2021/04/01/mean-squared-error-is-variance/
* motivation: golden ratio is better sampling. having it in sorted order is better for cache, and better for ray marching, to composite things properly.
* Show the right way to do it, show that it works.
* show your failed attempt and the things you found along the way
* log/log plots are in out.

1d function analysis...
* uniform white noise is not a great choice unsurprisingly
* regular sampling is ok, but has aliasing problems.
 * stratified converges at about the same rate but doesn't have the aliasing problems.
 * Random offset regular sampling isn't bad. it seems to beat stratified?? except in step
* random offset golden ratio is good. 
* step function makes golden ratio and regular sampling go nuts.
 * these functions are deterministic though so aren't getting anything from averaging a bunch of runs.

! regular with a random offset is pretty decent! motivation for using blue noise that way.
! removing "golden ratio" and "uniform" because they are cluttery and don't really belong, not a fair test.

*/