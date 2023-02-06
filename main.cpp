#include <stdio.h>
#include <vector>
#include <random>
#include <direct.h>

#include "pcg/pcg_basic.h"

#define DETERMINISTIC() true

static const size_t c_numSamples = 100;
static const size_t c_numtests = 10000;

static const float c_goldenRatio = 1.61803398875f;
static const float c_goldenRatioConjugate = 0.61803398875f; // 1/goldenRatio or goldenRatio-1
static const float c_smallGoldenRatioConjugate = 1.0f / (c_goldenRatio * c_goldenRatio);

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

// This will return >= <numSamples> samples.
std::vector<float> Sequence_FibonacciWordRandomOffset_Internal(size_t numSamples, float big)
{
    std::vector<float> ret;
    ret.reserve(numSamples); // Best guess for final sample count

    float small = big * c_goldenRatioConjugate;

    float run = c_goldenRatioConjugate;
    float sample = 0.0f;
    while (sample < 1.0f)
    {
        run += c_goldenRatio;
        float shift = std::floor(run);
        run -= shift;
        sample += (shift == 1.0f) ? small : big;
        ret.push_back(sample);
    }

    // apply a random offset
    pcg32_random_t rng = GetRNG();
    float offset = RandomFloat01(rng);
    for (float& f : ret)
        f = std::fmodf(f + offset, 1.0f);

    return ret;
}

std::vector<float> Sequence_FibonacciWordRandomOffset(size_t numSamples)
{
    // TODO: n(F) = Floor[ Log(F Sqrt(5) + 1/2)/Log(Phi)]
    // https://stackoverflow.com/a/5162856

    size_t numSamplesOrig = numSamples;
    int attemptCount = 1;
    std::vector<float> ret = Sequence_FibonacciWordRandomOffset_Internal(numSamples, 1.0f / float(numSamples));

    while (ret.size() > numSamplesOrig)
    {
        numSamples--;
        if (numSamples <= 0)
            break;
        std::vector<float> temp = Sequence_FibonacciWordRandomOffset_Internal(numSamples, 1.0f / float(numSamples));
        attemptCount++;
        if (temp.size() < numSamplesOrig)
            break;
        ret = temp;
    }

    return ret;
}

std::vector<float> Sequence_SmallGoldenRatioRandomOffset(size_t numSamples)
{
    pcg32_random_t rng = GetRNG();
    std::vector<float> ret(numSamples);
    ret[0] = RandomFloat01(rng);
    for (size_t index = 1; index < numSamples; ++index)
        ret[index] = fmodf(ret[index - 1] + c_smallGoldenRatioConjugate, 1.0f);
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

void IntegrationTests()
{
    _mkdir("out");

    Sequence sequences[] =
    {
        {"Uniform", Sequence_WhiteNoise},
        {"RegularRandomOffset", Sequence_RegularRandomoffset},
        {"Stratified", Sequence_Stratified},
        {"GoldenRatioRandomOffset", Sequence_GoldenRatioRandomOffset},
        {"FibonacciWordRandomOffset", Sequence_FibonacciWordRandomOffset},

        //{"Regular", Sequence_Regular}, // The same each run so unfair to include it, and doesn't add any new info anyways
        //{"GoldenRatio", Sequence_GoldenRatio}, // The same each run so unfair to include it, and doesn't add any new info anyways
        //{"SmallGoldenRatioRandomOffset", Sequence_SmallGoldenRatioRandomOffset},  // 99.9% identical to GoldenRatioRandomOffset, omitting it
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
            fprintf(file, "\"%i\"", (int)sampleCount + 1);
            for (size_t sequenceIndex = 0; sequenceIndex < _countof(sequences); ++sequenceIndex)
                fprintf(file, ",\"%f\"", std::sqrt(sequences[sequenceIndex].avgSquaredError[sampleCount]));
            fprintf(file, "\n");
        }

        // all done
        fclose(file);
    }
}

// this will return the <fibonacciIndex>th fibnonacci number - 1 numbers of samples
void SortedGoldenRatioSequenceTest(int fibonacciIndex)
{
    // Fibonacci index to Fibonacci number
    // https://r-knott.surrey.ac.uk/Fibonacci/fibFormula.html
    // This could (and show) be calculated or supplied before generating the sequence.
    // Note: The +2 is because there are differences in S_k and F_k. Example:
    // F_4 = 3 per this (https://en.wikipedia.org/wiki/Fibonacci_number)
    // S_4 = 01001010 = 8 (has 8 digits) per this https://en.wikipedia.org/wiki/Fibonacci_word#Other_properties
    // F_6 = 8, so S_k = F_(k+2)
    int fibonacciNumber = int(
        (
            std::pow(c_goldenRatio, float(fibonacciIndex + 2)) -
            std::pow(-c_goldenRatioConjugate, float(fibonacciIndex + 2))
        )
        / sqrt(5.0f)
    );

    printf("Fib(%i) = %i\n", fibonacciIndex, fibonacciNumber);

    // Calculate the large (A) and small (B) gaps between values
    // https://en.wikipedia.org/wiki/Fibonacci_word#Other_properties
    // in the talk about the unit circle and the golden angle
    float stepBig = std::pow(c_goldenRatio, float(-fibonacciIndex));
    float stepSmall = std::pow(c_goldenRatio, float(-fibonacciIndex - 1));

    printf("Big = %f\nSmall = %f\n", stepBig, stepSmall);

    // We will make <fibonacciNumber>-1 samples
    float sample = 0.0f;
    std::vector<float> samples;
    std::vector<bool> fibonacciWord;
    float run = c_goldenRatio - std::floor(c_goldenRatio);
    for (int sampleIndex = 1; sampleIndex <= fibonacciNumber + 1; ++sampleIndex)
    {     
        // Calculate the <sampleIndex>th digit of the infinite Fibonacci word
        // https://en.wikipedia.org/wiki/Fibonacci_word#Closed-form_expression_for_individual_digits
        bool useSmallStep = int(
            2.0f +
            std::floor(float(sampleIndex) * c_goldenRatio) -
            std::floor(float(sampleIndex + 1) * c_goldenRatio)
            ) == 1;

        // alternate method for calculating word
        run += c_goldenRatio;
        float shift = std::floor(run);
        run -= shift;

        // an alternate method to calculate the infinite fibonacci word
        bool useSmallStep2 = (shift == 1.0f);
        if (useSmallStep != useSmallStep2)
            printf("ERROR! useSmallStep != useSmallStep2!\n");

        samples.push_back(sample);
        sample += useSmallStep ? stepSmall : stepBig;
        fibonacciWord.push_back(useSmallStep);
    }

    // print out word
    printf("  Word: ");
    for (bool b : fibonacciWord)
        printf(b ? "1" : "0");
    printf("\n");

    // print out sequence
    printf("  Sorted Sequence:\n   ");
    for (float sample : samples)
        printf(" %0.3f", sample);
    printf("\n");

    // print out ground truth
    printf("  Unsorted Sequence, sorted:\n   ");
    std::vector<float> unsortedSamples;
    sample = 0.0f;
    for (size_t index = 0; index < samples.size(); ++index)
    {
        sample = std::fmod(sample + c_smallGoldenRatioConjugate, 1.0f);
        unsortedSamples.push_back(sample);
    }
    std::sort(unsortedSamples.begin(), unsortedSamples.end());
    for (size_t index = 0; index < unsortedSamples.size(); ++index)
        printf(" %0.3f", unsortedSamples[index]);
    printf("\n");

    printf("  This Word: ");
    for (size_t index = 0; index < unsortedSamples.size(); ++index)
    {
        float diff;
        if (index + 1 < unsortedSamples.size())
            diff = unsortedSamples[index + 1] - unsortedSamples[index];
        else
            diff = (1.0f + unsortedSamples[0]) - unsortedSamples[index];

        float diffBig = std::abs(diff - stepBig);
        float diffSmall = std::abs(diff - stepSmall);
        printf(diffSmall < diffBig ? "1" : "0");
    }
    printf("\n");

    // Next
    printf("\n");
}

int main(int argc, char** argv)
{
    IntegrationTests();

    //auto blah = Sequence_FibonacciWordRandomOffset_Internal(5, c_goldenRatioConjugate);

    // TODO: i think it's working now, but...
    // 2) The values are rotated i think (verify)
    // TODO: fib(5) = 13 may not match! look into it.

    // TODO: test hese samples in a convergence test for fibonacci counts of samples? with random offsets.
    // TODO: make a nice simple function to generate these samples, that you call from the testing function

    //for (int i = 1; i < 7; ++i)
    //    SortedGoldenRatioSequenceTest(i);

    // TODO: graph # of attempts for the fibonacci word thing? maybe need to find the nearest fibonacci number? seems like something to solve or improve.


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