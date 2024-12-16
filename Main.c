#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846
#define FRAME_SIZE 256
#define HOP_IN 128
#define HOP_OUT 121 // Adjust hop ratio for 0.95x speed

// Function to perform FFT (simplified)
void fft(float *input, float complex *output, int n) {
    if (n == 1) {
        output[0] = input[0];
        return;
    }

    float complex even[n / 2];
    float complex odd[n / 2];
    float input_even[n / 2];
    float input_odd[n / 2];

    for (int i = 0; i < n / 2; i++) {
        input_even[i] = input[2 * i];
        input_odd[i] = input[2 * i + 1];
    }

    fft(input_even, even, n / 2);
    fft(input_odd, odd, n / 2);

    for (int k = 0; k < n / 2; k++) {
        float complex t = cexp(-2.0 * I * PI * k / n) * odd[k];
        output[k] = even[k] + t;
        output[k + n / 2] = even[k] - t;
    }
}

// Function to perform IFFT (simplified)
void ifft(float complex *input, float *output, int n) {
    float complex conj_input[n];
    for (int i = 0; i < n; i++) {
        conj_input[i] = conj(input[i]);
    }

    float complex temp_output[n];
    fft((float *)conj_input, temp_output, n);

    for (int i = 0; i < n; i++) {
        output[i] = creal(temp_output[i]) / n;
    }
}

// Time-stretching function
void time_stretch(float *input, float *output, int input_size, float stretch_ratio) {
    int num_frames = (input_size - FRAME_SIZE) / HOP_IN + 1;
    float complex fft_frame[FRAME_SIZE];
    float complex fft_prev[FRAME_SIZE];
    float prev_phase[FRAME_SIZE] = {0};
    float phase_accum[FRAME_SIZE] = {0};
    float complex time_domain[FRAME_SIZE];
    float synthesis_frame[FRAME_SIZE];
    int output_index = 0;

    // Process each frame
    for (int i = 0; i < num_frames; i++) {
        int input_offset = i * HOP_IN;

        // Apply FFT to the current frame
        fft(&input[input_offset], fft_frame, FRAME_SIZE);

        // Calculate phase difference and adjust
        for (int k = 0; k < FRAME_SIZE / 2 + 1; k++) {
            float magnitude = cabs(fft_frame[k]);
            float phase = carg(fft_frame[k]);
            float delta_phase = phase - prev_phase[k];
            prev_phase[k] = phase;

            // Phase unwrapping
            delta_phase -= 2 * PI * k * HOP_IN / FRAME_SIZE;
            delta_phase = fmod(delta_phase + PI, 2 * PI) - PI;

            // Phase accumulation
            phase_accum[k] += delta_phase + 2 * PI * k * HOP_OUT / FRAME_SIZE;

            // Update FFT frame
            fft_frame[k] = magnitude * cexp(I * phase_accum[k]);
        }

        // Apply IFFT
        ifft(fft_frame, synthesis_frame, FRAME_SIZE);

        // Overlap-Add to output buffer
        for (int j = 0; j < FRAME_SIZE; j++) {
            if (output_index + j < input_size) {
                output[output_index + j] += synthesis_frame[j];
            }
        }
        output_index += HOP_OUT;
    }
}

int main() {
    // Example: Generate a sine wave as input
    int input_size = 1024;
    int output_size = (int)(input_size / 0.95);
    float input[input_size];
    float output[output_size];
    memset(output, 0, sizeof(output));

    for (int i = 0; i < input_size; i++) {
        input[i] = sin(2 * PI * i / input_size);
    }

    // Apply time-stretching
    time_stretch(input, output, input_size, 0.95);

    // Print the output
    for (int i = 0; i < output_size; i++) {
        printf("%.6f\n", output[i]);
    }

    return 0;
}
