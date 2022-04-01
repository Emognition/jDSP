package com.github.psambit9791.jdsp.signal;

import org.apache.commons.math3.util.MathArrays;
import org.jtransforms.fft.DoubleFFT_1D;

public class FastResample {

    private final int signalLength;

    private final double[] dft;

    public FastResample(double[] signal) {
        signalLength = signal.length;
        dft = new double[2 * signalLength];
        System.arraycopy(signal, 0, dft, 0, signalLength);
    }

    public double[] resample(int targetSamples) {
        calculateDft();
        double[] inverseDft = prepareInverseDftInput(targetSamples);
        calculateInverseDft(inverseDft);

        double scalingFactor = targetSamples / (double) signalLength;
        return MathArrays.scale(scalingFactor, getRealArray(inverseDft));
    }

    private void calculateDft() {
        DoubleFFT_1D fft = new DoubleFFT_1D(signalLength);
        fft.realForwardFull(dft);
    }

    private void calculateInverseDft(double[] inverseDft) {
        DoubleFFT_1D fft = new DoubleFFT_1D(inverseDft.length / 2);
        fft.complexInverse(inverseDft, true);
    }

    private double[] prepareInverseDftInput(int targetSamples) {
        double[] inverseDft = new double[2 * targetSamples];
        int N = Math.min(targetSamples, signalLength);

        int firstHalfLength = N % 2 == 0 ? N : N + 1;
        System.arraycopy(dft, 0, inverseDft, 0, firstHalfLength);

        int secondHalfLength = N % 2 == 0 ? N : N - 1;
        int srcOffset = dft.length - secondHalfLength;
        int dstOffset = inverseDft.length - secondHalfLength;
        System.arraycopy(dft, srcOffset, inverseDft, dstOffset, secondHalfLength);

        if (N % 2 == 0) {
            adjustInverseDftInput(inverseDft, N, targetSamples);
        }

        return inverseDft;
    }

    private void adjustInverseDftInput(double[] inverseDft, int N, int targetSamples) {
        if (N < signalLength) {
            // Downsampling adjustment
            inverseDft[N] += dft[N];
            inverseDft[N + 1] += dft[N + 1];
        } else if (N < targetSamples) {
            // Oversampling adjustment
            int componentIndex = 2 * targetSamples - N;
            inverseDft[componentIndex] /= 2;
            inverseDft[componentIndex + 1] /= 2;

            inverseDft[N] = inverseDft[componentIndex];
            inverseDft[N + 1] = inverseDft[componentIndex + 1];
        }
    }

    private double[] getRealArray(double[] complexArray) {
        double[] realComponents = new double[complexArray.length / 2];
        for (int i = 0; i < realComponents.length; i++) {
            realComponents[i] = complexArray[i * 2];
        }
        return realComponents;
    }
}
