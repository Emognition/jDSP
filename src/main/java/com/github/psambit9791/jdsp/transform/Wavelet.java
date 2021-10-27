/*
 *
 *  * Copyright (c) 2020 Sambit Paul
 *  *
 *  * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *  *
 *  * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *  *
 *  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

package com.github.psambit9791.jdsp.transform;

import com.github.psambit9791.jdsp.misc.UtilMethods;
import com.github.psambit9791.jdsp.signal.Convolution;
import com.github.psambit9791.jdsp.signal.Generate;
import org.apache.commons.math3.complex.Complex;

public class Wavelet {

    private double[] signal;
    private int[] widths;

    public Wavelet(double[] signal, int[] widths) {
        this.signal = signal;
        this.widths = widths;
    }

    private Complex[] ricker_cwt(double[] data, double[] wavelet) {
        Convolution c = new Convolution(data, wavelet);
        return UtilMethods.matToComplex(c.convolve("same"));
    }

    private Complex[] morlet_cwt(double[] data, Complex[] wavelet) {
        Complex[] wavelet_conjugate = new Complex[wavelet.length];
        for (int i=0; i<wavelet.length; i++) {
            wavelet_conjugate[i] = wavelet[i].conjugate();
        }
        double[][] decomp_wvlt = UtilMethods.transpose(UtilMethods.complexTo2D(wavelet_conjugate));

        Convolution c_real = new Convolution(data, decomp_wvlt[0]);
        Convolution c_imag = new Convolution(data, decomp_wvlt[1]);

        double[][] temp = {c_real.convolve("same"), c_imag.convolve("same")};
        temp = UtilMethods.transpose(temp);

        return UtilMethods.matToComplex(temp);
    }

    private Complex[] paul_cwt(double[] data, Complex[] wavelet) {
        double[][] decomp_wvlt = UtilMethods.transpose(UtilMethods.complexTo2D(wavelet));

        Convolution c_real = new Convolution(data, decomp_wvlt[0]);
        Convolution c_imag = new Convolution(data, decomp_wvlt[1]);

        double[][] temp = {c_real.convolve("same"), c_imag.convolve("same")};
        temp = UtilMethods.transpose(temp);

        return UtilMethods.matToComplex(temp);
    }

    // args value: Ignored for Ricker, omega0 for Morlet, order for Paul
    public Complex[][] cwt(String wavelet_type, double args) throws IllegalArgumentException{

        if (!wavelet_type.equals("ricker") && !wavelet_type.equals("morlet") && !wavelet_type.equals("paul")) {
            throw new ArithmeticException("wavelet_type must be 'ricker', 'morlet' or 'paul'");
        }

        Complex[][] output = new Complex[this.widths.length][this.signal.length];

        if (wavelet_type.equals("ricker")) {
            for (int i=0; i<this.widths.length; i++) {
                int N = Math.min(10*this.widths[i], this.signal.length);
                Generate gp = new Generate();
                double[] wavelet = gp.generateRicker(N, this.widths[i]);
                wavelet = UtilMethods.reverse(wavelet);
                output[i] = this.ricker_cwt(this.signal, wavelet);
            }
        }

        else if  (wavelet_type.equals("morlet")) {
            for (int i=0; i<this.widths.length; i++) {
                int N = Math.min(10*this.widths[i], this.signal.length);
                Generate gp = new Generate();
                Complex[] wavelet = gp.generateMorletCWTComplex(N, args, this.widths[i]);
                wavelet = UtilMethods.reverse(wavelet);
                output[i] = this.morlet_cwt(this.signal, wavelet);
            }
        }

        else if  (wavelet_type.equals("paul")) {
            for (int i=0; i<this.widths.length; i++) {
                Generate gp = new Generate();
                double norm = Math.sqrt(1.0/this.widths[i]);
                Complex[] wavelet = gp.generatePaulComplex((int)args, this.widths[i], (double)this.widths[i]);
                for (int w=0; w<wavelet.length; w++) {
                    wavelet[w] = wavelet[w].multiply(norm); //Normalization
                }
                output[i] = this.paul_cwt(this.signal, wavelet);
            }
        }

        return output;
    }
}
