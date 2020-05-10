package com.github.psambit9791.jdsp;

import com.github.psambit9791.jdsp.transform.Hilbert;
import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;

public class TestHilbert {

    private double[] signal = {1.0   ,  0.702, -0.168, -1.023, -1.148, -0.311,  0.884,  1.373,
            0.596, -0.816, -1.471, -0.607,  0.939,  1.461,  0.303, -1.208,
            -1.219,  0.301,  1.374,  0.576, -0.957, -1.051,  0.363,  1.12 ,
            0.127, -0.955, -0.418,  0.723,  0.534, -0.526, -0.542,  0.397,
            0.499, -0.337, -0.438,  0.332,  0.372, -0.373, -0.292,  0.449,
            0.172, -0.541,  0.026,  0.593, -0.322, -0.502,  0.666,  0.152,
            -0.866,  0.454,  0.633, -1.021,  0.16 ,  0.984, -1.083, -0.023,
            1.18 , -1.191,  0.035,  1.193, -1.382,  0.397,  0.909, -1.499,
            0.979,  0.202, -1.196,  1.392, -0.759, -0.245,  1.033, -1.238,
            0.852, -0.146, -0.529,  0.912, -0.921,  0.63 , -0.196, -0.221,
            0.512, -0.635,  0.602, -0.459,  0.258, -0.046, -0.147,  0.304,
            -0.418,  0.49 , -0.524,  0.525, -0.497,  0.446, -0.377,  0.296,
            -0.207,  0.116, -0.028, -0.053,  0.123, -0.176,  0.21 , -0.219,
            0.201, -0.15 ,  0.063,  0.061, -0.224,  0.422, -0.646,  0.884,
            -1.115,  1.311, -1.442,  1.475, -1.384,  1.154, -0.791,  0.328,
            0.177, -0.641,  0.978, -1.112,  1.007, -0.681,  0.214,  0.268,
            -0.626,  0.755, -0.62 ,  0.281,  0.13 , -0.449,  0.548, -0.388,
            0.045,  0.314, -0.498,  0.387, -0.017, -0.402,  0.573, -0.315,
            -0.248,  0.686, -0.566, -0.122,  0.804, -0.758};

    @Test
    public void testHilbertAmplitudeEnvelope() {
        double[] result = {1.199, 1.165, 0.985, 1.215, 1.263, 1.262, 1.362, 1.388, 1.41 ,
                1.462, 1.476, 1.484, 1.501, 1.5  , 1.49 , 1.478, 1.453, 1.418,
                1.385, 1.346, 1.293, 1.238, 1.186, 1.125, 1.059, 1.003, 0.942,
                0.871, 0.811, 0.767, 0.712, 0.648, 0.609, 0.589, 0.555, 0.512,
                0.499, 0.513, 0.515, 0.5  , 0.506, 0.551, 0.6  , 0.622, 0.636,
                0.687, 0.772, 0.845, 0.88 , 0.909, 0.975, 1.073, 1.16 , 1.206,
                1.223, 1.253, 1.323, 1.408, 1.468, 1.481, 1.459, 1.441, 1.457,
                1.499, 1.536, 1.536, 1.489, 1.41 , 1.328, 1.272, 1.251, 1.249,
                1.242, 1.209, 1.141, 1.041, 0.924, 0.807, 0.71 , 0.651, 0.629,
                0.637, 0.656, 0.676, 0.689, 0.689, 0.683, 0.665, 0.642, 0.611,
                0.578, 0.539, 0.497, 0.452, 0.403, 0.354, 0.3  , 0.25 , 0.2  ,
                0.166, 0.153, 0.176, 0.229, 0.297, 0.38 , 0.468, 0.567, 0.668,
                0.78 , 0.89 , 1.008, 1.123, 1.238, 1.344, 1.442, 1.525, 1.585,
                1.625, 1.629, 1.608, 1.542, 1.451, 1.321, 1.175, 1.017, 0.869,
                0.766, 0.71 , 0.723, 0.759, 0.788, 0.798, 0.75 , 0.678, 0.561,
                0.444, 0.387, 0.405, 0.502, 0.605, 0.64 , 0.652, 0.574, 0.489,
                0.548, 0.686, 0.855, 1.085, 0.956, 0.782};

        Hilbert h = new Hilbert(this.signal);
        h.hilbert_transform();
        double[] out = h.get_amplitude_envelope();
        assertArrayEquals(result, out, 0.01);
    }

    @Test
    public void testHilbertInstantaneousPhase() {
        double[] result = {-0.584,   0.924,   1.742,   2.571,   3.572,   4.463,   5.418,
                6.429,   7.418,   8.446,   9.506,  10.574,  11.672,  12.795,
                13.932,  15.094,  16.283,  17.493,  18.723,  19.978,  21.254,
                22.548,  23.873,  25.223,  26.583,  27.965,  29.385,  30.825,
                32.268,  33.743,  35.263,  36.788,  38.309,  39.879,  41.502,
                43.117,  44.712,  46.366,  48.092,  49.811,  51.489,  53.212,
                55.021,  56.854,  58.65 ,  60.442,  62.301,  64.222,  66.155,
                68.067,  69.98 ,  71.943,  73.966,  76.015,  78.057,  80.092,
                82.151,  84.26 ,  86.418,  88.599,  90.78 ,  92.956,  95.145,
                97.373,  99.651, 101.97 , 104.31 , 106.655, 108.993, 111.333,
                113.696, 116.105, 118.566, 121.072, 123.611, 126.167, 128.726,
                131.271, 133.797, 136.313, 138.849, 141.443, 144.103, 146.831,
                149.61 , 152.434, 155.292, 158.175, 161.083, 164.003, 166.939,
                169.877, 172.824, 175.765, 178.708, 181.633, 184.543, 187.407,
                190.207, 192.882, 195.411, 197.974, 200.651, 203.462, 206.332,
                209.242, 212.169, 215.108, 218.049, 220.988, 223.928, 226.859,
                229.786, 232.701, 235.611, 238.505, 241.393, 244.263, 247.122,
                249.962, 252.783, 255.582, 258.348, 261.081, 263.752, 266.365,
                268.889, 271.361, 273.843, 276.356, 278.937, 281.533, 284.14 ,
                286.732, 289.241, 291.661, 293.855, 295.992, 298.321, 300.716,
                303.19 , 305.641, 307.945, 310.147, 312.119, 314.132, 316.454,
                318.759, 321.014, 323.832};

        Hilbert h = new Hilbert(this.signal);
        h.hilbert_transform();
        double[] out = h.get_instantaneous_phase();
        assertArrayEquals(result, out, 0.01);
    }

    @Test
    public void testHilbertInstantaneousFrequency() {
        double Fs = 150.0;
        double[] result = {36.001, 19.528, 19.791, 23.897, 21.271, 22.799, 24.136, 23.611,
                24.542, 25.306, 25.497, 26.213, 26.81 , 27.144, 27.741, 28.385,
                28.887, 29.364, 29.961, 30.462, 30.892, 31.632, 32.229, 32.468,
                32.993, 33.9  , 34.377, 34.449, 35.213, 36.287, 36.407, 36.311,
                37.481, 38.746, 38.555, 38.078, 39.486, 41.205, 41.038, 40.059,
                41.134, 43.187, 43.76 , 42.876, 42.781, 44.38 , 45.86 , 46.147,
                45.646, 45.67 , 46.863, 48.296, 48.916, 48.749, 48.582, 49.155,
                50.349, 51.518, 52.068, 52.068, 51.948, 52.259, 53.19 , 54.383,
                55.362, 55.863, 55.983, 55.816, 55.863, 56.412, 57.511, 58.752,
                59.826, 60.614, 61.02 , 61.092, 60.757, 60.304, 60.065, 60.543,
                61.927, 63.503, 65.126, 66.344, 67.418, 68.23 , 68.827, 69.423,
                69.71 , 70.092, 70.14 , 70.354, 70.211, 70.259, 69.829, 69.471,
                68.373, 66.845, 63.861, 60.375, 61.187, 63.909, 67.108, 68.516,
                69.471, 69.877, 70.163, 70.211, 70.163, 70.187, 69.972, 69.877,
                69.59 , 69.471, 69.089, 68.946, 68.516, 68.254, 67.8  , 67.346,
                66.821, 66.033, 65.246, 63.765, 62.381, 60.256, 59.015, 59.253,
                59.993, 61.617, 61.975, 62.238, 61.879, 59.898, 57.773, 52.378,
                51.017, 55.601, 57.176, 59.062, 58.513, 55.004, 52.569, 47.078,
                48.057, 55.434, 55.028, 53.834, 67.275};

        Hilbert h = new Hilbert(this.signal);
        h.hilbert_transform();
        double[] out = h.get_instantaneous_frequency(Fs);
        assertArrayEquals(result, out, 0.5);

    }
}