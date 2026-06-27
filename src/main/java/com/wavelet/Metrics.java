package com.wavelet;

public class Metrics {

    // peak-to-peak SNR in dB (Eq. 18)
    public static double snrDb(double[] signal, double[] reference) {
        double pp = peakToPeak(reference);
        double sumSq = 0.0;
        for (int i = 0; i < signal.length; i++) {
            double d = signal[i] - reference[i];
            sumSq += d * d;
        }
        double rms = Math.sqrt(sumSq / signal.length);
        return 20.0 * Math.log10(pp / rms);
    }

    private static double peakToPeak(double[] x) {
        double min = x[0], max = x[0];
        for (double v : x) {
            if (v < min) min = v;
            if (v > max) max = v;
        }
        return max - min;
    }
}