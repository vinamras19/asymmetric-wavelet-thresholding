package com.wavelet;

public class SignalSynthesizer {

    // synthetic cw-ESR triplet (Tempol-like, paper Fig. 4)
    public static double[] synthesizeEsr(int n) {
        double[] field = linspace(3335.0, 3395.0, n);
        double[] peaks = {3348.0, 3365.0, 3382.0};
        double gamma = 1.2;
        double amplitude = 0.04;
        double baseline = -0.06;

        double[] main = new double[n];
        for (double b0 : peaks) {
            double[] line = lorentzianDerivative(field, b0, gamma);
            for (int i = 0; i < n; i++) main[i] += line[i];
        }

        double pp = peakToPeak(main);
        for (int i = 0; i < n; i++) main[i] *= (amplitude / pp);

        // satellites adjacent to each main peak
        double satOffset = 2.5;
        double satScale = 0.18;
        double satGamma = 0.6;
        double[] sat = new double[n];
        for (double b0 : peaks) {
            double[] line = lorentzianDerivative(field, b0 + satOffset, satGamma);
            for (int i = 0; i < n; i++) sat[i] += line[i];
        }
        double satPp = peakToPeak(sat);
        if (satPp > 0) {
            double scale = amplitude * satScale / satPp;
            for (int i = 0; i < n; i++) sat[i] *= scale;
        }

        double[] signal = new double[n];
        for (int i = 0; i < n; i++) signal[i] = main[i] + sat[i] + baseline;
        return signal;
    }

    public static double[] synthesizeEsr() {
        return synthesizeEsr(4096);
    }

    private static double[] lorentzianDerivative(double[] field, double b0, double gamma) {
        double[] y = new double[field.length];
        for (int i = 0; i < field.length; i++) {
            double u = (field[i] - b0) / gamma;
            y[i] = -2.0 * u / (gamma * Math.pow(1.0 + u * u, 2));
        }
        return y;
    }

    private static double[] linspace(double a, double b, int n) {
        double[] x = new double[n];
        double step = (b - a) / (n - 1);
        for (int i = 0; i < n; i++) x[i] = a + i * step;
        return x;
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