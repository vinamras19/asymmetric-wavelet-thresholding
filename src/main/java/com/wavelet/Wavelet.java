package com.wavelet;

public class Wavelet {

    // Daubechies-2 (db2): 4-tap orthogonal wavelet, 2 vanishing moments
    // Coefficients from Daubechies, "Ten Lectures on Wavelets" (1992), Table 6.1
    private static final double[] H = {
            0.4829629131445341,
            0.8365163037378079,
            0.2241438680420134,
            -0.1294095225512604
    };

    // High-pass derived from low-pass: g[n] = (-1)^n * h[L-1-n]
    private static final double[] G = {
            -0.1294095225512604,
            -0.2241438680420134,
            0.8365163037378079,
            -0.4829629131445341
    };

    public static double[] forward(double[] x, int k) {
        double[] out = x.clone();
        int n = x.length;
        for (int level = 0; level < k; level++) {
            forwardStep(out, n >> level);
        }
        return out;
    }

    public static double[] inverse(double[] coeffs, int k) {
        double[] out = coeffs.clone();
        int n = coeffs.length;
        for (int level = k - 1; level >= 0; level--) {
            inverseStep(out, n >> level);
        }
        return out;
    }

    private static void forwardStep(double[] x, int n) {
        int half = n / 2;
        double[] temp = new double[n];
        for (int k = 0; k < half; k++) {
            double a = 0.0, d = 0.0;
            for (int i = 0; i < H.length; i++) {
                int idx = (2 * k + i) % n;
                a += H[i] * x[idx];
                d += G[i] * x[idx];
            }
            temp[k] = a;
            temp[half + k] = d;
        }
        System.arraycopy(temp, 0, x, 0, n);
    }

    private static void inverseStep(double[] x, int n) {
        int half = n / 2;
        double[] temp = new double[n];
        for (int k = 0; k < half; k++) {
            double a = x[k];
            double d = x[half + k];
            for (int i = 0; i < H.length; i++) {
                int idx = (2 * k + i) % n;
                temp[idx] += H[i] * a + G[i] * d;
            }
        }
        System.arraycopy(temp, 0, x, 0, n);
    }

    // packed layout after k-level DWT: [A_k | D_k | D_{k-1} | ... | D_1]
    // A_k at [0, n/2^k), D_j at [n/2^j, n/2^(j-1))
    public static double[] extractDetail(double[] coeffs, int j, int n) {
        int start = n >> j;
        int end = n >> (j - 1);
        double[] d = new double[end - start];
        System.arraycopy(coeffs, start, d, 0, d.length);
        return d;
    }

    public static double[] extractApprox(double[] coeffs, int k, int n) {
        int len = n >> k;
        double[] a = new double[len];
        System.arraycopy(coeffs, 0, a, 0, len);
        return a;
    }

    public static void writeDetail(double[] coeffs, double[] detail, int j, int n) {
        int start = n >> j;
        System.arraycopy(detail, 0, coeffs, start, detail.length);
    }

    public static void writeApprox(double[] coeffs, double[] approx, int k, int n) {
        System.arraycopy(approx, 0, coeffs, 0, approx.length);
    }
}