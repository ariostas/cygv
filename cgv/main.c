/*
 * cgv driver: picks the number of prime lanes (2, 3 or 4) for the run.
 *
 * Each lane is a ~62-bit prime; the result needs enough primes to hold the
 * largest |GV| plus one to verify the CRT lift. With too few, cgv still gets
 * the right answer but spends extra passes; with too many, every step pays
 * for arithmetic it did not need. Unless -l is given, a cheap probe at a
 * lower max degree (2 lanes) measures the largest |GV| per degree, and its
 * growth is extrapolated to the requested degree.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cgv.h"

static int input_maxdeg(const char *path) {
    FILE *f = fopen(path, "r");
    int h, n, d;
    if (!f) return -1;
    if (fscanf(f, "%d %d %d", &h, &n, &d) != 3) d = -1;
    fclose(f);
    return d;
}

int main(int argc, char **argv) {
    int lanes = 0, quiet = 0;
    const char *in = NULL;
    char **args = malloc((argc + 1) * sizeof(char *));
    int na = 0;
    args[na++] = argv[0];
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-l") && i + 1 < argc) { lanes = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "-q")) quiet = 1;
        if (!strcmp(argv[i], "-t") && i + 1 < argc) { args[na++] = argv[i++]; args[na++] = argv[i]; continue; }
        if (!strcmp(argv[i], "-g")) { args[na++] = argv[i]; if (i + 1 < argc && strlen(argv[i + 1]) == 1 && argv[i + 1][0] >= '0' && argv[i + 1][0] <= '9') args[na++] = argv[++i]; continue; }
        if (argv[i][0] != '-') in = argv[i];
        args[na++] = argv[i];
    }
    args[na] = NULL;
    int D = in ? input_maxdeg(in) : -1;
    if (!lanes) {
        lanes = 3;
        /* probe at 3/4 of the degree: cheap for standard gradings (D-6..D-8) and for
         * gradings with large degree values (p-vectors) alike */
        int Dp = (int)(0.75 * D);
        if (D > 0 && D <= 12) lanes = 2;
        else if (D > 0 && Dp >= 8 && Dp < D) {
            /* probe quietly */
            char **pargs = malloc((na + 2) * sizeof(char *));
            int np = 0;
            for (int i = 0; i < na; i++) pargs[np++] = args[i];
            pargs[np++] = "-q"; pargs[np] = NULL;
            CgvProbe pr = {Dp, 0, NULL};
            if (cgv_entry_nl2(np, pargs, &pr) == 0 && pr.maxbits) {
                /* least-squares slope of max bits vs degree over the upper half */
                double sx = 0, sy = 0, sxx = 0, sxy = 0; int m = 0, dl = 0; double bl = 0;
                for (int d = pr.maxdeg / 2; d <= pr.maxdeg; d++) {
                    if (pr.maxbits[d] < 0) continue;
                    sx += d; sy += pr.maxbits[d]; sxx += (double)d * d; sxy += d * pr.maxbits[d]; m++;
                    dl = d; bl = pr.maxbits[d];
                }
                double slope = m >= 2 ? (m * sxy - sx * sy) / (m * sxx - sx * sx) : 1.0;
                if (slope < 0) slope = 0;
                /* scatter of the fit: erratic growth (e.g. p-vector gradings) gets a wider margin */
                double icpt = m ? (sy - slope * sx) / m : 0, ss = 0;
                for (int d = pr.maxdeg / 2; d <= pr.maxdeg; d++) if (pr.maxbits[d] >= 0) { double r = pr.maxbits[d] - (icpt + slope * d); ss += r * r; }
                double rmse = m ? sqrt(ss / m) : 0;
                double bits = bl + slope * (D - dl) + 4 + 2 * rmse; /* predicted log2 max|GV|, plus a margin */
                /* k primes (~2^62 each) lift |GV| < 2^(62k - 1); one more verifies */
                int k = 1;
                while (62.0 * k - 2 < bits) k++;
                lanes = k + 1;
                if (lanes < 2) lanes = 2;
                if (lanes > 4) lanes = 4;
                if (!quiet) fprintf(stderr, "lanes: probe to degree %d, max |GV| ~2^%.0f at degree %d, slope %.2f bits/degree (rmse %.1f) -> ~2^%.0f at %d: %d lanes\n",
                                    pr.maxdeg, bl, dl, slope, rmse, bits - 4 - 2 * rmse, D, lanes);
            }
            free(pargs);
        }
    }
    int rc;
    if (lanes <= 2) rc = cgv_entry_nl2(na, args, NULL);
    else if (lanes == 3) rc = cgv_entry_nl3(na, args, NULL);
    else rc = cgv_entry_nl4(na, args, NULL);
    free(args);
    return rc;
}
