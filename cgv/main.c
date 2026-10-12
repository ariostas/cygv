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
#define _POSIX_C_SOURCE 200809L   /* fileno, dup2, lseek */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#include <io.h>
#define dup2 _dup2
#define fileno _fileno
#define lseek _lseeki64
#else
#include <unistd.h>
#endif
#include "cgv.h"

static int input_maxdeg(const char *path) {
    FILE *f = fopen(path, "r");
    int h, n, d, ver;
    if (!f) return -1;
    if (fscanf(f, " cgv %d", &ver) != 1) { fclose(f); return -1; }   /* "cgv 2" header */
    if (fscanf(f, "%d %d %d", &h, &n, &d) != 3) d = -1;
    fclose(f);
    return d;
}

/* input on stdin is read by every entry called (probe, wider-key rerun): if fd 0 is not seekable (a pipe), spool it
 * into a temporary file on fd 0 (on Windows under %TEMP%: tmpfile() there needs a writable drive root) */
static int from_stdin;
static long long stdin_off;   /* where the input starts on fd 0 */
static void spool_stdin(void) {
    from_stdin = 1;
    if ((stdin_off = lseek(0, 0, SEEK_CUR)) >= 0) return;
    stdin_off = 0;
#ifdef _WIN32
    char *tn = _tempnam(getenv("TEMP"), "cgv");
    FILE *t = tn ? fopen(tn, "w+bTD") : NULL;   /* T: temporary, D: deleted when closed */
    free(tn);
#else
    FILE *t = tmpfile();
#endif
    char b[1 << 16]; size_t n;
    if (!t) { perror("cgv: temporary file for stdin"); exit(1); }
    while ((n = fread(b, 1, sizeof b, stdin)) > 0) if (fwrite(b, 1, n, t) != n) { perror("cgv: temporary file for stdin"); exit(1); }
    if (fflush(t) || dup2(fileno(t), 0) < 0) { perror("cgv: temporary file for stdin"); exit(1); }
}

static int wide;   /* key width in use: 0 = 128 bits, 1 = 256, 2 = 512 (raised when an entry asks for wider keys) */
static int entry(int lanes, int argc, char **argv, CgvProbe *pr) {
    static int (*const E[3][3])(int, char **, CgvProbe *) = {{cgv_entry_nl2, cgv_entry_nl3, cgv_entry_nl4},
        {cgv_entry_nl2w, cgv_entry_nl3w, cgv_entry_nl4w}, {cgv_entry_nl2x, cgv_entry_nl3x, cgv_entry_nl4x}};
    int l = lanes <= 2 ? 0 : lanes == 3 ? 1 : 2, rc;
    for (;;) {
        if (from_stdin) { clearerr(stdin); fseek(stdin, (long)stdin_off, SEEK_SET); }
        rc = E[wide][l](argc, argv, pr);
        if ((rc != CGV_RC_WIDE && rc != CGV_RC_WIDE + 1) || wide == 2) break;
        wide = rc == CGV_RC_WIDE + 1 ? 2 : wide + 1;
    }
    return rc;
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
    if (!in) spool_stdin();
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
            if (entry(2, np, pargs, &pr) == 0 && pr.maxbits) {
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
    rc = entry(lanes, na, args, NULL);
    free(args);
    if (rc == CGV_RC_WIDE || rc == CGV_RC_WIDE + 1) { fprintf(stderr, "cgv: keys wider than 512 bits asked for (CGV_FORCE_WIDE > 512?): not available\n"); return 1; }
    return rc;
}
