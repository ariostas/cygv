/* Entry points of the cgv core, one per number of prime lanes (see main.c). */
#ifndef CGV_H
#define CGV_H
typedef struct {
    int maxdeg_override;   /* > 0: run with this max degree instead of the input's */
    int maxdeg;            /* out: the max degree used */
    double *maxbits;       /* out: [0..maxdeg] largest log2|GV| seen at each degree (-1 if none) */
} CgvProbe;
int cgv_entry_nl2(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl3(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl4(int argc, char **argv, CgvProbe *probe);
/* the same with 256-bit lattice keys; narrower entries return CGV_RC_WIDE when the coordinates need wider keys
 * (after reading the input and reducing the curve lattice, before the passes) */
int cgv_entry_nl2w(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl3w(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl4w(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl2x(int argc, char **argv, CgvProbe *probe);   /* 512-bit keys */
int cgv_entry_nl3x(int argc, char **argv, CgvProbe *probe);
int cgv_entry_nl4x(int argc, char **argv, CgvProbe *probe);
#define CGV_RC_WIDE 75   /* 75: rerun with 256-bit keys, 76: with 512-bit keys */
#endif
