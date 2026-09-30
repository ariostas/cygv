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
#endif
