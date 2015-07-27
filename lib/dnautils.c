#include <Python.h>

#define BUF_SIZE 2048
#define DELETION 0
#define INSERTION 1
#define MATCH 2

static PyObject *DNAUtilError;

static PyObject *dnautils_equal(PyObject *self, PyObject *args) {
    PyObject *s1, *s2;
    char *str1, *str2;
    unsigned int i;

    if (!PyArg_ParseTuple(args, "SS", &s1, &s2)) {
        return NULL;
    }

    if (PyString_Size(s1) != PyString_Size(s2)) {
        PyErr_SetString(DNAUtilError, "Sequences have unequal lengths.");
        return NULL;
    }

    if ((str1 = PyString_AsString(s1)) == NULL) {
        return NULL;
    }
    if ((str2 = PyString_AsString(s2)) == NULL) {
        return NULL;
    }

    for (i = 0; i < PyString_Size(s1); i++) {
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N' &&
                str1[i] !='-' && str2[i] != '-') {
            return Py_BuildValue("O", Py_False);
        }
    }
    return Py_BuildValue("O", Py_True);
}

static PyObject *dnautils_hamming(PyObject *self, PyObject *args) {
    PyObject *s1, *s2;
    char *str1, *str2;
    unsigned int i;
    unsigned int distance = 0;

    if (!PyArg_ParseTuple(args, "SS", &s1, &s2)) {
        return NULL;
    }

    if (PyString_Size(s1) != PyString_Size(s2)) {
        PyErr_SetString(DNAUtilError, "Sequences have unequal lengths.");
        return NULL;
    }

    if ((str1 = PyString_AsString(s1)) == NULL) {
        return NULL;
    }
    if ((str2 = PyString_AsString(s2)) == NULL) {
        return NULL;
    }

    for (i = 0; i < PyString_Size(s1); i++) {
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N' &&
                str1[i] != '-' && str2[i] != '-') {
            distance += 1;
        }
    }
    return Py_BuildValue("I", distance);
}

struct alignment {
    char seq1[BUF_SIZE];
    char seq2[BUF_SIZE];
    int seq1_gaps;
    int seq2_gaps;
    int score;
};

int max(int del, int ins, int match) {
    int vals[3] = { del, ins, match };
    int max = vals[0];
    int max_i = 0;
    for (int i = 1; i < 3; i++) {
        if (vals[i] > max) {
            max = vals[i];
            max_i = i;
        }
    }
    return max_i;
}

int get_alignment(char *seq1, char *seq2, int insertion_penalty,
                  int deletion_penalty, int mismatch_penalty, int match_score,
                  int **p, int **q, struct alignment *a) {
    int offset = 0;
    int i = strlen(seq1);
    int j = strlen(seq2);

    a->score = 0;
    while (i > 0 || j > 0) {
        if (i <= 0 || j <= 0) {
            break;
        }
        if (q[i][j] == MATCH) {
            i -= 1;
            j -= 1;
            a->seq1[offset] = seq1[i];
            a->seq2[offset] = seq2[j];
            a->score += match_score;
        } else if (q[i][j] == INSERTION) {
            j -= 1;
            a->seq1[offset] = '-';
            a->seq2[offset] = seq2[j];
            a->score += insertion_penalty;
        } else if (q[i][j] == DELETION) {
            i -= 1;
            a->seq1[offset] = seq1[i];
            a->seq2[offset] = '-';
            a->score += deletion_penalty;
        } else {
            return -1;
        }
        offset++;
    }

    a->seq1[offset] = '\0';
    a->seq2[offset] = '\0';

    return 0;
}

int align_seqs(char *seq1, char *seq2, int insertion_penalty,
               int deletion_penalty, int mismatch_penalty, int match_score,
               struct alignment *a) {
    int l1 = strlen(seq1);
    int l2 = strlen(seq2);
    int **p = malloc(sizeof(int *) * (l1 + 1));
    int **q = malloc(sizeof(int *) * (l1 + 1));
    for (int i = 0; i < l1 + 1; i++) {
        p[i] = calloc(l2 + 1, sizeof(int));
        q[i] = calloc(l2 + 1, sizeof(int));
    }

    int del, ins, match;
    for (int i = 1; i < l1 + 1; i++) {
        for (int j = 1; j < l2 + 1; j++) {
            del = p[i - 1][j] + deletion_penalty;
            ins = p[i][j - 1] + insertion_penalty;
            if (seq1[i - 1] == seq2[j - 1] ||
                    seq1[i - 1] == 'N' ||
                    seq2[j - 1] == 'N' ||
                    seq1[i - 1] == '-' ||
                    seq2[j - 1] == '-') {
                match = p[i - 1][j - 1] + match_score;
            } else {
                match = p[i - 1][j - 1] + mismatch_penalty;
            }
            int max_i = max(del, ins, match);
            if (max_i == DELETION) {
                p[i][j] = del;
            } else if (max_i == INSERTION) {
                p[i][j] = ins;
            } else if (max_i == MATCH) {
                p[i][j] = match;
            }
            q[i][j] = max_i;
        }
    }

    int ret = get_alignment(seq1, seq2, insertion_penalty, deletion_penalty,
                            mismatch_penalty, match_score, p, q, a);

    for (int i = 0; i < l1 + 1; i++) {
        free(p[i]);
        free(q[i]);
    }
    free(p);
    free(q);

    return ret;
}

static PyObject *dnautils_align(PyObject *self, PyObject *args) {
    PyObject *s1, *s2;
    char *str1, *str2;
    int insertion_penalty = -20;
    int deletion_penalty = -20;
    int mismatch_penalty = -4;
    int match_score = 1;

    if (!PyArg_ParseTuple(args, "SS|iiii", &s1, &s2, &insertion_penalty,
                          &deletion_penalty, &mismatch_penalty,
                          &match_score)) {
        return NULL;
    }

    if ((str1 = PyString_AsString(s1)) == NULL) {
        return NULL;
    }
    if ((str2 = PyString_AsString(s2)) == NULL) {
        return NULL;
    }

    struct alignment a;
    if (align_seqs(str1, str2, insertion_penalty, deletion_penalty,
                   mismatch_penalty, match_score, &a) < 0) {
        PyErr_SetString(DNAUtilError, "Error aligning sequences.");
        return NULL;
    }

    return Py_BuildValue("SSi", PyString_FromString(a.seq1),
                         PyString_FromString(a.seq2), a.score);
}

static PyMethodDef DNAUtilsMethods[] = {
    {"equal", dnautils_equal, METH_VARARGS,
        "Checks if two sequences are equal."},
    {"hamming", dnautils_hamming, METH_VARARGS,
        "Gets the hamming distance between two sequences."},
    {"align", dnautils_align, METH_VARARGS, "Aligns two sequences."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initdnautils(void)
{
    PyObject *module;
    module = Py_InitModule("dnautils", DNAUtilsMethods);
    DNAUtilError = PyErr_NewException("dnautil.error", NULL, NULL);
    Py_INCREF(DNAUtilError);
    PyModule_AddObject(module, "error", DNAUtilError);
}
