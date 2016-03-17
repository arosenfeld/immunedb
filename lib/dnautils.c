#include <string.h>
#include <Python.h>

static PyObject *DNAUtilError;

static int process_input(PyObject *args, char *str1, char *str2,
                          unsigned int *length) {
    PyObject *s1, *s2;
    if (!PyArg_ParseTuple(args, "SS", &s1, &s2)) {
        return -1;
    }

    if (PyString_Size(s1) != PyString_Size(s2)) {
        PyErr_SetString(DNAUtilError, "Sequences have unequal lengths.");
        return -1;
    }

    if ((str1 = PyString_AsString(s1)) == NULL) {
        PyErr_SetString(DNAUtilError, "Could not cast first sequence to string.");
        return -1;
    }
    if ((str2 = PyString_AsString(s2)) == NULL) {
        PyErr_SetString(DNAUtilError, "Could not cast second sequence to string.");
        return -1;
    }

    *length = PyString_Size(s1);

    return 0;
}

static PyObject *dnautils_equal(PyObject *self, PyObject *args) {
    char *str1, *str2;
    unsigned int length;
    if (process_input(args, str1, str2, &length) < 0) {
        return NULL;
    }

    unsigned int i;
    for (i = 0; i < length; i++) {
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N') {
            return Py_BuildValue("O", Py_False);
        }
    }
    return Py_BuildValue("O", Py_True);
}

static PyObject *dnautils_hamming(PyObject *self, PyObject *args) {
    char *str1, *str2;
    unsigned int length;
    if (process_input(args, str1, str2, &length) < 0) {
        return NULL;
    }

    unsigned int i;
    unsigned int distance = 0;
    for (i = 0; i < length; i++) {
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N' &&
                str1[i] != '-' && str2[i] != '-') {
            distance += 1;
        }
    }
    return Py_BuildValue("I", distance);
}

static PyMethodDef DNAUtilsMethods[] = {
    {"equal", dnautils_equal, METH_VARARGS,
        "Checks if two sequences are equal."},
    {"hamming", dnautils_hamming, METH_VARARGS,
        "Calculates the hamming distance between two sequences."},
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
