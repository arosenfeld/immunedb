#include <string.h>
#include <Python.h>

#define BUF_SIZE 2048
#define DELETION 0
#define INSERTION 1
#define MATCH 2

static PyObject *DNAUtilError;

static PyObject*
dnautils_equal(PyObject *self, PyObject *args)
{
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
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N') {
            return Py_BuildValue("O", Py_False);
        }
    }
    return Py_BuildValue("O", Py_True);
}

static PyObject*
dnautils_hamming(PyObject *self, PyObject *args)
{
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

static PyMethodDef DNAUtilsMethods[] = {
    {"equal", dnautils_equal, METH_VARARGS,
        "Checks if two sequences are equal."},
    {"hamming", dnautils_hamming, METH_VARARGS,
        "Gets the hamming distance between two sequences."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initdnautils(void)
{
    PyObject *module;
    module = Py_InitModule("dnautils", DNAUtilsMethods);
    DNAUtilError = PyErr_NewException("dnautil.error", NULL, NULL);
    Py_INCREF(DNAUtilError);
    PyModule_AddObject(module, "error", DNAUtilError);
}
