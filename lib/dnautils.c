#include <Python.h>

static PyObject *DNAUtilError;

static PyObject *dnautils_equal(PyObject *self, PyObject *args) {
    PyObject *s1, *s2;
    char *str1, *str2;
    int i = 0;

    if (!PyArg_ParseTuple(args, "SS", &s1, &s2)) {
        return NULL;
    }
    if ((str1 = PyString_AsString(s1)) == NULL) {
        return NULL;
    }
    if ((str2 = PyString_AsString(s2)) == NULL) {
        return NULL;
    }

    if (PyString_Size(s1) != PyString_Size(s2)) {
        PyErr_SetString(DNAUtilError, "Sequences have unequal lengths.");
        return NULL;
    }

    for (i = 0; i < PyString_Size(s1); i++) {
        if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N') {
            return Py_BuildValue("O", Py_False);
        }
    }
    return Py_BuildValue("O", Py_True);
}

static PyMethodDef DNAUtilsMethods[] = {
    {"equal", dnautils_equal, METH_VARARGS, "Checks if two sequences are equal."},
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
