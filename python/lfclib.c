/**
 * @file lfclib.c
 * @author Pietro Bonfa
 * @date 9 Sep 2016
 * @brief Dipolar tensor calculator, Python extension
 *
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <numpy/arrayobject.h>
#include "dipolartensor.h"
#include "fastincommsum.h"
#include "rotatesum.h"
#include "simplesum.h"
#include "randomsample.h"

/* support numpy 1.6 - this macro got renamed and deprecated at once in 1.7 */
#ifndef NPY_ARRAY_IN_ARRAY
#define NPY_ARRAY_IN_ARRAY NPY_IN_ARRAY
#endif

#if NPY_API_VERSION <= NPY_1_6_API_VERSION
#define PyArray_SHAPE PyArray_DIMS
#endif

static char module_docstring[] = "This module provides two functions: Fields and DipolarTensor.";
static char py_lfclib_fields_docstring[] = "Calculate the Local Field components: dipolar, Lorentz and Contact\n"
        "\n"
        "    This function calculates the magnetic field (in Tesla) at the muon site.\n"
        "    Consider using the wrapper function locfield by importing the LFC module.\n"
        "\n"
        "    Parameters\n"
        "    ----------\n"
        "    calc_type : str \n"
        "        Type of calculation. Can be 's', 'r' or 'i'\n"
        "        's': simple sum of all the magnetic dipoles in the Lorentz sphere\n"
        "        'r': sum of all the magnetic dipoles in the Lorentz sphere rotated n_angles times by 360/n_angles degrees\n"
        "        'i': optimized version for incummensurate helical orders\n"
        "    positions : numpy.ndarray\n"
        "        Atomic positions in fractional coordinates.\n"
        "    FC : numpy.ndarray\n"
        "        Fourier components for the atoms in positions, in Cartesian coordinates\n"
        "    K  : numpy.ndarray\n"
        "        Propagation vector, in fractional coordinates of the reciprocal unit cell.\n"
        "    Phi: numpy.ndarray\n"
        "        Phases for the atoms in positions\n"
        "    Muon : numpy.ndarray\n"
        "        Muon positions, in fractional coordinates.\n"
        "    Supercell : numpy.ndarray (dtype=np.int32)\n"
        "        Number of replica along the a, b, and c lattice vectors.\n"
        "    Cell : numpy.ndarray\n"
        "        Lattice parameters (in cartesian axis):\n"
        "        v1(1)  v1(2)  v1(3)    ... 1st lattice vector\n"
        "        v2(1)  v2(2)  v2(3)    ... 2nd lattice vector\n"
        "        v3(1)  v3(2)  v3(3)    ... 3rd lattice vector\n"
        "    r : float\n"
        "        Lorentz sphere radius\n"
        "    nnn : int\n"
        "        Number of nearest neighbours of the muon to be included in the contact hyperfine field approximation.\n"
        "    rcont: float\n"
        "        Maximum distance from the muon for nearest neighbours of the muon to be included in the contact hyperfine field approximation.\n"
        "    nangles: int, optional\n"
        "        Number of divisions of the full turn in the 'r' and 'i' runs\n"
        "    rot_axis: numpy.ndarray, optional\n"
        "        axis for the rotations when using the 'r' option.\n"
        "\n"
        "    Returns\n"
        "    -------\n"
        "    Fields : list of 3 numpy.ndarray\n"
        "        the list contains (in order): Contact Field, Dipolar Field and Lorentz Field. Cartesian coordinates, units Tesla.\n";



static char py_lfclib_dt_docstring[] = "Dipolar tensor calculation.\n"
                                       "\n"
                                       "    This function calculates the dipolar tensor at the muon site.\n"
                                       "    Consider using the wrapper function dipten by importing the LFC module.\n"
                                       "\n"
                                       "    Parameters\n"
                                       "    ----------\n"
                                       "    positions : numpy.ndarray\n"
                                       "        Atomic positions in fractional coordinates.\n"
                                       "    Muon : numpy.ndarray\n"
                                       "        Muon position, in fractional coordinates.\n"
                                       "    Supercell : numpy.ndarray (dtype=np.int32)\n"
                                       "        Number of replica along the a, b, and c lattice vectors.\n"
                                       "    Cell : numpy.ndarray\n"
                                       "        Lattice parameters (in cartesian axis):\n"
                                       "        v1(1)  v1(2)  v1(3)    ... 1st lattice vector\n"
                                       "        v2(1)  v2(2)  v2(3)    ... 2nd lattice vector\n"
                                       "        v3(1)  v3(2)  v3(3)    ... 3rd lattice vector\n"
                                       "    r : float\n"
                                       "        Lorentz sphere radius\n"
                                       "\n"
                                       "    Returns\n"
                                       "    -------\n"
                                       "    Fields : numpy.ndarray\n"
                                       "        The array contains the dipolar tensor.\n";

typedef int bool;
enum { false, true };

typedef struct
{
    PyArrayObject *npy_obj;
    npy_int nd;
    npy_intp *shape;
    int typenum;
    void *p;
} npy_data_obj;

static void init_npy_data(npy_data_obj * data) {
    data->nd = 0;
    data->npy_obj = NULL;
    data->shape = NULL;
    data->p=NULL;
    data->typenum=-1;
}

static bool fill_npy_data(PyObject *o, int typenum, int min_size, int max_size, npy_data_obj * data, bool writable) {

    if (writable) {
        data->npy_obj = (PyArrayObject *) PyArray_FROMANY(o, typenum, min_size, max_size,
                        NPY_ARRAY_INOUT_ARRAY);
    } else {
        data->npy_obj = (PyArrayObject *) PyArray_FROMANY(o, typenum, min_size, max_size,
                        NPY_ARRAY_IN_ARRAY);
    }
    if (data->npy_obj == NULL) return false;

    data->typenum = typenum;
    data->nd = PyArray_NDIM(data->npy_obj);
    data->shape = PyArray_SHAPE(data->npy_obj);

    data->p = PyArray_DATA(data->npy_obj);

    return true;
}

static double *to_dbl(npy_data_obj data, bool * err) {
    if (data.typenum != NPY_DOUBLE && data.typenum != NPY_COMPLEX128) {
        *err=true;
        printf ("Error getting type");
    }
    return (double *) data.p;
}
static int *to_int(npy_data_obj data, bool * err) {
    if (data.typenum != NPY_INT32)  {
        *err=true;
        printf ("Error getting type");
    }
    return (int *) data.p;
}

static bool free_npy_obj(npy_data_obj data, bool err) {
#if NPY_API_VERSION >= 0x0000000c
    if (data.npy_obj) {
        if (PyArray_CHKFLAGS(data.npy_obj, NPY_ARRAY_WRITEBACKIFCOPY)) {
            if (!err) {
                PyArray_ResolveWritebackIfCopy(data.npy_obj);
            } else {
                PyArray_DiscardWritebackIfCopy(data.npy_obj);
            }
        }
    }
#endif
    Py_XDECREF(data.npy_obj);
    return err;
}

static bool free_many_npy_obj(int count, bool err, ...) {
    va_list list;
    int j = 0;
    va_start(list, err);
    for(j=0; j<count; j++)
    {
        err |= free_npy_obj(va_arg(list, npy_data_obj), err);
    }
    va_end(list);

    return err;
}

static PyObject * py_lfclib_fields(PyObject *self, PyObject *args, PyObject *kwargs) {
    /* input variables */
    char* calc_type = NULL;
    unsigned int nnn=0;
    unsigned int nangles=0;

    double r=0.0;
    double rcont=0.0;
    double ratms=-1.0;

    PyObject *opositions, *oFC, *oK, *oPhi;
    PyObject *omu, *osupercell, *ocell;
    PyObject *orot_axis = NULL;
    PyObject *oconstraints = NULL;
    PyObject *oconstraint_group = NULL;

    npy_data_obj positions, FC, K, Phi, mu, supercell, cell, rot_axis,
                 constraints, constraint_group;

    /* local variables */
    int num_atoms=0;
    int num_muons=0;
    int num_constraints=0;
    int icalc_type=0;
    bool err=false;
    int nd = 1;

    /* return variables */
    double * dip=NULL, *cont=NULL, *lor=NULL;
    npy_intp *out_dim=NULL;
    PyArrayObject *odip, *ocont, *olor;

    static char *kwlist[] = {"ctype", "positions", "FC", "K",
                             "phi", "muon", "supercell", "cell",
                             "r", "nnn", "rcont", "nangles", "axis",
                             "constraints", "constraint_group", "dist_from_atoms", NULL
                            };

    /* put arguments into variables */
    // if (!PyArg_ParseTuple(args, "sOOOOOOOdId|IO", &calc_type,
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOOOOOOOdId|IOOOd", kwlist,
                                     &calc_type,
                                     &opositions, &oFC, &oK, &oPhi,
                                     &omu, &osupercell, &ocell,
                                     &r,&nnn,&rcont,
                                     &nangles, &orot_axis,
                                     &oconstraints, &oconstraint_group,
                                     &ratms)) {
        return NULL;
    }

    /* turn inputs into numpy array types */
    init_npy_data(&positions);
    init_npy_data(&FC);
    init_npy_data(&K);
    init_npy_data(&Phi);
    init_npy_data(&mu);
    init_npy_data(&supercell);
    init_npy_data(&cell);
    init_npy_data(&rot_axis);
    init_npy_data(&constraints);
    init_npy_data(&constraint_group);

    err = false;
    err |= !fill_npy_data(opositions, NPY_DOUBLE, 2, 2, &positions, false);
    err |= !fill_npy_data(oFC, NPY_COMPLEX128, 2, 2, &FC, false);
    err |= !fill_npy_data(oK, NPY_DOUBLE, 1, 1, &K, false);
    err |= !fill_npy_data(oPhi, NPY_DOUBLE, 1, 1, &Phi, false);
    err |= !fill_npy_data(omu, NPY_DOUBLE, 1, 2, &mu, true);
    err |= !fill_npy_data(osupercell, NPY_INT32, 1, 1, &supercell, false);
    err |= !fill_npy_data(ocell, NPY_DOUBLE, 2, 2, &cell, false);


    /* Validate data */
    if (err) {
        PyErr_Format(PyExc_RuntimeError,
                     "Error parsing input data. Debug info: %s %p %p %p %p %p %p %p",
                     calc_type,
                     positions.npy_obj,
                     FC.npy_obj,
                     K.npy_obj,
                     Phi.npy_obj,
                     mu.npy_obj,
                     supercell.npy_obj,
                     cell.npy_obj);
        free_many_npy_obj(7, err, positions, FC, K, Phi, mu, supercell, cell);
        return NULL;
    }

    /* Select calculation type */
    if (strcmp(calc_type, "s")==0 || strcmp(calc_type, "sum")==0) {
        icalc_type = 1;
        nangles = 1;
    }  else if (strcmp(calc_type, "rnd")==0 || strcmp(calc_type, "random")==0) {
        icalc_type = 4;
        nangles = 1;
    }  else if (strcmp(calc_type, "r")==0 || strcmp(calc_type, "rotate")==0) {
        icalc_type = 2;
    }  else if (strcmp(calc_type, "i")==0 || strcmp(calc_type, "incommensurate")==0) {
        icalc_type = 3;
    }  else  {
        free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
        PyErr_Format(PyExc_ValueError,
                     "Valid calculations are 's', 'r', 'i', 'rnd'.  Unknown value %s", calc_type);
        return NULL;
    }

    /* General validity checks */

    /* Check optional arguments are present for calculation 'i' or 'r' */
    if ( icalc_type == 2 || icalc_type == 3 ) {
        if(nangles==0) {
            free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
            PyErr_Format(PyExc_ValueError,
                         "Number of angles required!");
            return NULL;
        }
    }


    /* Check optional arguments are present for calculation 'r' */
    if (icalc_type == 2) {
        if (!orot_axis) {
            free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
            PyErr_Format(PyExc_ValueError,
                         "Axis for rotation required!");
            return NULL;
        }
    }

    /* Check optional arguments are present for calculation 'rnd' */
    if (icalc_type == 4) {
        if ((!oconstraints) || (!oconstraint_group)) {
            free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
            PyErr_Format(PyExc_ValueError,
                         "Constraints and constraint_group for Random Sample required!");
            return NULL;
        }
    }

    if (nnn > 200) {
        free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
        PyErr_Format(PyExc_RuntimeError,
                     "Error, number of nearest neighbours exceedingly large.");
        return NULL;
    }
    if (nangles < 1) {
        free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
        PyErr_Format(PyExc_RuntimeError,
                     "Error, number of angles must be greater than 0.");
        return NULL;
    }


    /* Additional validation */
    if (r<=0.0) {
        PyErr_Warn( PyExc_UserWarning, "Radius for sums is <= 0!");
    }
    if (rcont<=0.0) {
        PyErr_Warn( PyExc_UserWarning, "Radius for contact hyperfine couplig is <= 0!");
    }


    /* Parse ptional argument for 'r' */
    if (icalc_type == 2) {
        err |= !fill_npy_data(orot_axis, NPY_DOUBLE, 1, 1, &rot_axis, false);
        if (err) {
            free_many_npy_obj(7, true, positions, FC, K, Phi, mu, supercell, cell);
            PyErr_Format(PyExc_ValueError,
                         "Axis for rotation required but not parsed!");
            return NULL;
        }
    }

    /* The first condition is impossible since it is implied by the
     *  PyArray_FROMANY function above. That function guaranties that the
     *  comparison of the shape is fine.
     */
    if (positions.nd != FC.nd ||
            positions.shape[0]!=FC.shape[0] ||
            positions.shape[1]!=FC.shape[1]) {
        free_many_npy_obj(8, true, positions, FC, K, Phi, mu, supercell, cell, rot_axis);
        PyErr_SetString(PyExc_RuntimeError, "positions and FC arrays must have "
                        "same shape.");
        return NULL;
    }

    if (Phi.shape[0] != positions.shape[0]) {
        free_many_npy_obj(8, true, positions, FC, K, Phi, mu, supercell, cell, rot_axis);
        PyErr_SetString(PyExc_RuntimeError, "positions and Phi arrays must have "
                        "same shape[0].");
        return NULL;
    }

    /* Get dimension of arrays for muons and atoms */
    num_atoms = positions.shape[0];
    num_muons = (mu.nd > 1) ? mu.shape[0] : 1;

    /* Parse optional argument for 'rnd' */
    if (icalc_type == 4) {

        err |= !fill_npy_data(oconstraints, NPY_DOUBLE, 2, 2, &constraints, false);
        err |= !fill_npy_data(oconstraint_group, NPY_INT32, 2, 2, &constraint_group, false);

        if (err) {
            free_many_npy_obj(10, err, positions, FC, K, Phi, mu, supercell, cell, rot_axis, constraints, constraint_group);
            PyErr_Format(PyExc_ValueError,
                         "Constrains and constraint group for random required but not parsed!");
            return NULL;
        } else {
            num_constraints = constraints.shape[0];
            if (constraint_group.shape[0] != num_constraints ||
                    constraint_group.shape[1] != num_atoms) {
                PyErr_Format(PyExc_ValueError, "Constraint group has invalid shape: %d %d, should be %d %d!",
                             constraint_group.shape[0], constraint_group.shape[1],
                             num_constraints, num_atoms);
                free_many_npy_obj(10, true, positions, FC, K, Phi, mu, supercell, cell, rot_axis, constraints, constraint_group);
                return NULL;
            }
        }
    }

    /* allocate output arrays */
    nd = 1;
    if (num_muons > 1) {
        nd = nd + 1;
    }
    if (icalc_type >= 2) {
        nd = nd + 1;
        /* finish creating return variables */
        out_dim = (npy_intp *) malloc(nd * sizeof(npy_intp));
        if (num_muons > 1) out_dim[0] = (npy_intp) num_muons;
        out_dim[nd-2] = (npy_intp) nangles;
        out_dim[nd-1] = (npy_intp) 3;
    } else {
        /* Simple sum case*/
        out_dim = (npy_intp *) malloc(nd * sizeof(npy_intp));
        if (num_muons > 1) out_dim[0] = (npy_intp) num_muons;
        out_dim[nd-1] = (npy_intp) 3;
    }

    odip = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
    ocont = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
    olor = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
    if (out_dim != NULL) {
        free(out_dim);
        out_dim = NULL;
    }

    if (!odip || !ocont || !olor) {

        Py_XDECREF(odip);
        Py_XDECREF(ocont);
        Py_XDECREF(olor);
        free_many_npy_obj(10, false, positions, FC, K, Phi, mu, supercell, cell, rot_axis, constraints, constraint_group);

        PyErr_SetString(PyExc_MemoryError, "Cannot create output arrays.");
        return NULL;
    }


    dip = (double *) PyArray_DATA(odip);
    cont = (double *) PyArray_DATA(ocont);
    lor = (double *) PyArray_DATA(olor);
    err = false;
    /* long computation starts here. No python object is touched so free thread execution */
    Py_BEGIN_ALLOW_THREADS
    switch (icalc_type)
    {
    case 1:
        SimpleSum(to_dbl(positions, &err),
                  to_dbl(FC, &err),
                  to_dbl(K, &err),
                  to_dbl(Phi, &err),
                  to_dbl(mu, &err),
                  to_int(supercell, &err),
                  to_dbl(cell, &err),
                  r, nnn, rcont, ratms, num_atoms,num_muons,cont,dip,lor);
        break;
    case 2:
        RotataSum(to_dbl(positions, &err),
                  to_dbl(FC, &err),
                  to_dbl(K, &err),
                  to_dbl(Phi, &err),
                  to_dbl(mu, &err),
                  to_int(supercell, &err),
                  to_dbl(cell, &err),
                  r, nnn,rcont,num_atoms,num_muons,to_dbl(rot_axis, &err),nangles,cont,dip,lor);
        break;
    case 3:
        FastIncommSum(to_dbl(positions, &err),
                      to_dbl(FC, &err),
                      to_dbl(K, &err),
                      to_dbl(Phi, &err),
                      to_dbl(mu, &err),
                      to_int(supercell, &err),
                      to_dbl(cell, &err),
                      r, nnn,rcont,num_atoms,num_muons,nangles,cont,dip,lor);
        break;
    case 4:
        RandomSample(to_dbl(positions, &err),
                     to_dbl(FC, &err),
                     to_dbl(K, &err),
                     to_dbl(Phi, &err),
                     to_int(supercell, &err),
                     to_dbl(cell, &err),
                     r, nnn,rcont, num_atoms, num_muons, to_dbl(constraints, &err), to_int(constraint_group, &err),
                     num_constraints, to_dbl(mu, &err), cont, dip, lor);
        break;
    }
    Py_END_ALLOW_THREADS

    free_many_npy_obj(10, false, positions, FC, K, Phi, mu, supercell, cell, rot_axis, constraints, constraint_group);

    return Py_BuildValue("NNN", ocont, odip, olor);
}


static PyObject * py_lfclib_dt(PyObject *self, PyObject *args) {

    double r=0.0;
    PyObject *opositions, *omu, *osupercell, *ocell;
    PyArrayObject *positions,  *mu, *supercell, *cell, *odt;

    int num_atoms=0;
    int * in_supercell;
    npy_intp * pShape, * out_dim;
    int nd = 2;


    /* put arguments into variables */
    if (!PyArg_ParseTuple(args, "OOOOd",
                          &opositions, &omu, &osupercell, &ocell,&r))
    {
        return NULL;
    }

    /* turn inputs into numpy array types */
    positions = (PyArrayObject *) PyArray_FROMANY(opositions, NPY_DOUBLE, 2, 2,
                NPY_ARRAY_IN_ARRAY);
    mu = (PyArrayObject *) PyArray_FROMANY(omu, NPY_DOUBLE, 1, 1,
                                           NPY_ARRAY_IN_ARRAY);

    supercell = (PyArrayObject *) PyArray_FROMANY(osupercell, NPY_INT32,
                1, 1, NPY_ARRAY_IN_ARRAY);
    cell = (PyArrayObject *) PyArray_FROMANY(ocell, NPY_DOUBLE, 2, 2,
            NPY_ARRAY_IN_ARRAY);

    /* Validate data */
    if ( !positions || !mu || !supercell || !cell) {
        Py_XDECREF(positions);
        Py_XDECREF(mu);
        Py_XDECREF(supercell);
        Py_XDECREF(cell);
        PyErr_Format(PyExc_RuntimeError,
                     "Error parsing numpy arrays.");
        return NULL;
    }

    /* General validity checks */

    pShape = PyArray_SHAPE(positions);

    num_atoms = pShape[0];

    if (num_atoms == 0) {
        PyErr_Format(PyExc_RuntimeError,
                     "No magnetic atoms specified.");
        return NULL;
    }


    /* to avoid problems with int types */
    in_supercell = (int *) malloc(3*sizeof(int));
    in_supercell[0] = *(npy_int64 *)PyArray_GETPTR1(supercell, 0);
    in_supercell[1] = *(npy_int64 *)PyArray_GETPTR1(supercell, 1);
    in_supercell[2] = *(npy_int64 *)PyArray_GETPTR1(supercell, 2);


    nd = 2;
    out_dim= (npy_intp *)malloc(2 * sizeof(npy_intp));
    out_dim[0] = (npy_intp) 3;
    out_dim[1] = (npy_intp) 3;

    odt = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);

    free(out_dim);

    if (!odt) {
        Py_XDECREF(positions);
        Py_XDECREF(mu);
        Py_XDECREF(supercell);
        Py_XDECREF(cell);
        if (in_supercell != NULL)
            free(in_supercell);
        PyErr_SetString(PyExc_MemoryError, "Cannot create output array.");
        return NULL;
    }

    /* long computation starts here. No python object is touched so free thread execution */
    Py_BEGIN_ALLOW_THREADS
    DipolarTensor( (double *) PyArray_DATA(positions),
                   (double *) PyArray_DATA(mu),
                   in_supercell,
                   (double *) PyArray_DATA(cell),
                   r, num_atoms,
                   (double *) PyArray_DATA(odt));

    Py_END_ALLOW_THREADS


    Py_DECREF(positions);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);

    free(in_supercell);
    return Py_BuildValue("N", odt);

}

static PyMethodDef lfclib_methods[] =
{
    {"Fields", (PyCFunction)py_lfclib_fields, METH_VARARGS | METH_KEYWORDS, py_lfclib_fields_docstring},
    {"DipolarTensor", (PyCFunction)py_lfclib_dt, METH_VARARGS, py_lfclib_dt_docstring},
    {NULL}  /* sentinel */
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "lfc",            /* m_name */
    module_docstring,    /* m_doc */
    -1,                  /* m_size */
    lfclib_methods,      /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};
#endif


PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_lfclib(void)
#else
initlfclib(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    (void) Py_InitModule("lfclib", lfclib_methods);
#endif

    import_array(); /* Must be present for NumPy */

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
