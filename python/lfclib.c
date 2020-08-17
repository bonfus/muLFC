/**
 * @file lfclib.c
 * @author Pietro Bonfa
 * @date 9 Sep 2016
 * @brief Dipolar tensor calculator, Python extension
 *     
 */

#define MULFC_VERSION "0.0.3"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <numpy/arrayobject.h>
#include "dipolartensor.h"
#include "fastincommsum.h"
#include "rotatesum.h"
#include "simplesum.h"

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



static PyObject * py_lfclib_fields(PyObject *self, PyObject *args) {
  /* input variables */
  char* calc_type = NULL;
  unsigned int nnn=0;
  unsigned int nangles=0;
  
  double r=0.0;
  double rcont=0.0;
    
  PyObject *opositions, *oFC, *oK, *oPhi;
  PyObject *omu, *osupercell, *ocell;
  PyObject *orot_axis = NULL;
  
  PyArrayObject *positions, *FC, *K, *Phi;
  PyArrayObject *mu, *supercell, *cell;
  PyArrayObject *rot_axis = NULL;

  /* local variables */
  int num_atoms=0;
  int icalc_type=0;
  
  
  double * in_positions=NULL;
  double * in_fc=NULL;
  double * in_K=NULL;
  double * in_phi=NULL;
  double * in_muonpos=NULL;
  int * in_supercell=NULL;
  double * in_cell=NULL;
  double * in_axis = NULL;

  /* return variables */
  double * dip=NULL, *cont=NULL, *lor=NULL;
  npy_intp *out_dim=NULL;
  PyArrayObject *odip, *ocont, *olor;
  /* General validity checks */

  npy_intp pndims=0; 
  npy_intp fcndims=0; 

  npy_intp * pShape=NULL; 
  npy_intp * fcShape=NULL; 

  npy_intp * phiSpahe=NULL; 

  int i;
  int nd = 1;
  npy_cdouble v;

  /* put arguments into variables */
  if (!PyArg_ParseTuple(args, "sOOOOOOOdId|IO", &calc_type, 
                                        &opositions, &oFC, &oK, &oPhi,
                                        &omu, &osupercell, &ocell,
                                        &r,&nnn,&rcont,
                                        &nangles,&orot_axis)) {
    return NULL;
  }
  
  /* turn inputs into numpy array types */
  positions = (PyArrayObject *) PyArray_FROMANY(opositions, NPY_DOUBLE, 2, 2,
                                              NPY_ARRAY_IN_ARRAY);
  FC = (PyArrayObject *) PyArray_FROMANY(oFC, NPY_COMPLEX128, 2, 2,
                                              NPY_ARRAY_IN_ARRAY);
  K = (PyArrayObject *) PyArray_FROMANY(oK, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
  Phi = (PyArrayObject *) PyArray_FROMANY(oPhi, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
  mu = (PyArrayObject *) PyArray_FROMANY(omu, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
                                              
  supercell = (PyArrayObject *) PyArray_FROMANY(osupercell, NPY_INT32,
                                                   1, 1, NPY_ARRAY_IN_ARRAY);
  cell = (PyArrayObject *) PyArray_FROMANY(ocell, NPY_DOUBLE, 2, 2,
                                             NPY_ARRAY_IN_ARRAY);

  /* Validate data */
  if (!calc_type || !positions || !FC || !K || !Phi || !mu || !supercell || !cell) {
    Py_XDECREF(positions);
    Py_XDECREF(FC);
    Py_XDECREF(K);
    Py_XDECREF(Phi);
    Py_XDECREF(mu);
    Py_XDECREF(supercell);
    Py_XDECREF(cell);
    PyErr_Format(PyExc_RuntimeError,
                    "Error parsing numpy arrays.");                 
    return NULL;
  }


  
  
  
  /* Select calculation type */
  if (strcmp(calc_type, "s")==0 || strcmp(calc_type, "sum")==0) {
    icalc_type = 1;
    nangles = 1;
  }  else if (strcmp(calc_type, "r")==0 || strcmp(calc_type, "rotate")==0) {
    icalc_type = 2;
  }  else if (strcmp(calc_type, "i")==0 || strcmp(calc_type, "incommensurate")==0) {
    icalc_type = 3;
  }  else  {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);    
    PyErr_Format(PyExc_ValueError,
                   "Valid calculations are 's', 'r', 'i'.  Unknown value %s", calc_type);
    return NULL;
  }

  /* Check optional arguments are present for calculation 'i' or 'r' */
  if (icalc_type >= 2) {
      if(nangles==0) {
        Py_DECREF(positions);
        Py_DECREF(FC);
        Py_DECREF(K);
        Py_DECREF(Phi);
        Py_DECREF(mu);
        Py_DECREF(supercell);
        Py_DECREF(cell);          
        PyErr_Format(PyExc_ValueError,
                        "Number of angles required!");
        return NULL;  
      }
  }

  
  /* Check optional arguments are present for calculation 'r' */
  if (icalc_type == 2) {
    if (!orot_axis) {
      Py_DECREF(positions);
      Py_DECREF(FC);
      Py_DECREF(K);
      Py_DECREF(Phi);
      Py_DECREF(mu);
      Py_DECREF(supercell);
      Py_DECREF(cell);          
      PyErr_Format(PyExc_ValueError,
                    "Axis for rotation required!");
      return NULL;        
    }
  }

  if (nnn > 200) {
    PyErr_Format(PyExc_RuntimeError,
                    "Error, number of nearest neighbours exceedingly large.");
    return NULL;
  }
  if (nangles < 1) {
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
    rot_axis = (PyArrayObject *) PyArray_FROMANY(orot_axis, NPY_DOUBLE, 1, 1,
                                             NPY_ARRAY_IN_ARRAY);
    if (!rot_axis) {
      Py_DECREF(positions);
      Py_DECREF(FC);
      Py_DECREF(K);
      Py_DECREF(Phi);
      Py_DECREF(mu);
      Py_DECREF(supercell);
      Py_DECREF(cell);          
      PyErr_Format(PyExc_ValueError,
                    "Axis for rotation required but not parsed!");
      return NULL;        
    }
  }
  
  pndims = PyArray_NDIM(positions);
  fcndims = PyArray_NDIM(FC);
  pShape = PyArray_SHAPE(positions);
  fcShape = PyArray_SHAPE(FC);
  /* The first condition is impossible since it is implied by the 
   *  PyArray_FROMANY function above. That function guaranties that the 
   *  comparison of the shape is fine.
   */ 
  
  if (pndims != fcndims || 
        pShape[0]!=fcShape[0] ||
        pShape[1]!=fcShape[1]) {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   /* Null in case it is optional */
	
    PyErr_SetString(PyExc_RuntimeError, "positions and FC arrays must have "
		    "same shape.");
    return NULL;
  }
  phiSpahe = PyArray_SHAPE(Phi);
  if (phiSpahe[0] != pShape[0]) {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   /* Null in case it is optional */
	
	PyErr_SetString(PyExc_RuntimeError, "positions and Phi arrays must have "
		    "same shape[0].");
    return NULL;
  }
  
  num_atoms = pShape[0]; 
  

  /* allocate variables needed for simulations */
  
  in_positions = (double *) malloc(3*num_atoms*sizeof(double));
  in_fc = (double *) malloc(6*num_atoms*sizeof(double));
  in_K = (double *) malloc(3*sizeof(double));
  in_phi = (double *) malloc(num_atoms*sizeof(double));
  in_muonpos = (double *) malloc(3*sizeof(double));
  in_supercell = malloc(3*sizeof(int));
  in_cell = (double *) malloc(9*sizeof(double));
  
  if (!in_positions || !in_fc || !in_K || !in_phi || !in_muonpos ||
      !in_supercell || !in_cell) {
    PyErr_SetString(PyExc_MemoryError, "Malloc failed.");
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   /* Null in case it is optional */
    if (in_positions != NULL){
        free(in_positions);
    }    
    if (in_fc != NULL){
        free(in_fc);
    }    
    if (in_K != NULL){
        free(in_K);
    }    
    if (in_phi != NULL){
        free(in_phi);
    }    
    if (in_muonpos != NULL){
        free(in_muonpos);
    }    
    if (in_supercell != NULL){
        free(in_supercell);
    }    
    if (in_cell != NULL){
        free(in_cell);
    }
    return NULL;          
  }
  
  for (i=0; i< num_atoms; i++){
    in_positions[3*i+0] = *(npy_float64 *) PyArray_GETPTR2(positions, i,0);
    in_positions[3*i+1] = *(npy_float64 *) PyArray_GETPTR2(positions, i,1);
    in_positions[3*i+2] = *(npy_float64 *) PyArray_GETPTR2(positions, i,2);
  }
  
  

  for (i=0; i< num_atoms; i++){
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,0);
    in_fc[6*i+0] = v.real;
    in_fc[6*i+1] = v.imag;
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,1);
    in_fc[6*i+2] = v.real;
    in_fc[6*i+3] = v.imag;    
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,2);
    in_fc[6*i+4] = v.real;
    in_fc[6*i+5] = v.imag;
  }

  
  in_K[0] = *(npy_float64 *) PyArray_GETPTR1(K,0);
  in_K[1] = *(npy_float64 *) PyArray_GETPTR1(K,1);
  in_K[2] = *(npy_float64 *) PyArray_GETPTR1(K,2);
  
  
  for ( i=0; i< num_atoms; i++){
    in_phi[i] = *(npy_float64 *) PyArray_GETPTR1(Phi,i);
  }
  
  
  in_muonpos[0] = *(npy_float64 *) PyArray_GETPTR1(mu,0);
  in_muonpos[1] = *(npy_float64 *) PyArray_GETPTR1(mu,1);
  in_muonpos[2] = *(npy_float64 *) PyArray_GETPTR1(mu,2);
  
  
  in_supercell[0] = *(npy_int64 *)PyArray_GETPTR1(supercell, 0);
  in_supercell[1] = *(npy_int64 *)PyArray_GETPTR1(supercell, 1);
  in_supercell[2] = *(npy_int64 *)PyArray_GETPTR1(supercell, 2);
  
  
  in_cell[0] = *(npy_float64 *) PyArray_GETPTR2(cell,0,0);
  in_cell[1] = *(npy_float64 *) PyArray_GETPTR2(cell,0,1);
  in_cell[2] = *(npy_float64 *) PyArray_GETPTR2(cell,0,2);
  in_cell[3] = *(npy_float64 *) PyArray_GETPTR2(cell,1,0);
  in_cell[4] = *(npy_float64 *) PyArray_GETPTR2(cell,1,1);
  in_cell[5] = *(npy_float64 *) PyArray_GETPTR2(cell,1,2);
  in_cell[6] = *(npy_float64 *) PyArray_GETPTR2(cell,2,0);
  in_cell[7] = *(npy_float64 *) PyArray_GETPTR2(cell,2,1);
  in_cell[8] = *(npy_float64 *) PyArray_GETPTR2(cell,2,2);
  
  
  if (icalc_type == 2) {
    in_axis = (double *) malloc(3*sizeof(double));
    in_axis[0] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 0);
    in_axis[1] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 1);
    in_axis[2] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 2);
  }
  
  

  /* allocate output arrays */
  nd = 1;
  if (icalc_type >= 2) {
    /* finish creating return variables */
    out_dim = (npy_intp *) malloc(2 * sizeof(npy_intp));
    out_dim[0] = (npy_intp) nangles;
    out_dim[1] = (npy_intp) 3;
    nd = 2;
  } else {
    out_dim = (npy_intp *) malloc(1 * sizeof(npy_intp));
    out_dim[0] = (npy_intp) 3;
    nd = 1;
  }
  
  odip = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  ocont = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  olor = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  if (out_dim != NULL) {
    free(out_dim); out_dim = NULL;
  }

  if (!odip || !ocont || !olor) {

    Py_XDECREF(odip);   
    Py_XDECREF(ocont);  
    Py_XDECREF(olor);
    if (in_positions != NULL)
      free(in_positions);
    if (in_fc != NULL)
      free(in_fc);
    if (in_K != NULL)
      free(in_K);
    if (in_phi != NULL)
      free(in_phi);
    if (in_muonpos != NULL)
      free(in_muonpos);
    if (in_supercell != NULL)
      free(in_supercell);
    if (in_cell != NULL)
      free(in_cell);          
    if (in_axis != NULL)
      free(in_axis);          
    PyErr_SetString(PyExc_MemoryError, "Cannot create output arrays.");
    return NULL;
  }


  dip = (double *) PyArray_DATA(odip);
  cont = (double *) PyArray_DATA(ocont);
  lor = (double *) PyArray_DATA(olor);
  
  /* long computation starts here. No python object is touched so free thread execution */
  Py_BEGIN_ALLOW_THREADS
  switch (icalc_type)
  {
    case 1:
      SimpleSum(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,cont,dip,lor);
      break;
    case 2:
      RotataSum(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,in_axis,nangles,cont,dip,lor);
      break;
    case 3:
      FastIncommSum(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,nangles,cont,dip,lor);
    
  }
  Py_END_ALLOW_THREADS

  if (in_positions != NULL){
      free(in_positions);
      in_positions = NULL;
  }
  if (in_fc != NULL)
    free(in_fc);
  if (in_K != NULL)
    free(in_K);
  if (in_phi != NULL)
    free(in_phi);
  if (in_muonpos != NULL)
    free(in_muonpos);
  if (in_supercell != NULL)
    free(in_supercell);
  if (in_cell != NULL)
    free(in_cell);          
  if (in_axis != NULL)            /* Null in case it is optional */
    free(in_axis);  
  Py_DECREF(positions);
  Py_DECREF(FC);
  Py_DECREF(K);
  Py_DECREF(Phi);
  Py_DECREF(mu);
  Py_DECREF(supercell);
  Py_DECREF(cell);
  Py_XDECREF(rot_axis);   /* Null in case it is optional */

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
  in_supercell = malloc(3*sizeof(int));
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
  {"Fields", (PyCFunction)py_lfclib_fields, METH_VARARGS, py_lfclib_fields_docstring},
  {"DipolarTensor", (PyCFunction)py_lfclib_dt, METH_VARARGS, py_lfclib_dt_docstring},
  {NULL}  /* sentinel */
};


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


PyMODINIT_FUNC PyInit_lfclib(void)
{
  PyObject *module = PyModule_Create(&moduledef);
  PyModule_AddStringConstant(module, "__version__", MULFC_VERSION);

  import_array(); /* Must be present for NumPy */

  return module;
}
