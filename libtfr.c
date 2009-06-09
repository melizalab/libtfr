/*
 * libtfr.c
 *
 * A python extension module interface to libtfrspec
 * Data are stored in numpy arrays.
 *
 */

#include <Python.h>

#include "numpy/arrayobject.h"
#include "tfr.h"
#include "version.h"

static double*
coerce_ndarray_double(PyArrayObject *in, PyArrayObject **out)
{
	int pcm_type = PyArray_TYPE(in);
	if (pcm_type==NPY_DOUBLE)
		return ((double *) PyArray_DATA(in));
	else if (PyArray_CanCastSafely(pcm_type, NPY_DOUBLE)) {
		*out = (PyArrayObject*)PyArray_Cast(in, NPY_DOUBLE);
		return ((double *) PyArray_DATA(*out));
	}
	else
		return NULL;
}

static void
hc2cmplx(const mfft *mtm, npy_cdouble *out)
{
	int nfft = mtm->nfft;
	int ntapers = mtm->ntapers;
	int real_count = nfft / 2 + 1;
	int imag_count = (nfft+1) / 2;  // not actually the count but the last index
	int t,n;
	double x;

	for (t = 0; t < ntapers; t++) {
		for (n = 0; n < real_count; n++)
			out[t*nfft+n].real = out[t*nfft+(nfft-n)].real = mtm->buf[t*nfft+n];
		for (n = 1; n < imag_count; n++) {
			x = mtm->buf[t*nfft+(nfft-n)];
			out[t*nfft+n].imag = x;
			out[t*nfft+(nfft-n)].imag = -x;
		}
	}
}


/* methods */
static PyObject*
libtfr_tfr_spec(PyObject *self, PyObject *args)
{
	/* arguments */
	PyObject *o = NULL;
	PyArrayObject *signal = NULL;
	PyArrayObject *signal_cast = NULL;
	int N;
	int step;
	int Np;
	int K = 6;
	double tm = 6.0;
	double flock = 0.01;
	int tlock = 5;

	/* output data */
	npy_intp out_shape[2];
	PyArrayObject *outdata = NULL;
	double *spec;

	/* internal stuff */
	mfft *mtmh;
	double *samples;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "Oiii|iddi", &o, &N, &step, &Np, &K, &tm, &flock, &tlock))
		return NULL;
	signal = (PyArrayObject*) PyArray_FromAny(o, NULL, 1, 1, NPY_CONTIGUOUS, NULL);
	if (signal==NULL) {
		PyErr_SetString(PyExc_TypeError, "Input signal must be an ndarray");
		return NULL;
	}

	/* allocate output array */
	out_shape[0] = N/2+1;
	out_shape[1] = PyArray_SIZE(signal) / step;
	outdata  = (PyArrayObject*) PyArray_ZEROS(2,out_shape,NPY_DOUBLE,1); // last arg give fortran-order
	spec = (double*) PyArray_DATA(outdata);

	/* coerce data to proper type */
	samples = coerce_ndarray_double(signal, &signal_cast);
	if (samples==NULL) {
		PyErr_SetString(PyExc_TypeError, "Unable to cast signal to supported data type");
		goto fail;
	}

	/* do the transform */
	mtmh = mtm_init_herm(N, Np, K, tm);
	tfr_spec(mtmh, spec, samples, PyArray_SIZE(signal),
		 -1, step, flock, tlock);
	mtm_destroy(mtmh);

	Py_DECREF(signal);
	Py_XDECREF(signal_cast);
	return PyArray_Return(outdata);
fail:
	Py_XDECREF(signal_cast);
	Py_XDECREF(signal);
	Py_XDECREF(outdata);
	return NULL;
}


static PyObject*
libtfr_hermf(PyObject *self, PyObject *args)
{
	/* arguments */
	int N,M;
	double tm;

	/* output data */
	npy_intp tapers_shape[2];
	PyArrayObject *h = NULL;
	PyArrayObject *Dh = NULL;
	PyArrayObject *Th = NULL;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "iid", &N, &M, &tm))
		return NULL;

	/* allocate output arrays */
	tapers_shape[0] = N;
	tapers_shape[1] = M;
	h  = (PyArrayObject*) PyArray_ZEROS(2,tapers_shape,NPY_DOUBLE,1);
	Dh  = (PyArrayObject*) PyArray_ZEROS(2,tapers_shape,NPY_DOUBLE,1);
	Th  = (PyArrayObject*) PyArray_ZEROS(2,tapers_shape,NPY_DOUBLE,1);

	double *hp = (double*) PyArray_DATA(h);
	double *Dhp = (double*) PyArray_DATA(Dh);
	double *Thp = (double*) PyArray_DATA(Th);

	hermf(N, M, tm, hp, Dhp, Thp);

	return Py_BuildValue("(OOO)", h, Dh, Th);
}


#ifndef NO_LAPACK
static PyObject*
libtfr_mtm_spec(PyObject *self, PyObject *args)
{
	/* arguments */
	PyObject *o = NULL;
	PyArrayObject *signal = NULL;
	PyArrayObject *signal_cast = NULL;
	int N;
	int step;
	double NW;
	int K = -1;
	int adapt = 1;

	/* output data */
	npy_intp out_shape[2];
	PyArrayObject *outdata = NULL;
	double *spec;

	/* internal stuff */
	mfft *mtmh;
	double *samples = NULL;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "Oiid|ii", &o, &N, &step, &NW, &K, &adapt))
		return NULL;
	signal = (PyArrayObject*) PyArray_FromAny(o, NULL, 1, 1, NPY_CONTIGUOUS, NULL);
	if (signal==NULL) {
		PyErr_SetString(PyExc_TypeError, "Input signal must be an ndarray");
		return NULL;
	}

	if (K < 1) {
		K = NW*2-1;
	}

	/* allocate output array */
	out_shape[0] = N/2+1;
	out_shape[1] = PyArray_SIZE(signal) / step;
	outdata  = (PyArrayObject*) PyArray_ZEROS(2,out_shape,NPY_DOUBLE,1); // last arg give fortran-order
	spec = (double*) PyArray_DATA(outdata);

	/* coerce data to proper type */
	samples = coerce_ndarray_double(signal, &signal_cast);
	if (samples==NULL) {
		PyErr_SetString(PyExc_TypeError, "Unable to cast signal to supported data type");
		goto fail;
	}

	/* do the transform */
	mtmh = mtm_init_dpss(N, NW, K);
	mtm_spec(mtmh, spec, samples, PyArray_SIZE(signal),
		 step, adapt);
	mtm_destroy(mtmh);

	Py_DECREF(signal);
	Py_XDECREF(signal_cast);
	return PyArray_Return(outdata);
fail:
	Py_XDECREF(signal_cast);
	Py_XDECREF(signal);
	Py_XDECREF(outdata);
	return NULL;
}

static PyObject*
libtfr_mtm_psd(PyObject *self, PyObject *args)
{
	/* arguments */
	PyObject *o = NULL;
	PyArrayObject *signal = NULL;
	PyArrayObject *signal_cast = NULL;
	int N;
	double NW;
	int K = -1;
	int adapt = 1;

	/* output data */
	npy_intp out_shape[1];
	PyArrayObject *outdata = NULL;
	double *spec;

	/* internal stuff */
	mfft *mtmh;
	double *samples = NULL;
	double sigpow;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "Od|ii", &o, &NW, &K, &adapt))
		return NULL;
	signal = (PyArrayObject*) PyArray_FromAny(o, NULL, 1, 1, NPY_CONTIGUOUS, NULL);
	if (signal==NULL) {
		PyErr_SetString(PyExc_TypeError, "Input signal must be an ndarray");
		return NULL;
	}

	if (K < 1) {
		K = NW*2-1;
	}
	N = PyArray_SIZE(signal);

	/* allocate output array */
	out_shape[0] = N/2+1;
	outdata  = (PyArrayObject*) PyArray_ZEROS(1,out_shape,NPY_DOUBLE,0);
	spec = (double*) PyArray_DATA(outdata);

	/* coerce data to proper type */
	samples = coerce_ndarray_double(signal, &signal_cast);
	if (samples==NULL) {
		PyErr_SetString(PyExc_TypeError, "Unable to cast signal to supported data type");
		goto fail;
	}

	/* do the transform */
	mtmh = mtm_init_dpss(N, NW, K);
	sigpow = mtfft(mtmh, samples, N);
	if (adapt)
		mtpower(mtmh, spec, sigpow);
	else
		mtpower(mtmh, spec, 0.0);
	mtm_destroy(mtmh);

	Py_DECREF(signal);
	Py_XDECREF(signal_cast);
	return PyArray_Return(outdata);
fail:
	Py_XDECREF(signal_cast);
	Py_XDECREF(signal);
	Py_XDECREF(outdata);
	return NULL;
}


static PyObject*
libtfr_mtfft(PyObject *self, PyObject *args)
{
	/* arguments */
	PyObject *o = NULL;
	PyArrayObject *signal = NULL;
	PyArrayObject *signal_cast = NULL;
	int N;
	double NW;
	int K = -1;

	/* output data */
	npy_intp out_shape[2];
	PyArrayObject *outdata = NULL;
	npy_cdouble *spec;

	/* internal stuff */
	mfft *mtmh;
	double *samples = NULL;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "Od|i", &o, &NW, &K))
		return NULL;
	signal = (PyArrayObject*) PyArray_FromAny(o, NULL, 1, 1, NPY_CONTIGUOUS, NULL);
	if (signal==NULL) {
		PyErr_SetString(PyExc_TypeError, "Input signal must be an ndarray");
		return NULL;
	}

	if (K < 1) {
		K = NW*2-1;
	}
	N = PyArray_SIZE(signal);

	/* allocate output array */
	out_shape[0] = N;
	out_shape[1] = K;
	outdata  = (PyArrayObject*) PyArray_ZEROS(2,out_shape,NPY_CDOUBLE,1);
	spec = (npy_cdouble*) PyArray_DATA(outdata);

	/* coerce data to proper type */
	samples = coerce_ndarray_double(signal, &signal_cast);
	if (samples==NULL) {
		PyErr_SetString(PyExc_TypeError, "Unable to cast signal to supported data type");
		goto fail;
	}

	/* do the transform */
	mtmh = mtm_init_dpss(N, NW, K);
	mtfft(mtmh, samples, N);
	/* data are stored in half-complex form */
	hc2cmplx(mtmh, spec);
	mtm_destroy(mtmh);

	Py_DECREF(signal);
	Py_XDECREF(signal_cast);
	return PyArray_Return(outdata);
fail:
	Py_XDECREF(signal_cast);
	Py_XDECREF(signal);
	Py_XDECREF(outdata);
	return NULL;
}


static PyObject*
libtfr_dpss(PyObject *self, PyObject *args)
{
	/* arguments */
	int N;
	double NW;
	int K;

	/* output data */
	npy_intp tapers_shape[2];
	PyArrayObject *tapers = NULL;
	PyArrayObject *lambdas = NULL;
	int rV;

	/* parse arguments */
	if (!PyArg_ParseTuple(args, "idi", &N, &NW, &K))
		return NULL;

	/* allocate output arrays */
	tapers_shape[0] = N;
	tapers_shape[1] = K;
	tapers  = (PyArrayObject*) PyArray_ZEROS(2,tapers_shape,NPY_DOUBLE,1);
	lambdas = (PyArrayObject*) PyArray_ZEROS(1,tapers_shape+1,NPY_DOUBLE,1);

	double *taperp = (double*) PyArray_DATA(tapers);
	double *lambdap = (double*) PyArray_DATA(lambdas);

	rV = dpss(taperp, lambdap, N, NW, K);

	if (rV==0)
		return Py_BuildValue("(OO)", tapers, lambdas);
	else {
		PyErr_SetString(PyExc_TypeError, "Invalid DPSS parameters");
		Py_XDECREF(tapers);
		Py_XDECREF(lambdas);
		return NULL;
	}
}


#endif

static PyMethodDef _libtfr_methods[] = {
	{"tfr_spec", libtfr_tfr_spec, METH_VARARGS,""},
	{"hermf",libtfr_hermf, METH_VARARGS,""},
#ifndef NO_LAPACK
	{"mtm_spec", libtfr_mtm_spec, METH_VARARGS,""},
	{"mtm_psd", libtfr_mtm_psd, METH_VARARGS,""},
	{"mtfft", libtfr_mtfft, METH_VARARGS,""},
	{"dpss", libtfr_dpss, METH_VARARGS,""},
#endif
	{NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_libtfr(void)
{
	import_array();
	PyObject* m;

	m = Py_InitModule3("_libtfr", _libtfr_methods,
			   "Compute time-frequency reassignment spectrograms");

}
