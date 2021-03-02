import numpy as np
import glassofast

def glassofast_f(S, rho, thr=1.0e-4, maxIt=1e4, start="cold", W=None, X=None):
	assert start in ['cold', 'warm']
	n = int(S.shape[0])
	if np.shape(rho) == ():
		L = np.ones((n, n)) * rho
	else:
		assert np.shape(rho)[0]==np.shape(rho)[1]
		L = rho
	if start =="cold":
		warm = 0
		W = X = np.zeros((n,n))
	else:
		warm = 1
		if (W==None or X==None):
			print("Warm start specified: W and X must not be None")
		else:
			W = W
			X = X
	S = np.array(S, order='f')
	X = np.array(X, order='f')
	W = np.array(W, order='f')
	L = np.array(L, order='f')
	thr = float(thr)
	maxIt = int(maxIt)
	warm = int(warm)

	glassofast.glasso(S, X, W, L, maxIt, warm, thr, n)
	cov_est = X
	pcs_est = W
	return cov_est, pcs_est