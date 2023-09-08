import numpy as np

def LyapSpec(Vfun,Jfun, T, x0, r,discard=100):
    # based on MatLAB code by Anton O. Belyakov
    # notable edits: remove transients
    n = len(x0)  # dimension of the system
    k = n  # by default calculate full Lyapunov spectrum
    t = T[0]  # initial time
    x = x0.copy()  # set initial state
    LE = np.zeros(k)  # Lyapunov exponents
    e = np.eye(n, k)  # make a matrix of variations (n x k)
    ee = np.eye(k)  # matrix for Gram-Schmidt orthonormalization
    trJ = 0  # integral of the trace of Jacobian matrix
    oldtraceJ = np.trace(Jfun(t, x))  # trace of Jacobian matrix in previous step
    warn = 0
    dtinv = np.sqrt(r)
    dt = 1 / r  # sampling interval
    Tstart_lya=  T[0]+discard*dt # time to start computing lya
    Tlist = np.arange(Tstart_lya, T[1], dt)
    # solve for discard steps
    for t in np.arange(T[0],Tstart_lya, dt):
        x = rk4(Vfun,t,x,dt)
        # Alternatively use Leapfrog:
        #J = Jfun(t,x)
        #v = Vfun(t,x)
        #x = x + dt * (v + 0.5 * dt * (J @ v))
    # compute lya
    for t in Tlist:
        # Alternatively use Leapfrog:
        #J = Jfun(t,x)
        #v = Vfun(t,x)
        #x = x + dt * (v + 0.5 * dt * (J @ v))
        J = Jfun(t,x)
        #rk4
        x = rk4(Vfun,t,x,dt)
        # Leapfrog prediction
        #v = Vfun(t,x)
        #x = x + dt * (v + 0.5 * dt * (J @ v))
        # Heun's step for e
        e_t = e + dt * J @ e  # Euler prediction of e
        J = Jfun(t + dt, x)  # next Jacobian matrix
        e = 0.5 * (e_t + e + dt * J @ e_t)  # Trapezoidal correction of e
        # Gram-Schmidt orthonormalization to calculate LEs
        for j in range(k):
            e[:, j] = e[:, :j+1] @ ee[:j+1, j]
            nrm = np.linalg.norm(e[:, j])
            e[:, j] = e[:, j] / nrm
            if j < k - 1:
                ee[j, j+1:k+1] = -(e[:, j+1:k+1].T @ e[:, j])
            LE[j] = LE[j] + np.log(nrm)  # accumulated LEs
        # trapezoidal integration of the sum of LEs
        traceJ = np.trace(J)  # derivative of the sum of LEs
        if np.abs(oldtraceJ - traceJ) - dtinv > 0 and warn == 0:
            warn = (oldtraceJ - traceJ) ** 2 - r
            print(f"Sampling rate r is too small. (oldtraceJ - traceJ)^2 - r = {warn}")
        trJ = trJ + 0.5 * dt * (oldtraceJ + traceJ)  # sum of LEs
        oldtraceJ = traceJ
        #print(LE)
    LE = LE / (Tlist[-1]-Tlist[0])  # LEs calculation
    trJ = trJ / (Tlist[-1]-Tlist[0])  # calculation of LEs' sum
    return LE, trJ, x , np.arange(Tstart_lya, T[1], dt,dtype='float64')