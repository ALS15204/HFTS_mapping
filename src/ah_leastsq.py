def leastsq(func, x0, dx, R=100., print_error_mesg=True,error_only=False, xtol=1.49012e-8,):
    """
    Rejig the inputs of scipy.optimize.leastsq so that they do what I
    want them to.
    \nInputs:\n
      func -- The same as for leastsq.
      x0 -- The same as for leastsq.
      dx -- A sequence of the same length as x0 containing the desired
      absolute stepsize to use when calculating the finite difference
      Jacobean.
      R -- The ratio of two step sizes: Dx/dx. Where Dx is the maximum
      stepsize taken at any time. Note that this is only valid for the
      first iteration, after which leastsq appears to approximately
      double the 'factor' parameter.
      print_error_mesg -- if True output error code and message if failure
    \nOutputs: (x,sigma)\n
    x -- array of fitted parameters
    sigma -- error of these
    The reason for doing this is that I found it difficult to tweak
    the epsfcn, diag and factor parametres of leastsq to do what I
    wanted, as far as I can determine these behave in the following
    way:
    dx = x*sqrt(epsfcn) ; x!=0,
    dx = 1*sqrt(epsfcn) ; x==0.
    Default epsfcn=2.2e-16 on scucomp2.
    Dx = abs(x*100)      ; x!=0, factor is not set,
    Dx = abs(x*factor)   ; x!=0, factor is set,
    Dx = abs(factor)     ; x==0, factor is set,
    Dx = 100             ; x==0, factor is not set, diag is not set,
    Dx = abs(100/diag)   ; x==0, factor is not set, diag is set,
    Dx = abs(factor/diag); x==0, factor is set, diag is set.
    Many confusing cases, particularly annoying when initial x==0 and
    it is not possible to control dx or Dx individually for each
    parameter.
    My solution was to add a large value to each parameter so that
    there is little or no chance it will change magnitude during the
    course of the optimisation. This value was calculated separately
    for each parameter giving individual control over dx. I did not
    think of a way to also control Dx individually, instead the ratio
    R=Dx/dx may be globally set.
    """
    from scipy import optimize
    ## limit the number of evaluation to a minimum number to compute
    ## the uncertainty from the second derivative - make R small to
    ## improve performance? - Doesn't work for very large number of
    ## parameters - errors are all nan, probably because of a bug in
    ## leastsq?
    if error_only:
        maxfev = len(x0)+1
        R = 1.
    else:
        maxfev = 0
    ## try and wangle actual inputs of numpy.leastsq to get the right
    ## step sizes
    x0=np.array(x0)
    dx=np.array(dx)
    epsfcn = 1e-15              # required that sqrt(epsfcn)<<dp/p
    xshift = x0+dx/np.sqrt(epsfcn)    # required that xshift>>p
    factor = R*np.sqrt(epsfcn)
    x = x0-xshift
    ## perform optimisation. try block is for the case where failure
    ## to calculte error
    try:
        (x,cov_x,info,mesg,success)=optimize.leastsq(
            lambda x:func(x+xshift),
            x,
            epsfcn=epsfcn,
            factor=factor,
            full_output=True,
            maxfev=maxfev,
            xtol = xtol,
            )
    except ValueError as err:
        if err.message=='array must not contain infs or NaNs':
            raise Exception('Bad covariance matrix in error calculation, residual independent of some variable?')
        else:
            raise
    ## warn on error if requested
    if (not success) & print_error_mesg:
        import warnings
        warnings.warn("leastsq exit code: "+str(success)+mesg)
    ## sometimes this is not an array
    if not np.iterable(x): x=[x]
    ## attempt to calculate covariance of parameters
    if cov_x is None:
        sigma_x = np.nan*np.ones(len(x))
    else:
        chisq=sum(info["fvec"]*info["fvec"])
        dof=len(info["fvec"])-len(x)+1        # degrees of freedom
        ## assumes unweighted data with experimental uncertainty
        ## deduced from fitted residual. ref gavin2011.
        std_y = np.sqrt(chisq/dof)
        sigma_x = np.sqrt(cov_x.diagonal())*std_y
    return(x+xshift,sigma_x)
