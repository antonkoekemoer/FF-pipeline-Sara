"""
Gaussian fit to histogram of distribution and related
functions.

These functions were ported from existing IDL packages.
They are used by ACS/WFC reference files generation and
statistics packages.

:Authors: Pey Lian Lim (Python), others (IDL)

:Organization: Space Telescope Science Institute

:History:
    * 2011/08/01 PLL changed documentation to Sphinx format.
    * 2010/03/04 PLL modified documentations to comply with Epydoc convention.
    * 2009/11/24 PLL converted these functions from IDL to Python.
    
"""
# External modules
import numpy, pylab, scipy
from scipy import optimize

def robust_sigma(in_y, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.

    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.

    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.

    :History:
        * H Freudenreich, STX, 8/90
        * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
        * Converted to Python by P. L. Lim, 11/2009

    Examples
    -------- 
    >>> result = robust_sigma(in_y, zero=1)

    Parameters
    ----------
    in_y: array_like
        Vector of quantity for which the dispersion is
        to be calculated

    zero: int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.

    Returns
    -------
    out_val: float
        Dispersion value. If failed, returns -1.
        
    """
    # Flatten array
    y = in_y.reshape(in_y.size, )

    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3
    
    if zero:
        y0 = 0.0
    else:
        y0 = numpy.median(y)

    dy    = y - y0
    del_y = abs( dy )

    # First, the median absolute deviation MAD about the median:

    mad = numpy.median( del_y ) / c1

    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = numpy.mean( del_y ) / c2
    if mad < eps:
        return 0.0

    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u*u
    q  = numpy.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        print 'ROBUST_SIGMA: This distribution is TOO WEIRD! Returning', c_err
        return c_err
 
    numerator = numpy.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
    n    = y.size
    den1 = numpy.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )

    if siggma > 0:
        out_val = numpy.sqrt( siggma )
    else:
        out_val = 0.0

    return out_val

def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=0):
    """
    Computes an iteratively sigma-clipped mean on a
    data set. Clipping is done about median, but mean
    is returned.

    .. note:: MYMEANCLIP routine from ACS library.

    :History:
        * 21/10/1998 Written by RSH, RITSS
        * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
        * 24/11/2009 Converted to Python. PLL.

    Examples
    --------
    >>> mean, sigma = meanclip(indata)

    Parameters
    ----------
    indata: array_like
        Input data.

    clipsig: float
        Number of sigma at which to clip.
    
    maxiter: int
        Ceiling on number of clipping iterations.
    
    converge_num: float
        If the proportion of rejected pixels is less than
        this fraction, the iterations stop.
    
    verbose: {0, 1}
        Print messages to screen?

    Returns
    -------
    mean: float
        N-sigma clipped mean.

    sigma: float
        Standard deviation of remaining pixels.

    """
    # Flatten array
    skpix = indata.reshape( indata.size, )

    ct = indata.size
    iter = 0; c1 = 1.0 ; c2 = 0.0

    while (c1 >= c2) and (iter < maxiter):
        lastct = ct
        medval = numpy.median(skpix)
        sig = numpy.float64( numpy.std(skpix) ) # Bug - Need to recast
        wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]

        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1
    # End of while loop

    mean  = numpy.mean( skpix )
    sigma = robust_sigma( skpix )

    if verbose:
        prf = 'MEANCLIP:'
        print '%s %.1f-sigma clipped mean' % (prf, clipsig)
        print '%s Mean computed in %i iterations' % (prf, iter)
        print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)

    return mean, sigma

def mytotal(inarray, axis, type='meanclip'):
    """
    Collapse 2-D array in one dimension.

    .. note:: MYTOTAL routine from ACS library.

    :History:
        * Obtained from M. Sirianni.
        * Modified and converted to Python by P. L. Lim in 2009.

    Examples
    --------
    >>> collapsed_array = mytotal(inarray, 1, type='median')

    Parameters
    ----------
    inarray: array_like
        Input 2-D array.

    axis: {1, 2}
        Axis to collapse.
            * 1: Return values along Y.
            * 2: Return values along X.

    type: {'median', 'meanclip', 'stdev'}
        Algorithm to use.

    Returns
    -------
    out_arr: array_like
        1-D array collapsed along desired axis with desired
        algorithm.

    """
    func_name = 'MYTOTAL'
    out_arr = 0.0

    # Check inarray
    if inarray.ndim != 2:
        print '%s: Input array must be 2D' % func_name
        return out_arr

    # Check axis
    if axis == 1:
        n_out = inarray.shape[0]
    elif axis == 2:
        n_out = inarray.shape[1]
    else:
        print func_name, ': Axis not supported -', axis
        return out_arr

    # Initialize output array
    out_arr = numpy.zeros(n_out)
    out_rng = range(0, n_out)

    # Check type
    if type == 'meanclip':
        for i in out_rng:
            if axis == 1:
                im_i = inarray[i,:]
            else:
                im_i = inarray[:,i]
            mmean, msigma = meanclip(im_i, maxiter=10, converge_num=0.001)
            out_arr[i] = mmean
    elif type == 'stdev':
        for i in out_rng:
            if axis == 1:
                im_i = inarray[i,:]
            else:
                im_i = inarray[:,i]
            mmean, msigma = meanclip(im_i, maxiter=10, converge_num=0.001)
            out_arr[i] = msigma
    elif type == 'median':
        for i in out_rng:
            if axis == 1:
                im_i = inarray[i,:]
            else:
                im_i = inarray[:,i]
            out_arr[i] = numpy.median(im_i)
    else:
        print func_name, ': Type not supported -', type
        out_arr = 0.0

    return out_arr

def gaussian(height, center_x, width_x):
    """
    Returns a gaussian function with the given parameters.
    This is used for least square fitting optimization.

    .. note:: This is used by `msky`.

    Parameters
    ----------
    height: float
        Peak amplitude.

    center_x: float
        Peak location.

    width_x: float
        Sigma of gaussian curve.

    Returns
    -------
    x: lambda function
        Function used for optimization.

    """
    return lambda x: height * numpy.exp(-(center_x-x)**2/(2.0*width_x**2))

def msky(inarray, do_plot=0, verbose=0, ptitle='', func=0):
    """   
    Find modal sky on an array.
    
    First step is determination of median value and sigma.
    Histogram of the data and fit parabola to the
    logaritmic histogram. The coefficient of the parabola
    are used to get mode and sigma of the sky on the
    assumption that it is well fitted by a gaussian or
    2nd-degree polynomial.

    .. note:: MYSKY5 routine from ACS library.

    :History:
        * 11/25/2009 Created by PLL based on IDL routine by MS.
        * 01/29/2010 PLL fixed broken where statement for histogram.

    Parameters
    ----------
    inarray: array_like
        Input data.
    
    do_plot: {0, 1}
        Do plot?
    
    verbose: {0, 1}
        Print info to screen?
    
    ptitle: string
        Title of plot. Only used if `do_plot`=1.
    
    func: {0, 1}
        Function for fitting:
            * 0: 2nd degree polynomial
            * 1: Gaussian

    Returns
    -------
    mmean: float
        Mode of fitted function.

    sigma: float
        Sigma of fitted function.

    """
    # Constants
    func_name = 'MSKY'
    nsig = 8.0
    c1 = 2.5 # 2.8
    c2 = 0.8 # 1.3

    # Min/max of input array
    arr_min = numpy.min(inarray)
    arr_max = numpy.max(inarray)

    # Get sigma
    mmean, sigma = meanclip(inarray, clipsig=5.0, maxiter=10, verbose=verbose)
    if sigma <= 0:
        print ' '
        print '%s: Weird distribution' % func_name
        print '%-6s: %.6f' % ('MEAN', mmean)
        print '%-6s: %.6f' % ('STDDEV', sigma)
        print '%-6s: %.6f' % ('MIN', arr_min)
        print '%-6s: %.6f' % ('MAX', arr_max)
        return mmean, sigma

    # Print info
    if verbose:
        print ' '
        print '%s input array info' % func_name
        print '%-6s: %.6f' % ('MIN', arr_min)
        print '%-6s: %.6f' % ('MAX', arr_max)

    # Flatten input array
    arr_1d = inarray.reshape( inarray.size, )

    # Define min and max for the histogram
    x = nsig * sigma
    mmean = numpy.median(arr_1d)
    minhist = mmean - x
    maxhist = mmean + x
    ufi = inarray[ numpy.where((inarray > minhist) & (inarray < maxhist)) ]

    # Calculate 25% and 75% percentile to get the interquartile range
    # IRQ = pc75-pc25
    # zenman, A. J. 1991. 
    sixd = numpy.argsort( ufi )
    ndata = ufi.size
    pc25 = ufi[ sixd[0.25*ndata] ]
    pc75 = ufi[ sixd[0.75*ndata] ] 
    irq = pc75 - pc25
    step = 2.0 * irq * ndata**(-1.0/3.0)

    # Calculate number of bins to use
    nbin = round(2*x / step - 1)

    # Histogram
    # http://www.scipy.org/Tentative_NumPy_Tutorial
    yhist, hbin = numpy.histogram(arr_1d, range=(minhist,maxhist), bins=nbin)
    xhist = 0.5 * ( hbin[1:] + hbin[:-1] )

    # Define xmin and xmax for the 2-0rder fit
    x1 = mmean - c1*sigma
    x2 = mmean + c2*sigma

    # Select the points beween x1 and x2 for the fit
    w = numpy.where((xhist > x1) & (xhist < x2) & (yhist > 0))
    count = len(w[0])
    xwg  = xhist[w]
    nywg = yhist[w]
    if count < 2:
        print ' '
        print '%s: Singular matrix' % func_name
        print '%-6s: %.6f' % ('X[W]', xwg)
        print '%-6s: %.6f' % ('NY[W]', nywg)
        print '%-6s: %.6f' % ('MEDIAN', mmean)
        return mmean, sigma

    # Change to log scale
    yhist = numpy.log10(yhist)
    iyh   = numpy.where( ~numpy.isinf(yhist) )
    xhist = xhist[iyh]
    yhist = yhist[iyh]
    nywg  = numpy.log10(nywg)

    # Calculate the fit coefficients
    ymax = nywg.max()

    # Gaussian
    # http://www.scipy.org/Cookbook/FittingData
    if func == 1:
        if verbose: print '%s: Fitting gaussian' % func_name
        
        # Initial guess
        ysum = numpy.sum(nywg)
        mmean = numpy.sum(xwg * nywg) / ysum
        sigma = numpy.sqrt(numpy.abs(numpy.sum((xwg-mmean)**2*nywg)/ ysum))
        params = (ymax, mmean, sigma)

        # Error function to minimize
        errorfunction = lambda p: scipy.ravel(gaussian(*p)(xhist) - yhist)

        # Linear least square fitting
        a_opt, a_success = optimize.leastsq(errorfunction, params)
        ymax  = a_opt[0]
        mmean = a_opt[1]
        sigma = a_opt[2]

        # Fit to entire range
        yall = ymax * numpy.exp(-(xhist-mmean)**2/(2.0*sigma**2))

    # 2nd degree polynomial
    else:
        if verbose: print '%s: Fitting 2nd deg polynomial' % func_name

        # Polynomial
        a_opt = numpy.polyfit(xwg, nywg, 2)
        mmean = -0.5*a_opt[1]/a_opt[0]
        sigma = numpy.sqrt(-0.5 / ( a_opt[0]*numpy.log(10) ) )

        # Fit to entire range
        yall = a_opt[0]*xhist**2 + a_opt[1]*xhist + a_opt[2]

    # Print results
    if verbose:
        print ' '
        print '%s: Results' % func_name
        print '%-6s: %.6f' % ('MODE', mmean)
        print '%-6s: %.6f' % ('SIGMA', sigma)
        print ' '

    plot_y1 = 0
    plot_y2 = ymax + 0.1

    # Plot results
    if do_plot:
        pylab.clf()
        # Data points
        pylab.plot(xhist, yhist, 'ko')
        # Mark fitted region
        pylab.plot(xwg, nywg, 'bo', hold=1)
        # Draw fitted function
        pylab.plot(xhist, yall, 'r--', hold=1)
        pylab.axvline(x=mmean, color='r')
        # Plot xmin xmax for the fit
        pylab.axvline(x=x1, ymin=0, ymax=0.1, color='b')
        pylab.axvline(x=x2, ymin=0, ymax=0.1, color='b')
        # Axis limits and labels
        pylab.axis([minhist, maxhist, plot_y1, plot_y2])
        pylab.xlabel('Pix value')
        pylab.ylabel('Log num of pix')
        pylab.title(ptitle)
        #pylab.show()

        # Pause for visual examination of plot
        x_user = raw_input('ENTER ANY CHAR TO CONTINUE: ')

    return mmean, sigma
