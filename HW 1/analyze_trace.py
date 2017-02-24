#!/usr/bin/env python

from __future__ import print_function
import numpy as np

def mean(trace):
    """ calculate the mean of a trace of scalar data
    results should be identical to np.mean(trace)
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the mean of this trace of scalars 
    """

    total = 0.0        # store sum of all data points
    num   = len(trace) # count total number of data points

    for i in range(num):
        total += trace[i]
    # end for i

    return total/float(num)

# end def mean

def std(trace):
    """ calculate the standard deviation of a trace of scalar data
    results should be identical to np.std(trace)
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the standard deviation of this trace of scalars 
    """

    stddev = 0.0 
    # calculate stadard deviation
    mu = mean(trace)
    N = len(trace)

    for i in range(N):
        stddev = stddev + (trace[i]-mu)**2

    stddev=np.sqrt(stddev/(N-1))

    return stddev

# end def std

def corr(trace):
    """ calculate the autocorrelation of a trace of scalar data
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the autocorrelation of this trace of scalars 
    """

    correlation = 1.0
    su =0.0
    n = len(trace)
    for k in range(1,n):
        cc = compcorr(trace,k)
        if cc<=0:
            break
        su = su + cc
        
        

    correlation = 1 + 2*su
    # calculate auto correlation
    # for t in range(n):
    #     for k in range(n):
    #         if t+k>=n:
    #             break
    #         su = su + (trace[t]-mu)*(trace[t+k]-mu)
    #         if su<=0:
    #             K=k
    #             break
    #     correlation = correlation+ 2*


    #print(correlation)
    return correlation

def compcorr(trace, k):
    R = 0.0
    n= len(trace)
    mu = mean(trace)
    sigma = std(trace)
    sigma = sigma**2
    for t in range(n-k):
        R = R + (trace[t] - mu)*(trace[t+k] - mu)
    
    R = R/((n-k)*sigma)
    return R

# end def corr

def error(trace):
    """ calculate the standard error of a trace of scalar data
    for uncorrelated data, this should match np.std(trace)/np.sqrt(len(trace))
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the standard error of this trace of scalars 
    """

    #stderr = 0.0
    N_eff = len(trace)/corr(trace)
    sigma = std(trace)
    stderr = sigma/np.sqrt(N_eff)
    # calculate standard error
    print(sigma/len(trace))
    return stderr

# end def error

def stats(trace):
    """ return statistics of a trace of scalar data
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return (mean,stddev,auto_corr,error)
    """
    # basically a composition of functions implemented above
    mymean = mean(trace)
    mystd  = std(trace)
    mycorr = corr(trace)
    myerr  = error(trace)

    return (mymean,mystd,mycorr,myerr)
# end def stats

if __name__ == '__main__':
    """ code protected by __main__ namespace will not be executed during import """
    import argparse

    # parse command line input for trace file name
    parser = argparse.ArgumentParser(description='analyze a trace')
    parser.add_argument('filename', type=str, help='filename containing a scalar trace')
    parser.add_argument('-p','--plot', action='store_true', help='plot data')
    parser.add_argument('-i','--initial_index', type=int, default=0, help='initial index of data to include')
    parser.add_argument('-f','--final_index', type=int, default=-1, help='final index of data to include, must be larger than initial_index')
    args = parser.parse_args()
    trace_file_name = args.filename

    # read trace file
    trace = np.loadtxt( trace_file_name )

    # determine final cutoff
    final_index = -1; # default
    if args.final_index == -1:
        final_index = len(trace)
    else:
        final_index = args.final_index
    # end if

    # guard against misuse
    if final_index > len(trace):
        raise RuntimeError("final_index exceeds length of trace")
    # end if

    # cut out interested portion
    trace = trace[args.initial_index:final_index]

    # calculate statistics
    mymean,mystd,mycorr,myerr = stats( trace )

    # formatted output of calculated statistics
    header = "%25s   mean   stddev   corr   err" % "filename"
    fmt    = "{filename:s}  {mean:1.4f}  {stddev:1.4f}   {corr:1.2f}  {err:1.4f}"
    output = fmt.format(
            filename = trace_file_name
          , mean     = mymean
          , stddev   = mystd
          , corr     = mycorr
          , err      = myerr )

    print( header )
    print( output )

    if (args.plot):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(
            1,2, gridspec_kw = {'width_ratios':[3, 1]})
        ax[0].set_xlabel("index", fontsize=14)
        ax[0].set_ylabel("data" , fontsize=14)
        ax[1].set_xlabel("freq.", fontsize=14)
        ax[1].get_yaxis().tick_right()

        # plot trace
        ax[0].plot(trace,c='black')

        # plot histogram
        wgt,bins,patches = ax[1].hist(trace, bins=30, normed=True
            , fc='gray', alpha=0.5, orientation='horizontal')
        # moving averge to obtain bin centers
        bins = np.array( [(bins[i-1]+bins[i])/2. for i in range(1,len(bins))] )
        def _gauss1d(x,mu,sig):
            norm  = 1./np.sqrt(2.*sig*sig*np.pi)
            gauss = np.exp(-(x-mu)*(x-mu)/(2*sig*sig)) 
            return norm*gauss
        # end def
        ax[1].plot(_gauss1d(bins,mymean,mystd),bins,lw=2,c="black")
        ax[1].set_xticks([0,0.5,1])

        # overlay statistics
        for myax in ax:
            myax.axhline( mymean, c='b', lw=2, label="mean = %1.4f" % mymean )
            myax.axhline( mymean+mystd, ls="--", c="gray", lw=2, label="std = %1.2f" % mystd )
            myax.axhline( mymean-mystd, ls="--", c="gray", lw=2 )
        # end for myax
        ax[0].legend(loc='best')

        plt.show()

    # end if plot

# end __main__
