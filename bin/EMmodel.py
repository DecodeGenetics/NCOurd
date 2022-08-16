#!/usr/bin/python3

import logging
import numpy as np
import numpy.random as nrand
import pandas as pd
import sys
from scipy import interpolate 
from scipy.special import xlogy,gamma
from scipy import stats
import scipy.optimize as optimize
import argparse
from grid import Grid
import re
from multiprocessing import Pool

np.set_printoptions(precision=4, suppress=False) 
tf=None

## Make sure the probability value negative binomial distribution doen't exceed 1-epsilon  
epsilon=1e-10

grid=Grid(type='linear')
kgrid=grid.k
wgrid=grid.w

## Detection probability function ##
def make_detection_probability_function(detdf):
  '''
  Compute an interpolant from the probability detection dataframe
  returned by the function compute_detection_probability.
  '''
  return interpolate.interp1d(detdf["k"],detdf["p"],kind="cubic")

def read_detection_probability_function(filepath):
  '''
  Read detection probability function data as a dataframe from
  file and return the corresponding interpolation function.
  '''
  detdf = pd.read_csv(filepath,sep="\t")
  return make_detection_probability_function(detdf)

def readparamfile(filepath):
    params={}
    df=pd.read_csv(filepath,sep='\t')
    
    keys  =['shape','k', 'scale','beta', 'lambda', 'part',     'weight']
    values=['Ks',   'Ks','betas','betas','lambdas','partition','weights']

    keyval=dict(zip(keys,values))
    for key in keys:
        if key in df.columns:
            params[keyval[key]]=list(df[key])

    return params

def strlist(string,n):
    return [string+str(i) for i in range(n)]

def weights2parts(ls,ks,weights,detdf):
    pdfs=[lambda x,k=k,l=l:  l**k*x**(k-1)*np.exp(-l*x)/gamma(k) for l,k in zip(ls,ks)]
    integrals=[np.sum(pdf(kgrid)*wgrid*detdf) for pdf in pdfs]
    observed=[w*i for w,i in zip(weights,integrals)]
    return [o/sum(observed) for o in observed] 



## EM functions ##
def nbinom_pmf(x,k,l):
    p=l
    n=k/(1-p)
    return  stats.nbinom.pmf(x,n,p,loc=args.nbinom_loc)

def nbinom_logpmf(x,k,l):
    p=l
    n=k/(1-p)
    return stats.nbinom.logpmf(x,n,p,loc=args.nbinom_loc)

geom_pmf=lambda x,k,l: nbinom_pmf(x,1,l)
geom_logpmf=lambda x,k,l: nbinom_logpmf(x,1,l)

def gamma_pdf(x,k,l):
    return stats.gamma(a=k,scale=1/l).pdf(kgrid)

def gamma_logpdf(x,k,l):
    return stats.gamma(a=k,scale=1/l).logpdf(kgrid)

def loglikelihod_mixed(k,l,gridweights,support_idx,w):
    gval=probfun(kgrid,k,l)
    lgval=logprobfun(kgrid,k,l)
    lgmaxidx=np.argmax(lgval)
    lsupport=lgval[np.array([max(support_idx[i],lgmaxidx) for i in range(len(support_idx))])]
    lsupport[support_idx==0]=0

    integrals=np.array([\
        np.sum(wgrid*gval*wts) if sidx==0 else np.sum(np.exp(lgval[sidx:]-lsup)*wts[sidx:]*wgrid[sidx:]) \
        for wts,lsup,sidx in zip(gridweights,lsupport,support_idx)]) 
    loglikelihoods=xlogy(w,integrals)+w*lsupport

    return loglikelihoods

def loadtwfun(twfunpaths,maxtracts,skips,rand,replace=True):
    if rand is None:
        return np.concatenate([np.loadtxt(twf,max_rows=mt,skiprows=skip) for twf,mt,skip in zip(twfunpaths,maxtracts,skips)])
    
    nrand.seed(rand)
    alltwfun=[np.loadtxt(twf) for twf in twfunpaths]
 
    for i in range(len(alltwfun)):
        if maxtracts[i] is None:
            maxtracts[i]=alltwfun[i].shape[0]
            if args.verbose:
                print("#sampling", maxtracts[i])
    return np.concatenate([twf[nrand.choice(twf.shape[0],mt,replace=replace),:] for twf,mt in zip(alltwfun,maxtracts)]) 

def weightedMaxLikelihoodK(params):
    k,l,wts=params
    w=np.append(wts,-sum(wts))
    F=lambda k,l=l,w=w: tf(k,l,w)
    res=optimize.minimize_scalar(F,method='bounded',bounds=[max(k/2,args.mink),2*k])
    return res.x

def weightedMaxLikelihoodL(params):
    k,l,wts=params
    w=np.append(wts,-sum(wts))
    F=lambda l,k=k,w=w: tf(k,l,w)
    res=optimize.minimize_scalar(F,method='bounded',bounds=[l/2,min(2*l,1-epsilon)])
    return res.x

# does not support multi process
def weightedMaxLikelihoodKL(f,k,l,wts,params): 
    w=np.append(wts,-sum(wts))
    F=lambda kl,w=w: f(*kl,w)
    res=optimize.minimize(F,(k,l),method=params['minmeth'],bounds=((max(args.mink,.5*k),max(args.mink,2*k)),(.5*l,2*l)))
    return res.x

# does not support multi process
def m_step_simul(tf,ks,ls,weights,params):
    kls=[weightedMaxLikelihoodKL(tf,k,l,wts,params) for k,l,wts in zip(ks,ls,weights)] # old values of k an l are used for bounds

    ks=[kl[0] for kl in kls]
    ls=[kl[1] for kl in kls]
    
    return ks,ls

def m_step_sep(tf,ks,ls,weights,params):
    if args.multiprocess:
        pool=params['pool']
        chunks=list(zip(ks,ls,weights))
        ls=pool.map(weightedMaxLikelihoodL,chunks)
    else:
        ls=[weightedMaxLikelihoodL([k,l,wts]) for k,l,wts in zip(ks,ls,weights)] 
    if args.multiprocess:
        pool=params['pool']
        chunks=list(zip(ks,ls,weights))
        ks=pool.map(weightedMaxLikelihoodK,chunks)
    else:
        ks=[weightedMaxLikelihoodK([k,l,wts]) for k,l,wts in zip(ks,ls,weights)]     

    return ks,ls

def EstimateLambdaK(twfun,detdf,n=2,maxIter=100,params=None,minK=0.000001,threshold=1e-7,method="simul"):
    global tf
    ks=params["Ks"] if params and "Ks" in params else [2]*n 
    if params and "betas" in params:
        bs=params["betas"]
    else:
        const=np.log(10)
        inc=(np.log(100000)-np.log(10))/(n+1)
        bs=[np.exp(const+(i+1)*inc)/ks[i] for i in range(n)]
    ls=params["lambdas"] if params and "lambdas" in params else [1/b for b in bs]
    partition=params["partition"] if params and "partition" in params else [1/n]*n
    if 'weights' in params:
        partition=weights2parts(ls,ks,params['weights'],detdf)

    ntracts=len(twfun)
    sumwts=[1]*ntracts+[-ntracts]
    ones=np.ones(ntracts+1)
    gridweights=np.append(twfun,[detdf],axis=0)
    nonzero=gridweights>0
    support_idx=np.array([np.arange(gridweights.shape[1])[nonzero[i]][0] for i in range(gridweights.shape[0])])
    support=kgrid[support_idx]

    tf=lambda k,l,w: -np.sum(loglikelihod_mixed(k,l,gridweights,support_idx,w))
    if args.verbose:
        print("#","tracts",twfun.shape[0],"starting parameters",*partition,*ks,*bs,*[k*b for k,b in zip(ks,bs)])

    if args.multiprocess:
        params['pool']=Pool(processes=len(ls))

    # print header
    print("LogLikelihood","Qfuct_after",*strlist("Qf_a",n),"Qfunct_before",*strlist("Qf_b",n),
            "iteration","delta",*strlist("partition",n),*strlist('k',n),*strlist("beta",n),*strlist("mean",n),"delta2",sep='\t')


    weights=np.array([[1/n]*ntracts  for i in range(n)])
    for i in range(maxIter):
        # E-step
        lastparams=partition+ks+ls
        oldks,oldls,oldpartition=ks,ls,partition
        oldwts=weights


        # Membership weights
        weights=None

        reajust=max(1,1 if not args.accwts else args.accwts-i)
        if reajust>1:
            if args.verbose:
                print ("# Reajusting weights:", reajust,"times")


        for j in range(reajust):
            logval=[]
            for k,l,part in zip(ks,ls,partition):
                logvalues=loglikelihod_mixed(k,l,gridweights,support_idx,ones)
                logval.append(logvalues[:-1]-logvalues[-1])
            logval=np.array(logval)
            if j==0:
                keeplogval=logval
            logval2=logval-np.max(logval,axis=0)
            preval=np.exp(logval2)
            val2=preval*np.array(partition).reshape(len(partition),1)
            sv2=np.sum(val2,axis=0)
            weights=val2/sv2

            partition=[sum(w) for w in weights]
            partition=[p/sum(partition) for p in partition]

        # M-step
        if method in ['sep','separated']:
            ks,ls=m_step_sep(tf,ks,ls,weights,params)
        elif method in ['simul','simultaneous']:
            ks,ls=m_step_simul(tf,ks,ls,weights,params)
        else:
            print("Invalid method",method)
            exit()

        ll=[]
        for k,l,wts in zip(ks,ls,weights):
            w=np.append(wts,-sum(wts))
            ll.append(tf(k,l,w))
        pll=[]
        for k,l,wts in zip(oldks,oldls,weights):
            w=np.append(wts,-sum(wts))
            pll.append(tf(k,l,w))


        # Calculate the value of the  loglikilehood function
        PartialLikelihood=np.exp(keeplogval)*np.array(oldpartition).reshape(len(partition),1)
        FullLikelihood=np.sum(PartialLikelihood,axis=0)

        for j in range(n): 
            if pll[j]<ll[j]:
                ks[j]=oldks[j]
                ls[j]=oldls[j]
                

        bs=[1/l for l in ls]

        currparams=partition+ks+ls
        delta=(sum([(1-x/y)**2 for x,y in zip(currparams,lastparams)])/(len(currparams)-1))**0.5
        dval=(weights-oldwts)/np.maximum(oldwts,1e-10)
        delta2=np.sqrt(np.nansum(dval*dval)/n/ntracts)
        print(np.sum(np.log(FullLikelihood)) ,sum(ll),*ll,sum(pll),*pll,i,delta,*np.array(partition),
                *np.array(ks),*np.array(bs),*np.array([k*b for k,b in zip(ks,bs)]),delta2,sep='\t')

        if delta<threshold:
            print("# Parameters have converged")
            break 

    if args.multiprocess:
        params['pool'].close()
        params['pool'].terminate()

def readIndex(idxfile):
    result=[]
    with open(idxfile) as fh:
        for line in fh:
            result.extend([int(item) for item in re.split('[\W,;]+',line) if item])
    return result

## Main ##
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run EM model on tractsets')
    parser.add_argument('--detfile', type=str, default='markerset/01.detdf.gz', help='file containing detection function')
    parser.add_argument('--ngroups', type=int,help='number mixture components to fit - Default is #tractfiles')
    parser.add_argument('--maxiter', type=int, default=1000, help='maximum number of iterations - Default is 1000')
    parser.add_argument('--verbose', action='store_true', help='more printing')

    # For simulated tracts
    parser.add_argument('--ntracts', type=str, help='number of tracts to take from each tractfile')
    parser.add_argument('--nskip', type=str, help='number of tracts to skip from each tractfile')


    parser.add_argument('--multiprocess', action='store_true', help='use multiprocessing')
    parser.add_argument('--paramfile',type=str,help='File containing initial parameters for the EM model')
    parser.add_argument('--method',type=str,default='sep',help='Simultaneously "simul" or separately "sep" evaluate shape and scale')
    parser.add_argument('--dclass',type=str,default='nbinom',help='Class of distributions excluding the default negative binomial distributions "nbinom" gamma "gamma" and geometric "geom" distributions are supported')
    parser.add_argument('--minmeth',type=str,help='Method to pass to scipy.optimize.minimize for simultaneous evaluation must support bounds')

    parser.add_argument('--maxspan',type=int,help='Discard tracts which are at least this long')
    parser.add_argument('--nbinom_loc',type=int,default=1,help='If nbiom distributions are used this sets location parameter of nbiom distributions - Default is 1') 
    parser.add_argument('--randomsample', type=int,help='Sample random tracts from input tracts with given seed')

    parser.add_argument('--accwts',type=int,help='Accelerate weights for the first N iterations')
    parser.add_argument('--cubicgrid',action='store_true',help='use cubic wgrid instead of linear for integration')
    parser.add_argument('--twidx',type=str,help='indexes to use in the tractfile')
    parser.add_argument('--mink', type=float, default=0.0, help='minimum value of k (minimum variance is mean**2*k)')
    parser.add_argument('tractfiles', nargs='+', help='tract weight function files to use')


    args = parser.parse_args()

    if args. verbose:
        print("#",args)
    nrand.seed(0)

    detdfpath=args.detfile
    twfunpaths=args.tractfiles
    if args.ngroups:
        n=args.ngroups
    else:
        n=len(twfunpaths)
    maxtracts=[int(s) for s in args.ntracts.split(",")] if args.ntracts is not None else [None]*n
    skips=[int(s) for s in args.nskip.split(",")] if args.nskip is not None else [0]*n
    if len(maxtracts)==1:
        maxtracts=[maxtracts[0]]*n
        skips=[skips[0]]*n

    logging.info("reading detection probability")
    dpfunc = read_detection_probability_function(detdfpath)
    detdf=dpfunc(kgrid)
    twfun = loadtwfun(twfunpaths,maxtracts,skips,rand=args.randomsample)
    if args.twidx:
        twfun=twfun[readIndex(args.twidx),:]

    params={}
    if args.paramfile:
        params=readparamfile(args.paramfile)

    if args.maxspan:
        twfun=twfun[(twfun[:,kgrid<args.maxspan]>0).any(axis=1),:]

    if args.cubicgrid:
        grid=Grid(type='cubic')
        kgrid=grid.k
        wgrid=grid.w

    # Family of distributions to fit
    if args.dclass=='nbinom':
        probfun=nbinom_pmf
        logprobfun=nbinom_logpmf
    elif args.dclass=='geom':
        probfun=geom_pmf
        logprobfun=geom_logpmf
        params['Ks']=[1]*n
    elif args.dclass=='gamma':
        probfun=gamma_pdf
        logprobfun=gamma_logpdf

    k=EstimateLambdaK(twfun,detdf,n=n,maxIter=args.maxiter,params=params,method=args.method)

