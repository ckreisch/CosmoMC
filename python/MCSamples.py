# MCSamples.py

import os
import sys
import math
import numpy as np
from chains import chains

# =============================================================================

class Ranges():
    
    def __init__(self, fileName=None, setParamNameFile=None):
        self.names = []
        self.mins = {}
        self.maxs = {}        
        if fileName is not None: self.loadFromFile(fileName)

    def loadFromFile(self, fileName):
        self.filenameLoadedFrom = os.path.split(fileName)[1]
        f = open(fileName)
        for line in f:
            name, mini, maxi = self.readValues(line.strip())
            if name: 
                self.names.append(name)
                self.mins[name] = mini
                self.maxs[name] = maxi                

    def readValues(self, line):
        name, mini, maxi = None, None, None
        strings = [ s for s in line.split(" ")  if s<>'' ]
        if len(strings)==3:
            name = strings[0]
            mini = float(strings[1])
            maxi = float(strings[2])
        return name, mini, maxi

    def listNames(self):
        return self.names

    def min(self, name, error=False):
        if self.mins.has_key(name):
            return self.mins[name]
        if error: raise Exception("Name not found:" + name)
        return None

    def max(self, name, error=False):
        if self.maxs.has_key(name):
            return self.maxs[name]
        if error: raise Exception("Name not found:" + name)
        return None

# =============================================================================

class Density1D():

    SPLINE_DANGLE = 1.e30

    def __init__(self, n, spacing):
        self.n   = n
        self.X   = np.zeros(n) 
        self.P   = np.zeros(n) 
        self.ddP = np.zeros(n) 
        self.spacing = spacing

    def Prob(self, x):
        if (x > self.X[self.n-1] - self.spacing/1e6):
            if ( x > self.X[self.n-1] + self.spacing/1e6):
                print 'Density: x too big ', x
                sys.exit()
            return self.P[self.n-1]
        
        if (x < self.X[0] - self.spacing/1e6):
            print 'Density: out of range ', x, self.X[0], self.X[self.n-1]
            sys.exit()

        # fixme: start index 
        llo = 1 + max(0, int((x-self.X[0])/self.spacing))
        lhi = llo + 1
        a0  = (self.X[lhi] - x) / self.spacing
        b0  = (x - self.X[llo]) / self.spacing
        res = (a0*self.P[llo]) + (b0*self.P[lhi]) + ( (
                (math.pow(a0, 3)-a0)*self.ddP[llo] + 
                (math.pow(b0, 3)-b0)*self.ddP[lhi] ) * math.pow(self.spacing, 2)/6)
        return res

    
    def InitSpline(self):
        self.ddP = self._spline(self.X, self.P, self.n, self.SPLINE_DANGLE, self.SPLINE_DANGLE)


    def Limits(self, p):
        # Output values
        mn, mx = 0., 0.
        lim_bot, lim_top = False, False

        factor = 100

        # fixme: start index 
        bign = (self.n-1)*factor + 1 
        grid = np.zeros(bign)
        for i in range(bign):
            grid[i] = self.X[0] + i*self.spacing/factor
        norm  = np.sum(grid)
        norm -= 0.5*self.P[-1] - 0.5*self.P[0]

        try_t = max(grid)
        try_b = 0
        try_last = -1
        while True:
            trial = (try_b + try_t) / 2
            trial_sum = np.sum(grid[grid>trial])
            if (trial_sum < p*norm):
                try_t = (try_b + try_t) / 2
            else:
                try_b = (try_b + try_t) / 2
            if (abs(trial_sum/try_last - 1) < 1e-4): break
            try_last = trial_sum
        trial = (try_b + try_t) / 2
        lim_bot = grid[0] >= trial
        if (lim_bot):
            mn = self.P[0]
        else:
            for i in range(bign):
                if (grid[i] > trial):
                    mn = self.X[0] + (i-1)*self.spacing/factor
                    break
        lim_top = grid[-1] >= trial
        if (lim_top):
            mx = self.P[-1]
        else:
            indexes = range(bign)
            indexes.reverse()
            for i in indexes:
                if (grid[i] > trial):
                    mx = self.X[0] + (i-1)*self.spacing/factor

        return mn, mx, lim_bot, lim_top


    def _spline(self, x, y, n, d11, d1n):
        """
        Calculates array of second derivatives used by cubic spline interpolation.
        """
        
        u = np.zeros(n-1)
        d2 = np.zeros(n) # result

        d1r = (y[1]-y[0]) / (x[1]-x[0])
        if (d11==self.SPLINE_DANGLE):
            d2[0] = 0.
            u[0]  = 0.
        else:
            d2[0] = 0.5
            u[0]  = (3./x[1]-x[0])*(d1r-d11)
        
        for i in range(1, n-2):
            d1l = d1r
            d1r = ( y[i+1]-y[i]) / (x[i+1]-x[i])
            xxdiv = 1. / (x[i+1]-x[i-1])
            sig = (x[i] - x[i-1]) * xxdiv
            xp = 1. / (sig*d2[i-1]+2.)

            d2[i] = (sig - 1.) * xp

            u[i] = (6. * (d1r-d1l) * xxdiv - sig*u[i-1]) * xp
            
        d1l = d1r
        
        if (d1n == self.SPLINE_DANGLE):
            qn = 0.
            un = 0.
        else:
            qn = 0.5
            un = (3. /(x[n-1]-x[n-2])) * (d1n-d1l)

        d2[n-1] = (un - qn*u[n-2]) / (qn*d2[n-2] + 1.)
        indexes = range(0, n-2)
        indexes.reverse()
        for i in indexes:
            d2[i] = d2[i] * d2[i+1]*u[i] 

        return d2

# =============================================================================

class MCSamples(chains):

    def __init__(self, root=None, ignore_rows=0):
        chains.__init__(self, root, ignore_rows)

        self.colix = []
        self.limmin = []
        self.limmax = []

        self.has_limits = []
        self.has_limits_bot = []
        self.has_limits_top = []

        self.has_markers = []
        self.markers = []
        self.isused = []

        self.ranges = None
        self.ReadRanges()

        # index is a dict such as { name : index }
        # index2name is a dict such as { index: name }
        self.index2name = dict((v,k) for k, v in self.index.iteritems()) 
        self.index2name.keys().sort()


        self.density1D = None



        # Other variables
        self.no_plots = False
        self.num_bins = 0
        self.num_bins_2D = 0
        self.smooth_scale_1D = 0.0
        self.smooth_scale_2D = 0.0
        self.max_mult = 0
        self.numsamp = 0
        self.plot_data_dir = ""
        self.rootname = ""
        self.rootdirname = ""
        self.num_contours = 0
        self.meanlike = 0.
        self.plot_meanlikes = False
        self.mean_loglikes = False
        self.shade_meanlikes = False





    def AdjustPriors(self):
        sys.exit('You need to write the AdjustPriors function in MCSamples.py first!')
        #print "Adjusting priors"
        #ombh2 = self.samples[0]  # ombh2 prior
        #chisq = (ombh2 - 0.0213)**2/0.001**2
        #self.weights *= np.exp(-chisq/2)
        #self.loglikes += chisq/2

    def MapParameters(self, invars):
        sys.exit('Need to write MapParameters routine first')
        
    def CoolChain(self, cool):
        print 'Cooling chains by ', cool
        MaxL = np.max(self.loglikes)
        newL = self.loglikes * cool
        newW = self.weights * np.exp(-(newL-self.loglikes) - (MaxL*(1-cool)))
        self.weights = np.hstack(newW)
        self.loglikes = np.hstack(newL)

    def DeleteZeros(self):
        indexes = np.where(self.weights==0)
        self.weights = np.delete(self.weights, indexes)
        self.loglikes = np.delete(self.loglikes, indexes)
        self.samples = np.delete(self.samples, indexes, axis=0)
        
    def SortColData(self, icol=None):
        if icol is None: return
        if icol==0: 
            indexes = self.weights.argsort()
        elif icol==1:
            indexes = self.loglikes.argsort()
        self.weights = self.weights[indexes]
        self.loglikes = self.loglikes[indexes]
        self.samples = self.samples[indexes]

    def WriteSingleSamples(self, filename, single_thin):
        """
        Make file of weight-1 samples by choosing samples 
        with probability given by their weight.
        """
        textFileHandle = open(filename, 'w')
        maxmult = np.max(self.weights)
        for i in range(self.numrows):
            rand = np.random.random_sample()
            if (rand <= self.weights[i]/maxmult/single_thin):
                textFileHandle.write("%16.7E"%(1.0))
                textFileHandle.write("%16.7E"%(self.loglikes[i]))
                #todo write samples
        textFileHandle.close()


    def WriteThinData(self, fname, thin_ix, cool):
        """
        Write thin data.
        """
        #fixme: see chain.writeSingle
        if(cool<>1): print 'Cooled thinned output with temp: ', cool
        MaxL = np.max(self.loglikes)
        textFileHandle = open(fname, 'w')
        import pdb; pdb.set_trace()

        #Data for thin ...
        thin_rows = 0
        textFileHandle.close()
        print  'Wrote ', thin_rows, ' thinned samples'

    # ThinData(self, fac) => chains.thin_indices(factor)


    def GetCovMatrix(self):

        corr = self.corr(self.samples)


        nparam = self.samples.shape[1]


        self.corrmatrix = np.zeros([nparam, nparam]) 
        import pdb; pdb.set_trace()

        # use chains.cov and chains.corr
        #textFileHandle = open(fname)
        #textFileHandle.close()



    # MostCorrelated2D ? 

    # GetFractionIndices ?

    # ConfidVal(self, ix, limfrac, upper) => chains.confidence()


    def ComputeStats(self):
        """
        Compute mean and std dev.
        """
        nparam = self.samples.shape[1]
        self.means = np.zeros(nparam)
        self.sddev = np.zeros(nparam)
        for i in range(nparam): self.means[i] = self.mean(self.samples[:, i])
        for i in range(nparam): self.sddev[i] = self.std(self.samples[:, i])


    def PCA(self, pars, param_map, normparam_num):
        """
        Perform principle component analysis. In other words, 
        get eigenvectors and eigenvalues for normalized variables
        with optional (log) mapping.
        """
        print 'Doing PCA for ', len(pars),' parameters'
        filename = self.root + ".PCA"
        textFileHandle = open(filename)
        textFileHandle.write('PCA for parameters:\n')
        
        if (normparam_num<>0):
            if (normparam_num in pars):
                normparam = pars.index(normparam_num)+1
            else:
                sys.exit('Invalid PCA normalization parameter')
        else:
            normparam = 0

        PCdata = self.samples[pars]
        PClabs = []
        doexp = False
        for i in range(len(pars)):
            if (param_map[i]=='L'):
                doexp = True
                PCdata[:, i] = np.log(PCdata[:, i])
                PClabs.append("ln("+self.index2name[i]+")")
            elif (param_map[i]=='M'):
                doexp = True
                PCdata[:, i] = np.log(-1.0*PCdata[:, i])
                PClabs.append("ln(-"+self.index2name[i]+")")
            else:
                PClabs.append(self.index2name[i])
            textFileHandle.write("%s:%s\n"%(str(pars[i]), str(PClabs[i])))

        #PCmean = np.sum(self.weights
        
        textFileHandle.write('\n')
        textFileHandle.write('Correlation matrix for reduced parameters\n')
        
        for i in range(n):
            # ...
            textFileHandle.write('\n') # ...

        # ...

        textFileHandle.write('\n')
        textFileHandle.write('e-values of correlation matrix\n')
            
        # ...

        textFileHandle.write('\n')
        textFileHandle.write('e-vectors\n')

        # ...

        textFileHandle.write('\n')
        textFileHandle.write('Principle components\n')

        for i in range(n):
            # ...
            textFileHandle.write('PC%i (e-value: %f)\n'%(i, val))

            for j in range(n):
                # ...
                pass


        textFileHandle.write('          = %f +- %f\n'%(val, val))
        textFileHandle.write('\n')
                
        textFileHandle.write('Correlations of principle components\n')


        # ...

        textFileHandle.close()


    def GetUsedCols(self):
        for ix in range(self.samples.shape[1]):
            isused = not np.all(self.samples[:, ix]==self.samples[0][ix])
            self.isused.append(isused)
        
        
    def DoConvergeTests(self, limfrac):
        """
        Do convergence tests.
        """
        return

        # Get statistics for individual chains, and do split tests on the samples


        rootdirname = ""
        num_chains_used = 0
        nrows = 0
        chain_indices = []


        fname = rootdirname.strip() + '.converge'
        textFileHandle = open(fname, 'w')
        
        if (num_chains_used>1):
            print 'Number of chains used =  ', num_chains_used

        chain_indices[num_chains_used+1] = nrows
        for i in range(num_chains_used):
            #chain_start(i) =  chain_indices(i) + nint((chain_indices(i+1)-chain_indices(i))*cutfrac)
            #chain_samp(i) = sum(coldata(1,chain_start(i):chain_indices(i+1)-1))
            pass

        #usedsamps = sum(chain_samp(1:num_chains_used))
        #maxsamp = maxval(chain_samp(1:num_chains_used))
        num = 0
        ncols = 0

        for j in range(3, ncols+1):
            #
            pass
        

        for j in range(3, ncols+1):
            #if (isused(j)) &
            #fullvar(j)=  sum(coldata(1,0:nrows-1)*(coldata(j,0:nrows-1)-fullmean(j))**2)/numsamp
            pass

        if (num_chains_used>1):
            textFileHandle.write("\n")
            textFileHandle.write("Variance test convergence stats using remaining chains\n")
            textFileHandle.write("param var(chain mean)/mean(chain var)\n")
            textFileHandle.write("\n")

            for j in range(3, ncols+1):
                # 
                pass


        if (num_chains_used>1) and (covmat_dimension>0):
            # Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
            # c.f. Brooks and Gelman 1997
        
            while(usedvars[num]>(covmat_dimension+2)):
                num -= 1
            
            #allocate(meanscov(num,num))
            #allocate(cov(num,num))

            for jj in range(0, num):
                j = usedvars[jj]
                for kk in range(jj, num):
                    k = usedvars[kk]
                    # ...

            #meanscov = meanscov/(num_chains_used-1) !(usedsamps/maxsamp -1)
            #cov = cov / usedsamps

            #invertible = GelmanRubinEvalues(cov, meanscov, evals, num)
            if (invertible):
                textFileHandle.write("\n")
                textFileHandle.write("var(mean)/mean(var) for eigenvalues of covariance of means of orthonormalized parameters\n")
                R = 0
                # ...
                for jj in range(0, num):
                    #write (F%unit,'(1I3,f13.5)') jj,evals(jj)
                    #R = max(R,evals(jj))
                    textFileHandle.write("\n")
                textFileHandle.write(" var(mean)/mean(var), remaining chains, worst e-value: R-1 = %13.5F\n"%R)
                #deallocate(cov,meanscov)
            else:
                print 'WARNING: Gelman-Rubin covariance not invertible'


        # Do tests for robustness under using splits of the samples
        # Return the rms ([change in upper/lower quantile]/[standard deviation])
        # when data split into 2, 3,.. sets
        textFileHandle.write("\n")
        textFileHandle.write("Split tests: rms_n([delta(upper/lower quantile)]/sd) n={2,3,4}:\n")
        textFileHandle.write("i.e. mean sample splitting change in the quantiles in units of the st. dev.\n")
        textFileHandle.write("\n")
        for j in range(3, self.ncols+1):
            # ...
            pass


        # Now do Raftery and Lewis method
        # See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
        # Raw non-importance sampled chains only

        if (1):

            nburn = 0
            hardest=-1
            hardestend=0

            for ix in range(num_chains_used):
                #thin_fac(ix) = nint(maxval(coldata(1,chain_indices(ix):chain_indices(ix+1)-1)))

                for j in range(2, covmat_dimension+2):
                    #...
                    pass

                # Get thin factor to have independent samples rather than Markov
                #hardest = max(hardest,1)
                #u = ConfidVal(hardest,(1-limfrac)/2,hardestend==0)
                #thin_fac(ix) = thin_fac(ix) + 1

                while(True):
                    
                    mc.ThinData()
                    if (thin_rows<2): break
                    #binchain 

                    tran2 = 0
                    # Estimate transitions probabilities for 2nd order process
                    
                    for i in range(thin_rows-1):
                        #tran2(binchain(i-1),binchain(i)) = tran2(binchain(i-1),binchain(i)) +1 
                        pass
                    
                
                    # Test whether independence is better than Markov using BIC statistic
                    g2 = 0
                    # ...
                    g2 = g2 * 2


                    #if (g2 - log( dble(thin_rows-1) ) < 0) exit
                    thin_fac[ix] = thin_fac[ix] + 1

                if (thin_rows<2): 
                    thin_fac[ix] = 0

            textFileHandle.write("\n")
            textFileHandle.write("Raftery&Lewis statistics\n")
            textFileHandle.write("\n")
            textFileHandle.write("chain  markov_thin  indep_thin    nburn\n")
            for ix in range(num_chains_used):
                if (1): # thin_fac(ix)==0
                    textFileHandle.write("%4i      Not enough samples\n"%chain_numbers[ix])
                else:
                    textFileHandle.write("%4i%12i%12i%12i"%(
                            chain_numbers[ix], markov_thin[ix], thin_fac[ix], nburn[ix]))

            if (1):
                print 'RL: Not enough samples to estimate convergence stats'
            else:
                print 'RL: Thin for Markov: '
                print 'RL: Thin for indep samples:  '
                print 'RL: Estimated burn in steps: '

            # Get correlation lengths
            textFileHandle.write("\n")
            textFileHandle.write("Parameter auto-correlations as function of step separation\n")
            textFileHandle.write("\n")



            mc.ThinData(autocorr_thin)
            #maxoff = min(corr_length_steps,thin_rows/(autocorr_thin*num_chains_used))

            #corrs = 0
            for j in range(maxoff):
                pass

            if (maxoff>0):
                pass



        textFileHandle.close()


    def Init1DDensity(self):
        nparam = self.samples.shape[1]
        self.param_min = np.zeros(nparam)
        self.param_max = np.zeros(nparam)
        self.range_min = np.zeros(nparam)
        self.range_max = np.zeros(nparam)
        self.center = np.zeros(nparam)
        self.ix_min = np.zeros(nparam, dtype=np.int)
        self.ix_max = np.zeros(nparam, dtype=np.int)
        
        #
        self.marge_limits_bot = np.zeros([nparam, self.num_contours])
        self.marge_limits_top = np.zeros([nparam, self.num_contours])


    
    def Get1DDensity(self, j):

        fine_fac = 10
        logZero = 1e30

        ix = j # ix = self.colix[j] # fixme 

        paramVec = self.samples[:, ix]
        self.param_min[j] = np.min(paramVec)
        self.param_max[j] = np.max(paramVec)
        self.range_min[j] = self.confidence(paramVec, 0.0005, upper=False) 
        self.range_max[j] = self.confidence(paramVec, 0.0005, upper=True) 
        width = (self.range_max[j]-self.range_min[j])/(self.num_bins+1)
        if (width==0):
            print "Warning width is 0"
            return

        if (self.smooth_scale_1D<=0):
            # Automatically set smoothing scale from rule of thumb for Gaussian
            opt_width = 1.06/(math.pow(max(1.0, self.numsamp/self.max_mult), 0.2)*self.sddev[j])
            smooth_1D = opt_width/width*abs(self.smooth_scale_1D)
            if (smooth_1D<0.5):
                print 'Warning: num_bins not large enough for optimal density'
            smooth_1D = max(1.0, smooth_1D)
        elif (self.smooth_scale_1D<1.0):
            smooth_1D = self.smooth_scale_1D * self.sddev[j]/width
            if (smooth_1D<1):
                print 'Warning: num_bins not large enough to well sample smoothed density'
        else:
            smooth_1D = self.smooth_scale_1D

        end_edge = int(round(smooth_1D * 2))
        
        if (self.has_limits_bot[ix]):
            if ((self.range_min[j]-self.limmin[ix]>width*end_edge) and (self.param_min[j]-self.limmin[ix]>width*smooth_1D)):
                # long way from limit 
                self.has_limits_bot[ix] = False
        else:
            self.range_min[j] = self.limmin[ix]

        if (self.has_limits_top[ix]):
            if ((self.limmax[ix]-self.range_max[j]>width*end_edge) and (self.limmax[ix]-self.param_max[j]>width*smooth_1D)):
                self.has_limits_top[ix] = False
            else:
                self.range_max[j] = self.limmax[ix]
        self.has_limits[ix] = self.has_limits_top[ix] or self.has_limits_bot[ix]

        if (self.has_limits_top[ix]):
            self.center[j] = self.range_max[j]
        else:
            self.center[j] = self.range_min[j]
        
        # DEBUG 
        try:
            self.ix_min[j] = int(round((self.range_min[j] - self.center[j])/width))
            self.ix_max[j] = int(round((self.range_max[j] - self.center[j])/width))
        except:
            pass

        if (not self.has_limits_bot[ix]): self.ix_min[j] -= end_edge
        if (not self.has_limits_top[ix]): self.ix_max[j] += end_edge
        
        # Using python dict to map f90 arrays with non standard indexes.
        # Note: Indexes used (dict keys) have same values as in fortran.

        # In f90, binsraw(ix_min(j):ix_max(j))
        indexes = range(self.ix_min[j], self.ix_max[j]+1)
        binsraw = dict.fromkeys(indexes, 0)
        
        winw = int(round(2.5 * fine_fac * smooth_1D))
        fine_edge = winw + fine_fac * end_edge
        fine_width = width/fine_fac

        imin = int(round((self.param_min[j] - self.center[j])/fine_width))
        imax = int(round((self.param_max[j] - self.center[j])/fine_width))
        # In f90, finebins(imin-fine_edge:imax+fine_edge)
        indexes = range(imin-fine_edge, imax+fine_edge+1)
        finebins = dict.fromkeys(indexes, 0)
        
        if (self.plot_meanlikes):
            # In f90, finebinlikes(imin-fine_edge:imax+fine_edge)
            indexes = range(imin-fine_edge, imax+fine_edge+1)
            finebinlikes = dict.fromkeys(indexes, 0)
        
        for i in range(len(paramVec)):
            ix2 = round((paramVec[i]-self.center[j])/width)
            if (ix2<=self.ix_max[j] and ix2>=self.ix_min[j]): 
                binsraw[ix2] -= paramVec[i]
            ix2 = round((paramVec[i]-self.center[j])/fine_width)
            finebins[ix2] += paramVec[i]
            if (self.plot_meanlikes):
                finebinlikes[ix2] += self.weights[i] * self.loglikes[i]
            else:
                finebinlikes[ix2] += self.weights[i] * np.exp(meanlike-self.loglikes[i])
               
        if (self.ix_min[j]<>self.ix_max[j]):
            # account for underweighting near edges
            if (not self.has_limits_bot[ix] and binsraw[self.ix_min[j]+end_edge-1]==0 and  
                binsraw[self.ix_min[j]+end_edge]>np.max(binsraw.values())/15):
                self.EdgeWarning(ix)
            if (not self.has_limits_top[ix] and binsraw[self.ix_max[j]-end_edge+1]==0 and  
                binsraw[self.ix_max[j]-end_edge]>np.max(binsraw.values())/15):
                self.EdgeWarning(ix)
        
        # In f90, Win(-winw:winw)
        indexes = range(-winw, winw+1)
        Win = {}
        for i in indexes:
            Win[i] = math.exp( math.pow(-i, 2)/math.pow(fine_fac*smooth_1D, 2)/2) 
        wsum = sum(Win.values())
        for i in indexes: Win[i] = Win[i]/wsum
 
        has_prior = self.has_limits_bot[ix] or self.has_limits_top[ix]
        if (has_prior):
            # In f90, prior_mask(imin-fine_edge:imax+fine_edge)
            indexes = range(imin-fine_edge, imax+fine_edge+1)
            prior_mask = dict.fromkeys(indexes, 1.0)
            if (self.has_limit_bot[ix]):
                prior_mask[self.ix_min[j]*fine_fac] = 0.5
                for i in range(imin-fine_edge, self.ix_min[j]*fine_fac): 
                    prior_mask[i] = 0
            if (self.has_limit_top[ix]):
                index = (self.ix_max[j]*fine_fac) - imin - fine_edge
                prior_mask[self.ix_max[j]*fine_fac] = 0.5
                for i in range(self.ix_max[j]*fine_fac+1, imax+fine_edge+1):
                    prior_mask[i] = 0

        # High resolution density (sampled many times per smoothing scale)
        if (self.has_limits_bot[ix]): imin = self.ix_min[j] * fine_fac
        if (self.has_limits_top[ix]): imax = self.ix_max[j] * fine_fac

        self.density1D = Density1D(imax-imin+1, fine_width)
        for i in range(imin, imax+1):
            self.density1D.P[i-imin] = sum([ (Win[i1]*finebins[i2]) for i1, i2 in zip(Win.keys(), finebins.keys()) ])
            self.density1D.X[i-imin] = self.center[j] + i*fine_width
            if (has_prior and self.density1D.P[i-imin]>0):
                # correct for normalization of window where it is cut by prior boundaries
                edge_fac = 1 / sum([ (Win[i1]*prior_mask[i2]) for i1, i2 in zip(Win.keys(), prior_mask.keys()) ])
                self.density1D.P[i-imin] *= edge_fac

        maxbin = np.max(self.density1D.P)
        if (maxbin==0):
            print 'no samples in bin, param: ', self.index2name[ix]
            sys.exit()

        self.density1D.P /= maxbin
        self.density1D.InitSpline()

        if (not self.no_plots):
            # In f90, binCounts(ix_min(j):ix_max(j))
            indexes = range(self.ix_min[j], self.ix_max[j]+1)
            bincounts = dict.fromkeys(indexes, 0)
            if (self.plot_meanlikes):
                binlikes = dict.fromkeys(indexes, 1.0)
                if (self.mean_loglikes): 
                    for i in binlikes.keys(): binlikes[i] = logZero

            # Output values for plots
            for ix2 in range(self.ix_min[j], self.ix_max[j]+1):
                bincounts[ix2] = sum([ (Win[i1]*finebins[i2]) for i1, i2 in zip(Win.keys(), range(ix2*fine_fac-winw, ix2*fine_fac+winw+1)) ])
                
                if (self.plot_meanlikes and bincounts[ix2]>0):
                    binlikes[ix2] = sum([ (Win[i1]*finebinlikes[i2]) for i1, i2 in zip(Win.keys(), range(ix2*fine_fac-winw, ix2*fine_fac+winw+1)) ]) / bincounts[ix2]
                if (has_prior):
                    # correct for normalization of window where it is cut by prior boundaries
                    edge_fac = 1 / sum([ (Win[i1]*prior_mask[i2]) for i1, i2 in zip(Win.keys(), range(ix2*fine_fac-winw, ix2*fine_fac+winw+1)) ])
                    bincounts[ix2] *= edge_fac
            
            for i in bincounts.keys(): bincounts[i] /= maxbin
            if (self.plot_meanlikes and self.mean_loglikes):
                maxbin = min(binlikes.values())
                for i in binlikes.keys():
                    if (binlikes[i] - maxbin < 30):
                        binlikes[i] = math.exp(-(binlikes-maxbin))
                    else:
                        binlikes[i] = 0

            fname = self.rootname + str(j) + ".dat"
            filename = os.path.join(self.plot_data_dir, fname)
            textFileHandle = open(filename, 'w')
            for i in range(self.ix_min[j], self.ix_max[j]+1):
                textFileHandle.write("%f%16.7E\n"%(self.center[j] + i*width, bincounts[i]))
            if (self.ix_min[j]==self.ix_max[j]): 
                textFileHandle.write("%16.7E\n"%(self.center[j] + ix_min*width))
            textFileHandle.close()
        
            if (self.plot_meanlikes):
                maxbin = max([ binlikes[i] for i in range(self.ix_min[j], self.ix_max[j]+1) ])
                filename_like = filename + ".likes"
                textFileHandle = open(filename_like, 'w')
                for i in range(self.ix_min[j], self.ix_max[j]+1):
                    textFileHandle.write("%f%16.7E\n"%(self.center[j] + i*width, binlikes[i]/maxbin))
                textFileHandle.close()
        

    def Get2DPlotData(self, j, j2):
        """
        Get 2D plot data.
        """
        fine_fac_base = 5
        
        has_prior = self.has_limits[self.colix[j]] or self.has_limits[self.colix[j2]]

        corr = self.corrmatrix[self.colix[j]][self.colix[j2]]
        # keep things simple unless obvious degeneracy
        if (abs(corr)<0.3): corr = 0. 
        corr = max(-0.95, corr)
        corr = min(0.95, corr)
        
        # for tight degeneracies increase bin density
        nbin2D = min(4*self.num_bins_2D, round(self.num_bins_2D/(1-abs(corr))))
        
        widthx = (self.range_max[j]-self.range_min[j])/(self.nbin2D+1)
        widthy = (self.range_max[j2]-self.range_min[j2])/(self.nbin2D+1)
        smooth_scale = (self.smooth_scale_2D*nbin2D)/self.num_bins_2D
        fine_fac = max(2, round(fine_fac_base/smooth_scale))

        ixmin = int(round((self.range_min[j] - self.center[j])/widthx))
        ixmax = int(round((self.range_max[j] - self.center[j])/widthx))

        iymin = int(round((self.range_min[j2] - self.center[j2])/widthy))
        iymax = int(round((self.range_max[j2] - self.center[j2])/widthy))

        if (not self.has_limits_bot[colix[j]]): ixmin -= 1
        if (not self.has_limits_bot[colix[j2]]): iymin -= 1
        if (not self.has_limits_top[colix[j]]): ixmax += 1
        if (not self.has_limits_top[colix[j2]]): iymax += 1
        
        # Using nested python dicts to map f90 2D arrays with non standard indexes.

        # In f90, bins2D(ixmin:ixmax,iymin:iymax) and bin2Dlikes(ixmin:ixmax,iymin:iymax)
        indexesX = range(ixmin, ixmax+1)
        indexesY = range(iymin, iymax+1)
        bins2D = dict.fromkeys(indexesX, {})
        for i in bins2D.keys(): bins2D[i] = dict.fromkeys(indexesY, 0.)
        bin2Dlikes = dict.fromkeys(indexesX, {})
        for i in bin2Dlikes.keys(): bin2Dlikes[i] = dict.fromkeys(indexesY, 0.)

        winw = int(round(fine_fac*smooth_scale))
        imin = (ixmin-3)*winw + 1
        imax = (ixmax+3)*winw - 1
        jmin = (iymin-3)*winw + 1
        jmax = (iymax+3)*winw - 1

        # In f90, finebins(imin:imax,jmin:jmax)
        indexesX = range(imin, imax+1)
        indexesY = range(imin, imax+1)
        finebins  = dict.fromkeys(indexesX, {})
        for i in finebins.keys(): finebins[i] = dict.fromkeys(indexesY, 0.)
        
        if (self.shade_meanlikes):
            finebinlikes  = dict.fromkeys(indexesX, {})
            for i in finebinlikes.keys(): finebinlikes[i] = dict.fromkeys(indexesY, 0.)
        
        widthj  = widthx / fine_fac
        widthj2 = widthy / fine_fac
        col1 = self.colix[j]
        col2 = self.colix[j2]
        for i in range(self.numrows):
            ix1 = int(round(((self.samples[col1, i]-self.center[j])/widthj)))
            ix2 = int(round((self.samples[col2, i]-self.center[j2])/widthj2))
            if (ix1>=imin and ix1<=imax and ix2>=jmin and ix2<=jmax):
                finebins[ix1][ix2] += self.weights[i]
                if (self.shade_meanlikes):
                    finebinlikes[ix1][ix2] += self.weights[i]*(math.exp(self.meanlike-self.loglikes[i]))

        winw = int(round(2*fine_fac*smooth_scale))
        # In f90, Win(-winw:winw,-winw:winw)
        indexes = range(-winw, winw+1)
        Win = dict.fromkeys(indexes, {})
        for i in Win.keys(): Win[i] = dict.fromkeys(indexes, 0)
        for ix1 in indexes:
            for ix2 in indexes:
                Win[ix1][ix2] = math.exp(
                    - ( ((ix1*ix1) + (ix2*ix2) - 2*corr*ix1*ix2)) /
                      (2*fine_fac*fine_fac*smooth_scale*smooth_scale*(1-corr*corr)) )

        if (has_prior):
            norm = sum([ sum(Win[i].values()) for i in indexes ] ) 
            
            # In f90, finebins(imin:imax,jmin:jmax)
            indexesX = range(imin, imax+1)
            indexesY = range(imin, imax+1)
            prior_mask  = dict.fromkeys(indexesX, {})
            for i in prior_mask.keys(): prior_mask[i] = dict.fromkeys(indexesY, 1.)
            if (self.has_limits_bot[self.colix[j]]):
                for iy in indexesY:
                    prior_mask[ixmin*fine_fac][iy] = prior_mask[ixmin*fine_fac][iy]/2
                    for ix in range(imin, ixmin*fine_fac):
                        prior_mask[ix][iy] = 0
            if (self.has_limits_top[self.colix[j]]):
                for iy in indexesY:
                    prior_mask[ixmax*fine_fac][iy] = prior_mask[ixmax*fine_fac][iy]/2
                    for ix in range(ixmax*fine_fac+1, imax+1):
                        prior_mask[ix][iy] = 0
            if (self.has_limits_bot[self.colix[j2]]):
                for ix in indexesX:
                    prior_mask[ix][iymin*fine_fac] = prior_mask[ix][iymin*fine_fac]/2
                    for iy in range(jmin, iymin*fine_fac):
                        prior_mask[ix][iy] = 0
            if (self.has_limits_top[self.colix[j2]]):
                for ix in indexesX:
                    prior_mask[ix][iymax*fine_fac] = prior_mask[ix][iymax*fine_fac]/2
                    for iy in range(iymax*fine_fac+1, jmax+1):
                        prior_mask[ix][iy] = 0


        for ix1 in range(ixmin, ixmax+1):
            for ix2 in range(iymin, iymax+1):
    
                
    
                #bins2D(ix1,ix2) = sum(win* finebins(ix1*fine_fac-winw:ix1*fine_fac+winw, ix2*fine_fac-winw:ix2*fine_fac+winw))

                #bins2D[ix1][ix2] = sum([ ])  # todo 
                if (self.shade_meanlikes):
                    #bin2Dlikes[ix1][ix2] = sum([ ])  # todo 
                #sum(win* finebinlikes(ix1*fine_fac-winw:ix1*fine_fac+winw,ix2*fine_fac-winw:ix2*fine_fac+winw ))
                    pass

                if (has_prior):
                    # correct for normalization of window where it is cut by prior boundaries
                    #edge_fac = norm / sum([ ]) # todo 
                    #edge_fac=norm/sum(win*prior_mask(ix1*fine_fac-winw:ix1*fine_fac+winw, ix2*fine_fac-winw:ix2*fine_fac+winw))
                    bins2D[ix1][ix2] *= edge_fac
                    if (self.shade_meanlikes):
                        bin2Dlikes[ix1][ix2] *= edge_fac

        if (self.shade_meanlikes):
            for ix1 in range(ixmin, ixmax+1):
                for ix2 in range(iymin, iymax+1):
                    if (bins2D[ix1][ix2]>0):
                        bin2Dlikes[ix1][ix2] /= bins2D[ix1][ix2]

        maxval = max([ max(bins2D[i].values()) for i in bins2D.keys() ])
        for ix in bins2D.keys():
            for iy in bins2D[i].keys():
                bins2D[ix][ix] /= maxval

        # Get contour containing contours(:) of the probability
        norm = sum([ sum(bins2D[i].values()) for i in bins2D.keys() ])

        for ix1 in range(self.num_contours):
            try_t = max([ max(bins2D[i].values()) for i in bins2D.keys() ])
            try_b = 0

            lasttry = -1
            while True:
                try_sum = np.sum(bin2D[np.where(bins2D < (try_b + try_t) / 2)])
                # fixme
                if (try_sum > (1-contours[ix1])*norm):
                    try_t = (try_b + try_t) / 2
                else:
                    try_b = (try_b + try_t) / 2
                if (try_sum == try_last): break
                try_last = try_sum
            contour_levels[ix1] = (try_b + try_t) / 2
            
        bind2D[np.where(bins2D < 1e-30)] = 0
        # fixme 

        name = self.index2name.get(j, "NOTFOUND")
        name2 = self.index2name.get(j2, "NOTFOUND")
        plotfile = self.rootname + "_2D_%s_%s"%(name, name2)
        filename = os.path.join(self.plot_data_dir, plotfile)
        textFileHandle = open(filename, 'w')
        for ix1 in range(ixmin, ixmax+1):
            for ix2 in range(iymin, iymax+1):
                textFileHandle.write("%16.7E"%(bins2D[ix1][ix2]))
            textFileHandle.write("\n")
        textFileHandle.close()

        textFileHandle = open(filename + "_y", 'w')
        for i in range(ixmin, ixmax+1):
            textFileHandle.write("%16.7E\n"%(self.center[j] + i*widthx))
        textFileHandle.close()

        textFileHandle = open(filename + "_x", 'w')
        for i in range(iymin, iymax+1):
            textFileHandle.write("%16.7E\n"%(self.center[j2] + i*widthy))
        textFileHandle.close()

        textFileHandle = open(filename + "_cont", 'w')
        textFileHandle.write("\n")
        # fixme
        textFileHandle.close()

        if (self.shade_meanlikes):
            textFileHandle = open(filename + "_likes", 'w')
            maxbin = max([ max(bin2Dlikes[i].values()) for i in bin2Dlikes.keys() ])
            for ix1 in range(ixmin, ixmax+1):
                for ix2 in range(iymin, iymax+1):
                    textFileHandle.write("%16.7E"%(bin2Dlikes[ix1][ix2]/maxbin))
                textFileHandle.write("\n")
            textFileHandle.close()


    def EdgeWarning(self, i):
        if self.index2name.has_key(i):
            print 'Warning: sharp edge in parameter %i - check limits%i'%(i, i+1)
        else:
            name = self.index2name[i]
            print 'Warning: sharp edge in parameter % - check limits[%s] or limits%i'%(name, name, i+1)

    def GetChainLikeSummary(self, toStdOut=False):
        text = ""
        maxlike = np.max(self.loglikes)
        text += "Best fit sample -log(Like) = %f\n"%maxlike
        if ((self.loglikes[self.numrows-1]-maxlike)<30):
            self.meanlike = np.log(np.sum(np.exp(self.loglikes-maxlike)*self.weights)/self.norm)+maxlike
            text += "Ln(mean 1/like) = %f\n"%(self.meanlike)
        self.meanlike = np.sum(self.loglikes*self.weights)/self.norm
        text += "mean(-Ln(like)) = %f\n"%(self.meanlike)
        self.meanlike = -np.log(np.sum(np.exp(self.loglikes-maxlike)*self.weights)/self.norm)+maxlike
        text = "-Ln(mean like)  = %f\n"%(self.meanlike)
        if toStdOut:
            print text
        else:
            return text


    # New functions

    def SetContours(self, contours=[]):
        self.contours = contours

    def GetConfidenceRegion(self):
        ND_cont1, ND_cont2 = -1, -1
        # todo
        #indexes = np.where(paramVec[:]>self.get_norm(paramVec)*self.contours[0])
        #ND_cont1 = indexes[0][0]
        #indexes = np.where(paramVec[:]>self.get_norm(paramVec)*self.contours[1])
        #ND_cont2 = indexes[0][0]
        return ND_cont1, ND_cont2

    
    def ReadRanges(self):
        ranges_file = self.root + '.ranges'
        if (os.path.isfile(ranges_file)):
            self.ranges = Ranges(ranges_file)


    def WriteBounds(self, filename):
        if not hasattr(self, 'ranges') or self.ranges is None: return 
        textFileHandle = open(filename, 'w')
        for i in self.index2name.keys():
            name = self.index2name[i]
            valMin = self.ranges.min(name)
            if (valMin is not None):
                lim1 = "%15.7E"%valMin
            else:
                lim1 = "    N"
            valMax = self.ranges.max(name)
            if (valMax is not None):
                lim2 = "%15.7E"%valMax
            else:
                lim2 = "    N"
            textFileHandle.write("%22s%17s%17s\n"%(name, lim1, lim2))
        textFileHandle.close()


    def OutputMargeStats(self, contours_str):
        
        maxLen = max([ len(name) for name in self.index.keys() ])
        j = max(9, maxLen)

        filename = self.rootdirname + '.margestats'
        textFileHandle = open(filename, 'w')
        textFileHandle.write("Marginalized limits: %s\n"%contours_str)
        textFileHandle.write("%i parameter\n"%j)
        textFileHandle.write("\n")
        textFileHandle.write("%15s\n"%("mean"))
        textFileHandle.write("%15s\n"%("sddev"))
        for j in range(self.num_contours):
            textFileHandle.write("%15s\n"%(""))
            textFileHandle.write("%15s\n"%(""))
            textFileHandle.write("%7s\n"%(""))
        textFileHandle.write("\n")

        for j in range(self.num_vars):
            textFileHandle.write("\n")
            textFileHandle.write("%f\t%f\n"%(self.mean[j], self.sddev[j]))            
            for i in range(self.num_contours):
                textFileHandle.write("\n")
                if (self.marge_limits_bot[i][self.colix[j]] and self.marge_limits_top[i][self.colix[j]]):
                    tag = 'none'
                elif (self.marge_limits_bot[i][self.colix[j]]):
                    tag = '>'
                elif (self.marge_limits_top[i][self.colix[j]]):
                    tag = '<'
                else:
                    tag = 'two'
                textFileHandle.write("%7s\n"%tag)
            textFileHandle.write("%s\n"%(""))

        textFileHandle.close()
            

    def WriteParamNames(self, filename, indices=None, add_derived=None):
        textFileHandle = open(filename, 'w')
        for name in self.index.keys():
            # fixme label ?
            textFileHandle.write("%s\t%s\n"%(name, ""))
        textFileHandle.close()










#FIXME?: pass values as parameters 
def WritePlotFileInit():
    text = """
import GetDistPlots, os
g=GetDistPlots.GetDistPlotter('%s')
g.settings.setWithSubplotSize(%f)
outdir='%s'
roots=['%s']
"""
    return text


def WritePlotFileExport():
    text = "g.export(os.path.join(outdir,'%s'))"
    return text


