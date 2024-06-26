import numpy,math,json
from matplotlib import pyplot
from matplotlib.collections import LineCollection
from .Sampling import Sampler
from .Sampling import Distribution

class GaussianProcess:
    def __init__(self):
        self.parameters = None
        self.kernel = None
        self.samples = None
        self.means = None
        self.queries = []
        self.weights = None
        self.cholesky = None

        self.constrained = False
        self.gaussian_constraints = False

        self.kernel_functions = {"rbf":lambda parameters,t1,t2 : (parameters[0]**2)*numpy.exp(-0.5*((numpy.abs(t1-t2)/parameters[1])**2))}
    @property
    def flat_query_locations(self):
        query_locations = []
        for query_group in self.queries:
            query_locations += [query.location for query in query_group]
        return numpy.array(query_locations)
    @property
    def query_locations(self):
        query_locations = []
        for query_group in self.queries:
            local_query = []
            for query in query_group:
                local_query += [query.location]
            query_locations += [local_query]            
        return query_locations
    
    @property
    def flat_means(self):
        to_concatenate = []
        for mean_group in self.means:
            if len(mean_group)==1:
                to_concatenate += [mean_group[0]]
            else:
                to_concatenate += [numpy.squeeze(mean_group)]
        return numpy.concatenate(to_concatenate)
        # return numpy.concatenate([numpy.squeeze(mean_group) for mean_group in self.means])

    def flatten(self,values):
        to_concatenate = []
        for group in values:
            if len(group)==1:
                to_concatenate += [group[0]]
            else:
                to_concatenate += [numpy.squeeze(group)]
        return numpy.concatenate(to_concatenate)
    def split(self,array):
        query_indices = [0]
        #array = numpy.transpose(array)

        for query_length in self._query_lengths:
            query_indices += [query_indices[-1]+query_length]
        output = []
        for query_length_index,query_length in enumerate(query_indices[1:]):
            output += [array[:,query_indices[query_length_index]:query_indices[query_length_index+1]]]
        output += [array[:,query_indices[-1]:-1]]
        return [array for array in output if array.size>0]
        
    def generateSeed(self):
        [group for group in self.query_locations]
        return [numpy.random.normal(loc=0.0,scale=1.0,size=len(group)) for group in self.query_locations]
    def perturbSeed(self,seed,standard_deviation):
        return [numpy.random.normal(loc=(seed)/(1+standard_deviation**2),scale=numpy.sqrt((standard_deviation**2)/(1+standard_deviation**2))) for seed,standard_deviation in zip(seed,standard_deviation,strict=True)]
    def constrain(self,constraints):
        output = GaussianProcess()
        output.constrained = True
        output.gaussian_constraints = numpy.all([constraint.type=="gaussian" or constraint.type=="loggaussian" for constraint in constraints])

        output.constraints = constraints
        return output
    def setQueryLocations(self,locations,bin_edges=None):
        if bin_edges is None:
            if self.constrained:
                bin_edges = self.constraints[0].bin_edges
            elif self.queries:
                bin_edges = self.queries[0][0].bin_edges
            else:
                raise(ValueError("Unable to determine bin edges"))
        self.queries = []
        self._query_lengths = [len(query_group) for query_group in locations]
        for query_location_group in locations:
            local_queries = []
            for query_location in query_location_group:
                local_queries += [Sampler(bin_edges,"Manual",None,"latin_hypercube",location=query_location)]
            self.queries += [local_queries]
        return self      
    def query(self,locations=None):
        if locations is None:
            if self.queries:
                locations = self.query_locations
            elif self.constrained:
                locations = [numpy.array([constraint.location for constraint in self.constraints])]
        elif not isinstance(locations[0],(list,tuple,numpy.ndarray)):
            locations = [locations]
        self.setQueryLocations(locations)
        if not self.constrained:
            self.means = [numpy.ones(len(query_group))*self.specified_mean for query_group in self.queries]
            self.covariances = self.kernelWrapper(self.parameters,self.flat_query_locations)
        else:
            if not self.gaussian_constraints:
                query_locations = self.query_locations
                constraint_approximations = [constraint.approximateGaussian() for constraint in self.constraints]
                approximate_gp = GaussianProcess().constrain(constraint_approximations).setKernel(self.kernel,self.parameters,self.inflation,self.specified_mean).query(query_locations)
                return approximate_gp
            else:
                query_locations = self.flat_query_locations

                constraints_mean = numpy.array([constraint.mean for constraint in self.constraints])
                constraint_locations = numpy.array([constraint.location for constraint in self.constraints])

                self.constraint_covariance = self.kernelWrapper(self.parameters,constraint_locations)
                self.cross_covariance = self.kernelWrapper(self.parameters,constraint_locations[:,numpy.newaxis],query_locations[numpy.newaxis,:])
                self.query_covariance = self.kernelWrapper(self.parameters,query_locations)
                inflation = self.inflation
                if len(inflation)==1:
                    inflation = [self.inflation[0] for constraint in self.constraints]

                constraint_uncertainty = numpy.array([inflation*constraint.standard_deviation**2 for constraint,inflation in zip(self.constraints,inflation)]) * numpy.eye(self.constraint_covariance.shape[0])

                matrix_inverse = numpy.linalg.inv(self.constraint_covariance+constraint_uncertainty)
                means = numpy.array(self.specified_mean + numpy.matmul(numpy.matmul(numpy.transpose(self.cross_covariance),matrix_inverse),constraints_mean-self.specified_mean))[numpy.newaxis]
                
                self.means = self.split(means)
                self.covariances = self.query_covariance - numpy.matmul(numpy.matmul(numpy.transpose(self.cross_covariance),matrix_inverse),self.cross_covariance)
        
        return self
    def getBounds(self,bounds,number,logspace):
        if not logspace:
            return numpy.linspace(bounds[0],bounds[1],number[0])
        else:
            return numpy.logspace(numpy.log10(bounds[0]),numpy.log10(bounds[1]),number[0])
    def quantile(self,value,group=None):
        if group is None:
            return numpy.array([query.quantile((value,)) for query_group in self.queries for query in query_group])
        else:
            return numpy.array([query.quantile((value,)) for query in self.queries[group]])

    def getCholesky(self):
        if self.cholesky is None:
            self.cholesky = numpy.linalg.cholesky(self.covariances+1e-6*numpy.eye(numpy.shape(self.covariances)[0]))
        return self.cholesky
    # Kernel
    def kernelWrapper(self,parameters,t1,t2=None):
        if t2 is None:
            if t1.ndim==1:
                t2 = t1[:,numpy.newaxis]
            elif t2 is None and t1.ndim==2:
                t1 = t1.reshape((t1.ndim[0],1,t1.ndim[1]))
                t2 = t1.reshape((1,t1.ndim[0],t1.ndim[1]))
            else:
                t1 = t1.reshape((numpy.product(t1.ndim[:-1]),1,t1.ndim[-1]))
                t2 = t1.reshape((1,numpy.product(t1.ndim[:-1]),t1.ndim[-1]))
        return self.kernel_function(parameters,t1,t2)
    def setKernel(self,kernel,parameters,inflation=None,specified_mean=False):
        output = GaussianProcess()
        if self.constrained:
            output = output.constrain(self.constraints)
        
        output.kernel = kernel
        output.kernel_function = output.kernel_functions[output.kernel]
        output.parameters = parameters
        output.inflation = inflation
        if output.inflation is None:
            if output.constrained:
                output.inflation = [1 for constraint in output.constraints]
            else:
                output.inflation = 1

        if not output.constrained and specified_mean is False:
            output.specified_mean = 0
        else:
            if specified_mean is False:
                output.specified_mean = numpy.mean(numpy.array([constraint.mean for constraint in output.constraints]))
            else:
                output.specified_mean = specified_mean
        return output

    # Sample Manipulation
    def getNextSampleByNumber(self,samples=None):
        if samples is None:
            samples = self.samples
        for sample in zip(*samples,strict=True):
            yield sample
    def getNextSampleByQuery(self,samples=None):
        if samples is None:
            samples = self.samples
        for sample_group in samples:
            for samples in numpy.transpose(sample_group):
                yield [sample for sample in samples]
    def getSamples(self,number_of_samples,seed=None):
        self.number_of_samples = number_of_samples
        if seed is not None:
            samples = numpy.transpose(numpy.transpose(self.flat_means[numpy.newaxis])+numpy.matmul(self.getCholesky(),numpy.transpose(self.flatten(seed)[numpy.newaxis])))
        else:
            samples = numpy.random.multivariate_normal(numpy.squeeze(self.flat_means),self.covariances,number_of_samples)
        split_samples = self.split(samples)
        self.assignSamples(split_samples)
        return self
    def redoMean(self):
        self.means = []
        for query_group in self.queries:
            self.means += [[query.mean for query in query_group]]
        return self
    def fromMCMCSamples(self,mcmc_samples):
        self.number_of_samples = len(mcmc_samples)
        samples = []
        for sample_group_index in range(len(mcmc_samples[0])):
            samples += [numpy.stack([mcmc_sample[sample_group_index] for mcmc_sample in mcmc_samples])]
        self.assignSamples(samples)
        self.means = [numpy.mean(sample,axis=0) for sample in self.samples]
        
        return self
    def assignSamples(self,samples,weights=None):
        self.samples = samples
        self.number_of_samples = numpy.shape(samples[0])[0]
        
        bin_edges = None
        if self.constrained:
            bin_edges = self.constraints[0].bin_edges
        elif self.queries[0][0].bin_edges is not None:
            bin_edges = self.queries[0][0].bin_edges


        sample_generator = self.getNextSampleByQuery(samples)

        global_queries = []
        for query_group in self.queries:
            local_queries = []
            for query in query_group:
                sample = next(sample_generator)

                if weights is None:
                    local_queries += [Distribution.fromSamples(sample,bin_edges,location=query.location).normalise()]
                else:
                    self.weights = weights
                    local_queries += [Distribution.fromSamples(sample,bin_edges,location=query.location,weights=weights).normalise()]
            global_queries += [local_queries]
        self.queries = global_queries
        return self
    def getSampleLikelihood(self,index=None,keep_separate=False,logspace=True):
        if index is None:
            samples = self.samples[0]
        else:
            samples = self.samples[index]
        if logspace:
            log_probabilities = numpy.empty(shape=samples.shape)
            for index,(query_samples,constraint) in enumerate(zip(numpy.transpose(samples),self.constraints,strict=True)):
                log_probabilities[:,index] = constraint.getLogProbability(query_samples)
            sum_log_probabilities = numpy.sum(log_probabilities,axis=1)
            if keep_separate:
                return sum_log_probabilities
            return numpy.sum(sum_log_probabilities)
        else:
            probabilities = numpy.empty(shape=samples.shape)
            for index,(query_samples,constraint) in enumerate(zip(numpy.transpose(samples),self.constraints,strict=True)):
                probabilities[:,index] = constraint.getProbability(query_samples)
            product_probabilities = numpy.prod(probabilities,axis=1)
            if keep_separate:
                return product_probabilities
            return numpy.prod(product_probabilities)
    def splitconstraints(self,number_of_folds):
        if number_of_folds>0:
            local_constraints = numpy.array(self.constraints)
            numpy.random.shuffle(local_constraints)

            number_of_constraints = len(self.constraints)
            
            number_in_fold = number_of_constraints//number_of_folds
            number_leftover = number_of_constraints%number_of_folds

            fold_indices = numpy.concatenate([numpy.repeat(numpy.arange(0,number_of_folds),number_in_fold),numpy.arange(0,number_leftover)])

            fold = [None]*number_of_folds
            for fold_index in fold_indices:
                fold[fold_index] = numpy.array(local_constraints[fold_indices==fold_index])
            return fold
        else:
            return self.constraints
    def log2Samples(self):
        samples = []
        for sample_group in self.samples:
            local_samples = []
            for sample in sample_group:
                local_samples += [numpy.log2(sample)]
            samples += [local_samples]
        
        self.assignSamples(samples)

        return self
    def pow2Samples(self):
        samples = []
        for sample_group in self.samples:
            local_samples = []
            for sample in sample_group:
                local_samples += [2**(sample)]
            samples += [local_samples]
        
        self.assignSamples(samples)

        return self
    def round(self,precision):
        samples = []
        for sample_group in self.samples:
            local_samples = []
            for sample in sample_group:
                local_samples += [numpy.round(sample,precision)]
            samples += [local_samples]
        self.assignSamples(samples)
        return self

    # Local Mean
    def removeLocalMean(self,fraction=(2,3),parameters=None):        
        constraints_mean = numpy.array([constraint.mean for constraint in self.constraints])
        constaints_location = numpy.array([constraint.location for constraint in self.constraints])

        self.yrange = numpy.max(constraints_mean)-numpy.min(constraints_mean)
        self.xrange = numpy.max(constaints_location)-numpy.min(constaints_location)

        mean_gp = GaussianProcess().constrain(self.constraints)
        if parameters is not None:
            mean_gp = mean_gp.setKernel("rbf",parameters,specified_mean=False).query()
        else:
            mean_gp = mean_gp.setKernel("rbf",(self.yrange/fraction[0],self.xrange/fraction[1]),specified_mean=False).query()        
        
        shifted_constraints = [constraint.shift(-current_mean) for constraint,current_mean in zip(self.constraints,mean_gp.flat_means)]
        self.constraints = shifted_constraints

        self.local_mean = mean_gp

        return (self,mean_gp)
    def addLocalMean(self,mean_gp=None,interpolation_locations=None):
        if mean_gp is None:
            mean_gp = self.local_mean
        if interpolation_locations is None:
            interpolation_locations = self.query_locations
        shifted_constraints = [constraint.shift(current_mean) for constraint,current_mean in zip(self.constraints,mean_gp.flat_means)]
        self.constraints = shifted_constraints

        if interpolation_locations is not None:
            mean_gp = mean_gp.query(interpolation_locations)
            for query_means,mean_means in zip(self.means,mean_gp.means):
                query_means += mean_means
                if self.samples is not None:
                    self.samples += mean_gp.flat_means
        
        return self
    
    # Fit
    def fit(self,method,kernel,bounds,kfold=0,number=(10,10),levels=5,specified_mean=False,logspace=False,start=None,iterations=10,tolerance=None):
        if method=="grid":
            parameter_1_values = self.getBounds(bounds[0],number,logspace)
            parameter_2_values = self.getBounds(bounds[1],number,logspace)
            
            if kfold>1:
                likelihood = numpy.empty((number[0],number[1],kfold))
            else:
                likelihood = numpy.empty((number[0],number[1]))

            split_constraints = self.splitconstraints(kfold)
            current_best_likelihood = None

            for index_1,parameter_1 in enumerate(parameter_1_values):
                for index_2,parameter_2 in enumerate(parameter_2_values):                    
                    if kfold>1:
                        for fold_index in numpy.arange(0,kfold):
                            current_constraints = numpy.concatenate([constraints for index,constraints in enumerate(split_constraints) if index!=fold_index])
                            current_queries = split_constraints[fold_index]
                            current_query_locations = numpy.array([query.location for query in current_queries])
                            gp = GaussianProcess().constrain(current_constraints)

                            gp = gp.setKernel(kernel,(parameter_1,parameter_2),specified_mean=specified_mean).query(current_query_locations).getSamples(1000)

                            query_gp = GaussianProcess().constrain(current_queries)
                            query_gp.samples = gp.samples
                            likelihood[index_1,index_2,fold_index] = query_gp.getSampleLikelihood(logspace=True,keep_separate=False)

                            a = 5

                        mean_current_likelihood = numpy.median(likelihood[index_1,index_2])
                        if current_best_likelihood is None or mean_current_likelihood>current_best_likelihood:
                            current_best_likelihood = mean_current_likelihood
                            best_parameters = (parameter_1,parameter_2)
                    else:
                        gp = GaussianProcess().constrain(split_constraints).setKernel(kernel,(parameter_1,parameter_2),specified_mean=specified_mean).getSamples(100)

                        likelihood[index_1,index_2] = gp.getSampleLikelihood()
                        
                        if likelihood[index_1,index_2] == numpy.min(likelihood):
                            best_parameters = (parameter_1,parameter_2)

            output_gp = GaussianProcess().constrain(self.constraints).setKernel(kernel,best_parameters,specified_mean=specified_mean)
            return output_gp
        elif method=="multigrid":
            for level in range(levels):
                if level==0:
                    new_bounds = bounds
                parameters_estimated = GaussianProcess().constrain(self.constraints).fit("grid",kernel,new_bounds,kfold,number,logspace=logspace)
                
                old_bounds_1 = self.getBounds(new_bounds[0],number,logspace)
                old_bounds_2 = self.getBounds(new_bounds[1],number,logspace)

                old_bounds_1_index = numpy.where(old_bounds_1==parameters_estimated.parameters[0])[0]
                old_bounds_2_index = numpy.where(old_bounds_2==parameters_estimated.parameters[1])[0]

                if old_bounds_1_index==0:
                    new_bounds_1 = (old_bounds_1[0],old_bounds_1[1])
                elif old_bounds_1_index==len(old_bounds_1)-1:
                    new_bounds_1 = (old_bounds_1[-2],old_bounds_1[-1])
                else:
                    new_bounds_1 = (old_bounds_1[old_bounds_1_index-1][0],old_bounds_1[old_bounds_1_index+1][0])
                
                if old_bounds_2_index==0:
                    new_bounds_2 = (old_bounds_2[0],old_bounds_2[1])
                elif old_bounds_2_index==len(old_bounds_2)-1:
                    new_bounds_2 = (old_bounds_2[-2],old_bounds_2[-1])
                else:
                    new_bounds_2 = (old_bounds_2[old_bounds_2_index-1][0],old_bounds_2[old_bounds_2_index+1][0])
                
                new_bounds = (new_bounds_1,new_bounds_2)            
            parameters_estimated = GaussianProcess().constrain(self.constraints).fit("grid",kernel,new_bounds,kfold,number,logspace=logspace,levels=levels)   
            output_gp = GaussianProcess().constrain(self.constraints).setKernel(kernel,parameters_estimated.parameters,specified_mean=specified_mean)
            return output_gp
        elif method=="gradient":
            if start is None and bounds is not None:
                start = ((bounds[0][0]+bounds[0][1])/2,(bounds[0][1]+bounds[1][1])/2)
            for iteration in range(iterations):
                if iteration==0:
                    centre = start
                
                for fold in range(kfold):
                    current_constraints = numpy.concatenate([constraints for index,constraints in enumerate(split_constraints) if index!=fold])
                    current_queries = split_constraints[fold]

                    gp =  GaussianProcess().constrain(self.observations)
                    gp.setKernel("rbf",centre)

                    gp.getSamples(1000)

                    query_gp = GaussianProcess().constrain(current_queries)
                    query_gp.samples = gp.samples
                    central_likelihood = query_gp.getSampleLikelihood()

                    side_gp = GaussianProcess().constrain(self.observations)
                    side_gp.setKernel("rbf",(1.1*centre[0],centre[1]))

                    side_gp.getSamples(1000)

                    query_gp = GaussianProcess().constrain(current_queries)
                    query_gp.samples = side_gp.samples
                    side_likelihood = query_gp.getSampleLikelihood()

                    top_gp = GaussianProcess().constrain(self.observations)
                    top_gp.setKernel("rbf",(centre[0],1.1*centre[1]))

                    top_gp.getSamples(1000)

                    query_gp = GaussianProcess().constrain(current_queries)
                    query_gp.samples = top_gp.samples
                    top_likelihood = query_gp.getSampleLikelihood()

                    a=5

    # Display
    def plotMean(self,axis=None,group=None,**kwargs):        
        if axis is None:
            axis = pyplot.gca()
        if group is None:
            query_locations = self.flat_query_locations
            means = self.flat_means
        else:
            query_locations = self.query_locations[group]
            means = numpy.squeeze(self.means[group])
        axis.plot(query_locations,means,**kwargs)
    def plotMedian(self,axis=None,group=None,**kwargs):        
        if axis is None:
            axis = pyplot.gca()
        if group is None:
            query_locations = self.flat_query_locations
            medians = numpy.array([query.quantile([0.5]) for query_group in self.queries for query in query_group])
        else:
            query_locations = self.query_locations[group]
            medians = numpy.array([query.quantile([0.5]) for query in self.queries[group]])
        axis.plot(query_locations,medians,**kwargs)
    def plotSamples(self,group=None,indices=None,axis=None,scatter=False,**kwargs):
        if self.queries:
            if group is not None:
                locations = self.query_locations[group]
            else:
                locations = self.flat_query_locations
        elif self.constrained:
            locations = numpy.array([constraint.location for constraint in self.constraints])

        if indices is None:
            indices = range(0,self.number_of_samples,1)
        if axis is None:
            axis = pyplot.gca()

        if group is not None:
            sample_generator = self.getNextSampleByNumber([self.samples[group]])
        else:
            sample_generator = self.getNextSampleByNumber()
        count = -1
        for sample_index in indices:
            while count<sample_index:
                count += 1
                sample = next(sample_generator)
            if not scatter:
                axis.plot(locations,numpy.squeeze(sample),**kwargs)
            elif scatter:
                axis.scatter(locations,numpy.squeeze(sample),**kwargs)
    def plotConstraints(self,indices=None,axis=None,flipX=False,fmt='err',quantiles=(0.025,0.975),**kwargs):
        if indices is None:
            indices = range(0,len(self.constraints))
        if axis is None:
            axis = pyplot.gca()
        for index in indices:
            if fmt=="err":
                axis.plot([self.constraints[index].location,self.constraints[index].location],self.constraints[index].quantile(quantiles),**kwargs)
            if fmt=="segment":
                edges_one = numpy.vstack((numpy.repeat(self.constraints[index].location,100),self.constraints[index].quantile(numpy.linspace(0.01,0.98,100)))).T
                edges_two = numpy.vstack((numpy.repeat(self.constraints[index].location,100),self.constraints[index].quantile(numpy.linspace(0.02,0.99,100)))).T
                segments = [(edge_one,edge_two) for edge_one,edge_two in zip(edges_one,edges_two,strict=True)]
                alphas = 0+1.0*numpy.sin(numpy.linspace(0.01,0.98,100)*numpy.pi)
                line_segments = LineCollection(segments,alpha=alphas,**kwargs)
                axis.add_collection(line_segments)
            else:
                axis.scatter(self.constraints[index].location,self.constraints[index].mean,**kwargs)
        if flipX:
            axis.invert_xaxis()
    def plotQueries(self,axis=None,**kwargs):
        if axis is None:
            axis = pyplot.gca()
        for index in range(0,len(self.queries)):
            axis.errorbar(self.queries[index].location,self.queries[index].mean,yerr=self.queries[index].standard_deviation*2,fmt="o",**kwargs)
    def plotArea(self,axis=None,group=0,set_axes=True,**kwargs):
        import matplotlib.patches as patches
        if axis is None:
            axis = pyplot.gca()
        query_locations = self.query_locations[group]
        boundaries = numpy.transpose(numpy.array([query.quantile((0.025,0.975)) for query in self.queries[group]]))
        boundaries[1] = numpy.flip(boundaries[1])
        boundaries_flat = boundaries.flatten()
        queries_flat = numpy.concatenate((query_locations,numpy.flip(query_locations)))
        patch = patches.Polygon(numpy.concatenate((queries_flat[:,numpy.newaxis],boundaries_flat[:,numpy.newaxis]),axis=1),**kwargs)
        axis.add_patch(patch)
        if set_axes:
            axis.set_xlim((numpy.min(query_locations),numpy.max(query_locations)))
            axis.set_ylim((numpy.min(boundaries),numpy.max(boundaries)))
    def plotPcolor(self,axis=None,colourbar=True,invert_x=False,map=None,mask=False,**kwargs):
        if axis is None:
            axis = pyplot.gca()
        x = self.queries[-1][0].bin_midpoints
        locations = numpy.array([query.location for query in self.queries[-1]])
        probabilities = numpy.transpose(numpy.array([query.probabilities for query in self.queries[-1]]))

        if map is not None:
            pyplot.set_cmap(map)
        mappable = axis.pcolormesh(locations,x,probabilities,**kwargs)
        if invert_x:
            axis.invert_xaxis()
        if colourbar:
            pyplot.colorbar(mappable)
        if mask and "vmin" in kwargs:
            mappable.cmap.set_under(alpha=0)

        a = 5

    # Output
    def toJSON(self,filename):
        # Want kernel function, parameters, means, samples
        json_data = json.dumps(self,cls=GPEncoder,indent=4)
        json_data_stripped = json_data.replace('"xxx',"").replace('xxx"',"").replace('xxx',"")
        with open(filename,"w") as file:
            file.write(json_data_stripped)
    def fromJSON(self,filename):
        with open(filename,"r") as file:
            raw_data = file.read()
        json_data = json.loads(raw_data)
        self.setQueryLocations(json_data["locations"],bin_edges=numpy.array(json_data["edges"]))
        if "samples" in json_data:
            self.fromMCMCSamples(json_data["samples"])

        return self

class GPEncoder(json.JSONEncoder):
    def listOfArraysToString(self,list_of_arrays):
        return "xxx["+", ".join([self.arrayToString(group) for group in list_of_arrays])+"]xxx"
    def arrayToString(self,array):
        if numpy.squeeze(array).size>1:
            return "["+", ".join([str(value) for value in numpy.squeeze(array)])+"]"
        else:
            return "["+str(numpy.squeeze(array))+"]"
    def toStr(self,array):
        if isinstance(array[0],int) or isinstance(array[0],float):
            return "xxx["+", ".join(str(x) for x in array)+"]xxx"
        else:
            return "xxx["+", ".join(str(x) for group in array for x in group)+"]xxx"
    def default(self,obj):
        if isinstance(obj,GaussianProcess):
            output = {}
            if obj.kernel is not None:
                output["kernel"] = obj.kernel
                output["parameters"] = self.toStr(obj.parameters)
            if obj.means is not None:
                output["means"] = self.listOfArraysToString(obj.means)
            if obj.queries is not None:
                output["locations"] = self.listOfArraysToString(obj.query_locations)
                output["edges"] = self.toStr(obj.queries[0][0].bin_edges)
            if obj.means is not None:
                output["means"] = self.listOfArraysToString(obj.means)
            if obj.samples is not None:
                output["samples"] = []
                sample_generator = obj.getNextSampleByNumber()

                for sample in sample_generator:
                    output["samples"] += [self.listOfArraysToString(sample)]
            
            if obj.weights is not None:
                output["weights"] = self.toStr(obj.weights)

            return output
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self,obj)

