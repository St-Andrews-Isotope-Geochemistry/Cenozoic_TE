import numpy,copy,math,warnings,json,numbers
from matplotlib import pyplot

class Distribution:
    @property
    def bin_midpoints(self):
        return numpy.array((self.bin_edges[0:-1]+self.bin_edges[1:])/2)
    @property
    def mean(self):
        self.checkNormalised()
        return numpy.sum(self.bin_midpoints*self.probabilities)
    @property
    def standard_deviation(self):
        self.checkNormalised()
        return numpy.sqrt(numpy.sum((self.bin_midpoints-self.mean)**2 * self.probabilities))

    def approximateGaussian(self,inflation=1):
        gaussian = Distribution(self.bin_edges,"Gaussian",(self.mean,self.standard_deviation*inflation),location=self.location).normalise()
        return gaussian
    @property
    def equally_spaced(self):
        first_difference = numpy.diff(self.bin_edges)
        if numpy.all(numpy.abs(first_difference-first_difference[0])<1e-10):
            return True
        return False

    def __init__(self,bin_edges,type,values,location=None):
        self.type = type.lower()
        self.location = location
        if self.type!="manual":
            self.values = values
        self.bin_edges = numpy.array(bin_edges)
        if self.type=="flat":
            assert len(values)==2,"Number of values must be 2"
            assert values[0]<values[1],"First value must be less than second"

            probabilities = numpy.zeros(len(self.bin_edges)-1)
            probabilities[numpy.logical_and(self.bin_midpoints>=values[0],self.bin_midpoints<=values[1])] = 1
            self.probabilities = probabilities
        elif self.type=="gaussian" or self.type=="normal":
            assert len(values)==2,"Number of values must be 2"

            mu = values[0]
            sigma = values[1]

            gaussian = (1/(sigma*numpy.sqrt(2*numpy.pi))) * numpy.exp(-0.5*(((self.bin_midpoints-mu)/sigma)**2))
            self.probabilities = gaussian
        elif self.type=="loggaussian":
            mu = values[0]
            sigma = values[1]

            loggaussian = (1/(self.bin_midpoints*sigma*numpy.sqrt(2*numpy.pi))) * numpy.exp(-(((numpy.log(self.bin_midpoints)-mu)**2)/((2*sigma)**2)))
            self.probabilities = loggaussian
        elif self.type=="subsample":
            pass
        elif self.type=="manual":
            if values is not None:
                assert len(bin_edges)==len(values)+1
            self.bin_edges = bin_edges
            self.probabilities = values
        else:
            ValueError("Unknown distribution type")
    
    def normalise(self,value=None,by=None):
        assert numpy.sum(self.probabilities)!=0,"Sum of probabilities must not be 0"
        
        if value is None:
            self.probabilities = (self.probabilities/numpy.sum(self.probabilities))
        else:
            if by is None or by=="Area":
                self.probabilities = value*(self.probabilities/numpy.sum(self.probabilities))
            elif by=="Maximum":
                self.probabilities = value*(self.probabilities/numpy.max(self.probabilities))
            else:
                raise ValueError("Unknown 'by' parameter")
        return self
    def checkNormalised(self,tolerance=1e-6):
        difference = numpy.abs(numpy.sum(self.probabilities)-1)
        if difference>tolerance:
            warnings.warn("Distribution not normalised - off by "+str(difference))
    
    def quantile(self,input):
        output = ()
        for value in input:
            cumulative_probabilities = numpy.concatenate(([0],numpy.cumsum(self.probabilities)))
            values = cumulative_probabilities-value
            #if any(values==0)
            #    output(self_index,value_index) = mean(self(self_index).bin_midpoints(values==0));
            if any(numpy.isnan(values)):
                output += (numpy.nan,)
            else:
                output += (self.piecewiseInterpolate(cumulative_probabilities,self.bin_edges,value),)
        return output

    def shift(self,amount):
        output = copy.copy(self)
        output.bin_edges = numpy.copy(self.bin_edges) + amount
        return output

    def getProbability(self,value):
        if self.type=="gaussian":
            return (1/(self.standard_deviation*numpy.sqrt(2*numpy.pi))) * numpy.exp(-0.5*(((value-self.mean)/self.standard_deviation)**2)) 
        elif self.type=="loggaussian":
            return (1/(value*self.standard_deviation*numpy.sqrt(2*numpy.pi))) * numpy.exp(-(((numpy.log(value)-self.mean)**2)/((2*self.standard_deviation)**2)))
        else:
            if len(value)>1:
                output = []
                for val in value:
                    if val>numpy.min(self.bin_midpoints) and val<numpy.max(self.bin_midpoints):
                        output += [self.piecewiseInterpolate(self.bin_midpoints,self.probabilities,val)]
                    else:
                        output += [0]
                return numpy.array(output)
            else:
                if value>numpy.min(self.bin_midpoints) and value<numpy.max(self.bin_midpoints):
                    return self.piecewiseInterpolate(self.bin_midpoints,self.probabilities,value)
                else:
                    return 0
    def getLogProbability(self,value):
        if self.type=="gaussian":
            return -numpy.log(self.standard_deviation)-(0.5*numpy.log(2*numpy.pi))-(0.5*((value-self.mean)/self.standard_deviation)**2)
        elif self.type=="loggaussian":
            return -numpy.log(value*self.standard_deviation)-(0.5*numpy.log(2*numpy.pi))-(((numpy.log(value)-self.mean)**2)/(2*self.standard_deviation**2))
        else:
            if len(value)>1:
                output = []
                for val in value:
                    output += [numpy.log(self.piecewiseInterpolate(self.bin_midpoints,self.probabilities,val))]
                return numpy.array(output)
            else:
                return numpy.log(self.piecewiseInterpolate(self.bin_midpoints,self.probabilities,value))
    
    @staticmethod
    def fromSamples(samples,bin_edges=None,location=None,weights=None):
        if bin_edges is None:
            histogram,bin_edges = numpy.histogram(samples)
        elif weights is None:
            histogram = numpy.histogram(samples,bins=bin_edges)[0]
        else:
            histogram = numpy.histogram(samples,bins=bin_edges,weights=weights)[0]
        return Distribution(bin_edges,"Manual",histogram,location=location)
    
    
    def plot(self,axis=None,**kwargs):
        if axis is None:
            axis = pyplot.gca()
        axis.plot(self.bin_midpoints,self.probabilities,**kwargs)
        axis.set_ylim(bottom=0)
    def area(self,axis=None,**kwargs):
        if axis is None:
            axis = pyplot.gca()
        axis.fill_between(self.bin_midpoints,self.probabilities,**kwargs)
        axis.set_ylim(bottom=0)
    def __repr__(self):
        return "Distribution("+str(self.bin_midpoints)+str(self.probabilities)+")"
    
    def toJSON(self,filename):
        json_data = json.dumps(self,cls=MCEncoder,indent=4)
        json_data_stripped = json_data.replace('"xxx',"").replace('xxx"',"").replace('xxx',"")
        with open(filename,"w") as file:
            file.write(json_data_stripped)
    

    @staticmethod
    def piecewiseInterpolate(x,y,xq):
        signed_distance = xq-x
        distance_sign = numpy.sign(signed_distance)
        distance_sign[distance_sign<0] = 0
        crossover = (distance_sign[0:-1]-distance_sign[1:]).astype("bool")
        distances = numpy.abs(numpy.concatenate((signed_distance[numpy.append(crossover,[False])],signed_distance[numpy.append([False],crossover)])))
        weights = 1-(distances*(1/numpy.sum(distances)))
        values =  [y[numpy.append(crossover,[False])],y[numpy.append([False],crossover)]]
        output = numpy.dot(weights,values)
        return output[0]
    @staticmethod
    def collapseArrayByLocation(distributions,tolerance):
        change_made = False
        distribution_locations = numpy.array([distribution.location for distribution in distributions])
        distributions_available = set(range(0,len(distribution_locations)))
        output_distributions = []
        for first_index in range(0,len(distribution_locations)):
            for second_index in range(0,len(distribution_locations)):
                if first_index!=second_index and first_index in distributions_available and second_index in distributions_available:
                    distance = numpy.abs(distributions[first_index].location-distributions[second_index].location)
                    if distance<tolerance:
                        change_made = True
                        distributions_available.remove(first_index)
                        distributions_available.remove(second_index)

                        new_location = (distributions[first_index].location+distributions[second_index].location)/2
                        new_probabilities = distributions[first_index].probabilities*distributions[second_index].probabilities
                        output_distributions += [Distribution(distributions[first_index].bin_edges,"Manual",new_probabilities,location=new_location).normalise()]
        for index in distributions_available:
            output_distributions += [distributions[index]]
        if change_made:
            return Distribution.collapseArrayByLocation(output_distributions,tolerance)
        else:
            return output_distributions
# class Distributions:
#     def __init__(self,distributions):
#         self.distributions = distributions
#     def 

class Sampler(Distribution):
    def __init__(self,bin_edges,type,values,method,location=None):
        super().__init__(bin_edges,type,values,location=location)
        self.method = method
    def getSamples(self,number_of_samples):
        if self.method.lower()=="monte_carlo":
            self.getMonteCarloSamples(number_of_samples)
            return self
        else:
            raise(ValueError("Unknown strategy for generating samples"))
    def getMonteCarloSamples(self,number_of_samples):
        r = numpy.random.rand(number_of_samples)
        samples = numpy.full([number_of_samples],numpy.nan)        
        
        cumulative_distribution = numpy.concatenate(([0],numpy.cumsum(self.probabilities)))
        if self.equally_spaced:
            a = 5
                
        for sample_index in range(number_of_samples):
            output_bins = numpy.zeros(len(cumulative_distribution))
            
            # Find which bin the random number falls in
            for bin_index in range(len(output_bins)):
                if r[sample_index]>=cumulative_distribution[bin_index] and r[sample_index]<cumulative_distribution[bin_index+1]:
                    output_bins[bin_index:bin_index+2] = 1
                    #output_bins = bool(output_bins)
                    break
            
            distance_from_left = r[sample_index] % 1/len(self.bin_midpoints)
            distance_from_right = 1/len(self.bin_midpoints)-distance_from_left
            
            values = self.bin_edges[output_bins.astype(bool)]
            distances = numpy.array((distance_from_left,distance_from_right))
            normalised_distance = distances/sum(distances)            
            
            samples[sample_index] = sum(normalised_distance*values)
        self.samples = samples
        return self

class MarkovChain:
    def __init__(self):
        self.samples = []
        self.burnt = 0
    def __len__(self):
        return len(self.samples)
    
    @property
    def burnt_index(self):
        if self.burnt>0:
            return math.floor(self.burnt/100*len(self))
        else:
            return 0
    def addSample(self,sample):
        output = MarkovChain()
        output.samples = self.samples+[sample]
        return output
    def accumulate(self,name):
        return tuple(sample.__dict__[name] for sample in self.samples[self.burnt_index:])
    def final(self,name):
        return self.samples[-1].__dict__[name]
    def burn(self,percentage):
        self.burnt = percentage
        return self
    def fromSamples(self,samples):
        for sample in samples:
            markov_chain_sample = MarkovChainSample()
            for (name,value) in sample.items():
                markov_chain_sample = markov_chain_sample.addField(name,value)
            self = self.addSample(markov_chain_sample)
        return self
    def round(self,precision):
        output = MarkovChain()
        for sample in self.samples:
            markov_chain_sample = MarkovChainSample()
            for (name,value) in sample.items():
                if isinstance(value,numpy.ndarray) or isinstance(value,numbers.Number):
                    markov_chain_sample = markov_chain_sample.addField(name,numpy.round(value,precision))            
                else:
                    markov_chain_sample = markov_chain_sample.addField(name,[numpy.round(array,precision) for array in value])
            output = output.addSample(markov_chain_sample)
        return output
    def toJSON(self,filename):
        json_data = json.dumps(self.samples,cls=MCEncoder,indent=4)
        json_data_stripped = json_data.replace('"xxx',"").replace('xxx"',"").replace('xxx',"")
        with open(filename,"w") as file:
            file.write(json_data_stripped)
    def fromJSON(self,filename):
        with open(filename,"r") as file:
            raw_data = file.read()
        json_data = json.loads(raw_data)
        self = self.fromSamples(json_data)

        return self
class MarkovChainSample:
    def __init__(self):
        pass
    def addField(self,name,value,precision=None):
        output = MarkovChainSample()
        # Copy over any previous values
        for item in self.__dict__.items():
            output.__dict__[item[0]] = item[1]
        # Add new value
        if name not in output.__dict__:
            if precision is None:
                output.__dict__[name] = value
            else:
                output.__dict__[name] = self.iterround(value,precision)
        else:
            raise ValueError(name+" already in MarkovChainSample")
        return output
    def items(self):
        return self.__dict__.items()
    def iterround(self,value,precision):
        if isinstance(value,numbers.Number):
            return round(value,precision)
        elif isinstance(value,numpy.ndarray):
            return numpy.round(value,precision)
        else:
            return [self.iterround(value,precision) for value in value]

class MCEncoder(json.JSONEncoder):
    def listOfArraysToString(self,list_of_arrays):
        return "xxx["+", ".join([self.arrayToString(group) for group in list_of_arrays])+"]xxx"
    def arrayToString(self,array):
        return "xxx["+", ".join([str(value) for value in numpy.squeeze(array)])+"]xxx"
    def toStr(self,array):
        if isinstance(array[0],int) or isinstance(array[0],float):
            return "xxx["+", ".join(str(x) for x in array)+"]xxx"
        else:
            return "xxx["+", ".join(str(x) for group in array for x in group)+"]xxx"
    def default(self,obj):
        output = {}
        if isinstance(obj,MarkovChainSample):
            for name,value in obj.__dict__.items():
                if isinstance(value,list) and isinstance(value[0],numpy.ndarray):
                    output[name] = self.listOfArraysToString(value)
                elif isinstance(value,numpy.ndarray) and len(value)==1:
                    output[name] = value[0]
                elif isinstance(value,numpy.ndarray):
                    output[name] = self.arrayToString(value)
                else:
                    output[name] = value
            return output
        elif isinstance(obj,Distribution):
            output["type"] = obj.type
            output["location"] = obj.location
            if obj.type!="manual":
                output["values"] = obj.values
            output["bin_edges"] =  self.arrayToString(obj.bin_edges)
            output["probabilities"] =  self.arrayToString(obj.probabilities)
            return output
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self,obj)
