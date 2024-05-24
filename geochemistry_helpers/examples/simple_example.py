
# Normally this wouldn't be necessary as you should install the geochemistry_helpers library
# But this will work whether or not it is installed by adding the files to the path
import sys
sys.path.append(".")

# Import the necessary pieces
from geochemistry_helpers import GaussianProcess,Sampling
import numpy # For matrices and arrays
from matplotlib import pyplot # For plotting


# Assume we're looking at something that can be anywhere between 0 and 10
x_values = numpy.arange(0,10,0.01) 
# And time is anywhere between 0 and 100
t_values = numpy.arange(0,100,0.1)

# Create some random constraints
constraints = [] # Create an empty list
constraints += [Sampling.Distribution(x_values,"Gaussian",(4,0.2),location=10).normalise()] # Create a distribution (using inputs: axis, type of distribution, (mean,standard_deviation), t_value)
constraints += [Sampling.Distribution(x_values,"Gaussian",(2,0.5),location=50).normalise()]
constraints += [Sampling.Distribution(x_values,"Gaussian",(5,0.1),location=80).normalise()]

# Create the Gaussian Process
# Uses a 'builder pattern'
# Create an empty gaussian process, then constrain it with a list of distributions (ideally Gaussian), set the kernel (type and hyperparameters), then query the fit (at the t_values), and get samples
gaussian_process = GaussianProcess().constrain(constraints).setKernel("rbf",(1,20)).query(t_values).getSamples(1000)

# Create some subplots
figure,axes = pyplot.subplots(nrows=1)
axes = [axes] # (if you create more than one subplot, they come out as a list, but if you create one it is not in a list, so I've put it into a list so that I can use the same syntax for one vs more than one subplot)

gaussian_process.plotArea(axis=axes[0],alpha=0.5)
gaussian_process.plotConstraints(axis=axes[0],color="black")
gaussian_process.plotSamples(axis=axes[0],color="grey",indices=[0,1,2])

axes[0].set_xlabel("t")
axes[0].set_ylabel("y")

pyplot.show()


