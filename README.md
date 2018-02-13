# multi-dimensional-hits
Matlab codes used to implement and test the algorithm for the multi-dimensional HITS centrality introduced in 

		"Multi-dimensional HITS: An always computable ranking 	
	       	for temporal multi-layer directed networks" 

by Francesca Arrigo and Francesco Tudisco

To run, the codes require the Sandia Lab Matlab Tensor Toolbox available at:
http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html

The folder contains:
data		         - folder containing the dataset used in manuscript
multi_dimensional_HITS.m - function implementing the algorithm described in the manuscript.
rand_alpha.m             - function that generates random parameters alpha fulfilling the 				   convergence criterion
main.m                   - script that describes the guidelines for using the method
	

>> Reference:
>> "Multi-dimensional HITS: An always computable ranking for temporal multi-layer directed networks", 
>>  F. Arrigo (http://arrigofrancesca.wixsite.com/farrigo), 
>>  F. Tudisco (http://personal.strath.ac.uk/f.tudisco) 
