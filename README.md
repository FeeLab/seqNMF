# seqNMF 
Unsupervised discovery of temporal sequences in high-dimensional datasets, with applications to neuroscience
 
Emily Mackevicius and Andrew Bahle - FeeLab :notes: :bird: 2018

### Description
SeqNMF is an algorithm which uses regularized convolutional non-negative matrix factorization to extract repeated sequential patterns from high-dimensional data. It has been validated using neural calcium imaging, spike data, and spectrograms, and allows the discovery of patterns directly from timeseries data without reference to external markers.

For more information see our [**paper**](https://elifesciences.org/articles/38471); [**COSYNE talk**](https://www.youtube.com/watch?reload=9&v=XyWtCtZ_m-8); tutorial [**video**](https://cbmm.mit.edu/video/unsupervised-discovery-temporal-sequences-high-dimensional-datasets) and [**materials**](https://stellar.mit.edu/S/project/bcs-comp-tut/materials.html); and Simons foundation [**article**](https://www.simonsfoundation.org/2018/05/04/finding-neural-patterns-in-the-din/).



### Usage
The main function is seqNMF.m and it can be called 
```matlab
[W,H,cost,loadings,power] = seqNMF(X,'K',K,'L',L,'lambda',0.01)
```
Where X is the data matrix, K and L are the factorization parameters and lambda is a parameter controling the strength of regularization.

Specifically seqNMF factorizes the NxT data matrix X into K factors. Factor exemplars are returned in the NxKxL tensor W. Factor timecourses are returned in the KxT matrix H

                                    ----------    
                                L  /         /|
                                  /         / |
        ----------------         /---------/  |          ----------------
        |              |         |         |  |          |              |
      N |      X       |   =   N |    W    |  /   (*)  K |      H       |           
        |              |         |         | /           |              |
        ----------------         /----------/            ----------------
               T                      K                         T

### Demo
See the [demo script](demo.m), for a demonstration of the seqNMF algorithm on synthetic data and songbird imaging data. This demo also gives examples of how to cross validate, test for significance and select parameters.




