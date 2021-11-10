# GRASS

GRAph-based Spatial Smoothing (GRASS) of fMRI data.

This toolbox is an implementation of the method proposed in the following [paper](https://doi.org/10.1101/2021.05.04.442605):

> Behjat, H., Westin, C.F. and Aganj, I., 2021. Cortical Surface-Informed Volumetric Spatial Smoothing of fMRI Data via Graph Signal Processing. bioRxiv 2021.05.04.442605.

A brief overview of the method and some results are presented in this [video](https://www.youtube.com/watch?v=gUtLZBCto-E).

### Prerequisites
- [Matlab](https://se.mathworks.com/products/matlab.html)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- [Human Connectome Project](http://www.humanconnectomeproject.org/) (HCP) data, which can be downloaded from [here](https://db.humanconnectome.org/). 

The current implementation is tailored for HCP data. An implementation that is suited for data other than the HCP will be released in the near future. 

### Usage and a Brief Description
The main function is `demo_hcpdata.m`. Open it and adjust the initial few settings, further description of which is given in the following. 

There are three settings that are related to the HCP data: specify the HCP root directory, select a desired HCP subject ID, and select a desired fMRI volume (out of the 19 available volumes for each HCP subject). The code relies on the assumption that preprocessed HCP data are properly extracted and reside in their original directory structure as obtained by extracting the downloaded data, as described [here](https://www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf), under section *Directory structure for preprocessed MR data*. 

There is a setting related to the graph design type `gtype`. A specific graph is created for either the left or the right hemisphere, and you also specify the spatial resolution of the graph. The graph represents the cerebral hemisphere cortex (CHC) of the individual, with graph vertices representing voxels that fall within the cortical ribbon and graph edges defined based on two principles: (i) adjacency of voxels in 3D space, (ii) pruning out of annatomically invalid connections that results from (i), for instance, at opposite banks of narrow sulci. You can create a graph that matches the resolution of voxels in the fMRI volume, e.g. 2 mm cubic, or a graph that has higher resolution, e.g. 1 mm cubic, to better benefit from the higher resolution provided by the structural data. In the latter case, fMRI data are upsampled to match the resolution of the strutural data. 

The smoothing parameter is called `tau`, which can be seen as an equivalent to the FWHM parameter used for Gaussian smoothing. This parameter controls the spatial size of the filter: by using a larger `tau`, wider spatial filter are created, and thus, more smoothing is applied to the data. Try few different `tau`s and visually inspect the smoothed volumes. 

Note that, unlike FWHM, which is usually specified in units of milimeter, `tau` does not have a unit. As such, its interpretation varies depending on the resolution of fMRI voxels, i.e., resolution of the graph design. In particular, a larger `tau` is needed on higher resolution graphs to obtain the same level of smoothing as obtained on lower resolution graphs; as a rule of thumb, the ratio can be assigned relative to the inverse ratio of voxel volumes, e.g., for a 1 mm cubic design, you would use a      `tau` that is 8 times larger than that used on a 2 mm cubic design, to obtain an approximately similar extent of smoothing. The rational behind this approximation is related to how heat kernels defined on the graph spectrum can be approximated by polynomials, and in turn, how the order of these polynomials is linked to the spatial width of the resulting graph filters in terms of number of *hops* on the graph.   

      
For further details on the method, please refer to the aforementioned [paper](https://doi.org/10.1101/2021.05.04.442605) and [video](https://www.youtube.com/watch?v=gUtLZBCto-E).
    