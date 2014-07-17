Code for Bright Field electron tomographic reconstruction for crystalline
or amorphous samples 

This package contains source codes to reconstruct STEM HAADF tomographic data using the 
Model Based Iterative Reconstruction Method. The source codes are open source
under a BSD License. See the License.txt file for more information.

The current set of dependencies are:

Boost version 1.44 or newer
Qt version 4.7.4 or newer
libTiff

In order to build the sources you will also need the following:

CMake version 2.8.6 or newer.

The codes are known to compile on OS X (10.6, 10.7, 10.8), Windows (Visual Studio 2008 and 2010)
and RedHat Enterprise Linux 5.x, CentOS 5.x and OpenSUSE 12.1

To compile code on Unix based OS: 

1) cd to the Folder with the full package

2) mkdir Build 

3) cd Build 

4) cmake ../
 
   If you wish to adjust the parameters, try using ccmake instead and switch off things like the GUI 

5) make 


Note : All arguments within a [] are optional parameters.The <> indicate the default values of the parameters.

Running the command line code : 

./MbirReconstruction  -s <> : Full path of the input MRC file. Expects the header to be in the FEI format to read tilt 
                              angles and pixel sizes
                     --subvolume : Sub-volume to choose of the form x_start,y_start,tilt_start,x_end,y_end,tilt_end
                                   The code assumes the data is centered i.e. the center of rot is at (x_end-x_start)/2
                      --outputfile : Full path of the output .rec/mrc file
                      --thickness : Sample thickness in nano-meters 
                      --sigma_x  : The scaling parameter in the q-GGMRF prior in units of nm^{-1}
                      --target_gain <> : Average value of the direct beam (normalizing) measurement in counts; This effects 
                                         the choice of sigma_x. The GUI version can automatically estimate the sigma_x given 
                                         this value
                      [--default_offset] :  The value in counts when the material is not present (dark counts). If not present 
                                            the program automatically initializes it  
                      [--use_default_offset]    :  A flag to tell the code to use the default input offset value. MUST BE 
                                                  used if default_offset is used
                      [--diffuseness <0.2>]  : The value used to compute the prior model parameter "p" (p = diffuseness + 1)
                      [--default_variance]    : The scaling value for the variance (\sigma^{2} in the paper)                      
                      [--tilt_selection <0>] : Which was the tilt axis : 0 is default (y-axis)                                          
                      [--final_resolution <1>] : The resolution of the reconstructed pixels in integer multiples 
                                                of the detector size 
                      [--stop_threshold <.005>] : Stopping criteria as a percentage change of voxels relative 
                                                  to the previous iteration
                      [--num_resolutions <3>]  : Number of resolutions of the multi-resolution reconstruction
                                                 before log(.) is applied
                      [ --outer_iterations <30>] : Maximum number of outer iterations
                      [ --inner_iterations <10>] : Maximum number of inner iterations 
                                                    at the coarsest resolution. 
                                                    This is done to provide a reasonable initial condition to the algo 
                      [--default_recon_value <0>] : At the coarsest resolution the object is initialized to this value
                      [--extend_object]       : Typical microscope samples extend out on the sides. An accurate 
                                                reconstruction requires using this flag to reconstruct 
                                                a large volume
                      [--delete_tmp_files]   : This flag is used to clear the temporary files created
                      [--exclude_views]      : Used to exclude certain views. Indicate the views to exclude 
                                               separated by "," (Ex: --exclude_views 5,10,30)
                      [--num_threads]        : Specify the number of threads. If not specified the default is to use the number 
                                               of cores present in the machine. (Only present in CryoTomo branch)
                      [--tilts_file]         : Full path name of text file containing a list of tilts. 
                                               If this is present the header tilts will be OVER-WRITTEN. (Only present in CryoTomo branch)

***********************
Running the GUI 
***********************

   Requires Qt and cmake to correctly point to the GUI software. 
   Go through same steps as before. 
   After "make"
   cd Bin/
   open TEMBIR.app
