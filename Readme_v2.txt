Code for Bright Field electron tomographic reconstruction for crystalline
or amorphous samples 

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
                      --default_dosage :  The dosage value in counts when the material is not present
                      --use_default_dosage    :  A flag to tell the code to use the default input dosage value
 
                      [--diffuseness <0.2>]  : The value used to compute the prior model parameter "p" (p = diffuseness + 1)
                      [--bragg_T <3>]   : Forward model Bragg threshold parameter  
                      [--bragg_delta <0.5>]  : Forward model Bragg delta (weighting) parameter 
                      [--default_variance]    : The scaling value for the variance (\sigma^{2} in the paper)                      
                      [--bf_offset <0>]        : Offset in the measured counts. This value is added to the input file 
                      [--tilt_selection <0>] : Which was the tilt axis : 0 is default (y-axis)                    
                      
                      [--final_resolution <1>] : The resolution of the reconstructed pixels in integer multiples 
                                                of the detector size 
                      [--stop_threshold <.001>] : Stopping criteria as a percentage change of voxels relative 
                                                  to the previous iteration
                      [--num_resolutions <3> ]  : Number of resolutions of the multi-resolution reconstruction
                                                 before log(.) is applied
                      [ --outer_iterations <600>] : Maximum number of outer non-homgenous subiterations
                      [ --inner_iterations <100>] : Maximum number of inner non-homgenous subiterations 
                                                    at the coarsest resolution. 
                                                    This is done to provide a reasonable initial condition to the algo 
                      [--default_recon_value <0>] : At the coarsest resolution the object is initialized to this value
                      [--extend_object]       : Typical microscope samples extend out on the sides. An accurate 
                                                reconstruction requires using this flag to reconstruct 
                                                a large volume
                      [--delete_tmp_files]   : This flag is used to clear the temporary files created
                      [--exclude_views]      : Used to exclude certain views. Indicate the views to exclude 
                                               separated by "," (Ex: --exclude_views 5,10,30)

* Running the GUI 

   Requires Qt and cmake to correctly point to the GUI software. 
   Go throug same steps as before. 
   After "make"
   cd Bin/
   open TEMBIR.app
