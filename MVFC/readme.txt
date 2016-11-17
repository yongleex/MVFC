This is the implementation package of MVFC method for the paper "Outlier detection for particle image velocimetry data using a locally estimated noise variance";
We submit it to the Journal "Measurement Science & Technology" in september 2016;

*********************************************************************************
Software Requirement: We successfully tested this package in Matlab 2015b and Matlab 2016a, and therefore we recommend you use the latest version of Matlab.

Installation: 
            1.Unzip the files, and run "Initialization.m" to add the paths to the Matlab temporarily.
            2.Enjoy the examples in the folder "..\Example".
            3.The syntax of each method please refer to the "*.m" files and the examples.


Flow Cases:
            1. A synthetic cellular vortex flow.
               |- case 1: Homogeneous Gaussian white noise + 15% Scattered outlier +  4*4 cluster;
               |- case 2: Non-Homogeneous Gaussian noise   + 15% Scattered outlier +   0  cluster;
               |- Stat 1: Find the best parameters H and L in a parameter grid net; Note that it is very time comsuming.
               |- Stat 2: 100 Montel-Carlo simulation with addition of scattered outlers. 
               |- Stat 3: 100 Montel-Carlo simulation with addition of clustered outlers.
               |- Stat 4: 100 Montel-Carlo simulation with variable cluster sizes.
               |- Stat 5: 100 Montel-Carlo simulation with variable vortex number (vortex resolution).
               |- Stat 6: 100 Montel-Carlo simulation with variable level of Non-Homogeneous Gaussian noise + scattered outlers.
               |- Stat 7: 100 Montel-Carlo simulation with variable level of Non-Homogeneous Gaussian noise + clustered outlers
               
            2. Experimental turbulent jet flow.               
            3. A synthetic source flow.
            4. Experimental vortex pair flow
         Note: The results of these experiments are restored in the floder "..\figs\".

Other Methods:
            0. MVFC(in the folder "..\MVFC_dir\MVFC.m"): Our proposed modified vector field correction method.
            1. NMT/CON(in the folder "..\convl\convl2.m"): a classical Normalized Median Test(Outlier Validation) + replaced by 5*5 (it is 3*3 previous) Median value + smoothed with 3*3 average kernel.
            2. VTM(in the folder "..\VTM\vtmedian.m"): variable threshold median test.
            3. FADV(in the folder "..\FADV\FADV.m"): flow adaptive data validation method.
            4. VFC(in the folder "..\VFC\VFC.m"): Vector Field Correction method in Exp Fluid.
            5. DCT-PLS(in the fold "..\pppiv\pppiv.m"): Garcia's work, A robust smooth method combining Discrete Cosine Transform and Penalized Least Square.

Evaluation criteria £¨in the folder "..\tools"£©:
            1. Visualized comparison: This is the direct method to see the difference. (coutour-plot of the velocity magnitude)
            2. ODC(L_odc.m): over-detected vector count. 'L' means Lee Yong, the author of this package
            3. UDC(L_udc.m): Un-detected vector count.
            % following criteria are not mentioned in the paper
            4. NRMSE: normalized root mean square error. Ref£ºGarcia D (2011) A fast all-in-one method for automated post-processing of PIV data. Exp Fluids.
            5. SSIM: structure similarity.  Ref£ºZhou Wang, Alan C Bovik, Hamid R Sheikh, and Eero P Simoncelli. Image quality assessment: from error visibility to structural similarity. IEEE Trans.Image Processing¡£
            6. TES:  turbulent energy spectrum. Ref£ºPope, Turbulent Flows, 1st Edition, Page 219~221.
               
Reference: 
            1. Yong Lee and Hua Yang "Outlier detection for particle image velocimetry data using a locally estimated noise variance", submitted to Meas Sci & Tech
            2. Westerweel J, Scarano F (2005) Universal outlier detection for PIV data. Exp Fluids.
            3. Vennemann P (2008) Particle image velocimetry for microscale blood flow measurement. TU Delft, Delft University of Technology
            4. Garcia D (2011) A fast all-in-one method for automated post-processing of PIV data. Exp Fluids.
            5. Garcia D (2010) Robust smoothing of gridded data in one and higher dimensions with missing values. Comput Stat Data Anal 54(4):1167¨C1178
            6. Shinneeb, A. M., Bugg, J. D., & Balachandar, R. (2004). Variable threshold outlier identification in PIV data. Meas Sci Technol, 15(9), 1722-1732.
            7. Liu, Z. L., Jia, L. F., Zheng, Y., & Zhang, Q. K. (2008). Flow-adaptive data validation scheme in PIV. Chem Eng Sci 63(1), 1-11.
            8. Lee Y, Yang H, Yin Z. A robust vector field correction method via a mixture statistical model of PIV signal[J]. Experiments in Fluids, 2016, 57(3): 1-20.
            9. Wang H, Gao Q, Feng L, Wei R, Wang J (2015) Proper orthogonal decomposition based outlier correction for piv data. Exp Fluids 56(2):1¨C15

This is the end. 
Yong Lee (Leeyong@hust.edu.cn)
2016.09.04 @ Huazhong University of Science and technology.