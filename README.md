### About the paper
This is a code package for the paper: 
R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.

@ARTICLE{10364735,
  author={Liu, Rang and Li, Ming and Liu, Qian and Swindlehurst, A. Lee},
  journal={IEEE Transactions on Wireless Communications}, 
  title={SNR/CRB-Constrained Joint Beamforming and Reflection Designs for RIS-ISAC Systems}, 
  year={2023},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TWC.2023.3341429}}

- If you use this simulation code package in any way, please cite the original paper above.
- All codes are contributed by Rang Liu (email: rangl2@uci.edu; website: https://rangliu0706.github.io/). 
   Please feel free to contact with her if you have any suggestions. 
- The link of this paper is: https://ieeexplore.ieee.org/document/10364735
- More information can be found at: https://www.minglabdut.com/resource.html
- Copyright Notice: This code is licensed for personal, non-commercial use only, specifically for academic purposes. Copyright reserved by the MingLab (led by Prof. Ming Li), School of Information and Communication Engineering, Dalian University of Technology, Dalian 116024, China. 


### Software platform
- Please note that the MATLAB2022b is used for this simulation code package, and there may be some imcompatibility problems among different sofrware versions. 
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)

### Content of this simulation code package
- The folder "code_SI" contains the simulations for the scenario with self-interference, and the folder "code" for the scenario without self-interference.
- In the folder "code", run the file "main_plot_BP" to obtain Fig. 2 and run the file "main_plot_convergence" to obtain Fig. 3.
- In these two folders, the files "main_SR_P" are used for Fig. 4(a), the files "main_SR_P_CRB" for Fig. 4(b), the files "main_SR_N" for Fig. 5(a), the files "main_SR_N_CRB" for Fig. 5(b), the files "main_SR_SNR" for Fig. 6(a), and the files "main_SR_CRB" for Fig. 6(b). 

Abstract of the paper: 
In this paper, we investigate the integration of integrated sensing and communication (ISAC) and reconfigurable intelligent surfaces (RIS) for providing wide-coverage and ultrareliable communication and high-accuracy sensing functions. In particular, we consider an RIS-assisted ISAC system in which a multi-antenna base station (BS) simultaneously performs multiuser multi-input single-output (MU-MISO) communications and radar sensing with the assistance of an RIS. We focus on both target detection and parameter estimation performance in terms of the signal-to-noise ratio (SNR) and Cramer-Rao bound (CRB), respectively. Two optimization problems are formulated for maximizing the achievable sum-rate of the multi-user communications under an SNR constraint for target detection or a CRB constraint for parameter estimation, the transmit power budget, and the unit-modulus constraint of the RIS reflection coefficients. Efficient algorithms are developed to solve these two complicated non-convex problems. We then extend the proposed joint design algorithms to the scenario with imperfect self-interference cancellation. Extensive simulation results demonstrate the advantages of the proposed joint beamforming and reflection designs compared with other schemes. In addition, it is shown that more RIS reflection elements bring larger performance gains for directof-arrival (DoA) estimation than for target detection.





