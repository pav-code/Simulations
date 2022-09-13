# Simulations
Digital modulations, Error Control Coding schemes, Channel Estimation Simulations, Medical Images reconstructions and FEM simulations

+ _Error Control Coding_:  
  + An implementation of a Low Density Parity Check (LDPC) Stochastic decoder based on research done by Qualcomm
+ _finite_difference_electric_dipole_: 
  + A solution of a spherical version of the Laplace equation to analyze the effect of a dipole in the human brain. Solved via two boundary conditions: Neumann and Dirichlet.
+ _wireless_link_:
  + A Complete 100 meter link between a receiver and a transmitter with a total bandwidth of 0.5 MB/s
  + A high frequency PCB was designed in Keysight ADS and manufactured. Electronic components include:
    + LNA
    + Modulator/Demodulator
    + LO/IF Amplifers
    + Up and Downconversion Mixers
  + Software, implemented in MATLAB, system components include:
    + Course/Fine frequency and frame synchronizations and corrections
    + Phase estimations and corrections
    + Constellation normalization
    + Optimal Maximum-Likelihood (ML) Detector
    + Message display (PC)
+ _communication_system_over_a_fading_channel_:
  + Simulation of an over the air transmission channel
  + Based on a Rayleigh distribution for a non-line of sight, urban city environment
  + Orthogonal frequency division multiplexing (OFDM) transmission over the channel
  + Coded and uncoded transmission (error control coding scheme: repetition)

+ _90_nm_differential_amplifier_:
  + Four circuits make up the design: push-pull, differential pair, high swing cascode current sink and bootstrapped current source
  + Design is compared versus a production 90 nm differential amplifier

+ _image_processing_:
  + _computed_tomography_:
    + Backprojection algorithm to reconstruct a medical CT image
  + _diffusion_tensor_imaging_:
    + Analysis of a MRI mamalogram producing FA (Fractional Anisotropy), VR (Volume Ratio), MD (Mean Diffusivity) 
    + Diffusion Signal Curve Fitting
    + PDD (Principal Direction of Diffusion) on major cerebral arteries 
  + _image_segmentation_and_reconstruction_:
    + ART (Argebraic Reconstruction Techniques) Kaczmarz algorithm
