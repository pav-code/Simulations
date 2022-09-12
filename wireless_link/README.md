# Physical Wireless System  
+ 100 meter range
+ 16 QAM supported
+ Handles phase offset, frequency offset, and frequency drift
+ Physically tested on BPSK, QPSK, 16 QAM
+ 64 QAM supported as well, above this noise cannot be handled
+ 2.4 GHz transmission
+ 500kHz bandwidth


# Hardware
+ Keysight ADS
+ Receiver Gain: 23.8 dB 
+ Transmitter Gain: 10 dB
+ Components:  
   + Low Noise Amplifer: MGA-86563
   + LO/IF Amplifers: HMC311ST89E
   + Up/down conversion mixers: HMC272AMS8E

# Software
+ Transmitter:  
   + Character Input
   + Bits to Symbol conversion
   + Baker Sequence addition
   + Upsample and filter
   + Add pilot sine wave
   + Modulate to IF
   
+ Receiver:  
   + Demodulate to baseband
   + Course frequency synchronization
   + Signal detection and frame synchronization
   + Downsample
   + Fine frequency synchronization
   + Phase Estimation and Correrction
   + Normalizating Constellation
   + Optimal Maximum-Likelihood detection
   + Bits to Ascii
   + Display the message

