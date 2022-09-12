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

Eye diagram at 100 meters, no Line-Of-Sight between transmitter and receiver  
![image](https://user-images.githubusercontent.com/32593957/189559268-f4a6f7b3-d11d-47e9-8c01-ad4c7c481fdc.png)

Constellation for the 16 QAM of the above  
![image](https://user-images.githubusercontent.com/32593957/189559293-b99512f7-1356-4bdc-b85e-3d41edb47d82.png)



