# A simulation of a dipole in the brain  
+ Useful for prediction of ECG voltage on the surface of the scalp
+ Modelling of epileptic seizures

# Finite Difference Model
+ Project is a numerical solution of the Laplace equation
+ Spherical coordinate system  
<img src="https://user-images.githubusercontent.com/32593957/189561597-ddc5dd1a-9fc9-4ec5-a58f-275f55001c06.png" width="250"/> 

+ Group terms and simplify  
<img src="https://user-images.githubusercontent.com/32593957/189561656-25a431e7-14a2-47c4-af68-e603abf59228.png" width="500"/>

+ Each term will be named A,B,C,D and E, respectively
  + Dirichlet boundary: A * X = b. Substituted:  
<img src="https://user-images.githubusercontent.com/32593957/189561795-807e82b1-eaa0-4975-9b41-49cd4a686bfd.png" width="300"/>

  + Neumann boundary:  
<img src="https://user-images.githubusercontent.com/32593957/189561827-e554fb45-c6be-4261-9f68-4c6d5afed9c9.png" width="100"/>
  
  + Substituted:  
<img src="https://user-images.githubusercontent.com/32593957/189561974-aed48b7f-4a54-41ac-9c57-9588dec51335.png" width="300"/>

+ FD Stencil:  
<img src="https://user-images.githubusercontent.com/32593957/189562203-d46d799f-c4bd-412c-b6f1-8bb49d82efbd.png" width="200"/>

# Results
<img src="https://user-images.githubusercontent.com/32593957/189562419-4ad2d600-3eb8-4e09-a18e-a022e1983e81.png" width="350" />
