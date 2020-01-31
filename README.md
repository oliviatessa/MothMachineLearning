# MothMachineLearning


Insect flight is a highly non-linear dynamical system.  As such, strategies for understanding control have typically relied on either simulation methods (e.g., Model Predictive Control (MPC), genetic algorithms) or linearization of the dynamical system. Here we develop a new framework that combines MPC and deep learning to create an efficient method for solving the inverse problem of flight control. We used a feedforward, fully-connected neural network to answer the question, “What is the temporal pattern of forces required to follow a complex trajectory?” Combining neural networks with simulations based on dynamical systems models yields a data-driven controller where the data are derived from a non-linear physical model. We first trained a deep neural network (4 hidden layers, with hundreds of nodes) on ~8 million simulated 2D insect trajectories. Our network accurately predicted the force, force angle, abdomen angle, and tangential and angular velocities (7 outputs), when it was provided with initial conditions and a goal location (12 inputs). The coefficient of determination (r^2) for all predictions was > 0.999 on a validation dataset (1 million additional trajectories). Next, we evaluated the neural network’s ability to control a simulated insect.  We used the aforementioned predictions and compared the final conditions generated to simulations. Again, we found that network-prescribed final conditions were nearly identical to numerically solved conditions (r^2 > 0.999). Overall, this work shows that machine-learning may be an efficient approach for controlling nonlinear dynamical systems.


This is a work in progress (code is still messy). 


# Here's how I set up the Anaconda environment on my computer

## create conda environment with tools for generating data and training network
## note that this environment requires a GPU and the correct NVIDIA drivers
```conda env create -f environment_ODE_DL.yml```


## Create kernelspec
```
conda activate deeplearn_V4
python -m ipykernel install --user --name deepLearn_V4
conda deactivate
```
