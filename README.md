# complexRenewableEnergyNetworks-PropagationSimulation

Generate electric fields using diffusion.
Required libraries:

numpy, matplotlib.pyplot, time

Required libraries for cuda:
cupy



# Comparison off speeds using cuda
# CPU: i3 8350k
# GPU: Nvidea RTX 3070
# at 500x500 grid and 8 charges

# Using CPU single thread
------------------
Starting iterating
------------------
0.0 Calculated 92821.36667124118 Mpixel*steps/sec
0.04 Calculated 493.2773321025157 Mpixel*steps/sec
0.08 Calculated 499.2477159777741 Mpixel*steps/sec


# Using Cuda
------------------
Starting iterating
------------------
0.0 Calculated 2680.7579363176774 Mpixel*steps/sec
0.04 Calculated 3870.1454481690603 Mpixel*steps/sec
0.08 Calculated 3950.6346402739905 Mpixel*steps/sec
0.12 Calculated 4091.2581464373316 Mpixel*steps/sec
0.16 Calculated 3894.164340951664 Mpixel*steps/sec
0.2 Calculated 3891.23234473775 Mpixel*steps/sec
0.24 Calculated 4057.247124250416 Mpixel*steps/sec