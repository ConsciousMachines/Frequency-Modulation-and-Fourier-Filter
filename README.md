# Frequency-Modulation-and-Fourier-Filter
Applications of Fourier Transforms to Frequency Modulation and filtering the frequency space 

This is a generalization of the beautiful work at:
http://www.mrcolson.com/2016/04/21/Simple-Python-FM-Synthesis.html

The famous FM8 synthesizer that is the basis of famous musical works by artists like Skrillex is for the most part various combinations of Frequency Modulation operators. Anyone who made their own FM syntehsizer, whether analog or digital like in Reaktor or using python, will realize that very quickly the sound starts to become noisy and difficult to listen to as you add larger parameters and operators. This work is an attempt to filter out the noise out of a many-operator FM combination. 

It does so by converting the original FM signal into its Fourier space, and then taking out chunks of it. This is completely subjective and open to experimentation- the width of the frequency window you take out, the amount (i used linear but it should be logarithmic) of frequency "points" you filter at (this is basically a comb filter). But making these variations makes the resulting sound much more clean, pleasant, mysterious, and outer-space-like while maintaining the metallic properties so commonly associated with FM. 

Example 1: 110 Hz (warning: lower volume for the pre-filter output) the file name is based on the notation by the original python-FM8 article, and I try to keep it descrptive by showing the power and frequency of each modulator. 

![alt tag](https://github.com/ConsciousMachines/Frequency-Modulation-and-Fourier-Filter/blob/master/ex1.png)

Example 2: 300 Hz

![alt tag](https://github.com/ConsciousMachines/Frequency-Modulation-and-Fourier-Filter/blob/master/ex2.png)

There is a lot of fun to be had with this tool especially for DSP enthusiasts. 
