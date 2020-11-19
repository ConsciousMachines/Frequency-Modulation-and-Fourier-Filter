import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import write, read
out_dir = r'C:\Users\pwnag\Desktop'

# Part 1 - emulating FM8 growl 

# well i got really close with :
# E = np.sin(2 * np.pi * time * E_f)
# F = np.sin(2 * np.pi * time * F_f + E * 13.0)
# but once I added D->E, the whole thing got messed. using the same scale, it's too noisy. so 13/5 and 13 don't work as 100;20 should.
# all I can guess is that it has some crazy post processing to sound good. I mean, why not, FM sounds like shit most of the time. 
# whatever. I think I got close enough, logic dictates there is something else going on with the sound, since using 13/5 and 13 
# gives me crazy artefacts while their 100 and 20 are fine. But my 10 and 2 sound similar but off. 
# when using 
# D = np.sin(2 * np.pi * time * D_f)
# E = np.sin(2 * np.pi * time * E_f + D * 0.5)
# F = np.sin(2 * np.pi * time * F_f + E * 13.0 * lfo1)
# it feels like there is a traveling filter, as in a filter that moves along some frequency band, taking out some nasty artefacts. 
# this might be what makes FM8 sound good. 

# based on doing a 440hz sine -100%-> sine, i think FM8 amplitude goes to 2.0
# Combos 
# car 30 mod 69.420 / 1 
# car 30 mod 30 / 1 
# car 20 mod 40 / 1

sps = 44100             # Samples per second
duration_s = 5.0        # Duration in seconds 

fund = 99.0
D_f = fund * 14.0
E_f = fund * 0.5
F_f = fund * 2.0        # Frequency / pitch of the sine wave

time = np.arange(duration_s * sps) / sps
lfo1 = 0.5 * (np.sin(2 * np.pi * time * 4.0) + 1.0)

D = np.sin(2 * np.pi * time * D_f)
E = np.sin(2 * np.pi * time * E_f + D * 0.4)
F = np.sin(2 * np.pi * time * F_f + E * 13.0 * lfo1)
wave_ints = np.int32(F * 0.2 * 2147483647) # because waveform is float so * 2**16 - also make quiet using *0.8 
write(out_dir + '\\first_sine_wave.wav', sps, wave_ints)


#A = np.sin(2 * np.pi * time * 440.0)
#B = np.sin(2 * np.pi * (440.0 + A * 2.0))
#wave_ints = np.int16(B * 0.8 * 32767) # because waveform is float so * 2**16 - also make quiet using *0.8 
#write(out_dir + '\\first_sine_wave.wav', sps, wave_ints)



GRAF_LEN = 50000
plt.subplot(4, 1, 1)
plt.title('Frequency Modulation')
plt.plot(lfo1)
plt.ylabel('LFO')
plt.subplot(4, 1, 2)
plt.plot(D[:GRAF_LEN])
plt.ylabel('D')
plt.subplot(4, 1, 3)
plt.plot(E[:GRAF_LEN])
plt.ylabel('E')
plt.subplot(4, 1, 4)
plt.plot(F[:GRAF_LEN])
plt.ylabel('F')
plt.show()









# Part 2 - cleaning up the old code 

def go(): # this is a working cleanup of the old FM code based on Mr Colson. 
    # parameters & setup 
    if True:
        length = 5
        sampleRate = 44100
        N = sampleRate * length
        specttrum_divisions = 40 # how many top frequencies we want (has to divide 20000 for now)
        interval_size = 20 # how fat of an interval we want around main freq's
        window_len = 20000//specttrum_divisions
        car = 300 
        carrier = 110
        mod1 = 121
        mod2 = 61
        mod3 = car 
        mod4 = 962
        p1 = 5
        p2 = 5
        p3 = 0 
        p4 = 5
        r1=1
        r2=1
        r3=1
        r4=1

        name = '%dHz_%sp1_%sr1_%dm1_%sp2_%sr2_%dm2_%sp3_%sr3_%dm3_%sp4_%sr4_%dm4.wav' % (carrier, str(p1), str(r1),mod1,str(p2), str(r2),mod2,str(p3),str(r3),mod3,str(p4),str(r4),mod4)
        sampleInc = 1.0/sampleRate

        x = np.arange(0,length, sampleInc)
        x0 = np.arange(0,length, sampleInc) # E N V E L O P E   0
        ramp0 = [rmp/(length) for rmp in x0]
        ramp0 = ramp0[::-1]
        ramp01= ramp0[::-1]

        mod4 = p4*np.sin(2*np.pi*r4*mod4*x)
        mod3 = p3*np.sin(2*np.pi*r3*mod3*x)
        mod2 = p2*np.sin(2*np.pi*r2*mod2*x)
        mod1 = p1*np.sin(2*np.pi*r1*mod1*x + mod2 )
        y = np.sin(2*np.pi*carrier*x*ramp0 + ramp0*mod1 - ramp01*mod4)
                    
        mx = 1.059*(max(abs(y))) # scale to max pk of -.5 dB
        y = y/mx
        wavData = np.asarray(32000*y, dtype = np.int16)
        #write(music_dir + str('Orig'+name), 44100, wavData)
            
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone start
        _yfft = np.fft.rfft(y)
        yfft = np.abs(_yfft) * 2./N # scale and filter out negative frequency component
        scale_factor = max(yfft)
        yfft = yfft / scale_factor # magnitude of strongest partial

        new_fft = np.zeros(20000)
        freq_indices = []
        for i in range(specttrum_divisions):
            interval = yfft[i*window_len:(i+1)*window_len]
            # TODO: CREATE LOGARITHMIC SLICES
            top = np.max(interval)
            top_idx = np.where( interval == top )[0][0] # top index is relative to the 0-500 window. the others are relative to 0-20k
            freq_indices.append( top_idx + i * window_len)
            bot_idx = max(0,top_idx - interval_size) + i*window_len
            next_idx = min(window_len, top_idx + interval_size) + i*window_len
            # TODO: CREATE ADDITIONAL HARMONICS
            new_fft[bot_idx:next_idx] = _yfft[bot_idx:next_idx]
    new_fft_s = np.abs(new_fft) * 2./N # rescale the sev fft 
    new_fft_s = new_fft_s / max(new_fft_s) 

    #~~~~~~~~~~~~~~~~~~~~ R E S O N A N T   F I L T E R ~~~~~~~~~~~~~~
    y_log = np.where(abs(new_fft_s)>0, 20*np.log10(abs(new_fft_s)), -101) # log transform into dB

    cutoff = 5010 # cutoff = 5000; res_interval = 20; cutoff + res_interval//2
    lin_slope = 0.001
    filter_x = np.arange(cutoff,20000, 1)
    filter_y = np.clip((filter_x-cutoff+1)*lin_slope,a_min=0,a_max=101)*(-1) # filter is a negative sloped line
    filtered_ylog = np.minimum(y_log[cutoff:] , filter_y)       # frequencies above cutoff are ceilinged by filter_y
    diff_ratios = y_log[cutoff:] / filtered_ylog                # power of filtered frequencies relative to original 

    # [original under cutoff] - [filtered part] - [rest, padded w zeros for IFFT]
    new_fft_f = np.concatenate([new_fft[:cutoff],new_fft[cutoff:20000] * diff_ratios, np.zeros(((N//2)-20000))],axis=0) # apply these ratios to the new_fft freqs 

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone end
    if True:
        y = np.fft.irfft(new_fft_f) # movement comes from phase interference. Since ifft removes any time information, and throws all freqs together. 
        y = y/(1.059*max(abs(y)))
        wavData = np.asarray(32000*y, dtype = np.int16)
        write(out_dir + '\\' + str(name), sampleRate, wavData)
        a = read(out_dir + '\\' + str(name))
        b = read(r'C:\Users\pwnag\Desktop\ORIG.wav')
    return np.sum(np.square(a[1] - b[1]))

go()






# Test Zone 
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


specttrum_divisions = 50 # how many top frequencies we want (has to divide 20000 for now)
interval_size = 10 # how fat of an interval we want around main freq's
# parameters & setup 
if True:
    length = 5
    sampleRate = 44100
    N = sampleRate * length
    window_len = 20000//specttrum_divisions
    car = 50
    carrier = car
    mod1 = car + 0.01
    mod2 = car + 0.02
    mod3 = car + 0.03
    mod4 = car + 0.04
    p1 = 5
    p2 = 5
    p3 = 0 
    p4 = 5
    r1=1
    r2=1
    r3=1
    r4=0.5

    name = '%dHz_%sp1_%sr1_%dm1_%sp2_%sr2_%dm2_%sp3_%sr3_%dm3_%sp4_%sr4_%dm4.wav' % (carrier, str(p1), str(r1),mod1,str(p2), str(r2),mod2,str(p3),str(r3),mod3,str(p4),str(r4),mod4)
    sampleInc = 1.0/sampleRate

    x = np.arange(0,length, sampleInc)
    x0 = np.arange(0,length, sampleInc) # E N V E L O P E   0
    ramp0 = [rmp/(length) for rmp in x0]
    ramp0 = ramp0[::-1]
    ramp01= ramp0[::-1]

    mod4 = p4*np.sin(2*np.pi*r4*mod4*x)
    mod3 = p3*np.sin(2*np.pi*r3*mod3*x)
    mod2 = p2*np.sin(2*np.pi*r2*mod2*x)
    mod1 = p1*np.sin(2*np.pi*r1*mod1*x + mod2 )
    y = np.sin(2*np.pi*carrier*x*ramp0 + ramp0*mod1 - ramp01*mod4)
    '''
    fund = 99.0
    D_f = fund * 14.0
    E_f = fund * 0.5
    F_f = fund * 2.0        # Frequency / pitch of the sine wave

    time = np.arange(length * sampleRate) / sampleRate
    lfo1 = 0.5 * (np.sin(2 * np.pi * time * 4.0) + 1.0) 

    D = np.sin(2 * np.pi * time * D_f)
    E = np.sin(2 * np.pi * time * E_f + D * 0.4)
    y = np.sin(2 * np.pi * time * F_f + E * 13.0 * lfo1)
    '''
    mx = 1.059*(max(abs(y))) # scale to max pk of -.5 dB
    y = y/mx
    wavData = np.asarray(32000*y, dtype = np.int16)
    #write(music_dir + str('Orig'+name), 44100, wavData)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone start
    _yfft = np.fft.rfft(y)
    yfft = np.abs(_yfft) * 2./N # scale and filter out negative frequency component
    scale_factor = max(yfft)
    yfft = yfft / scale_factor # magnitude of strongest partial

    new_fft = np.zeros(20000)
    freq_indices = []
    for i in range(specttrum_divisions):
        interval = yfft[i*window_len:(i+1)*window_len]
        # TODO: CREATE LOGARITHMIC SLICES
        top = np.max(interval)
        top_idx = np.where( interval == top )[0][0] # top index is relative to the 0-500 window. the others are relative to 0-20k
        freq_indices.append( top_idx + i * window_len)
        bot_idx = max(0,top_idx - interval_size) + i*window_len
        next_idx = min(window_len, top_idx + interval_size) + i*window_len
        # TODO: CREATE ADDITIONAL HARMONICS
        new_fft[bot_idx:next_idx] = _yfft[bot_idx:next_idx]

new_fft_s = np.abs(new_fft) * 2./N # rescale the sev fft 
new_fft_s = new_fft_s / max(new_fft_s) 

#~~~~~~~~~~~~~~~~~~~~ R E S O N A N T   F I L T E R ~~~~~~~~~~~~~~
y_log = np.where(abs(new_fft_s)>0, 20*np.log10(abs(new_fft_s)), -101) # log transform into dB

cutoff = 5010 # cutoff = 5000; res_interval = 20; cutoff + res_interval//2
lin_slope = 0.001
filter_x = np.arange(cutoff,20000, 1)
filter_y = np.clip((filter_x-cutoff+1)*lin_slope,a_min=0,a_max=101)*(-1) # filter is a negative sloped line
filtered_ylog = np.minimum(y_log[cutoff:] , filter_y)       # frequencies above cutoff are ceilinged by filter_y
diff_ratios = y_log[cutoff:] / filtered_ylog                # power of filtered frequencies relative to original 



#   F I L T E R   A N A L Y S I S
y_log_inv = y_log[cutoff:] + 101
y_log_inv = np.power(y_log_inv / max(y_log_inv), 10)

filter_inv = filter_y + 101
filter_inv = np.power(filter_inv / max(filter_inv), 10)

fig, axs = plt.subplots(3)
axs[0].plot(y_log_inv)
axs[0].plot(filter_inv)
axs[1].plot(new_fft_s[cutoff:])
plt.show()

#np.sum(np.square(y_log_inv - new_fft_s[cutoff:])) # there is error since i used np.where to only do certain points
# I D E A : next filter : use sigmoid to avoid the log 

# I D E A : do this for every 1024 size window 

# [original under cutoff] - [filtered part] - [rest, padded w zeros for IFFT]
new_fft_f = np.concatenate([new_fft[:cutoff],new_fft[cutoff:20000] * diff_ratios, np.zeros(((N//2)-20000))],axis=0) # apply these ratios to the new_fft freqs 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone end
if True:
    y = np.fft.irfft(new_fft_f) # movement comes from phase interference, from big window size. ifft removes any time information, and throws all freqs together. 
    y = y/(1.059*max(abs(y))) 
    wavData = np.asarray(32000*y, dtype = np.int16)
    write(out_dir + '\\' + str(name), sampleRate, wavData)
    a = read(out_dir + '\\' + str(name))
    b = read(r'C:\Users\pwnag\Desktop\ORIG.wav')
np.sum(np.square(a[1] - b[1]))
