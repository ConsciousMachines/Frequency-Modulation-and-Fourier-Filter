import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write
from scipy.signal import argrelextrema
import warnings
from scipy.misc import imsave
warnings.filterwarnings('ignore')
'''
mod2 ---> mod1 ---> y
         e3 * mod3--^
         e4 * mod4--^
IDEAS: make a saw and square wave! then a fractal FM synthesizer            
'''
dir0 = '/path/' # path to create music file

scale = np.array([2**(i/12) for i in range(12)])
scale3 = [0.99,1.0293022366434921, 1.0905077326652577, 1.1553526968722729, 1.2240535433046553, 1.2968395546510096, 1.3739536474580891, 1.4556531828421873, 1.5422108254079407, 1.6339154532411, 1.7310731220122859, 1.8340080864093424, 1.9430638823072117]
notes = np.array(['C','C#','D','D#','E','F','F#','G','G#','A','A#','B'])

def f2n(freq):#converts a frequency to the note
    base=16.35#tuning 
    freq=np.asarray(freq)
    a=np.floor(np.subtract(np.log2(freq),np.log2(base)))
    b=np.divide(freq,np.multiply(base,np.power(2,a)))
    if freq.shape!=():
        out = np.array([]) 
        for f in range(len(freq)): 
            if b[f]>=1.9430638823072117:out=np.append(out,['A',np.int8(a[f]+1)])
            else:out=np.append(out,[notes[np.searchsorted(scale3,b[f])-1],np.int8(a[f])])
        return np.squeeze(out).reshape(len(freq),2)
    elif freq.shape==():
        if b>=1.9430638823072117:return['A',a+1]
        else:return[np.squeeze(notes[np.searchsorted(scale3,b)-1]),a]

music_dir = '/path2/'
plotpath = '/path3/'
length = 5
car = 300 # main carrier frequency



class FM2(object):
    def __init__(self, carrier,  \
                 #mod1,mod2,mod3,mod4, \
                 mod1=car, mod2=car, mod3=car, mod4=car, \
                 p1=1,p2=0,p3=0,p4=0, r1=1,r2=1,r3=1,r4=1, \
                 length=length, sampleRate=44100,waveFile=True):
        specttrum_divisions = 40 # how many top frequencies we want (has to divide 20000 for now)
        interval_size = 20 # how fat of an interval we want around main freq's

        self.increment = .01 # p = power index, r = ratio
        self.carrier = carrier
        self.mod1 = mod1
        self.mod2 = mod2
        self.mod3 = mod3
        self.mod4 = mod4
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4
        self.rate = sampleRate
        self.ident = id(self)
        self.name = '%dHz_%sp1_%sr1_%dm1_%sp2_%sr2_%dm2_%sp3_%sr3_%dm3_%sp4_%sr4_%dm4.wav' % (self.carrier, str(self.p1), str(self.r1),self.mod1,str(self.p2), str(self.r2),self.mod2,str(self.p3),str(self.r3),self.mod3,str(self.p4),str(self.r4),self.mod4)
        sampleInc = 1.0/self.rate

        x = np.arange(0,length, sampleInc)
        x0 = np.arange(0,length, sampleInc) # E N V E L O P E   0
        ramp0 = [rmp/(length) for rmp in x0]
        ramp0 = ramp0[::-1]
        ramp01= ramp0[::-1]
        
        mod4 = p4*np.sin(2*np.pi*self.r4*self.mod4*x)
        mod3 = p3*np.sin(2*np.pi*self.r3*self.mod3*x)
        mod2 = p2*np.sin(2*np.pi*self.r2*self.mod2*x)
        mod1 = p1*np.sin(2*np.pi*self.r1*self.mod1*x + mod2 )
        y = np.sin(2*np.pi*self.carrier*x*ramp0 + ramp0*mod1 - ramp01*mod4)
                   
        mx = 1.059*(max(abs(y))) # scale to max pk of -.5 dB
        y = y/mx
        signal0 = y
        wavData = np.asarray(32000*y, dtype = np.int16)
        self.wavData = wavData
        write(music_dir + str('Orig'+self.name), 44100, self.wavData)
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone start
        N = self.rate * length
        y = np.fft.rfft(y)
        yfft = np.abs(y[0:N/2])

        yfft = 2./N*yfft # scale and filter out negative frequency component
        scale_factor = max(yfft)
        yfft = yfft / scale_factor # magnitude of strongest partial
        plt.figure(figsize=(7,8))

        plt.subplot(6, 1, 1) # Original wave
        inc = self.increment
        T = 1./self.carrier
        x_ax = np.arange(0,T*1.01, T/100.0)
        plt.plot(signal0, label = 'Original FM Waveform')
        plt.axis([0,400,-1,1])
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.title(self.name)
        plt.grid(True)
        
        plt.subplot(6, 1, 2) # Original wave FFT
        plt.xlim(0,20000)
        y_fftLOG = np.where(abs(yfft)>0, 20*np.log10(abs(yfft)), -101)
        plt.plot(y_fftLOG,label='Original FFT')
        plt.legend()
        y_new = []
        new_fft = []
        tops = [] # these are gona be frequencies of interest
        freq_indices = []
        sev_window = 20000//specttrum_divisions
        for i in range(specttrum_divisions):
            segment = np.zeros(sev_window)
            y_new_segment = np.zeros(sev_window)
            interval = yfft[i*sev_window:(i+1)*sev_window]
            y_interval = y[i*sev_window:(i+1)*sev_window]
            # TODO: CREATE LOGARITHMIC SLICES
            top = np.max(interval)
            freq_indices.append(np.where(interval==top)[0][0]+i*sev_window)
            top_index = np.where(interval==top)[0]
            bot_idx = max(0,top_index-interval_size)
            next_idx = top_index+interval_size
            # TODO: CREATE ADDITIONAL HARMONICS

            # TODO: CREATE THYCC BASS
            
            segment[bot_idx:next_idx] = interval[bot_idx:next_idx]
            y_new_segment[bot_idx:next_idx] = y_interval[bot_idx:next_idx]
            new_fft = np.concatenate([new_fft,y_new_segment],axis=0)
            y_new = np.concatenate([y_new,y_new_segment],axis=0)
        y_new0 = y_new
        yfft2 = np.abs(y_new[0:N/2])
        yfft2 = 2./N*yfft2  # rescale the sev fft 
        yfft2 = yfft2 / max(yfft2) 

        #~~~~~~~~~~~~~~~~~~~~ R E S O N A N T   F I L T E R ~~~~~~~~~~~~~~
        # Parameters
        cutoff = 5000
        res_interval = 20 # size of hanning window
        cutoff = cutoff + res_interval//2
        resonance = 0.1 # resonance boost STILL NEED TO MAKE RESONANCE
        lin_slope = 0.001
        slope = 10000
        
        # Use x to modulate over time
        #x = np.arange(0,length, sampleInc)

        y_log = np.where(abs(yfft2)>0, 20*np.log10(abs(yfft2)), -101)

        filter_x = np.arange(cutoff,20000, 1)
        # L O G A R I T H M I C   F I L T E R 
        #filter_y = np.clip(np.log2((filter_x-cutoff+1)*slope),a_min=0,a_max=101)*(-1)
        #print(filter_y[:100])
        # L I N E A R   F I L T E R 
        filter_y = np.clip((filter_x-cutoff+1)*lin_slope,a_min=0,a_max=101)*(-1)
        filtered_ylog = np.minimum(y_log[cutoff:] , filter_y)
        filt_graph = np.concatenate([y_log[:cutoff],filtered_ylog],axis=0)

        diff_ratios = filtered_ylog / y_log[cutoff:] # needs work
        #print(diff_ratios[:1000])
        newer_y = y_new0[cutoff:] / diff_ratios
        y_new = np.concatenate([y_new[:cutoff],newer_y],axis=0)
        
        # Virtical slope filter
        '''
        cutoff = 7500
        y_new[cutoff:] = 0
        '''
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        y_new = np.concatenate([y_new,np.zeros((110250-20000))],axis=0)
        y_new0 = np.concatenate([y_new0,np.zeros((110250-20000))],axis=0)

        plt.subplot(6, 1, 3) # Sev Styl FFT graph
        plt.cla()
        plt.xlim(0,20000)
        plt.plot(y_log,label = 'Sev Style Sliced FFT')
        plt.legend()

        plt.subplot(6, 1, 4) # Filtered FFT 
        plt.cla()
        plt.xlim(0,20000)
        plt.plot(filt_graph,label='Filtered FFT')
        plt.legend()

        plt.subplot(6, 1, 5) # Frequencies and Notes
        plt.cla()
        plt.axis('off')
        note_values = f2n(freq_indices)
        notes = []
        for p in note_values:
            notes.append(p[0]+' '+p[1])
        plt.text(0,0.6,'Frequencies: \n' + str(freq_indices),fontsize=7)
        plt.text(0,0.2,'Chromatic Notes: \n' + str(notes),fontsize=7)

        y = np.fft.irfft(y_new) 
        mx = 1.059*(max(abs(y))) 
        y = y/mx


        plt.subplot(6, 1, 6)
        plt.cla()
        x_ax = np.arange(0,T*1.01, T/100.0)
        plt.plot(y, label = 'Sev Style Modified Waveform')
        plt.axis([0,800,-1,1])
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.grid(True)
        plt.draw()
        plt.pause(0.0001)
        plt.clf()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone end
        wavData = np.asarray(32000*y, dtype = np.int16)
        self.wavData = wavData
        write(music_dir + str(self.name), 44100, self.wavData)

#FM2(carrier=61,mod1=240,mod2=150,mod4=300,p1=30,p2=0,p4=0)
FM2(carrier=110,mod1=121,mod2=61,mod4=962,p1=5,p2=5,p4=5)





class FM(object):
    def __init__(self, carrier,  \
                 #mod1,mod2,mod3,mod4, \
                 mod1=car, mod2=car, mod3=car, mod4=car, \
                 p1=1,p2=0,p3=0,p4=0, r1=1,r2=1,r3=1,r4=1, \
                 length=length, sampleRate=44100,waveFile=True):
        specttrum_divisions = 20 # how many top frequencies we want (has to divide 20000 for now)
        interval_size = 40 # how fat of an interval we want around main freq's

        self.increment = .01 # p = power index, r = ratio
        self.carrier = carrier
        self.mod1 = mod1
        self.mod2 = mod2
        self.mod3 = mod3
        self.mod4 = mod4
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4
        self.rate = sampleRate
        self.ident = id(self)
        self.name = '%dHz_%sp1_%sr1_%dm1_%sp2_%sr2_%dm2_%sp3_%sr3_%dm3_%sp4_%sr4_%dm4.wav' % (self.carrier, str(self.p1), str(self.r1),self.mod1,str(self.p2), str(self.r2),self.mod2,str(self.p3),str(self.r3),self.mod3,str(self.p4),str(self.r4),self.mod4)
        sampleInc = 1.0/self.rate

        x = np.arange(0,length, sampleInc)
        x0 = np.arange(0,length, sampleInc) # E N V E L O P E   0
        ramp0 = [rmp/(length) for rmp in x0]
        ramp0 = ramp0[::-1]
        ramp01= ramp0[::-1]
        
        mod4 = p4*np.sin(2*np.pi*self.r4*self.mod4*x)
        mod3 = p3*np.sin(2*np.pi*self.r3*self.mod3*x)
        mod2 = p2*np.sin(2*np.pi*self.r2*self.mod2*x)
        mod1 = p1*np.sin(2*np.pi*self.r1*self.mod1*x + mod2 )
        y = np.sin(2*np.pi*self.carrier*x*ramp0 + ramp0*mod1 - ramp01*mod4)
                   
        mx = 1.059*(max(abs(y))) # scale to max pk of -.5 dB
        y = y/mx
        signal0 = y
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone start
        N = self.rate * length
        y = np.fft.rfft(y)
        yfft = np.abs(y[0:N/2])

        yfft = 2./N*yfft # scale and filter out negative frequency component
        scale_factor = max(yfft)
        yfft = yfft / scale_factor # magnitude of strongest partial
        plt.figure(figsize=(7,8))

        plt.subplot(6, 1, 1) # Original wave
        inc = self.increment
        T = 1./self.carrier
        x_ax = np.arange(0,T*1.01, T/100.0)
        plt.plot(signal0, label = 'Original FM Waveform')
        plt.axis([0,400,-1,1])
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.title(self.name)
        plt.grid(True)
        
        plt.subplot(6, 1, 2) # Original wave FFT
        plt.xlim(0,20000)
        y_fftLOG = np.where(abs(yfft)>0, 20*np.log10(abs(yfft)), -101)
        plt.plot(y_fftLOG,label='Original FFT')
        plt.legend()
        y_new = []
        new_fft = []
        tops = [] # these are gona be frequencies of interest
        freq_indices = []
        sev_window = 20000//specttrum_divisions
        for i in range(specttrum_divisions):
            segment = np.zeros(sev_window)
            y_new_segment = np.zeros(sev_window)
            interval = yfft[i*sev_window:(i+1)*sev_window]
            y_interval = y[i*sev_window:(i+1)*sev_window]
            
            top = np.max(interval)
            freq_indices.append(np.where(interval==top)[0][0]+i*sev_window)
            top_index = np.where(interval==top)[0]
            bot_idx = max(0,top_index-interval_size)
            next_idx = top_index+interval_size
            
            segment[bot_idx:next_idx] = interval[bot_idx:next_idx]
            y_new_segment[bot_idx:next_idx] = y_interval[bot_idx:next_idx]
            new_fft = np.concatenate([new_fft,y_new_segment],axis=0)
            y_new = np.concatenate([y_new,y_new_segment],axis=0)
        y_new0 = y_new
        yfft2 = np.abs(y_new[0:N/2])
        yfft2 = 2./N*yfft2  # rescale the sev fft 
        yfft2 = yfft2 / max(yfft2) 

        #~~~~~~~~~~~~~~~~~~~~ R E S O N A N T   F I L T E R ~~~~~~~~~~~~~~
        # Parameters
        cutoff = 500
        res_interval = 20 # size of hanning window
        cutoff = cutoff + res_interval//2
        resonance = 0.1 # resonance boost STILL NEED TO MAKE RESONANCE
        lin_slope = 0.005
        slope = 10000
        
        # Use x to modulate over time
        #x = np.arange(0,length, sampleInc)

        y_log = np.where(abs(yfft2)>0, 20*np.log10(abs(yfft2)), -101)

        filter_x = np.arange(cutoff,20000, 1)
        # L O G A R I T H M I C   F I L T E R 
        #filter_y = np.clip(np.log2((filter_x-cutoff+1)*slope),a_min=0,a_max=101)*(-1)
        #print(filter_y[:100])
        # L I N E A R   F I L T E R 
        filter_y = np.clip((filter_x-cutoff+1)*lin_slope,a_min=0,a_max=101)*(-1)
        filtered_ylog = np.minimum(y_log[cutoff:] , filter_y)
        filt_graph = np.concatenate([y_log[:cutoff],filtered_ylog],axis=0)

        diff_ratios = filtered_ylog / y_log[cutoff:] # needs work
        #print(diff_ratios[:1000])
        newer_y = y_new0[cutoff:] / diff_ratios
        y_new = np.concatenate([y_new[:cutoff],newer_y],axis=0)
        
        # Viertical slope filter
        '''
        cutoff = 7500
        y_new[cutoff:] = 0
        '''
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        y_new = np.concatenate([y_new,np.zeros((110250-20000))],axis=0)
        y_new0 = np.concatenate([y_new0,np.zeros((110250-20000))],axis=0)

        plt.subplot(6, 1, 3) # Sev Styl FFT graph
        plt.cla()
        plt.xlim(0,20000)
        plt.plot(y_log,label = 'Sev Style Sliced FFT')
        plt.legend()

        plt.subplot(6, 1, 4) # Filtered FFT 
        plt.cla()
        plt.xlim(0,20000)
        plt.plot(filt_graph,label='Filtered FFT')
        plt.legend()

        plt.subplot(6, 1, 5) # Frequencies and Notes
        plt.cla()
        plt.axis('off')
        note_values = f2n(freq_indices)
        notes = []
        for p in note_values:
            notes.append(p[0]+' '+p[1])
        plt.text(0,0.6,'Frequencies: \n' + str(freq_indices),fontsize=7)
        plt.text(0,0.2,'Chromatic Notes: \n' + str(notes),fontsize=7)

        y = np.fft.irfft(y_new) 
        mx = 1.059*(max(abs(y))) 
        y = y/mx
        
        plt.subplot(6, 1, 6)
        plt.cla()
        x_ax = np.arange(0,T*1.01, T/100.0)
        plt.plot(y, label = 'Sev Style Modified Waveform')
        plt.axis([0,800,-1,1])
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.grid(True)
        plt.draw()
        plt.pause(0.0001)
        plt.clf()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sev filter zone end
        wavData = np.asarray(32000*y, dtype = np.int16)
        self.wavData = wavData
        write(music_dir + str(self.name), 44100, self.wavData)

#FM(carrier=450,mod1=76,mod4=10,p1=40,p4=0)




class FMs(object):
    def __init__(self, carrier,  \
                 #mod1,mod2,mod3,mod4, \
                 mod1=car, mod2=car, mod3=car, mod4=car, \
                 p1=1,p2=0,p3=0,p4=0, r1=1,r2=1,r3=1,r4=1, \
                 length=length, sampleRate=44100,waveFile=True,plot=False,movie=False):
        self.increment = .01 # p = power index, r = ratio
        self.carrier = carrier
        self.mod1 = mod1
        self.mod2 = mod2
        self.mod3 = mod3
        self.mod4 = mod4
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4
        self.rate = sampleRate
        self.ident = id(self)
        #self.name = '%dHz_car%sp1_%sr1_%sp2_%sr2_%sp3_%sr3_.wav' % (self.carrier, str(self.p1), str(self.r1),     ,str(self.p2), str(self.r2),str(self.p3), str(self.r3))
        self.name = '%dHz_%sp1_%sr1_%dm1_%sp2_%sr2_%dm2_%sp3_%sr3_%dm3_%sp4_%sr4_%dm4.wav' % (self.carrier, str(self.p1), str(self.r1),self.mod1,str(self.p2), str(self.r2),self.mod2,str(self.p3),str(self.r3),self.mod3,str(self.p4),str(self.r4),self.mod4)

        sampleInc = 1.0/self.rate

        x = np.arange(0,length, sampleInc)
        x0 = np.arange(0,length, sampleInc) # E N V E L O P E   0
        ramp0 = [rmp/(length) for rmp in x0]
        ramp0 = ramp0[::-1]
        '''
        x1 = np.arange(0,length/2.0, sampleInc) # E N V E L O P E   1
        ramp = [rmp/(length/2.0) for rmp in x1]
        ramp2 = ramp # /
        ramp = ramp[::-1] # \
        #ramp = ramp + ramp2 # \ / over length
        ramp = ramp2+ramp
        ramp3 = np.sin(x*3)
        ramp3_5 = np.sin(x*5)

        x2 = np.arange(0,length/8.0, sampleInc)
        ramp4 = [rmp/(length/8.0) for rmp in x2]
        ramp5 = ramp4 # /
        ramp4 = ramp4[::-1] # \
        #full_ramp = ramp4+ramp5+ramp4+ramp4+ramp2# \ / \ \ then slow /
        full_ramp = ramp2+ramp4+ramp5+ramp4+ramp4
        full_ramp = full_ramp[:-2]
        
        full_ramp = np.array(full_ramp[:110250])
        full_ramp = np.append(full_ramp,np.zeros(110250))
        #print(len(full_ramp))
        
        ramp3_6 = np.array(ramp3_5[110250:])
        ramp3_6 = np.append(np.zeros(110250),ramp3_6)
        #print(len(ramp3_6))
        '''
        
        mod4 = p4*np.sin(2*np.pi*self.r4*self.mod4*x)
        mod3 = p3*np.sin(2*np.pi*self.r3*self.mod3*x)
        mod2 = p2*np.sin(2*np.pi*self.r2*self.mod2*x)
        mod1 = p1*np.sin(2*np.pi*self.r1*self.mod1*x + mod2)
        #y = np.sin(2*np.pi*self.carrier*x*ramp + full_ramp*mod1 -full_ramp*mod3 + mod4)
        y = np.sin(2*np.pi*self.carrier*x*ramp0 + ramp0*mod1 )
                   
        mx = 1.059*(max(abs(y))) # scale to max pk of -.5 dB
        y = y/mx
        wavData = np.asarray(32000*y, dtype = np.int16)
        self.wavData = wavData
        '''
        sev_yfft = np.fft.rfft(self.wavData)
        sev_inv = np.fft.irfft(sev_yfft)
        '''
        #write(sev_dir + str(self.name), 44100, sev_inv)
        write(sev_dir + str(self.name), 44100, self.wavData)
        
        if movie == True:
            data = []
            for i in range(1,len(x)):
                mod4 = p4*np.sin(2*np.pi*self.r4*self.mod4*i)
                mod3 = p3*np.sin(2*np.pi*self.r3*self.mod3*i)
                mod2 = p2*np.sin(2*np.pi*self.r2*self.mod2*i)
                mod1 = p1*np.sin(2*np.pi*self.r1*self.mod1*i + mod2)
                y = np.sin(2*np.pi*self.carrier*i*ramp[i] + full_ramp[i]*mod1 + ramp3_6[i]*mod3)
                data.append(y)
                if i >= 4005 :
                    if i % 500 == 0: 
                        rate = self.rate
                        data = data[-4000:]
                        T = 1.0/rate # sample period
                        N = len(data) # N samples
                        xfft = np.linspace(0.0, rate/2.0, N/2) # x axis
                        #print((np.min(xfft),np.max(xfft),'x'))
                        yfft = np.fft.rfft(data) # perform fft
                        yfft = 2./N*np.abs(yfft[0:N/2])# scale and filter out negative frequency component
                        #print((np.min(yfft),np.max(yfft),'y'))
                        yfft = yfft / max(yfft) # magnitude of strongest partial
                        yfftLog = np.where(abs(yfft)>0, 20*np.log10(abs(yfft)), -101)#where magnutude of yfft>0, it gets transformed, and the rest is filled with -101dB
                        plt.plot(xfft,yfftLog)
                        plt.ylim(-100,0)
                        plt.xlim(0,5000)
                        plt.grid(True)
                        plt.title('FFT of %s' % self.name)
                        plt.savefig(plotpath+str(i), bbox_inches="tight")
                        plt.clf()
   
        if plot == True:
            '''
            inc = self.increment
            T = 1./self.carrier
            x = np.arange(0,T*1.01, T/100.0)
            plt.plot(y, label = self.name)
            plt.axis([0,100,-1.1,1.1])
            labels = [str(round(1000*n,2)) for n in np.arange(0,T*1.01,T/10.0)]
            plt.xticks([n for n in np.arange(0, 101, 10)], labels)
            leg = plt.legend(fancybox=True)
            leg.get_frame().set_alpha(0.5)
            plt.xlabel('Milliseconds')
            plt.title('Source Waveform')
            plt.grid(True)
            plt.show()'''

            rate = self.rate
            data = self.wavData
            T = 1.0/rate # sample period
            N = len(data) # N samples
            xfft = np.linspace(0.0, rate/2.0, N/2) # x axis
            yfft = np.fft.rfft(data) # perform fft
            yfft = 2./N*np.abs(yfft[0:N/2])# scale and filter out negative frequency component
##            if fft2 == True:
##                plt.plot(xfft,yfft)
##                plt.ylim(-100,0)
##                plt.xlim(0,5000)
##                plt.show()
##                plt.close('all')
            yfft = yfft / max(yfft) # magnitude of strongest partial
            yfftLog = np.where(abs(yfft)>0, 20*np.log10(abs(yfft)), -101)#where magnutude of yfft>0, it gets transformed, and the rest is filled with -101dB
            plt.plot(xfft,yfftLog)
            plt.ylim(-100,0)
            plt.xlim(0,5000)
            plt.grid(True)
            plt.title('FFT/%s' % self.name)
            plt.show()
            
                
                
FM0 = FMs(carrier=300,mod1=302,p1=1,p2=0,p3=0,p4=0,r1=1,r2=1,r3=7,plot=False ,length=1)

FM0 = FMs(carrier=300,mod1=151,p1=3,p2=0,p3=0,p4=0,r1=1,r2=1,r3=1,plot=True ,length=5)



sampleInc = 1.0/44100.
x = np.arange(0,length, sampleInc)
x1 = np.arange(0,length/2.0, sampleInc) # E N V E L O P E   Z O N E
ramp = [rmp/(length/2.0) for rmp in x1]
ramp2 = ramp # /
ramp = ramp[::-1] # \
ramp = ramp + ramp2 # \ / over length
ramp3 = np.sin(x*3)
ramp3_5 = np.sin(x*5)
#plt.plot(ramp)
#plt.show()
#plt.plot(ramp3)
#plt.show()
x2 = np.arange(0,length/8.0, sampleInc)
ramp4 = [rmp/(length/8.0) for rmp in x2]
ramp5 = ramp4 # /
ramp4 = ramp4[::-1] # \
full_ramp = ramp4+ramp5+ramp4+ramp4+ramp2 # \ / \ \ then slow /
full_ramp = full_ramp[:-2]
#plt.plot(full_ramp)
#plt.show()
ramp6 = np.sin(full_ramp)






























