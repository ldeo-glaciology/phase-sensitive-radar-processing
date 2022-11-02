# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 20:06:36 2020

Class definitions for ApRES processing code

DataFileObject.
==============
Instantiated with one required string, giving the name of the data file
eg fileDescriptor = DataFileObject('DATAFILENAME.DAT')

Methods:
    ExtractBurst(BurstNumber (integer))
        Output is an instance of a BurstObject
eg Burst = fileDescriptor.ExtractBurst(3)

Instance variables:
    Filename           : Name of data file
    BurstLocationList  : Python list of byte offset of each burst in file
    NoBurstsInFile     : Number of bursts in the file (len(BurstLocationList))
    
BurstObject.
============
Typically instantiated with a call to the ExtractBurst method on a DataFileObject
eg Burst = fileDescriptor.ExtractBurst(3)

Methods:
    ExtractChirp(ChirpList (Python list))
        Output is an instance of a ChirpObject, in which all chirps in the ChirpList
        have been averaged
    PlotBurst()
        Plots the full raw burst as a time series

Instance variables:
    v         : Array containing burst data
    Filename  : Name of data file
    Header    : Burst header (Python dictionary), with additional entries:
    BurstNo   : Burst number in data file

ChirpObject
===========
Typically instantiated with a call to the ExtractChirp method on a BurstBbject
eg Chirp = Burst.ExtractChirp([1,3])

Methods:
    FormProfile(StartFreq, StopFreq, padfactor, ref)
        StartFreq, StopFreq: start and end frequencies to use (eg 2e8 and 4e8)
        padfactor:           zero padding for the fft (eg. 2)
        ref:                 whether or not to apply Paul Brennan's reference
                             phase (1 or 0, for yes or no)
        Returns and instance of a ProfileObject
    PlotChirp()
        Plots the chirp as function of time

Instance variables:
    vdat:       Array containing chirp data
    t:          Array containing time for chirp samples
    ChirpList:  List of chirps averaged to make vdat
    Filename:   Name of data file
    BurstNo:    Number of burst within data file
    Header:     Burst header, as created by ExtractBurst method on FileDataObject

ProfileObject.
==============
Typically instantiated with a call to the FormProfile method on a ChirpObject
eg Profile = Chirp.FormProfile(StartFreq, StopFreq, padfactor, ref)

Methods:
    PlotProfile(MaxDepth (double))
        MaxDepth:  Maximum depth (in metres) to which to plot profile
        
Instance variables:
    Range:     Array with depth in metres each profile depth bin
    Profile:   Array containing profile (complex double)
    F0:        Start frequency used to form profile
    F1:        End frequency used to form profile
    pad:       Pad factor used when zeropadding
    ChirpList: List of chirps averaged to form profile
    Filename:  Name of original data file 
    BurstNo:   Number of burst in data file
    Header:    Burst header, as produced using ExtractBurst method on DataFileObject 
    rad2m:     radians to metres of range conversion factor
    bin2m:     bin to metres of range conversion factor
"""
import gcsfs
import numpy as np
import matplotlib.pyplot as plt
import math
import warnings
import copy
class DataFileObject:
    
    def __init__(self, Filename, remote_load=False):
        self.BurstLocationList = []
        self.Filename = Filename
        self.remote_load = remote_load
        

        if remote_load:
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode = 'rb')
        else: 
            datafile = open(self.Filename, "rb")
            

        inbuff = datafile.read()
        datafile.close()     
        
        a = "*** Burst Header ***"
        b = a.encode()
        locn = inbuff.find(b)
        while locn != -1:
            self.BurstLocationList.append(locn)
            locn = inbuff.find(b, locn + len(b))
       
        self.NoBurstsInFile = len(self.BurstLocationList)

    def ExtractBurst(self, BurstNo):
        Burst = BurstObject()
        Burst.Filename = self.Filename
        Burst.BurstNo = BurstNo
        Burst.Header = {"BurstNo":BurstNo}
        
        if self.remote_load:
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode = 'rb')
        else: 
            datafile = open(self.Filename, "rb")
 
        
        datafile.seek(self.BurstLocationList[BurstNo])
        inbuff = datafile.read(2000)
        locn = 0
        locn1 = inbuff.find(b'\x0D\x0A', locn, locn + 80)
        inline = inbuff[locn:locn1].decode()
        while inline.count("End Header") == 0:
            tmp = inline.split("=")
            if len(tmp) == 2:
                if tmp[0].lower() == "rxant" or \
                   tmp[0].lower() == "txant" or \
                   tmp[0].lower() == "afgain":
                   Burst.Header[tmp[0]] = \
                        [int(x) for x in tmp[1].split(',') if x]
                        
                elif tmp[0].lower() == "triples" or \
                     tmp[0].lower() == "attenuator1" or \
                     tmp[0].lower() == "batterycheck":
                     Burst.Header[tmp[0]] = \
                        [float(x) for x in tmp[1].split(',') if x]
                                      
                elif tmp[0].lower() == "latitude" or \
                     tmp[0].lower() == "longitude" or \
                     tmp[0].lower() == "temp1" or \
                     tmp[0].lower() == "temp2" or \
                     tmp[0].lower() == "batteryvoltage" or \
                     tmp[0].lower() == "tstepup" or \
                     tmp[0].lower() == "tstepdn" or \
                     tmp[0].lower() == "fsc" or \
                     tmp[0].lower() == "sw_issue" or \
                     tmp[0].lower() == "er_ice" or \
                     tmp[0].lower() == "position_depth_conversion" or \
                     tmp[0].lower() == "maxdepthtograph":
                    Burst.Header[tmp[0]] = float(tmp[1])
                    
                elif tmp[0].lower() == "rmb_issue" or \
                     tmp[0].lower() == "vab_issue" or \
                     tmp[0].lower() == "reg00" or \
                     tmp[0].lower() == "reg01" or \
                     tmp[0].lower() == "reg02" or \
                     tmp[0].lower() == "reg03" or \
                     tmp[0].lower() == "reg0b" or \
                     tmp[0].lower() == "reg0c" or \
                     tmp[0].lower() == "reg0d" or \
                     tmp[0].lower() == "reg0e":
                    Burst.Header[tmp[0]] = tmp[1]
                elif tmp[0].lower() == "time stamp":
                    Burst.Header[tmp[0]] = tmp[1]
                else:
                    Burst.Header[tmp[0]] = int(tmp[1])
                    
            locn = locn1 + 2
            locn1 = inbuff.find(b'\x0D\x0A', locn,locn + 80)
            inline = inbuff[locn:locn1].decode()

        # Re-open the file for binary read to get burst data
        
        if self.remote_load:
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode = 'rb')
        else: 
            datafile = open(self.Filename, "rb")      
        
        datafile.seek(self.BurstLocationList[BurstNo] + locn1 + 2)
        NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"] * Burst.Header["NSubBursts"]
        if Burst.Header["Average"] == 1:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"]
            inbuff = datafile.read(NsamplesInBurst * 4)
            Burst.v = np.frombuffer(inbuff, dtype=np.float32)/2**16*2.5-1.25
        elif Burst.Header["Average"] == 2:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"]
            inbuff = datafile.read(NsamplesInBurst * 4)
            Burst.v = np.frombuffer(inbuff, dtype=np.uint32)
        else:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"] * Burst.Header["NSubBursts"] *  Burst.Header["nAttenuators"]  # because each subburst conatins a certain number of chirps given by Burst.Header["nAttenuators"]
            inbuff = datafile.read(NsamplesInBurst * 2)
            Burst.v = (np.frombuffer(inbuff, dtype=np.uint16))/2**16*2.5-1.25
        datafile.close()
        
        if "SamplingFreqMode" in Burst.Header:
            if Burst.Header["SamplingFreqMode"] == 0:
                Burst.Header["dt"] = 1/40000
            else:
                Burst.Header["dt"] = 1/80000
        else:
            Burst.Header["dt"] = 1/40000
            
        Burst.Header["NChirps"] = Burst.Header["NSubBursts"] * Burst.Header["nAttenuators"] *\
            Burst.Header["TxAnt"].count(1) * Burst.Header["RxAnt"].count(1)
        
        if "FreqStepUp" in Burst.Header:
            Burst.Header["K"] = Burst.Header["FreqStepUp"] / Burst.Header["TStepUp"]
        else:
            Burst.Header["K"] = 2e8
            Burst.Header["StartFreq"] = 2e8
            Burst.Header["StopFreq"] = 4e8
            
        Burst.Header["c0"] = 3e8
        if  not ("ER_ICE" in Burst.Header):
            Burst.Header["ER_ICE"] = 3.18
        Burst.Header["CentreFreq"] = (Burst.Header["StartFreq"] + Burst.Header["StopFreq"])/2
        Burst.Header["B"] = (Burst.Header["StopFreq"] - Burst.Header["StartFreq"])
        
        # deal out attenuator settings to the chirps
        setting_counter = 0
        Burst.Header['Attenuator1_allChirps'] = []
        Burst.Header['AFGain_allChirps'] =[]
        for chirp in range(Burst.Header['NChirps']):
            Burst.Header['Attenuator1_allChirps'].append(Burst.Header['Attenuator1'][setting_counter])
            Burst.Header['AFGain_allChirps'].append(Burst.Header['AFGain'][setting_counter])
            
            # keep track of which setting to use
            setting_counter += 1
            if setting_counter >= Burst.Header['nAttenuators']: # if the counter reaches nAttenuators, reset it to zero 
                setting_counter = 0


        return(Burst)
        
class BurstObject:
    
    def __init__(self):
        self.v = 0
        self.Filename = ''
        self.Header = 0
        self.BurstNo = 0
        
    def ExtractChirp(self, ChirpList):
        Chirp = ChirpObject()
        Chirp.t = np.array(range(self.Header["N_ADC_SAMPLES"])) * self.Header["dt"]
        Chirp.Header = copy.deepcopy(self.Header)  # we dont want chirp-specific new entries in the chirp.header updating the burst.header
        Chirp.ChirpList = ChirpList
        
        Chirp.Header['Attenuator1_thisChirp'] = [self.Header['Attenuator1_allChirps'][i] for i in ChirpList]
        Chirp.Header['AFGain_thisChirp'] = [self.Header['AFGain_allChirps'][i] for i in ChirpList]

        if any(i != Chirp.Header['Attenuator1_thisChirp'][0] for i in Chirp.Header['Attenuator1_thisChirp'])\
            or any(i != Chirp.Header['AFGain_thisChirp'][0] for i in Chirp.Header['AFGain_thisChirp']):
            warnings.warn('This is stacking over chirps with different attenuator settings.')


        if self.Header["Average"] == 1:
            Chirp.vdat = self.v          
        elif self.Header["Average"] == 2:
            Chirp.vdat = self.v          
        else:
            Chirp.vdat = np.zeros(self.Header["N_ADC_SAMPLES"])
            no = 0
            for ind in ChirpList:
                if ind < self.Header["NChirps"]:
                    no += 1
                    chirpoffset = ind * self.Header["N_ADC_SAMPLES"]
                    Chirp.vdat = Chirp.vdat + self.v[chirpoffset:chirpoffset + self.Header["N_ADC_SAMPLES"]]
                else:
                    print('chirp index > number of chirps.')
            Chirp.vdat = Chirp.vdat/no

        return(Chirp)
        
    def PlotBurst(self):
        t = np.array(range(len(self.v))) * self.Header["dt"]
        plt.plot(t, self.v)
        plt.axis([0,t[-1],-1.25,1.25])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (V)")
        plt.grid()
        return(0)
    
class ChirpObject:
    def __init__(self):
        self.vdat = 0
        self.t = 0
        self.ChirpList = 0
        self.Filename = ''
        self.BurstNo = 0
        self.Header = 0
        
    def PlotChirp(self):
        plt.plot(self.t, self.vdat)
        plt.axis([0,self.t[-1],-1.25,1.25])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (V)")
        plt.grid("on")
        return(0)

    def FormProfile(self, F0=200000000, F1=400000000, pad=2, ref=1):
        Profile = ProfileObject()
        if self.Header["StartFreq"] > F0:
            F0 = self.Header["StartFreq"]
        if self.Header["StopFreq"] < F1:
            F1 = self.Header["StopFreq"]
        T0 = (F0-self.Header["StartFreq"])/self.Header["K"]
        T1 = (F1-self.Header["StartFreq"])/self.Header["K"]
        chirp = self.vdat[math.ceil(T0/self.Header["dt"]):\
                          math.floor(T1/self.Header["dt"])]

        Nt = len(chirp)
        Nt = math.floor(Nt/2) * 2
        winchirp = np.multiply(chirp[0:Nt],np.blackman(Nt))
        Nfft = math.floor(Nt*pad)

        padchirp = np.zeros(Nfft)
        padchirp[0:math.floor(Nt/2)-1] = winchirp[math.floor(Nt/2):-1]
        padchirp[-math.floor(Nt/2):-1] = winchirp[0:math.floor(Nt/2)-1]
        Profile.Profile = np.fft.fft(padchirp)/Nfft * math.sqrt(2*pad) 
        Profile.bin2m = self.Header["c0"]/(2.*(T1-T0)*pad*math.sqrt(self.Header["ER_ICE"])*self.Header["K"])
        Profile.Range = np.asarray([i for i in range(Nfft)]) * Profile.bin2m       
        Profile.Profile = Profile.Profile[0:math.floor(Nfft/2)-1]
        Profile.Range = Profile.Range[0:math.floor(Nfft/2)-1]
        if ref == 1:
          m = np.asarray([i for i in range(len(Profile.Profile))])/pad
          phiref = 2*math.pi*self.Header["CentreFreq"]*m/self.Header["B"] -\
             m * m * 2*math.pi * self.Header["K"]/2/self.Header["B"]**2
          Profile.Profile = Profile.Profile * np.exp(phiref*(-1j));
        
        Profile.BurstNo = self.BurstNo
        Profile.Header = self.Header
        Profile.Filename = self.Filename
        Profile.ChirpList = self.ChirpList
        Profile.pad = pad
        Profile.rad2m = self.Header["CentreFreq"]*math.sqrt(self.Header["c0"])/ \
            (4.*math.pi*self.Header["ER_ICE"])
        return(Profile)  
       
class ProfileObject:
    
    def __init__(self):
        self.Range = 0
        self.Profile = 0
        self.F0 = 2e8
        self.F1 = 4e8
        self.pad = 2
        self.ChirpList = 0
        self.Filename = ''
        self.BurstNo = 0
        self.Header = 0
        self.rad2m = 0
        self.bin2m = 0

    def PlotProfile(self,dmax):
        plt.plot(self.Range, 20*np.log10(np.abs(self.Profile)))
        plt.xlim(0, dmax)
        plt.xlabel("Range (m)")
        plt.ylabel("Amplitude (dB)")
        plt.grid("on")

        return(0)
        
    
