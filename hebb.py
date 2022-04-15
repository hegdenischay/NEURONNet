from neuron import h
from neuron.units import ms,mV

h.load_file('stdrun.hoc')

# M Cell for M1, G1 in 4-coupled

class MCELL:
    
    def MCell(self, gid, M):
        self._gid: int = gid
        self.M: int = M+1
        
        # Set the morphology by setting soma, dendrite and axon
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.axon = h.Section(name='axon', cell=self)

        self.all = [self.axon, self.dend, self.soma]

        # explicitly connect cells in the way we intend to

        self.dend.connect(self.soma, 1, 0)
        self.axon.connect(self.soma, 0, 0)

        self._spike_detector = h.NetCon(self.axon(0.5)._ref_v, None, sec=self.axon)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        self.axon_v = h.Vector().record(self.axon(0.5)._ref_v)
        self._ncs = []

        #Defining geometry of soma
        self.soma.L = 18.8
        self.soma.diam = 18.8 #in microns
        self.soma.nseg = 1                  #No. of segments

        #Defining geometry of dend
        self.dend.nseg = 1                  #No. of segments
        self.dend.L = 701.9                 #in microns
        self.dend.diam = 3.18                #in microns
        self.dendexcisyn = h.ExpSyn(self.dend(0.5))
        self.dendexcisyn.tau = 1 *ms   # tau is decay time constant
        self.dendexcisyn.e = 0    # reversal potential
        self.dendinhisyn = h.ExpSyn(self.dend(0.1))
        self.dendinhisyn.tau = 8
        self.dendinhisyn.e = -70

        #Defining geometry for axon
        self.axon.nseg = 1
        self.axon.L = 152
        self.axon.diam = 3.18
        self.axonexcisyn = h.ExpSyn(self.axon(0.8))
        self.axonexcisyn.tau = 2 #Decay time constant
        self.axonexcisyn.e = 0 #Reversal potential
        self.axoninhisyn = h.ExpSyn(self.axon(0.1))
        self.axoninhisyn.tau = 8
        self.axoninhisyn.e = -70

        #Setting biophysics
        for sec in self.all:
            sec.Ra = 123    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2

        self.soma.insert('hh')            #Inserting HH neurons
        self.axon.insert('hh')
        self.dend.insert('pas')
        for seg in self.dend:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -78    # Leak reversal potential mV

    def __init__(self, gid, M):
        self.MCell(gid, M)
            
    def __repr__(self):
        return 'Set [{}]_Mcell [{}]'.format(self.M,self._gid)
        #This shows how to represent each part when called upon

#This is the ORN class
class ORN:
    def __init__(self,gid,M):
        self._gid = gid
        self.M = M+1

        # Set morphology
        self.soma = h.Section(name='soma', cell=self)
        self.axon = h.Section(name='axon', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dendriticknob = h.Section(name='dendriticknob', cell=self)
        self.ciliumArr = [h.Section(name="cilium%d" % i, cell=self) for i in range(4)]
        self.all = [self.soma, self.axon, self.dend, self.dendriticknob]
        self.all.extend(self.ciliumArr)
        self.axon.connect(self.soma(0),0)
        self.dend.connect(self.soma(1),0)
        
        self.dendriticknob.connect(self.dend(0),1)
        
        self.dendriticknob.nseg = 1
        self.dendriticknob.diam = 2
        self.dendriticknob.L = 2
        self.dendriticknob.insert('ciliaProp')

        for i in range(4):
            self.ciliumArr[i].connect(self.dendriticknob(1),0)

        self._spike_detector = h.NetCon(self.axon(0.5)._ref_v, None, sec=self.axon)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        self.axon_v = h.Vector().record(self.axon(0.5)._ref_v)
        self._ncs = [] 

        # anatomical and biophysical properties
        self.soma.nseg = 1 
        self.soma.L = 9 # micrometer
        self.soma.diam = 6
        self.soma.insert('hh1')

        self.axon.nseg = 1
        self.axon.L = 100
        self.axon.diam = 1
        self.axon.insert('hh')
    
        self.dend.nseg = 1
        self.dend.L = 50
        self.dend.diam = 1.5
        self.dend.insert('dendProp')
        self.dendexcisyn = h.ExpSyn(self.dend(0.5))
        self.dendexcisyn.tau = 1 *ms   # tau is decay time constant
        self.dendexcisyn.e = 0    # reversal potential
        # self.dendArr[0].e_pas = -65
        # self.dendArr[0].g_pas = 0.001

        for i in range(4):
            self.ciliumArr[i].nseg = 1
            self.ciliumArr[i].diam = 0.18
            self.ciliumArr[i].L = 200
            self.ciliumArr[i].insert("blr300%d" % i)
        # print(dir(self))
        for i in range(4):
            self.stim = h.IClamp(self.ciliumArr[i](0.5)) 
            # put it in middle of all cilia
            self.stim.delay = 1    # [ms] delay
            self.stim.dur = 60   # [ms] duration
            self.stim.amp = 100
            
        for sec in self.all:
            sec.Ra = 123    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
#         for seg in self.dend:
#             seg.pas.g = 0.001  # Passive conductance in S/cm2
#             seg.pas.e = -78    # Leak reversal potential mV

        self.dendriticknob(0.5).ciliaProp._ref_cilMemPot1 = self.ciliumArr[0](0.5).blr3000._ref_memPot
        self.dendriticknob(0.5).ciliaProp._ref_cilMemPot2 = self.ciliumArr[1](0.5).blr3001._ref_memPot
        self.dendriticknob(0.5).ciliaProp._ref_cilMemPot3 = self.ciliumArr[2](0.5).blr3002._ref_memPot
        self.dendriticknob(0.5).ciliaProp._ref_cilMemPot4 = self.ciliumArr[3](0.5).blr3003._ref_memPot
        self.dend(0.5).dendProp._ref_ciliaMemPoten = self.dendriticknob(0.5).ciliaProp._ref_ciliaMemPot
        self.soma(0.5).hh1._ref_dMemPot = self.dend(0.5).dendProp._ref_dmemPot

        # for sect in self.ciliumArr:
            # odStim = h.IClamp(0.5, sec=sect)
            # odStim.delay = 1
            # odStim.dur = 60
            # odStim.amp = 0.9

        self.tstop = 6
        
    def __repr__(self):
        return 'Set [{}]_ORNcell [{}]'.format(self.M,self._gid)
        #This shows how to represent each part when called upon

class GCELL:
    def __init__(self, gid, M):
        self._gid = gid #Neuron no.
        self.M = M+1
        #Setting morphology
        #Creating soma,dend and axon
        self.soma=h.Section(name='soma',cell=self)
        self.dend=h.Section(name='dend',cell=self)
        self.axon=h.Section(name='axon',cell=self)
        
        self.all = [self.axon, self.soma, self.dend] 
        #list of all the sections in the cell.
        #We could explicitly specify the connection location  by self.dend.connect(self.soma(0.5))
        
        self.dend.connect(self.soma,1,0)
        self.axon.connect(self.soma,0,0)
        
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        self.soma_v = h.Vector().record(self.soma(0.5)._ref_v)
        self._ncs = [] 
        
        #Defining geometry of soma
        self.soma.L = self.soma.diam = 30 #in microns
        self.soma.nseg = 1                  #No. of segments
        
        #Defining geometry of dend
        self.dend.nseg = 1                  #No. of segments
        self.dend.L = 100                 #in microns
        self.dend.diam = 3.18                #in microns
        self.dendexcisyn = h.ExpSyn(self.dend(0.8))
        self.dendexcisyn.tau = 2 *ms
        self.dendexcisyn.e = 0 
        self.dendinhisyn = h.ExpSyn(self.dend(0.1))
        self.dendinhisyn.tau = 8
        self.dendinhisyn.e = -70
        
        #Defining geometry for axon
        self.axon.nseg = 1
        self.axon.L = 100
        self.axon.diam = 2.18
        self.axonexcisyn = h.ExpSyn(self.axon(0.8))
        self.axonexcisyn.tau = 2 #Decay time constant
        self.axonexcisyn.e = 0 #Reversal potential
        self.axoninhisyn = h.ExpSyn(self.axon(0.1))
        self.axoninhisyn.tau = 8
        self.axoninhisyn.e = -70
        
        #Setting biophysics
        for sec in self.all:                                                  
            sec.Ra = 100    # Axial resistance in Ohm * cm                    
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        
        self.soma.insert('hh')            #Inserting HH neurons
        self.axon.insert('hh')
        self.dend.insert('pas')
        for seg in self.dend:                               
            seg.pas.g = 0.001  # Passive conductance in S/cm2 
            seg.pas.e = -78    # Leak reversal potential mV
        
    def __repr__(self):
        return 'Set [{}]_Gcell [{}]'.format(self.M,self._gid)
        #This shows how to represent each part when called upon
        
