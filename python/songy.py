"""
.. module:: songy
    :synopsis: Python module to read binary output files from SONG.
.. moduleauthor:: Thomas Tram <thomas.tram@port.ac.uk>

This module defines a parent class called SongBinary as well as one child class for each
type of binary file: {FixedTauFile, Fixedk1k2File}. SongBinary holds the general functions
for manipulating the binary file, while the child class holds filetype specific methods
and the initialisation procedure.
"""
import numpy as np

class SongBinary:
    """
    This class holds the basic routines for manipulating the
    binary files created by SONG. Normally this class should
    not be initialised directly, but you should use the
    appropriate child class instead.
    Example:
    dat = songy.SongBinary(filename)
    """
    def __init__(self,filename):
        self.filename = filename
        self.filetype = ''

    def get_binary_mapping(self):
        """
        Extracts the filemap from the header and store it in self.mapping.
        Also performs a test that self.filetype matches the actual file.
        """
        f = open(self.filename)
        line = f.readline()
        #Check that the file type matches first line
        if self.filetype not in line:
            raise Exception('songy.'+self.__class__.__name__+' expects filetype <'+self.filetype+'>, but '+self.filename+' reports:\n'+line)

        while line.find('BLOCK')==-1:
            line = f.readline()
        line = f.readline()
        self.mapping = []
        while line[0]=='#':
            self.mapping.append(line)
            line = f.readline()
        f.close()
        
    def read_block(self,blocknumber):
        """
        Read block [blocknumber] from the binary file using the binary map
        computed by self.get_binary_mapping(). If no mapping is found,
        the function will call self.get_binary_mapping().
        """
        if not hasattr(self, 'mapping'):
            self.get_binary_mapping()
        mapelement = self.mapping[blocknumber];
        f = open(self.filename,'rb')
        mapelement=mapelement[1:].split()
        f.seek(int(mapelement[4]))
        if mapelement[1][0]=='c':
            val = np.fromfile(f,dtype='uint8', count=int(mapelement[2]), sep="")
        else:
            val = np.fromfile(f,dtype=mapelement[1][0], count=int(mapelement[2]), sep="")
        f.close()
        return val

    def get_name_array(self,blockstartnum):
        """
        Returns a list of names from three consecutive blocks in the binary file.
        The initial block is assumed to be the input parameter <blockstartnum>
        Block 1: Number of names <numsources>
        Block 2: Length of a name string <lennamestring>
        Block 3: Character array of size (<numsources>*<lennamestring>)
        """
        numsources = self.read_block(blockstartnum)[0]
        lennamestring = self.read_block(blockstartnum+1)[0]      
        names = self.read_block(blockstartnum+2)
        names = names.reshape(numsources,lennamestring)
        name_array = []
        for name in names:
            name = name[name.nonzero()].tostring()
            name_array.append(name)
        return name_array

class Fixedk1k2File(SongBinary):
    """
    This class is used for accessing SONG binary files for fixed k1 and k2.
    These files are usually called <sources_song_kXXX.dat>.
    Example:
    data = songy.Fixedk1k2File('output/sources_song_k001.dat')
    All information from the file is extracted as part of the
    object initialisation, except for the actual sources. These
    must be subsequent calls to
    data.get_source()
    """
    def __init__(self,filename):
        self.filetype = 'fixed k1 and k2'
        self.filename = filename
        self.get_binary_mapping()
        tmp = self.read_block(0)
        self.header = tmp[tmp.nonzero()].tostring()

        #Case of fized k1 and k2
        self.k1 = self.read_block(1)[0]
        self.k2 = self.read_block(2)[0]
        self.tau = self.read_block(4)
        self.a = self.read_block(5)
        self.H = self.read_block(6)
        self.k3 = self.read_block(8)
        self.mu = self.read_block(9)
        self.sourceblockstart = 19
        # Read first order perturbations into dictionaries
        name_array_1st = self.get_name_array(10)
        self.first_order_sources_k1 = {}
        self.first_order_sources_k2 = {}
        sourcek1 = self.read_block(13).reshape(len(self.tau),len(name_array_1st))
        sourcek2 = self.read_block(14).reshape(len(self.tau),len(name_array_1st))
        for i in range(len(name_array_1st)):
            name = name_array_1st[i]
            self.first_order_sources_k1[name] = sourcek1[:,i]
            self.first_order_sources_k2[name] = sourcek2[:,i]
        #Store names of second order perturbations:
        self.sourcenames = self.get_name_array(15)

    def get_source(self,sourcename):
        """
        Returns the source with name <sourcename>. For a list of all
        possible source names, see self.sourcenames.
        The output will be a [Ntau x Nk] matrix.
        """
        for i in range(len(self.sourcenames)):
            if sourcename in self.sourcenames[i]:
                val = self.read_block(19+i)
                #Case of fized k1 and k2. Reshape matrix:
                return val.reshape(len(self.tau),len(self.k3))
    

class FixedTauFile(SongBinary):
    """
    This class is used for accessing SONG binary files for fixed tau or z.
    These files are usually called <sources_song_tauXXX.dat> or
    <sources_song_zXXX.dat>.
    Example:
    data = songy.FixedTauFile('output/sources_song_z001.dat')
    All information from the file is extracted as part of the
    object initialisation, except for the actual sources. These
    must be subsequent calls to
    data.get_source()
    The format of self.k3 is a list of length N*(N+1)/2 where N is
    the number of k1 values. Each element of the list will contain
    a 1 dimensional array of k3 values for a specific (k1,k2)-pair.
    The flattened index in the list can be recovered from index_k1
    and index_k2 by self.flatidx[index_k1,index_k2].
    """
    def __init__(self,filename):
        self.filetype = 'fixed tau'
        self.filename = filename
        self.get_binary_mapping()
        tmp = self.read_block(0)
        self.header = tmp[tmp.nonzero()].tostring()
        
        #Case of fixed tau
        self.k1 = self.read_block(7)
        self.k2 = [self.k1[0:i+1] for i in range(len(self.k1))]
        self.tau = self.read_block(2)[0]
        self.z = self.read_block(1)[0]
        self.a = self.read_block(4)[0]
        self.H = self.read_block(5)[0]
        #k3 is naturally a list of lists of numpy arrays. However, it is easier to deal with a single list,
        #so we will use a 2d array to access the list.
        #Create lower triangular matrix such that self.flatidx[index_k1][index_k2] is the index in the flattened list.
        N = len(self.k1)
        self.flatidx = -np.ones((N,N),dtype='int')
        self.flatidx[np.tril_indices(N)] = range(0,N*(N+1)/2)
        #Read the flattened lower triangular matrix of k3 sizes and form the cumulative sum.
        #We use [:-1] not to get the total number, which would result in an empty list after np.split()
        self.k3sizes_cumsum = self.read_block(8)[:-1].cumsum()
        #Read packed k3 array and split it into a (flattened) lower triangular list:
        self.k3 = np.split(self.read_block(9),self.k3sizes_cumsum)
        # Read first order perturbations into dictionaries
        name_array_1st = self.get_name_array(10)
        self.first_order_sources = {'k':self.read_block(7)}
        source = self.read_block(13).reshape(len(self.first_order_sources['k']),len(name_array_1st))
        for i in range(len(name_array_1st)):
            self.first_order_sources[name_array_1st[i]] = source[:,i]
        #Store names of second order perturbations:
        self.sourcenames = self.get_name_array(14)

    def get_source(self,sourcename):
        """
        Returns the source with name <sourcename>. For a list of all
        possible source names, see self.sourcenames.
        The output will be the same format as self.k3, i.e. a list of
        length N*(N+1)/2 where N is the number of k1 values. Each
        element of the list will be the source as a function of k3.
        The flattened index in the list can be recovered from index_k1
        and index_k2 by self.flatidx[index_k1,index_k2].
        """
        for i in range(len(self.sourcenames)):
            if sourcename in self.sourcenames[i]:
                val = self.read_block(17+i)
                #Case of fixed tau: do same split as k3:
                return np.split(val,self.k3sizes_cumsum)
                

    
