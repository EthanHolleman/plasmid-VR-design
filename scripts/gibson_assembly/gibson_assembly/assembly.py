from pydna.assembly import Assembly

class GibsonAssembler():

    def __init__(self, construct):
        self.construct = construct
    
    # this object should take over all assembly functions
    # also should implement reverse complmentation of inserted
    # regions if required as that is currently not implemented
    # additionally should search for regions of existing homology
    # before converting seqrecords to amplicon instances.
