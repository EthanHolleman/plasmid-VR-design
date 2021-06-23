import os

class Primer3File():

    def __init__(self, filepath):
        self.filepath = filepath
    
    @property
    def filepath:
        return self._filepath
    
    @filepath.setter
    def filepath(self, newpath):
        if os.path.isfile(newpath):
            self._filepath = newpath
        else:
            raise TypeError('Most be valid file path!')
    
    def parse(self):
        pass