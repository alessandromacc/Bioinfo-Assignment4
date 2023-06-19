import random
import pandas as pd

class Constants:
    nucleotides = ['A', 'C', 'G', 'T']

class Generator:
    '''Utility class that contains useful methods for generating random sequences or retrieving entire ones from files.'''
    @staticmethod
    def randomSequence(l: int) -> str:
        '''Random generation of a nucleotides string, length "l" has to be parsed as parameter.'''
        return ''.join(random.choices(Constants.nucleotides, k=l))
    
    @staticmethod
    def randomSequenceFromFile(path: str) -> str:
        '''By providing a file path as argument, the method retrieves all the sequences stored by rows in the file and returns a random one among them.'''
        file = open(path, 'r')
        lines =[i.strip('\n') for i in file.readlines()]
        file.close()
        return random.choice(lines)
    
    @staticmethod
    def randomGenome(l: int = 1000) -> str:
        '''Random generation of a genome, useful in scanning since the length parameter is by default set to 1000.'''
        return ''.join(random.choices(Constants.nucleotides, k=l))
    
    @staticmethod
    def randomGenomeFromFile(path: str, l: int = 1000) -> str:
        '''Generates a genome from the whole concatenated sequence content of the file "path", starting at a random position
        and with length "l", by default set to 1000. Takes care not to have any unassigned nucleotides in the chosen sequence.'''
        sf = open(path, 'r')
        sequence = ''.join([i.strip('\n') for i in sf.readlines()]).upper()
        sf.close()
        rsp = random.randint(0, len(sequence)-l-1)
        while 'N' in sequence[rsp:rsp+l]:
            rsp = random.randint(0, len(sequence)-l-1)
        return sequence[rsp:rsp+l]
    
    @staticmethod
    def genomeFromFile(path: str) -> str:
        '''Generates a genome from the whole sequence content of the file "path", concatenated. Mind the size of the èarsed file!'''
        sf = open(path, 'r')
        return ''.join([i.strip('\n') for i in sf.readlines()]).upper()

class Model:
    '''Base class for creation of more specific models; contains the attribute isModel, used for compatibility with the MarkovModel class
    if provided as a model. The model system can easily be extended as shown below by the custom classes present.'''
    isModel = True

class NoModel:
    '''Counterpart of the Model class, used as default in MarkovModel, but can be replaced by parsing a Model.'''
    isModel = False

class StandardModel(Model):
    '''Subclass of Model, inherits the isModel attribute, and can be used to create a custom model by the user,
    which can use this class as interface to host data representing the computed model, ensuring compatibility
    with the MarkovModel class when loading such model on an active spot for query evaluation or scanning.'''
    def __init__(self, mod: pd.DataFrame = None, avl: int = None) -> None:
        self.model = None
        self.average_source_length = None


class CpGInModel(Model):
    '''Inside model for CpG island, pre-computed and stored in a StandardModel-like object. Fast access for repetitive executions.'''
    model = pd.DataFrame([[0.19, 0.28, 0.40, 0.14],[0.19, 0.36, 0.25, 0.20],[0.17, 0.33, 0.36, 0.14],[0.09, 0.34, 0.38, 0.19]], columns=Constants.nucleotides, index=Constants.nucleotides)
    average_source_length = 566


class CpGOutModel(Model):
    '''Inside model for CpG island, pre-computed and stored in a StandardModel-like object. Fast access for repetitive executions.'''
    model = pd.DataFrame([[0.29, 0.20, 0.29, 0.23],[0.32, 0.29, 0.07, 0.31],[0.26, 0.23, 0.29, 0.21],[0.18, 0.23, 0.29, 0.29]], columns=Constants.nucleotides, index=Constants.nucleotides)
    average_source_length = 566
