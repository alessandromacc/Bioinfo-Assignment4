import pandas as pd
import numpy
from utils import Constants, Model, StandardModel, NoModel
from matplotlib import pyplot as plt

class MarkovModel(Model):
    '''The class MarkovModel implements a markov chain of order 1, in a structurate object. It can generate
    the model from sequences in a file, provided as argument path, or from a single sequence provided as the 
    argument s, or a model compliant with the standards of the Model class in utils.py can be provided, from which
    the data about the model can be quickly retrieved. Providing a model as argument will override other arguments parsed.'''
    def __init__(self, s: str = None, path: str = None, model: StandardModel = NoModel) -> None:
        if not model.isModel:
            if s == None and path == None:
                raise ValueError('ArgumentError: no input string or filepath provided')
            elif path == None:
                s = s.upper()
                self.source = s
                self.average_source_length = len(self.source)
            elif s == None:
                self.source = []
                totrows = 0
                filez = open(path, 'r')
                for line in filez:
                    line = line[:-1].upper()
                    self.source.append(line)
                    totrows += 1
                filez.close()
                self.source = ''.join(self.source)
                self.average_source_length = round(len(self.source)/totrows)
            else:
                raise ValueError('ArgumentError: both input string and filepath provided, only 1 required')
            counter = {}
            mod = pd.DataFrame([[0 for i in range(4)] for j in range(4)], columns=Constants.nucleotides, index=Constants.nucleotides)
            for nucl in range(1, len(self.source)):
                if self.source[nucl-1] in Constants.nucleotides and self.source[nucl] in Constants.nucleotides:
                    mod.loc[self.source[nucl-1], self.source[nucl]] += 1
                    counter[self.source[nucl]] = counter.get(self.source[nucl], 0) + 1
            for i in Constants.nucleotides:
                for j in mod.index:
                    mod.loc[j, i] = round(mod.loc[j,i]/counter[j], 2)
            self.model = mod
        else:
            self.model = model.model
            self.average_source_length = model.average_source_length
    
    def scoreQuery(self, q: str, log: bool = False) -> float:
        '''Outputs the score of a query sequence q, evaluated on the model represented by the class instance.
        By parsing bool = True one makes the method switch to calculating such score by sum of logarithms,
        instead of by product of probabilities from the model.'''
        q = q.upper()
        if not log:
            score = 0.25
            for i in range(1, len(q)):
                score *= self.model[q[i]][q[i-1]]
            return score
        else:
            score = numpy.log(0.25)
            for i in range(1, len(q)):
                score += numpy.log(self.model[q[i]][q[i-1]])
            return score

class GenomeInOutWindow:
        '''Class designed to extend the evaluation of query sequences to larger sequences, which can be assumed to represent genomes.
        Provides the necessary methods for scanning such genomes with a sliding window, evaluating each position and representing
        graphically the results.'''
        @staticmethod
        def evaluate(genome: str, inmod: MarkovModel, outmod: MarkovModel, wsize: int = None, log: bool = True) -> tuple:
            '''Receives as input the reference genome to scan as a string, the inside model and the outside model, biologically significant, both required to be
            an instance of the MarkovModel class for compatibility, the size of the sliding window, which can affect the scores produced and a boolean "log", which
            will turn off unnecessary logging if set to False.
            The method returns to the user a tuple containing as first element a list with the score corresponding to the window, of length "wsize",
            starting from each legal position of the genome, and the window size itself, since if not parsed, it needs to be set before proceding.
            Mind that in case the window size was greater than the genome size, the method would automatically raise an error and notify the user.'''
            if wsize == None:
                wsize = inmod.average_source_length
            data = []
            if wsize > len(genome):
                raise ValueError('GenomeLengthError: window size exceeds genome length')
            for s in range(0, len(genome)-wsize+1):
                score = round(logRatioEvaluate(inmod.scoreQuery(genome[s:s+wsize], True), outmod.scoreQuery(genome[s:s+wsize], True), False), 1)
                if log:
                    if score > 0:
                        print(f'Genome position: {s} : {s+wsize}')
                        print(f'Window score: {score}')
                        print('-'*40)
                data.append(score)
            return data, wsize
        
        @staticmethod
        def quickPlot(dataArray: list|tuple, ws: int, stringency: int = 20, peaks = None) -> None:
            '''Useful when the user is simply in need of plotting data already known. Takes an array of scores,
            the window size of the genome scanning, the stringency of the peak calling and the location of the peaks.
            Devised to quicjly plot the informations from a file written by FileHandler.writeEvaluation'''
            plt.plot([i for i in range(len(dataArray))], dataArray, linewidth = 1)
            plt.xlabel('Genome starting position')
            plt.ylabel('CpG island probability')
            plt.title(f'Genome scanning, width: {ws}, peak sharpness: {round(numpy.log(ws)*stringency)}'+'\n'+f'peak calling threshold: {round(numpy.log2(ws), 1)}')
            plt.axhline(0, 0, color = 'black', linewidth = 2)
            for i in peaks:
                plt.axvline(i, color = 'red', linewidth = 0.5)
            plt.axhline(numpy.log2(ws), color = 'green', linewidth=0.8)
            plt.yticks(list(plt.yticks()[0]) + [numpy.log2(ws)])
            plt.grid(True)
            plt.show()
        
        @staticmethod
        def plotScore(dataArray: list|tuple, ws: int, stringency: int = 20, peaks: bool = False) -> None:
            '''The method is used for plotting the data one produced with genome scanning. The parameters include "dataArray", the list of scores
            for each position of the genome, "ws", the window size, employed in peak calling, and the stringency, again in peak calling.'''
            plt.plot([i for i in range(len(dataArray))], dataArray, linewidth = 1)
            plt.xlabel('Genome starting position')
            plt.ylabel('CpG island probability')
            plt.title(f'Genome scanning, width: {ws}, peak sharpness: {round(numpy.log(ws)*stringency)}'+'\n'+f'peak calling threshold: {round(numpy.log2(ws), 1)}')
            plt.axhline(0, 0, color = 'black', linewidth = 2)
            if peaks:
                putativeStartSites = GenomeInOutWindow.callPeaks(dataArray, ws, stringency)[0]
                for j in putativeStartSites:
                    plt.axvline(j, color = 'red', linewidth=0.5)
            plt.axhline(numpy.log2(ws), color = 'green', linewidth=0.8)
            plt.yticks(list(plt.yticks()[0]) + [numpy.log2(ws)])
            plt.grid(True)
            plt.show()
        
        @staticmethod
        def callPeaks(dataArray: list|tuple, ws: int, stringency: int = 20) -> tuple:
            '''Algorithm for peak calling; receives the scores data as parameter, as well as the window size and the set stringency of
            the operation, which are used to compute a probability threshold for peaks and height. Returns a tuple containing the list of the possible
            start sites of the considered genomic feature identified and the list of their scores. Fine tuning of the stringency
            value is recommended.'''
            lastMax = 0
            lastMin = 0
            putativeStartSites = []
            potentialMax = 0
            potentialMin = 0
            for i in range(len(dataArray)):
                if dataArray[i] > dataArray[potentialMax]:
                    potentialMax = i
                    if dataArray[i] - dataArray[potentialMin] >= numpy.log(ws)*stringency:
                        lastMin, potentialMin = potentialMin, i
                if dataArray[i] < dataArray[potentialMin]:
                    potentialMin = i
                    if dataArray[i] - dataArray[potentialMax] <= -numpy.log(ws)*stringency and dataArray[potentialMax] >= numpy.log2(ws):
                        lastMax, potentialMax = potentialMax, i
                        if lastMax not in putativeStartSites:
                            putativeStartSites.append(lastMax)
            #if dataArray[potentialMax] - dataArray[lastMin] >= numpy.log(ws)*stringency/2 and dataArray[potentialMax] >= numpy.log2(ws):
            #    putativeStartSites.append(potentialMax)
            return putativeStartSites, [round(dataArray[i], 1) for i in putativeStartSites]


def logRatioEvaluate(s1: float, s2: float, log: bool = True) -> float:
    '''Computes the log ratio evaluation from the scores of a sequence tested on two alternative models.'''
    if log:
        return numpy.log2(s1) - numpy.log2(s2)
    else:
        return round(s1 - s2, 2)