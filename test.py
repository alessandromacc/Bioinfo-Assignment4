import numpy
from markov import MarkovChain, GenomeInOutWindow, logRatioEvaluate
from utils import Generator, CpGInModel, CpGOutModel, FileHandler
import sys

'''
Markov chains bioinformatics assignment:
The file markov.py contains the code for implementing a Markov chain and evaluate query sequences with inside/outside models that can be trained starting from
appropriate sequence files, as well as the class for scanning a genome to test for loci score relative to inside/outside model, with an included plotting feature
for such genome-wide test.

The file utils.py contains the code for some useful purposes, such as constants, the Generator for random sequences or for sequences from files and the Model
system of the software, developed to be easily extended in custom features.

The test.py file contains the code for testing the Markov chain construction and performance in evaluating queries and scanning genomes. The CLI interface is sys
based, so the user can easily provide the software the desired paramenters and instructions for the actions to be performed. The interface only tests with
the human chromosome 22 CpG islands in/out models In detail, the user can declare
the following flags:

-q <str>: set the query string;
-f <str>: set the file path for sequence extraction;
-r : generate randomly the query sequence;
-l <int>: set the length of the random query to generate (default 1000);
-L , --log : use sum of log probabilities to evaluate a query (default False, non-changeable in scan mode);
-s : set the mode to genome scanning;
-P , --plot : enable data plotting when in scan mode (default False);
-w <int>: set the window size in scan mode;
-F , --fast: use pre-computed models for CpG islands instead of runtime train them (default False);
-S <int>: set the stringency for peak calling when plotting scan mode data (default 20);
-M , --mute : turns off unnecessary logging, which does not include the outcome of the requested operation (default True);
-k , --peak : enables peak calling on the positional scores in scan mode, including graphical representation;
--save <str>: when in scan mode, writes the scores and some parameters from the scannning in a file whose path has to be provided.

Mind that declaring scan mode overrides query declaration, however can be combined with the randomness flag to generate a completely random genome
or with the path flag to provide the file the user wants to extract the genome from. Furthermore, only declaring the path will result in
reading the whole file and using it as the genome, whereas declaring also the randomness flag will result in only extracting
a string from such genome.

The files used to pre-compute the Markov models are also made available, with chr22.fa containing the fasta sequence of hg19 chromosome 22, 
chr22_annot.txt containing the annotation of all the CpG islands found in hg19 chromosome 22, CpG.txt containing row by row the sequence of
all the CpG islands in chromosome 22 and outside.txt containing row by row random sequences extracted from chr22.fa to train the outside model.
The code for preparing these files from the original annotation (annot.txt) and sequence (chr22.fa) is in the file cpg_data_setup.py.

In general, the hierarchy of provided inputs for normal evaluation is:
1. Path - entire row extraction from file
2. Query
3. Random generation of modifiable length

Whereas for the scan mode:
1. Path - random fragment of modifiable length
2. Path - whole sequence from file
3. Random generation of modifiable length
x. Any declared query sequence here is ignored
'''

args = sys.argv[1:]
query = None
path = None
log_prob = False
scan = False
random = False
plot = False
l = 1000
wsize = None
fast = False
stringency = 20
logging = True
callPeaks = False
save = False
savename = None
read = False
readpath = None

for i in range(len(args)):
    if args[i] == '-q':
        query = args[i+1]
    elif args[i] == '-f':
        path = args[i+1]
    elif args[i] == '-r':
        random = True
    elif args[i] == '-l':
        l = int(args[i+1])
    elif args[i] == '-L' or args[i] == '--log':
        log_prob = True
    elif args[i] == '-s':
        scan = True
        read = False
    elif args[i] == '-P' or args[i] == '--plot':
        plot = True
    elif args[i] == '-w':
        wsize = int(args[i+1])
    elif args[i] == '-F' or args[i] == '--fast':
        fast = True
    elif args[i] == '-S':
        stringency = int(args[i+1])
    elif args[i] == '-M' or args[i] == '--mute':
        logging = False
    elif args[i] == '-k' or args[i] == '--peak':
        callPeaks = True
    elif args[i] == '--save':
        save = True
        savename = args[i+1]
    elif args[i] == '--read':
        read = True
        scan = False
        readpath = args[i+1]

if read:
    scores, l, wsize, stringency, *peaks = FileHandler.evaluationFromFile(readpath)
elif scan:
    if path != None and random:
        query = Generator.randomGenomeFromFile(path, l)
    elif path != None:
        query = Generator.genomeFromFile(path)
    elif random:
        query = Generator.randomGenome(l)
    else:
        raise ValueError('FlagError: scan mode declared, either random generation of filepath must be declared as well')
elif path != None:
    query = Generator.randomSequenceFromFile(path)
elif query != None:
    pass
elif random:
    query = Generator.randomSequence(l)
else:
    if query == None:
        raise ValueError('ArgumentError: no query sequence provided')

if not read:
    if fast:
        insidemod = MarkovChain(model = CpGInModel)
        outsidemod = MarkovChain(model = CpGOutModel)
    else:
        insidemod = MarkovChain(path = 'CpG.txt')
        outsidemod = MarkovChain(path = 'outside.txt')

    if not scan:
        if logging:
            print(f'Parameters set: scan: {scan}; filepath: {path}; random: {random}; length: {l}; fast: {fast}')
            print(f'Evaluating query sequence: "{query}"')
            print(f'Reference inside/outside model derived from hg19 chromosome 22 CpG islands sequences and other random chromosome 22 sequences')
        if log_prob:
            if logging:
                print('Using sum of the logarithms to evaluate query score')

        if not log_prob:
            insidescore = insidemod.scoreQuery(query)
            outsidescore = outsidemod.scoreQuery(query)
            print(f'Inside score:', insidescore)
            print(f'Outside score:', outsidescore)
            logratio = logRatioEvaluate(insidescore, outsidescore, True)
        else:
            insidescore = insidemod.scoreQuery(query, True)
            outsidescore = outsidemod.scoreQuery(query, True)
            print(f'Inside score:', insidescore)
            print(f'Outside score:', outsidescore)
            logratio = logRatioEvaluate(insidescore, outsidescore, False)

        print(f'Final log ratio evaluation: {logratio}')

    else:
        if wsize == None:
            if logging:
                print(f'Parameters set: scan: {scan}; filepath: {path}; random: {random}; length: {l}; fast: {fast}; plot: {plot}; peak call: {callPeaks}; stringency: {stringency}; window size: {insidemod.average_source_length}')
                print(f'Scanning {l} bases long genome, window size: {insidemod.average_source_length}, peak sharpness: {round(numpy.log(insidemod.average_source_length)*stringency)}, peak calling threshold: {round(numpy.log2(insidemod.average_source_length), 1)}', '\n')
        else:
            if logging:
                print(f'Parameters set: scan: {scan}; filepath: {path}; random: {random}; length: {l}; fast: {fast}; plot: {plot}; peak call: {callPeaks}; stringency: {stringency}; window size: {wsize}')
                print(f'Scanning {l} bases long genome, window size: {wsize}, peak sharpness: {round(numpy.log(wsize)*stringency)}, peak calling threshold: {round(numpy.log2(wsize), 1)}', '\n')
        data, wsize = GenomeInOutWindow.evaluate(query, insidemod, outsidemod, wsize, logging)
        call, scores = GenomeInOutWindow.callPeaks(data, wsize, stringency)
        print(f'Potential start sites identified [position, score]: {[[call[i], scores[i]] for i in range(len(call))]}')
        if save:
            FileHandler.writeEvaluation(savename, data, wsize, stringency, GenomeInOutWindow.callPeaks(data, wsize, stringency))
        if plot:
            GenomeInOutWindow.plotScore(data, wsize, stringency, callPeaks)
else:
    GenomeInOutWindow.quickPlot(scores, wsize, stringency, peaks)
