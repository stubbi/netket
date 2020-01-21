import matplotlib.pyplot as plt
import collections
import pickle
import sys
import json

# evaluation of a specific configuration

experimentFolder = sys.argv[1]
systemSizes = sys.argv[2].split(',')
listOMPNodes = sys.argv[3].split(',')
listOMPTasks = sys.argv[4].split(',')
listOMPThreads = sys.argv[5].split(',')
listSamples = sys.argv[6].split(',')
listIterations = sys.argv[7].split(',')
listInitialHidden = sys.argv[8].split(',')
listSampleSteps = sys.argv[9].split(',')
numRuns = int(sys.argv[10])

class Evaluation:
    def __init__(self, experimentFolder, systemSizes, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns):
        self.experimentFolder=experimentFolder
        self.systemSizes=systemSizes
        self.listOMPNodes=listOMPNodes
        self.listOMPTasks=listOMPTasks
        self.listOMPThreads=listOMPThreads
        self.listSamples=listSamples
        self.listIterations=listIterations
        self.listInitialHidden=listInitialHidden
        self.listSampleSteps=listSampleSteps
        self.numRuns=numRuns

    def plotTVDIterations(self):
        pass

    def plotTVDSamples(self):
        pass

    def plotTVDSamplesIterations(self):
        pass

    def loadHistograms(self, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps):
        histograms=[]
        for r in range(self.numRuns):
            with open("{directory}/histogram.json".format(directory=self.directory(nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, r)), 'r') as f:
                histograms.append(json.load(f))
        return histograms

    def loadExact(self):
        with open("{directory}/exact.json".format(directory=self.experimentFolder), "rb") as f:
             return pickle.load(f)

    def tvd(self, exact, histogram):
        tvd = 0
        for i in range(len(exact)):
            exact_prob = abs(exact[i])**2
            nqs_prob = histogram.get(i, 0.0)
            tvd += abs(exact_prob-nqs_prob)
        return tvd
        
    def generateAll(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        with open(results_file, 'w') as f:
            f.write('system_size,nodes,tasks,threads,numSamples,numIterations,numInitialHidden,numSampleSteps,tvd')
        
        for size in self.systemSizes:
            for nodes in self.listOMPNodes:
                for tasks in self.listOMPTasks:
                    for threads in self.listOMPThreads:
                        for numSamples in self.listSamples:
                            for numIterations in self.listIterations:
                                for numInitialHidden in self.listInitialHidden:
                                    for numSampleSteps in self.listSampleSteps:

                                        histograms = self.loadHistograms(nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps)
                                        tvd = self.tvd(self.loadExact(), self.mergeAndNormalise(histograms))

                                        line = "{},{},{},{},{},{},{},{},{}".format(size, nodes,tasks,threads,numSamples,numIterations,numInitialHidden,numSampleSteps,tvd)
                                        with open(results_file, 'a') as f:
                                            f.write(line)

    def directory(self, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        return "{experimentFolder}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations/{initial_hidden}initialHidden/{sample_steps}sampleSteps/run{run}".format(
                experimentFolder=self.experimentFolder,
                nodes=nodes,
                tasks=tasks,
                threads=threads,
                samples=numSamples,
                iterations=numIterations,
                initial_hidden=numInitialHidden,
                sample_steps=numSampleSteps,
                run=run
            )

    def mergeAndNormalise(self, histograms):
        c = collections.Counter()
        for h in histograms:
            c = c+h
        total = sum(c.values)
        for key in c:
            c[key] /= total
        return c

ev = Evaluation(experimentFolder, systemSizes, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns)
ev.generateAll()