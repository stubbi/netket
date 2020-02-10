import matplotlib.pyplot as plt
import collections
import pickle
import sys
import json
import pandas
import os
import shutil

# evaluation of a specific configuration

experimentFolder = sys.argv[1]
listSystemSizes = sys.argv[2].split(',')
listCycles = sys.argv[3].split(',')
numCircuits = int(sys.argv[4])
listOMPNodes = sys.argv[5].split(',')
listOMPTasks = sys.argv[6].split(',')
listOMPThreads = sys.argv[7].split(',')
listSamples = sys.argv[8].split(',')
listIterations = sys.argv[9].split(',')
listInitialHidden = sys.argv[10].split(',')
listSampleSteps = sys.argv[11].split(',')
numRuns = int(sys.argv[12])

class Evaluation:
    def __init__(self, experimentFolder, listSystemSizes, listCycles, numCircuits, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns):
        self.experimentFolder=experimentFolder
        self.listSystemSizes=listSystemSizes
        self.listCycles=listCycles
        self.numCircuits=numCircuits
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

    def loadHistogram(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        with open("{directory}/histogram.json".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)), 'r') as f:
            return json.load(f)

    def loadDuration(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        with open("{directory}/duration.time".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)), 'r') as f:
            return float(f.readline())

    def loadExact(self,qubits, cycles, circuit):
        with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/exact.json".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit), "rb") as f:
            return pickle.load(f, encoding='latin1')

    def numberOfHadamards(self, qubits, cycles, circuit):
        hadamards = 0
        with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/in.qc".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit), "r") as f:
            for line in f:
                if(line[0] == 'H'):
                    hadamards = hadamards + 1
        return hadamards

    def tvd(self, exact, histogram):
        tvd = 0.0
        for i in range(len(exact)):
            exact_prob = abs(exact[i])**2
            nqs_prob = histogram.get(str(i), 0.0)
            tvd += abs(exact_prob-nqs_prob)
        return tvd/2.0

    def plot(self, df, x, grouped, fixed):
        fig, ax = plt.subplots()
        groupDFBy = grouped + [x]
        df = df.groupby(groupDFBy, as_index = False).mean()
        filterBy = (df.groupby(grouped, as_index = False).mean())[grouped]
        df = df.loc[(df[list(fixed)] == pandas.Series(fixed)).all(axis=1)]

        for y in ['tvd', 'duration']:
            title = ''.join(['{} {} '.format(key, val) for key, val in fixed.items()])
            name = '{}_{}_{}'.format(x, y, title.replace(' ', '_'))
            for index, row in filterBy.iterrows(): 
                print(df)
                print(row.to_frame().T)
                print()
                toPlot = df.merge(row.to_frame().T, 'left')
                l = ''.join(['{} {} '.format(toPlot[[i]], i) for i in grouped])
                ax.plot(toPlot[[x]], toPlot[[y]], label = l)
            plt.legend()
            plt.title(title)
            plt.xlabel(x)
            plt.ylabel(y)
            plt.savefig('plots/{}.pdf'.format(name))
            plt.close()


    def generatePlots(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        if not os.path.exists('plots'):
            os.makedirs('plots')
        df = pandas.read_csv(results_file)
        df = df[df['success'] == True]
        df = df.astype({'tvd': 'float64', 'duration': 'float64'})

        for i in pandas.unique(df['#iterations']):
            for s in pandas.unique(df['#samples']):
                for t in pandas.unique(df['#sampleSteps']):
                    for x in ['#hadamards', '#qubits', '#cycles']:
                        grouped = []
                        if(x == '#qubits'):
                            grouped = ['#cycles']
                        if(x == '#cycles'):
                            grouped = ['#qubits']
                        if(x == '#hadamards'):
                            grouped = ['#cycles', '#qubits']
                        
                        fixed = {'#iterations': i,
                                 '#samples': s,
                                 '#sampleSteps': t}

                        d = df.copy()
                        self.plot(d, x, grouped, fixed)


        shutil.make_archive('plots', 'zip', 'plots')

          
    def generateCSV(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        with open(results_file, 'w') as f:
            f.write('#qubits,#cycles,circuit,#nodes,#tasks,#threads,#samples,#iterations,#initialHidden,#sampleSteps,run,#hadamards,tvd,duration,success\n')
        
        for size in self.listSystemSizes:
            for cycles in self.listCycles:
                for circuits in range(self.numCircuits):
                    hadamards = self.numberOfHadamards(size, cycles, circuits)
                    for nodes in self.listOMPNodes:
                        for tasks in self.listOMPTasks:
                            for threads in self.listOMPThreads:
                                for numSamples in self.listSamples:
                                    for numIterations in self.listIterations:
                                        for numInitialHidden in self.listInitialHidden:
                                            for numSampleSteps in self.listSampleSteps:
                                                for run in range(self.numRuns):

                                                    try:
                                                        histogram = self.loadHistogram(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)
                                                        tvd = '{:f}'.format(self.tvd(self.loadExact(size, cycles, circuits), self.normalise(histogram)))
                                                        duration = '{:f}'.format(self.loadDuration(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run))
                                                        success = True
                                                    except:
                                                        tvd = '-'
                                                        duration = '-'
                                                        success = False

                                                    line = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(size,cycles,circuits,nodes,tasks,threads,numSamples,numIterations,numInitialHidden,numSampleSteps,run,hadamards,tvd,duration,success)
                                                    with open(results_file, 'a') as f:
                                                        f.write(line)

    def directory(self, qubits, cycles, circuit, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        return "{experimentFolder}/{qubits}qubits/{cycles}cycles/circuit{circuit}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations/{initial_hidden}initialHidden/{sample_steps}sampleSteps/run{run}".format(
                experimentFolder=self.experimentFolder,
                nodes=nodes,
                tasks=tasks,
                threads=threads,
                samples=numSamples,
                iterations=numIterations,
                initial_hidden=numInitialHidden,
                sample_steps=numSampleSteps,
                run=run,
                qubits=qubits,
                cycles=cycles,
                circuit=circuit
            )

    def normalise(self, histogram):
        total = sum(histogram.values())
        for key in histogram:
            histogram[key] = float(histogram[key])/float(total)
        return histogram

ev = Evaluation(experimentFolder, listSystemSizes, listCycles, numCircuits, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns)
ev.generateCSV()
ev.generatePlots()