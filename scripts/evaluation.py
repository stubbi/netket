import matplotlib.pyplot as plt
import collections
import pickle
import sys
import json
import pandas
import os
import shutil
import numpy as np
import math
import nqs as nq

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
            return pickle.load(f)
    
    def loadRBM(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        f = "{directory}/parameters.json".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run))
        nqs = nq.nqs.NQS(size,numInitialHidden,numSamples)
        print('start load')
        nqs.load(f)
        print('finished load')
        return nqs

    def numberOfHadamards(self, qubits, cycles, circuit):
        hadamards = 0
        with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/in.qc".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit), "r") as f:
            for line in f:
                if(line[0] == 'H'):
                    hadamards = hadamards + 1
        return hadamards

    def tvd(self, exact, rbm):
        def toBinaryArray(i):
            return [int(b) for b in format(i, '0{}b'.format(math.log(len(exact),2)))]

        tvd = 0.0
        for i in range(len(exact)):
            exact_prob = abs(exact[i])**2
            print(toBinaryArray(i))
            nqs_prob = abs(rbm.psi(toBinaryArray(i)))**2
            print(nqs_prob)
            tvd += abs(exact_prob-nqs_prob)
        return tvd/2.0

    def f_xeb(self, exact, histogram, qubits):
        f_xeb = 0.0
        shots = 0
        for h in histogram:
            shots += histogram[h]
            f_xeb += histogram[h] * (2**qubits * abs(exact[int(h)])**2 + 1)
        return f_xeb/shots

    def plotPDF(self, qubits, cycles, circuit):
        exact = self.loadExact(qubits, cycles, circuit)
        exact_probs_sorted = sorted([len(exact) * abs(e)**2 for e in exact])

        fig, ax = plt.subplots()
        ax.plot(range(len(exact_probs_sorted)), exact_probs_sorted, label = 'no errors')
        plt.legend()
        plt.suptitle(self.experiment(), fontsize=14, fontweight='bold')
        plt.title('{} qubits {} cycles circuit {} entropy: {} Porter-Thomas: {}'.format(qubits, cycles, circuit, self.circuitEntropy(qubits, cycles, circuit), self.porterThomasEntropy(qubits)), fontdict={'size':10})
        plt.ylabel('p(j)')
        plt.xlabel('Bit string index j (ordered)')
        plt.savefig('plots/circuits/pdf_{}qubits_{}cycles_circuit{}.pdf'.format(qubits, cycles, circuit))
        plt.close()

    def plotDepthEntropy(self, qubits):
        entropies = [np.average([self.circuitEntropy(qubits, cycles, circuit) for circuit in range(self.numCircuits)]) for cycles in self.listCycles]

        fig, ax = plt.subplots()
        ax.plot(self.listCycles, entropies)
        ax.axhline(self.porterThomasEntropy(qubits))
        plt.suptitle(self.experiment(), fontsize=14, fontweight='bold')
        plt.title('{} qubits'.format(qubits), fontdict={'size':10})
        plt.ylabel('Entropy')
        plt.xlabel('Depht')
        plt.savefig('plots/circuits/entropy_{}qubits.pdf'.format(qubits))
        plt.close()

    def circuitEntropy(self, qubits, cycles, circuit):
        exact = self.loadExact(qubits, cycles, circuit)
        exact_probs = [abs(e)**2 for e in exact]
        return np.sum([e * math.log(e) if e != 0 else 0 for e in exact_probs])

    def porterThomasEntropy(self, qubits):
        return math.log(2**int(qubits)) - 1.0 + 0.577

    def experiment(self):
        return '{}\n'.format(self.experimentFolder.split('/')[-1])

    def plot(self, df, x, grouped, fixed):
        all_keys = ['#qubits', '#cycles', '#iterations', '#samples', '#sampleSteps', 'initialHidden', '#hadamards'] 
        experiment = self.experiment()
        title = ', '.join(['{}:{}'.format(key, fixed.get(key, 'all')) for key in all_keys if key not in grouped + [x]])

        filterBy = df[grouped].drop_duplicates()

        groupDFBy = grouped + [x]
        df = df.loc[(df[list(fixed)] == pandas.Series(fixed)).all(axis=1)]
        df = df.groupby(groupDFBy, as_index = False).mean() 

        for y in ['tvd', 'duration', 'f_xeb']:
            fig, ax = plt.subplots()
            name = '{}_{}_{}'.format(x.replace('#', ''), y.replace('#', ''), title.replace(', ', '_').replace('#', ''))
            for index, row in filterBy.iterrows(): 
                toPlot = pandas.merge(df, row.to_frame().T, how='inner')
                if(toPlot.shape[0] > 0):
                    l = ''.join(['{}:{} '.format(i, toPlot[i].tolist()[0]) for i in grouped])
                    ax.plot(toPlot[x].tolist(), toPlot[y].tolist(), label = l)
            plt.legend()
            plt.suptitle(experiment, fontsize=14, fontweight='bold')
            plt.title(title, fontdict={'size':10})
            plt.xlabel(x)
            plt.ylabel(y)
            plt.savefig('plots/{}.pdf'.format(name))
            plt.close()

    def generateReport(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        df = pandas.read_csv(results_file)

        report_file = "{directory}/report.txt".format(directory=self.experimentFolder)
        with open(report_file, 'w') as f:
            f.write('experiment: {}\n\n'.format(self.experimentFolder.split('/')[-1]))
            f.write('total number of experiments: {}\n'.format(df.shape[0]))
            f.write('number of succeeded experiments: {}\n'.format(df[df['success'] == True].shape[0]))
            f.write('number of failed experiments: {}\n\n'.format(df[df['success'] == False].shape[0]))
            f.write('tested number of qubits: {}\n'.format(self.listSystemSizes))
            f.write('tested number of cycles: {}\n'.format(self.listCycles))
            f.write('tested number of samples: {}\n'.format(self.listSamples))
            f.write('tested number of iterations: {}\n'.format(self.listIterations))
            f.write('tested number of sample steps: {}\n'.format(self.listSampleSteps))
            f.write('number of runs: {}\n\n'.format(self.numRuns))
            f.write('failed combinations:\n')
            failed = df[df['success'] == False].copy().drop(columns=['tvd', 'duration', 'f_xeb']).drop_duplicates()
            for index, row in failed.iterrows():
                f.write('{} \n'.format(row))


    def generatePlots(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        if not os.path.exists('plots'):
            os.makedirs('plots')
        if not os.path.exists('plots/circuits'):
            os.makedirs('plots/circuits')
        df = pandas.read_csv(results_file)
        df = df[df['success'] == True]
        df = df.astype({'tvd': 'float64', 'duration': 'float64', 'f_xeb': 'float64'})

        for i in pandas.unique(df['#iterations']):
            for s in pandas.unique(df['#samples']):
                for x in ['#qubits', '#cycles']:
                    grouped = []
                    if(x == '#qubits'):
                        grouped = ['#cycles']
                    if(x == '#cycles'):
                        grouped = ['#qubits']
                    
                    fixed = {'#iterations': i,
                                '#samples': s}

                    self.plot(df.copy(), x, grouped, fixed)

        for i in pandas.unique(df['#iterations']):
            self.plot(df.copy(), '#samples', ['#qubits', '#cycles'], {'#iterations': i})

        for s in pandas.unique(df['#samples']):
            self.plot(df.copy(), '#iterations', ['#qubits', '#cycles'], {'#samples': s})

        self.plot(df.copy(), '#sampleSteps', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#iterations', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#iterations', ['#qubits', '#sampleSteps'], {})
        self.plot(df.copy(), '#sampleSteps', ['#iterations'], {})
        self.plot(df.copy(), '#iterations', ['#sampleSteps'], {})
        self.plot(df.copy(), '#samples', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#qubits', ['#cycles'], {})
        self.plot(df.copy(), '#cycles', ['#qubits'], {})
        for i in [1,2]:
            for j in [1,2]:
                self.plot(df.copy(), '#nodes', ['#qubits', '#cycles'], {'#tasks': i, '#threads': j})
                self.plot(df.copy(), '#tasks', ['#qubits', '#cycles'], {'#nodes': i, '#threads': j})
                self.plot(df.copy(), '#threads', ['#qubits', '#cycles'], {'#nodes': i, '#tasks': j})

        for q in self.listSystemSizes:
            self.plotDepthEntropy(q)
            for c in self.listCycles:
                for i in range(self.numCircuits):
                    self.plotPDF(q, c, i)

        shutil.make_archive('plots', 'zip', 'plots')

          
    def generateCSV(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        with open(results_file, 'w') as f:
            f.write('#qubits,#cycles,circuit,#nodes,#tasks,#threads,#samples,#iterations,#initialHidden,#sampleSteps,run,#hadamards,tvd,f_xeb,duration,success\n')
        
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
                                            step_size = int(len(self.listSampleSteps)/len(self.listSystemSizes))
                                            index_sample_steps = self.listSystemSizes.index(size) * step_size
                                            for numSampleSteps in self.listSampleSteps[index_sample_steps:index_sample_steps+step_size]:
                                                for run in range(self.numRuns):

                                                    try:
                                                        print('load histogram')
                                                        histogram = self.loadHistogram(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)
                                                        print('load rbm')
                                                        rbm = self.loadRBM(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)
                                                        print('load tvd')
                                                        tvd = '{:f}'.format(self.tvd(self.loadExact(size, cycles, circuits), rbm))
                                                        print('load f_xeb')
                                                        f_xeb = '{:f}'.format(self.f_xeb(self.loadExact(size, cycles, circuits), self.normalise(histogram), int(size)))
                                                        print('load duration')
                                                        duration = '{:f}'.format(self.loadDuration(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run))
                                                        print('load success')
                                                        success = True
                                                    except:
                                                        print("Unexpected error:", sys.exc_info()[0])
                                                        tvd = '-'
                                                        f_xeb = '-'
                                                        duration = '-'
                                                        success = False

                                                    line = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(size,cycles,circuits,nodes,tasks,threads,numSamples,numIterations,numInitialHidden,numSampleSteps,run,hadamards,tvd,f_xeb,duration,success)
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
ev.generateReport()
ev.generatePlots()