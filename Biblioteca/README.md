VLCDA
## Virtual Logic Circuits Design Automation

Python module dedicated to the evolution of logic circuits through genetic algorithms

- Evolving Combinational Logic Circuits through genetic algorithm
- A set of Odd Parity Generetors evolveds until the 100% Fitness 


> To exploit the true potential of emerging nanotechnologies, automation of the
> electronic design (EDA) of digital systems must be rethought. The digital projects
> for the nanoelectronic age it must be immune to defects and accept variability. curiosity-
> Mind you, these are characteristics of biological systems shaped by evolution. algorithms
> genetics can be used for the synthesis of digital systems when hard-
> fully integrated ware, such as those imposed by programmable devices.
> It has even been pointed out that as digital circuits evolve, fault handling strategies
> recognized by evolution reappear naturally.


## Installation

VEDA module requires [Python 3](https://www.python.org/) and some basics modules to run. 

- Random
- Datetime
- Bisect

If you don't have theses modules, run the following commands:

```sh
pip install random
pip install datetime
pip install bisect
```

To install the module VEDA run the following commands...

```sh
pip install VLCDA
```
Be sure to maintain pip updated ...
```sh
pip install --upgrade VLCDA
```

## Usage
The VLCDA module provides classes for evolving logic circuits. Here's an overview of the main classes and their functionality:

### Class `Genoma`

The `Genoma` class represents an individual genome in the genetic algorithm.

#### Method `__init__(self, numberOfGenes, nInputs, nOutputs)`

- `numberOfGenes`: The number of genes in the genome.
- `nInputs`: The number of inputs.
- `nOutputs`: The number of outputs.

#### Method `generate_parent(self)`

Generates a random parent genome.

#### Method `mutate(self)`

Applies mutation to the genome.

#### Method `copyGene(self, childGenome)`

Copies the genes from the current genome to another genome.

#### Method `calculateFitness(self, logicFunction)`

Calculates the fitness of the genome based on a specified logic function.

#### Method `calculateNoiseFitness(self)`

Calculates the noise fitness of the genome.

#### Method `setFaultChance(self)`

Sets the fault chance for the genome.

#### Method `calculateIndexOfActiveGenes(self)`

Calculates the index of active genes in the genome.

#### Method `calculateInputSignals(self, input)`

Calculates the input signals of the genome based on the given input.

#### Method `showInputs(self)`

Displays the input signals of the genome.

#### Method `showOutputs(self)`

Displays the output signals of the genome.

#### Method `calculateFaultRate(self)`

Calculates the fault rate of the genome.

#### Method `display(self)`

Displays the genome's details.

#### Method `__str__(self)`

Overrides the default string representation of the genome.

#### Method `addNoise(self)`

Adds noise to the genome.

#### Method `setNumberOfActiveGenes(self)`



### Class `GeneticAlgorithm`

The `GeneticAlgorithm` class implements a genetic algorithm to evolve genomes towards a desired solution.

#### Method `__init__(self, y=10, maxGeneration=4000000)`

- `y`: The number of children generated per generation.
- `maxGeneration`: The maximum number of generations allowed.

#### Method `display(self, guess, fitness, noiseFitness, totalGeneration)`

Displays information about the current generation, including genotype, fitness, noise fitness, and the total number of generations.

#### Method `addToFitnessList(self, childFitness)`

Adds a child's fitness to the fitness list, ensuring sorting.

#### Method `showBarPlot(self, genome_id, samplingLen)`

Displays a bar plot representing the frequency distribution of fitnesses in the population.

#### Method `getBestGenomeWithSize(self, listChild)`

Returns the genome with the best fitness and shortest genome size in a list of children.

#### Method `getBestGenome(self, listChild)`

Returns the genome with the best fitness in a list of children.

#### Method `evolution(self, genome, logicFunction)`

Executes the process of genetic algorithm evolution. This includes generating children, mutation, selecting the best genomes, and displaying progress information.

This class is used to evolve genomes towards an optimal solution, applying mutations and natural selection over multiple generations until a stopping criterion is reached or the desired solution is found.


## License

MIT

