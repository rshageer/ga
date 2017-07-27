/* @file
 * Provide simple javascript Genetic Algorithm
 *
 * This file will eventually be split into a single JS file for the main
 * library functionality.
 * I've included some execution code below the library code
 */


/* @type {DNA}
 * Contains a chromosome and all that is necessary to evaluate it
 *
 * @prop {DNA~constructor}
 *   @param {integer} chromoLength
 *     the length of the required chromosome
 *
 *   @param {callback} randomGeneFunc
 *     callback function to return a random gene in the chromosome
 *     
 *   @param {callback} fitnessFunc
 *     callback function to return a fitness of a chromosome
 *     
 * @prop {integer} fitness
 *   gives the fitness score of the chromosomes
 *
 * @prop {Array} genes
 *   array to hold the chromosome, each item is a gene
 *
 * @prop {integer} chromoLength
 *   the length of a chromosome
 * 
 * @prop {variable} randomGene
 *   function to return a random gene 
 *
 * @prop {function} customFitnessFunc
 *   method to evaluate fitness of the chromosome
 *    
 * @prop {DNA~mutate} mutate
 *   @param {integer} mutationRate
 *     mutation rate of each gene per chromosome
 *
 *   @param {callback} randomGeneFunc
 *     function to return a random gene. Default is DNA~randomGene
 *     
 * @prop {DNA~generateRandomChromosome} generateRandomChromosome
 *   generate a random chromosome to fill the DNA~genes array
 *
 * @prop {DNA~calculateFitness} calculateFitness
 *   uses the custom fitness function provided in order to calculate fitness and
 *   store it in the DNA~fitness property
 */
class DNA {
    constructor(chromoLength, randomGeneFunc, fitnessFunction) {
        this.fitness = 0.0;
        this.genes = [];
        this.chromoLength = chromoLength;
        this.randomGene = randomGeneFunc;
        this.customFitnessFunc = fitnessFunction;
    }

    mutate(mutationRate, randomGeneFunc = this.randomGene) {
		// loop through each gene and randomly replace it based on the 
		// mutation rate.
        
         if(Math.random() < mutationRate){
            let randomGene = Math.floor(Math.random() * this.genes.length);
            this.genes[randomGene] = randomGeneFunc();
        }
  
    }

    generateRandomChromosome() {
		// generate a random gene and fill the DNA~genes array. 
        for (let i = 0; i < this.chromoLength; i++) {
            this.genes.push(this.randomGene());
        }
    }

    calculateFitness() {
        this.fitness = this.customFitnessFunc(this.genes);
    }
}

/* @type {Population}
 *   contains a population of chromosomes that are all being evaluated in a generation.
 *
 * @prop {int} size
 *   Amount of DNA objects to hold in the Population
 *
 * @prop {Array} items
 *   Array of DNA objects in the population
 *
 * @prop {Population~generateRandom} generateRandom
 *   method to generate a random population and populate the Population~items Array
 *
 * @prop {Population~getFittest} getFittest
 *   return the DNA item with the highest fitness score
 *
 * @prop {Population~getTotalFitness} getTotalFitness
 *   @return total fitness value of the population.
 *
 * @prop {Population~getAverageFitness} getAverageFitness
 *   @return average fitness value of the population.
 */
class Population {
	
    constructor() {
        this.size = 30;
        this.items = [];
    }
	
    generateRandom(chromoLength, randomGeneFunc, fitnessFunc) {
		// fill the Population~items array with new DNA objects 
		// Each object should have a random chromosome and it's fitness calculated
        for (let i = 0; i < this.size; i++) {
            let dna = new DNA(chromoLength, randomGeneFunc, fitnessFunc);
            dna.generateRandomChromosome();
            dna.calculateFitness();
            this.items.push(dna);
        }
    }

    getFittest() {
		//Choose the item with the highest fitness in the population
        this.items.sort((a, b) => {
            return (b.fitness - a.fitness);
        });
        return this.items[0];
    }
    
    getTotalFitness() {
        let totalFitness = 0;
        this.items.forEach((v, i) => {
            totalFitness += v.fitness;
        });

        return totalFitness;
    }
    
    getAverageFitness() {
        return Math.floor(this.getTotalFitness() / this.items.length);
    }
};

/* @type {Pool} 
 * The pool object is where all the computation happens, it holds a Population object
 * which in turn holds a host of DNA objects. It then manages generates a new generation
 * based on
 *
 * @prop {Population} currentGeneration
 *   the current population being evaluated by the pool.
 *
 * @prop {float} mutationRate
 *   the rate at which each gene in a chromosome is mutated.
 *
 * @prop {integer} generation
 *   the count of the generation being evaluated currently
 *
 * @prop {Pool~generatePopulation} generatePopulation
 *   generates a new population. if the population is undefined, it generates
 *   it randomly. and if it's not undefined it generates based on the current
 *   generation.
 *
 * @prop {Pool~crossover} crossover
 *   takes two random parents from the current population and mates them to 
 *   produce offsprings that would eventually be fed into the general population.
 *
 *   @param {integer} type = 1
 *     integer indicating the type of crossover function to perform. 
 *
 *   @return {DNA[]}
 *     returns an array of newly generated child DNA Objects.
 *
 * @prop {Pool~selectParents} selectParents
 *   randomly selects parents from the current population based on their probability
 *   of selection. 
 *
 *   @param {int} amt = 2
 *     the amount of parents to select
 *
 *   @return {DNA[]}
 *     returns an array of DNA objects that were selected from the current population
 *
 * @prop {Pool~selection_roulette} selection_roulette
 *   uses roulette selection to get a parent (DNA Object)
 * 
 *   @param {DNA[]} pool
 *     an Array of DNA items to be selected from
 *
 *   @param {integer} totalFitness
 *     total of all the fitness scores of the pool array
 *
 *   @return {DNA}
 *
 * @prop {Pool~selection_bag} selection_bag
 *   uses bag selection to get a parent (DNA Object). items are treated as though they
 *   are all put in a bag with the probability of selection based on fitness
 * 
 *   @param {DNA[]} pool
 *     an Array of DNA items to be selected from
 *
 *   @return {DNA}
 *
 * @prop {Pool~selection_fittestWins} selection_fittestWins
 *   always pairs the fittest item with another.
 * 
 *   @param {DNA[]} pool
 *     an Array of DNA items to be selected from
 *
 *   @param {integer} totalFitness
 *     total of all the fitness scores of the pool array
 *
 *   @param {integer} index
 *     zero index of which parent is being chosen currently
 *
 *   @return {DNA}
 * */
class Pool {
    constructor(mutationRate = .01, crossoverProbability = 0.8) {
        this.currentGeneration = new Population();
        this.mutationRate = mutationRate;
        this.crossoverProbability = crossoverProbability;
        this.generation = 0; // which generation are we currently in
    }

    generatePopulation() {
        //Generate a new population if it does not exist else create a new one 
        if (typeof (this.currentGeneration) == 'undefined') {
			this.currentGeneration = new Population();
            this.currentGeneration.generateRandom();

        } else {
            let newGeneration = [];
			
            // uses the crossover to populate newly born DNA objects
            for (let i = 0; i < this.currentGeneration.size; i) {
                let parents = [];
                
                //select parents and crossover as probable
                if(Math.random() < this.crossoverProbability) {
                    parents = this.crossover().slice(0);
                } else {
                    parents = this.selectParents().slice(0);
                }
                
                parents.forEach((v, i) => {
                    newGeneration.push(v);
                });
                i = newGeneration.length;
            }

            //mutate the new generation and calculate its new fitness
            for (let i = 0; i < newGeneration.length; i++) {
                newGeneration[i].mutate(this.mutationRate);
                newGeneration[i].calculateFitness();
            }
			
			// clear the current generation and copy the newly generated DNA's to it's population.
            this.currentGeneration.items = [];
            this.currentGeneration.items = newGeneration.slice(0);
        }
        this.generation++ // keep count of the generation.
    };

    // TODO: include other crossover types. Even consider randomising the process.
    crossover(type = 1) {
        //Choose differnt types of Crossover (To be completed)  
        switch (type) {
            case 1:
                //simple single point crossover
                let parents = this.selectParents().slice(0);
                let crossoverPoint = Math.floor(Math.random() * parents[0].genes.length);
                let parAFront = parents[0].genes.slice(0, crossoverPoint);
                let parBFront = parents[1].genes.slice(0, crossoverPoint);
                let parABack = parents[0].genes.slice(crossoverPoint, parents[0].genes.length);
                let parBBack = parents[1].genes.slice(crossoverPoint, parents[1].genes.length);
                parents[0].genes = [];
                parents[1].genes = [];
                parents[0].genes = parAFront.concat(parBBack.slice(0)).slice(0);
                parents[1].genes = parBFront.concat(parABack.slice(0)).slice(0);
                return parents.slice(0);
                break;

                //two point crossover
            case 2:

                break;


                //Uniform and Half/uniform crossover
            case 3:

                break;
        }

    };

    selectParents(amt = 2, selectionType = 2) {
        let parents = [];
        let pool = this.currentGeneration.items.slice(0); //all the DNA objects in the Population 
        
        let totalFitness = this.currentGeneration.getTotalFitness();
        for (let i = 0; i < amt; i++) {

			//select a parent randomly if no items are fit
            if(totalFitness == 0) {
                let index = Math.floor(Math.random() * pool.length);
                parents.push(pool[index]);
            } else {
                selectionType = selectionType > 0? selectionType : Math.floor(Math.random() * (3 - 1)) + 1;
                switch(selectionType){
                    case 1: //roulette selection
                        parents.push(this.selection_roulette(pool.slice(0), totalFitness));
                        break;
                        
                    case 2: // fittest chromosome is always selected and paired based on roulette
                        parents.push(this.selection_fittestWins(pool.slice(0), totalFitness, i));
                        break;
                        
                    case 3: //bag selection
                        parents.push(this.selection_bag(pool.slice(0)));
                        break;
                }
                //
                
                //
            }
 
        }

        return parents.slice(0);
    }
    
    selection_roulette(pool, totalFitness) {
        let fitItems = [];
        
        //remove unfit items from pool
        pool.forEach((v, i) => {
            if(v.fitness > 0)
                fitItems.push(v);
        });
        
        let r = Math.floor(Math.random() * totalFitness); //random pointer to the DNA item
        let accFitness = 0;

        for (let i = 0; i < fitItems.length; i++) {
            accFitness += fitItems[i].fitness;
            if (accFitness >= r) {
                return fitItems[i];
            }
        }
        
    }

    selection_bag(pool) {
        let indices = [];
        
        pool.forEach((v, i, a) => {
            let s = Math.ceil( v.fitness * 100 );
            for(let c = 0; c < s; c++){ indices.push(i);}
        });
        
        let grabber = Math.floor(Math.random() * indices.length);
        
        return pool[indices[grabber]];
    }

    selection_fittestWins(pool, totalFitness, index) {
        
        if(index == 0) {
            pool.sort((a, b) => {return ( b - a );});
            return pool[0];
        } else {
            return this.selection_roulette(pool, totalFitness);
        }
    }
}
