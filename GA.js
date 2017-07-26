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
 *   the lenght of a chromosome
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
        this.genes = this.genes.map(function (v, i) {
            if (Math.random() <= mutationRate) {
                return randomGeneFunc();
            } else {
                return v;
			}
        });

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
        let fittestItem = this.items[0];
		
        for (let i = 1; i < this.items.length; i++) {
            if (this.items[i].fitness > fittestItem.fitness) {
                fittestItem = this.items[i];
            }
        }

        return fittestItem;
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
 */
class Pool {
    constructor() {
        this.currentGeneration = new Population(); 
        this.mutationRate = 0.005;
        this.crossoverProbability = 0.8;
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

    selectParents(amt = 2) {
        let parents = [];
        let pool = this.currentGeneration.items.slice(0); //all the DNA objects in the Population 
        
        for (let i = 0; i < amt; i++) {
			//roulette wheel selection.
			
			//get the total fitness of the population
            let totalFitness = 0;
            pool.forEach((val, index) => {
                return totalFitness += val.fitness;
            });

			//select a parent
            if(totalFitness == 0) {
                let index = Math.floor(Math.random() * pool.length);
                parents.push(pool[index]);
            } else {
                let r = Math.floor(Math.random() * totalFitness); //random pointer to the DNA item
                let accFitness = 0;

                for (let c = 0; c < pool.length; c++) {
                    accFitness += pool[c].fitness;
                    if (accFitness >= r) {
                        parents.push(pool[c]);
                        break;
                    }
                }
            }
 
        }

        return parents.slice(0);
    }
};

//___________________________________________________________________________________//

// end of library code: 
//Implementation of above lib
//run the program
(function () {

    //random goal of this genetic algorithm
    var goal = "RiazSHage";
    let p = new Pool()
    p.currentGeneration.size = 1000;
    p.mutationRate = 0.001;
    p.crossoverProbability = .8;
    let timeout = 300;
    let b = 0;

    //custom fitness func
    var calculateFitness = function (genes) {
        let strDNA = [];
        let score = 0;
        //convert ascii to string
        for (let i = 0; i < genes.length; i++) {
            strDNA[i] = String.fromCharCode(genes[i]);
        }

        let strGoal = goal.split("");

        for (let i = 0; i < strGoal.length; i++) {
            if (strDNA[i] == strGoal[i]) {
                score++;
            }
        }

        return score;
    };

    //custom randome gene func
    var randomGene = function() {
        return Math.floor(Math.random() * 127);
    };

    //func to translate chromosome into a string representation 
    var translate = function(dna) {
        strArr = [];
        for (let i = 0; i < dna.genes.length; i++) {
            strArr.push(String.fromCharCode(dna.genes[i]));
        }
        let str = strArr.join("");
        return str;
    };


    p.currentGeneration.generateRandom(goal.length, randomGene, calculateFitness);

    while (p.currentGeneration.getFittest().fitness <= goal.length) {
        let d = p.currentGeneration.getFittest();
        console.log("Generation: " + p.generation + "\t\tFitness: " + Math.floor((d.fitness / goal.length) * 100) + "%\t\tResult: " + translate(d));
        if (d.fitness == goal.length) break;

        p.generatePopulation();

        b++;
        if (b >= timeout) break;
    }

}());