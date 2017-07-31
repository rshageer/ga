

/*
let data = {
  labels: ["Mon", "Tue", "Wed", "Thur", "Fri"],
  series: [[1, 3, 5, 3, 7]],
};

let options = {
  width:900,
  height:300,
};

new Chartist.Bar('.ct-chart', data, options);
*/
(function(){

    let goal                    = "Riaz Shageer";
    let mutationRate            = 0.01;
    let crossoverProbability    = 0.7;
    let p                       = new Pool(mutationRate, crossoverProbability);
    p.currentGeneration.size    = 1000;
    let timeout                 = 1000;
    let b                       = 0;
    let completeCriteria        = goal.length;

    
    let randomGene = function(gene){
      let allChars = "abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ";
      let randomPointer = Math.floor(Math.random() * allChars.length);
      return allChars[randomPointer];
    };
    
    let getFitness = function(genes) {
      let score = 0;

      for(let i = 0; i <= genes.length; i++){
        if(genes[i] == goal[i]) score++;
      }

      return score;
    }
    
    let translate = function(dna){
       return dna.genes.join("");
    };

    p.currentGeneration.generateRandom(goal.length, randomGene, getFitness);
    while (p.currentGeneration.getFittest().fitness <= completeCriteria) {
        let d = p.currentGeneration.getFittest();
       
        console.log("Generation: " + p.generation
            + "\t\tFitness: " + Math.floor((d.fitness / completeCriteria) * 100) 
            + "%\t\tResult: " + translate(d)
            + "\tTotal Fitness: " + p.currentGeneration.getTotalFitness()
            + "\tAvg Fitness: " + p.currentGeneration.getAverageFitness()
            );
        
        
        if (d.fitness == completeCriteria) break;
        p.generatePopulation();

        b++;
        if (b >= timeout) break;
    }

})();
