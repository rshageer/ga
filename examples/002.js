

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

    let goal = "Riaz Shageer was here";
    let mutationRate            = 0.001;
    let crossoverProbability    = 0.6;
    let p                       = new Pool(mutationRate, crossoverProbability);
    p.currentGeneration.size    = 500;
    let timeout                 = 30;
    let b                       = 0;
    let completeCriteria        = goal.length;

    
    let randomGene = function(gene){
        //randomize
        if (typeof(gene) == "undefined"){ 
            return (Math.floor(Math.random()*2-1) + 1).toString();
        }
        //flip bit if it's already defined
        let newGene = (gene == "1")? "0":"1";
        return newGene;
    };
    
    let getFitness = function(genes) {
        let score = 0;
        let ch = [];
        let goalArr = goal.split("");
        for(let i = 0; i < genes.length; i++){
            ch.push(genes[i]);

            if (ch.length == 7) {
                if (String.fromCharCode(parseInt(ch.join(""), 2)) == (goalArr[((i + 1)/7)-1])) score++;
                ch = [];
                    
            }
        }
        return score;
    }
    
    let translate = function(dna){
        let result = "";
        let ch = [];
        for(let i = 0; i < dna.genes.length; i++){
            ch.push(dna.genes[i]);
            
            if (ch.length == 7) {
                result += String.fromCharCode(parseInt(ch.join(""), 2))
                ch = [];
            } else {
            }
        }
        
        return result;
        
    }
    p.currentGeneration.generateRandom(goal.length * 7, randomGene, getFitness);
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
