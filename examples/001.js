(function () {

    //random goal of this genetic algorithm
    let goal                    = "Riaz Shageer";
    let mutationRate            = 0.015;
    let crossoverProbability    = 0.65;
    let p                       = new Pool(mutationRate, crossoverProbability);
    p.currentGeneration.size    = 750;
    let timeout                 = 2000;
    let b                       = 0;
    let completeCriteria        = goal.length;

    //custom fitness func
    let calculateFitness = function (genes) {
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
        //return Math.pow(score, 2);
    };

    //custom randome gene func
    let randomGene = function() {
        return Math.floor(Math.random() * (127 - 32)) + 32;
    };

    //func to translate chromosome into a string representation 
    let translate = function(dna) {
        strArr = [];
        for (let i = 0; i < dna.genes.length; i++) {
            strArr.push(String.fromCharCode(dna.genes[i]));
        }
        let str = strArr.join("");
        return str;
    };


    p.currentGeneration.generateRandom(goal.length, randomGene, calculateFitness);
    let labels      = []
    let totalValues = []
    let avgValues   = []
    
    while (p.currentGeneration.getFittest().fitness <= completeCriteria) {
        let d = p.currentGeneration.getFittest();
       
        console.log("Generation: " + p.generation
            + "\t\tFitness: " + Math.floor((d.fitness / completeCriteria) * 100) 
            + "%\t\tResult: " + translate(d)
            + "\tTotal Fitness: " + p.currentGeneration.getTotalFitness()
            + "\tAvg Fitness: " + p.currentGeneration.getAverageFitness()
            );
        
/*         labels.push(p.generation);
        totalValues.push(p.currentGeneration.getTotalFitness());
        avgValues.push(p.currentGeneration.getAverageFitness()); */
        
        if (d.fitness == completeCriteria) break;
        p.generatePopulation();

        b++;
        if (b >= timeout) break;
    }

    /*
    //draw chart
    let ctx = document.getElementById("chart").getContext("2d");
    
    let chart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [
            {
                label:"Total Fitness",
                backgroundColor: 'rgb(255, 99, 132)',
                borderColor: 'rgb(255, 99, 132)',
                data: totalValues,
            }
            ]
        },    
    });
    */
}());
