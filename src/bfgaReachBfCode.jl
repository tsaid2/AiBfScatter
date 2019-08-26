
module bfgaReachBfCode

    import Base.isless
    #using ResumableFunctions

    include("Types.jl")
    using .Types

    #using GeneticAlgorithms
    #include("GeneticAlgorithms.jl")
    #using .GeneticAlgorithms

    include("bfga.jl")
    using .bfga

    using Distributed



    function fitness(input, goal)
        input.fitness = goal.m_length - sum((input.dna-goal.dna).^2)
        input.fitness
    end


    function simulate_entity(input, goal)
        #bft = bfType(ent.dna)
        #target_goal = length(goal)
        _res = "goal code : $(goal.program )  \n"
        _res *= "goal code : $(input.program ) "
        _res
    end

    #=function getTargetFitness()
        256* length(join(words, ""))
    end=#


    #=function getParams()
        #=
        populationSize ::  Int
        generations ::  Int
        genomeSize ::  Int
        maxGenomeSize :: Int
        crossoverRate :: Float64
        mutationRate :: Float64
        elitism  :: Bool # Keep previous generation's fittest individual in place of worst in current
        historyPath :: String  # Path to save log of fitness history at each generation. Can be used to plot on Excel chart etc.

        totalFitness :: Float64
        targetFitness :: Float64
        targetFitnessCount ::  Int
        currentGeneration ::  Int

        #thisGeneration :: Array
        #nextGeneration :: Array
        #fitnessTable :: Array{Float64,1}
        =#

        logfile = open("../Results/logRepeat.log", "w")

        tgFitness =  getTargetFitness()
        println("targetFitness = $tgFitness ")
        write(logfile, "targetFitness = $tgFitness \n")
        return Main.GeneticAlgorithms.Types.GAParams(60, 10000 , 36, 150, 0.7, 0.01, true, logfile ,  0.0 , tgFitness, 0.0 , 0 )
    end=#

    function getBfCode(ent)
        ent.program
    end


end
