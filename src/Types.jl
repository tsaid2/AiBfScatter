
module Types

    using RandomExtensions
    import Dates

    export Entity, EntityData, GAmodel, targetParams, GAParams

    abstract type Entity end
    # -------

    mutable struct EntityData
        entity
        generation::Int
        scores :: Array

        EntityData(entity, generation::Int) = new(entity, generation)
        EntityData(entity, model) = new(entity, model.params.currentGeneration)
    end

    mutable struct GAParams
        populationSize ::  Int
        generations ::  Int
        genomeSize ::  Int
        maxGenomeSize :: Int
        crossoverRate :: Float64
        mutationRate :: Float64
        elitism  :: Bool # Keep previous generation's fittest individual in place of worst in current
        historyPath :: IOStream  # Path to save log of fitness history at each generation. Can be used to plot on Excel chart etc.

        totalFitness :: Float64
        targetFitness :: Float64
        targetFitnessCount ::  Int
        currentGeneration ::  Int

        #thisGeneration :: Array
        #nextGeneration :: Array
        #fitnessTable :: Array{Float64,1}
    end

    # -------

    mutable struct GAmodel
        initial_pop_size::Int
        gen_num::Int

        population::Array
        refSet ::Array
        subSets :: Array
        #pop_data::Array{EntityData}
        #freezer::Array{EntityData}

        rng::AbstractRNG

        ga
        specific_fitness

        scores :: Array
        params :: GAParams

        instructionsSet #:: AbstractDict

        GAmodel() = new(0, 1, Any[], Any[], Any[], MersenneTwister(time_ns()), nothing, nothing, [],
                GAParams(136, 1000 , 150, 150, 0.7, 0.01, true, open("output.txt", "w") ,  0.0 , 0, 0.0 , 0 ),
                Dict() )
        GAmodel(params :: GAParams) = new(0, 1, Any[], Any[], Any[], MersenneTwister(time_ns()),
                nothing, nothing, [], params, Dict())
    end

end  # module Types
