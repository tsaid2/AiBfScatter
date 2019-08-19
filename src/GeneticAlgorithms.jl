
module GeneticAlgorithms

    # -------
    using RandomExtensions
    using Distributed

    include("Types.jl")
    using .Types

    import Dates

    import Base, Base.show, Base.isless

    export  runga,
            isless,
            freeze,
            defrost,
            generation_num,
            population

    # -------

    #isless(lhs::Entity, rhs::Entity) = lhs.fitness < rhs.fitness

    fitness!(ent, fitness_score) = ent.fitness = fitness_score

    distance(ind1 , ind2 ) = abs(ind1.fitness - ind2.fitness)

    global _g_model

    # -------

    function freeze(model::GAmodel, entity::EntityData)
        push!(model.freezer, entity)
        println("Freezing: ", entity)
    end

    function freeze(model::GAmodel, entity)
        entitydata = EntityData(entity, model.params.currentGeneration)
        freeze(model, entitydata)
    end

    freeze(entity) = freeze(_g_model, entity)


    function defrost(model::GAmodel, generation::Int)
        filter(model.freezer) do entitydata
            entitydata.generation == generation
        end
    end

    defrost(generation::Int) = defrost(_g_model, generation)


    generation_num(model::GAmodel = _g_model) = model.params.currentGeneration


    population(model::GAmodel = _g_model) = model.population

    function show_simulation(model :: GAmodel, ent)
        printed = model.specific_fitness.simulate_entity(ent, model.instructionsSet)
        @show ent.fitness
        #@show model.specific_fitness.fitness(ent)
        println(printed)
        return "fitness : $(ent.fitness) \n $printed"
    end

    function runga(mdl::Module, fit_mdl :: Module )
        # get the parameters of the specific fitness
        model = GAmodel(fit_mdl.getParams())
        model.specific_fitness = fit_mdl

        # save the ga (mutation and crossover units function)
        model.ga = mdl
        # save once the dictionary of instructions Set
        model.instructionsSet = mdl.getInstructionsSet()
        #run the ga
        runga(model; resume = false)
    end

    function addNew(refSet, population)
        l = []
        k = length(population[1].dna)
        for e in population
            mini = k
            for s in refSet
                mini = min(mini, sum(e.dna-s.dna).^2)
            end
            push!(l, mini)
        end
        index = argmax(l)
    end

    function findFirstNotIn(ens, i)
        if i in ens
            return findFirstNotIn(ens, i+1)
        else
            return i
        end
    end

    function subSetsGeneration(n)
        subSets = []
        currentTypeSet = []
        oldTypeSet = []
        #Type2
        for i in 1:(n-1)
            map( a -> (push!(subSets, [a,i]); push!(currentTypeSet, [i, a]) ), i+1:n )
        end
        # Type3
        oldTypeSet = currentTypeSet
        oldTypeSet = map( a -> (i = findFirstNotIn(a, 1); push!(subSets, [a;i]); [a;i]), oldTypeSet )
        #Type 4
        oldTypeSet = map( a -> (i = findFirstNotIn(a, 1); push!(subSets, [a;i]); [a;i]), oldTypeSet )
        # Type 5
        map(a -> push!(subSets, 1:a), 5:n)
        subSets
    end

    function searchInCurrentRefSet(model)
        newSolutions = true
        j = 0
        while newSolutions
            # 3. Generate NewSubsets with the subset generation method. Make NewSolutions = FALSE.
            lrefSet = length(model.refSet)
            model.subSets = subSetsGeneration(lrefSet)
            newSolutions = false
            while !isempty(model.subSets)
                # 4. Select the next subset s in NewSubsets
                s = pop!(model.subSets)
                # 5. Apply the solution combination method to s to obtain one or more new trial solutions x.
                # Apply the improvement method to the trial solutions.
                ls = length(s)
                pool = []
                for i in 1:ls-1
                    map( a -> ((child1, child2) = model.ga.crossover(Any[model.refSet[a],model.refSet[i]]);
                            push!(pool, child1); push!(pool, child2) ) , i+1:ls )
                end

                pmap(
                    ent -> (model.ga.entityToBfInstructions!(ent); model.specific_fitness.fitness(ent, model.instructionsSet))
                    , pool)
                sort!(pool, rev=true)
                # 6. Apply the reference set update method.
                newEntAdded = []
                refSetHasChanged = false
                #@show length(model.refSet)
                for e in pool
                    if !isempty(model.refSet) && e.fitness > (model.refSet[1].fitness)
                        push!(newEntAdded, e)
                        deleteat!(model.refSet, 1)
                        refSetHasChanged = true
                    end
                end
                append!(model.refSet, newEntAdded)
                evaluate_refSet(model)

                lastIdx = length(model.refSet)
                _fitness = model.refSet[lastIdx].fitness
                # displays
                j += 1
                if j % 100 == 0
                    #@show refSetHasChanged
                    #_log = ""
                    _log = "    $(Dates.now()) , "
                    _log *= "Gen : $(model.params.currentGeneration) , "
                    _log *= "BEST: $_fitness , \n"
                    if j%500 == 0
                        _log *= show_simulation(model, model.refSet[lastIdx]) * "\n"
                    end
                    print(_log)
                    if model.params.historyPath != nothing
                        write(model.params.historyPath, _log)
                    end
                end

                if (model.params.targetFitness > 0 && _fitness >= model.params.targetFitness)
                    model.params.targetFitnessCount = model.params.targetFitnessCount +1
                    if (model.params.targetFitnessCount > 500) # default : set to 1000
                        # Stop the algo
                        newSolutions = false
                        break;
                    end
                else
                    # in case we expand and lose the best one
                    model.params.targetFitnessCount = 0;
                end

                if refSetHasChanged
                    # 7. Make NewSolutions = TRUE.
                    newSolutions = true
                end
                # 8. Delete s from NewSubsets.
                # --> already done with pop!
            end
        end
    end




    function runga(model::GAmodel; resume = false)
        stop = false

        #set _expandAmount & _expandRate, TODO I think they are not placed very well in case we want specifi paramters for each fitness
        _expandAmount = 0
        _expandRate = 5000

        if (!resume)
            #  initialize params.
            model.params.totalFitness = 0;
            model.params.targetFitness = model.specific_fitness.getTargetFitness();
            model.params.targetFitnessCount = 0;
            model.params.currentGeneration = 0;
            #stop = false;

            #reset_model(model)
            #create_initial_population(model)
            #evaluate_population(model)
        end

        if model.params.historyPath != nothing
            #to save params in the history file
            write(model.params.historyPath, "params : " * string(model.params) * " \n")
        end

        # beginning of the scatter implementation
        #=
        1. Start with P = Ø. Use the diversification generation method to construct a solution
        and apply the improvement method. Let x be the resulting solution.
        If x ∉ P then add x to P (i.e., P = P ∪ x ), otherwise, discard x.
        Repeat this step until |P| = PSize.
        =#
        lastIdx = model.params.populationSize #length(model.population)
        reset_model(model)
        #=pSize = 0
        popu = []
        while pSize < lastIdx
            create_initial_population(model)
            evaluate_population(model)
            sol = model.ga.create_entity(model.population[lastIdx].dna)
            #model.ga.mutate(sol) # Apply the improvement method TODO
            if !model.ga.isPresent(sol, popu)
                push!(popu, sol)
                pSize += 1
            end
            reset_model(model)
        end
        model.population = popu=#
        create_initial_population(model)
        evaluate_population(model)
        for i in 1:7
            push!(model.refSet, pop!(model.population))
        end
        #=
        2. Use the reference set update method to build RefSet = { x 1 , ..., x b } with the “best” b solutions in P.
         Order the solutions in RefSet according to their objective function value such that x 1 is
         the best solution and x b the worst. Make NewSolutions = TRUE.
        =#

        found = false
        tour = 0
        while !found
            tour += 1
            _log = "** Tour n. $tour at $(Dates.now()) \n"
            print(_log)
            if model.params.historyPath != nothing
                write(model.params.historyPath, _log)
            end

            for i in 1:7
                index = addNew(model.refSet, model.population)
                push!(model.refSet, model.population[index])
                deleteat!(model.population, index)
            end

            evaluate_refSet(model)
            #@show length(model.refSet)
            searchInCurrentRefSet(model)
            best = model.refSet[length(model.refSet)]
            found = best.fitness >= model.params.targetFitness
            if !found
                #append!(model.population, model.refSet)
                create_initial_population(model)
                evaluate_population(model)
                model.refSet = model.refSet[8:14]
            end
        end


        print("\n **********DONE*******")
        if model.params.historyPath != nothing
            best = model.refSet[length(model.refSet)] #[length(model.population)]
            _log = show_simulation(model, best)
            #_log *= "\n Generation : $(model.params.currentGeneration) \n"
            _log *= "fitness $(best.fitness) \n"
            _log *= "dna : $(best.dna) \n"
            write(model.params.historyPath, _log)
            close(model.params.historyPath)
        end
        model

    end

    # -------

    function reset_model(model::GAmodel)
        global _g_model = model

        model.params.currentGeneration = 1
        empty!(model.population)
    end



    function create_initial_population(model::GAmodel)
        for i = 1:model.params.populationSize
            entity = model.ga.create_entity(i, model.params.genomeSize)
            push!(model.population, entity)
        end
    end




    function evaluate_population(model::GAmodel)
        #pmap(model.ga.entityToBfInstructions!, model.population)
        #scores = [ model.specific_fitness.fitness(ent, model.instructionsSet) for ent in model.population ]

        # for each entity, translate dna to bf instructions, then perform the specific fitness,
        # it saves its fitness & bonus in itself and return the fitness, see the fitness function
        scores = pmap(
            ent -> (model.ga.entityToBfInstructions!(ent); model.specific_fitness.fitness(ent, model.instructionsSet))
            , model.population)

        model.params.totalFitness = sum(scores)

        #pmap(fitness!, model.population, scores) # TODO decomment for normal use

        # the isless for the entities is defined in bfga
        sort!(model.population; rev = false)
        model.scores = sort!(scores; rev = false)
    end

    function evaluate_refSet(model::GAmodel)
        #pmap(model.ga.entityToBfInstructions!, model.population)
        #scores = [ model.specific_fitness.fitness(ent, model.instructionsSet) for ent in model.population ]

        # for each entity, translate dna to bf instructions, then perform the specific fitness,
        # it saves its fitness & bonus in itself and return the fitness, see the fitness function
        scores = pmap(
            ent -> (model.ga.entityToBfInstructions!(ent); model.specific_fitness.fitness(ent, model.instructionsSet))
            , model.refSet)

        #model.params.totalFitness = sum(scores)

        #pmap(fitness!, model.population, scores) # TODO decomment for normal use

        # the isless for the entities is defined in bfga
        sort!(model.refSet; rev = false)
        #model.scores = sort!(scores; rev = false)
    end




    function crossover_population(model::GAmodel, groupings)
        thisGeneration = (population(model))
        _length = 1

        model.population = Any[]

        model.params.currentGeneration += 1
        g = model.ga.create_entity(model.params.currentGeneration, model.params.genomeSize)
        g2 = model.ga.create_entity(model.params.currentGeneration, model.params.genomeSize)
        #println("Generation n' $(model.params.currentGeneration) ")
        if model.params.elitism
            l = model.params.populationSize #length(thisGeneration)
            g = thisGeneration[l]
            l_1 = (model.params.populationSize) -1
            g2 = thisGeneration[l-1]
            #push!(model.population, g)
            #push!(model.pop_data, EntityData(g, model.params.currentGeneration))
            #push!(model.population, g2)
            #push!(model.pop_data, EntityData(g2, model.params.currentGeneration))

            _length += 2;
        end
        #println("_______________________________________________________________________________")
        #@show model.scores

        for i in _length:2:model.params.populationSize
            pidx1 = rouletteSelection(model)
            pidx2 = rouletteSelection(model)
            #println("$pidx1    $pidx2") # for debugging
            child1, child2 = nothing, nothing
            parent1 = thisGeneration[pidx1]
            parent2 = thisGeneration[pidx2]

            if (rand() < model.params.crossoverRate)
                child1, child2 = model.ga.crossover(Any[parent1, parent2])
            else
                child1, child2 = parent1, parent2
                parent1.bonus = 0
                parent2.bonus = 0
            end

            push!(model.population, child1)
            push!(model.population, child2)
        end

        if model.params.elitism
            push!(model.population, g)
            #push!(model.pop_data, EntityData(g, model.params.currentGeneration))
            push!(model.population, g2)
            #push!(model.pop_data, EntityData(g2, model.params.currentGeneration))
        end

        model.params.populationSize = length(model.population)
        #pmap(model.ga.clearCode!, model.population)

        # Expand genomes.
        if model.population[1].m_length != model.params.genomeSize
            #println("expannndd")
            newGenomeSize = model.params.genomeSize
            for m in model.population
                if m.m_length != newGenomeSize
                    model.ga.expand(m, newGenomeSize)
                end
            end
        end
    end

    function mutate_population(model::GAmodel) # deprecated, the mutation is performed in the crossover
        pmap(model.ga.mutate, model.population)
        #=for entity in model.population
            model.ga.mutate(entity)
        end=#
    end

    # We can add strangers in the population
    function add_stranger(model :: GAmodel, num)
        #println("-_-")
        for i in 1:num
            entity = model.ga.create_entity(model.params.currentGeneration, model.population[1].m_length)
            push!(model.population, entity)
            #@async push!(model.pop_data, EntityData(entity, model.params.currentGeneration))
        end
    end

    function rouletteSelection(model :: GAmodel)
        #idx = trunc(Int, rand()*length(model.scores))
        #_idx = idx ==0 ? 1 : idx
        #rand(1:trunc(Int, _idx))

        n = model.params.populationSize #length(model.scores)
        randomFitness = rand() * (model.scores[n] == 0 ? 1 : model.scores[n])

        idx = -1
        first = 1
        last = model.params.populationSize
        mid = trunc(Int , (last - first)/2)

        #  ArrayList's BinarySearch is for exact values only
        #  so do this by hand.
        while (idx == -1 && first <= last)
            if (randomFitness < model.scores[mid])
                last = mid;
            elseif (randomFitness > model.scores[mid])
                first = mid;
            end
            mid = trunc(Int, (first + last)/2)
            #  lies between i and i+1
            if ((last - first) == 1)
                idx = last;
            end
        end

        return idx
    end
end
