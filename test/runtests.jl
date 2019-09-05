using Test
using Random
#using LinearAlgebra
#Random.seed!(9874984737484)

for tests in [
			#"runtestString.jl",		#
            #"runtestAdd.jl", 			#
			#"runtestLogicalOr.jl",		#
			#"runtestLogicalXor.jl", 	#
			#"runtestReverseString.jl",	#
			#"runtestLengthString.jl",	#
			#"runtestTimesTwo.jl",		#
			"runtestTimesThree.jl"
			#"runtestLogicalAnd.jl",	#
			#"runtestRepeat.jl",		#
			#"runtestWarningCountdown.jl" # BROKEN ! Doesn't work, NENI
			#"runtestCountdown.jl"
			#"runtestExtractInQuote.jl",
			#"runtestExtractInQuoteInside.jl" # Never performed
			#"runtestFibonacci.jl"
			#"runtestTrimLefToFQuote.jl"
			#"runtestSubtract.jl"
			]
    include(tests)
end
