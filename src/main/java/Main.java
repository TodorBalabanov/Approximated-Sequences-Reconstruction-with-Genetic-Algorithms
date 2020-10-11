import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Application single entry point class.
 * 
 * @author Todor Balabanov
 */
public class Main {
	/**
	 * Pseudo-random number generator.
	 */
	private static final Random PRNG = new Random();

	/**
	 * Chunk size is related with the size of the visible part of the sequence.
	 */
	private static final int CHUNKS_SIZE = 3;

	/**
	 * Genetic algorithm population size.
	 */
	private static final int POPULATION_SIZE = 137;

	/**
	 * The mutation rate shows the probability for each gene to mutate. Because
	 * it is a probability it has values from 0 to 1.
	 */
	private static final double MUTATION_RATE = 0.005;

	/**
	 * Number of generations to be created as depth of the recursion.
	 */
	private static final int RECURSION_DEPTH = 7;
	
	/**
	 * The histogram threshold is the minimum number of the appearance of the
	 * less probable unique chunk. When the size of the sequence is unknown by
	 * estimating the less probable unique chunk appearance the amount of chunks
	 * sample can be estimated.This estimation should be done with confidence
	 * level of about 95% and according to the rules of the normal probability
	 * distribution.
	 */
	private static final int HISTOGRAM_THRESHOLD = 100;

	/**
	 * How many genetic algorithm generations to be evolved.
	 */
	private static final long EVOLUTION_EPOCHS = 10000;

	/**
	 * Original sequences which should be reconstructed.
	 */
	private static final int[][][] ORIGINAL_STRIPS = {{
			{9, 7, 12, 9, 5, 10, 12, 3, 11, 10, 4, 12, 10, 6, 12, 11, 6, 10, 12,
					8, 9, 5, 12, 8},
			{3, 12, 11, 7, 10, 4, 11, 12, 7, 11, 12, 5, 11, 9, 6, 10, 1, 1, 1,
					3, 9, 4, 11, 9, 5, 8, 11, 7, 8, 6, 11, 8},
			{3, 11, 10, 7, 11, 8, 3, 12, 10, 7, 11, 10, 4, 12, 10, 4, 9, 6, 12,
					10, 6, 11, 12, 6, 8, 10, 5, 11, 12, 7, 8, 5, 10, 12, 7, 11,
					9, 3, 9, 8, 7, 12, 4, 9, 5, 12, 8, 6, 12, 11, 7, 1, 1, 1,
					12, 7, 11, 8, 4, 11, 10},
			{8, 9, 5, 11, 12, 3, 11, 12, 5, 10, 12, 6, 3, 10, 6, 3, 10, 6, 12,
					11, 4, 10, 12, 6, 11, 8, 7, 9, 11, 4, 8, 9, 5, 8, 9, 4, 8,
					12, 5, 11, 12, 5, 10, 12, 6, 11, 12, 3, 10, 11, 4, 6, 12, 5,
					11, 10, 3, 12, 11, 7, 8, 9, 4, 8, 9, 5, 8, 9, 4, 8, 11, 7,
					12, 10, 5, 11, 8, 6, 4, 12, 6, 11, 12, 5, 10, 3, 12, 11, 3,
					9, 8, 4, 9, 8, 5, 10, 8, 4, 11, 12, 5, 12, 3, 7, 11, 6, 12,
					10, 6, 12, 10, 3, 11, 7, 12, 11, 3, 9, 7, 4, 8, 9, 5, 8, 11,
					12, 6, 10, 12, 5, 10, 6, 11, 10, 6, 12, 10, 6, 11, 12, 7, 3,
					1, 1, 1, 3, 8, 9, 4, 8, 9, 5, 8, 9, 4},
			{10, 4, 6, 12, 8, 5, 12, 10, 3, 6, 11, 4, 12, 9, 5, 11, 8, 6, 11,
					10, 7, 11, 5, 4, 9, 8, 7, 6, 8, 3, 10, 9, 5, 10, 4, 9, 3,
					7},},
			{{9, 7, 12, 11, 9, 10, 3, 12, 11, 8, 10, 4, 12, 6, 10, 5, 12, 11, 7,
					10, 6, 12, 8, 9, 5, 11, 8},
					{3, 12, 11, 7, 10, 4, 11, 12, 7, 11, 12, 5, 11, 9, 6, 10, 3,
							9, 4, 11, 9, 5, 8, 11, 7, 8, 6, 10, 8},
					{3, 11, 10, 7, 11, 9, 3, 10, 7, 11, 9, 4, 12, 10, 4, 9, 6,
							5, 10, 6, 5, 12, 4, 5, 11, 6, 7, 8, 10, 9, 3, 9, 8,
							7, 12, 4, 9, 5, 12, 8, 6, 12, 7, 1, 1, 1, 12, 11, 8,
							4, 11, 10},
					{8, 9, 5, 11, 12, 3, 11, 12, 5, 10, 12, 6, 3, 10, 6, 3, 10,
							6, 12, 11, 4, 10, 12, 6, 11, 8, 7, 9, 11, 4, 8, 9,
							5, 8, 9, 4, 8, 12, 5, 11, 12, 5, 10, 12, 6, 11, 12,
							3, 10, 11, 6, 4, 12, 6, 11, 10, 3, 12, 11, 7, 8, 9,
							4, 8, 9, 5, 8, 9, 4, 8, 11, 7, 12, 10, 5, 11, 8, 6,
							4, 12, 6, 11, 12, 5, 10, 3, 12, 11, 3, 9, 8, 4, 9,
							8, 5, 10, 8, 4, 11, 12, 5, 12, 3, 7, 11, 6, 12, 10,
							6, 12, 10, 3, 11, 7, 12, 11, 3, 9, 7, 4, 8, 9, 5, 8,
							11, 12, 5, 10, 12, 5, 10, 6, 11, 10, 6, 12, 10, 6,
							11, 12, 7, 3, 1, 1, 1, 3, 8, 9, 4, 8, 9, 5, 8, 9,
							4},
					{10, 4, 6, 12, 8, 5, 12, 10, 3, 6, 11, 4, 12, 9, 5, 11, 8,
							6, 11, 10, 7, 11, 5, 4, 9, 7, 8, 6, 3, 8, 10, 9, 5,
							10, 4, 9, 3, 7},},};

	/**
	 * Creates a population from an original sequence pattern.
	 * 
	 * @param reel
	 *            Virtual reel as numbers array.
	 * 
	 * @param original
	 *            The chromosome of the original sequence.
	 * 
	 * @param size
	 *            Size of the population.
	 * 
	 * @return Randomly generated population.
	 */
	private static List<Chromosome> initializeRandomPopulation(int[] reel,
			Chromosome original, int size) {
		List<Chromosome> result = new ArrayList<Chromosome>();

		/* Create random initial chromosomes. */
		for (int i = 0; i < size; i++) {
			Chromosome candidate = Chromosome.initializeRandom(original);

			/* Evaluate randomly generated chromosome. */
			candidate.sampling(original);
			candidate.fitness(-candidate.distance(original));

			/* Add randomly generated chromosome to the population. */
			result.add(candidate);
		}

		return result;
	}

	/**
	 * Finds the best-found solution.
	 * 
	 * @param population
	 *            Current generation as population of individuals.
	 * 
	 * @return A reference to the best-found solution into the population.
	 */
	private static Chromosome bestFound(List<Chromosome> population) {
		/*
		 * There is no way to have a best-found solution if the population is
		 * empty.
		 */
		if (population.size() <= 0) {
			throw new RuntimeException(
					"Population size should be greater than zero!");
		}

		/* The best-found solution is the one with the highest fitness value. */
		Chromosome result = population.get(0);
		for (Chromosome candidate : population) {
			if (candidate.fitness() > result.fitness()) {
				result = candidate;
			}
		}

		return result;
	}

	/**
	 * Do selection of parents and a child place into the population.
	 * 
	 * @param population
	 *            Current generation as population of individuals.
	 * 
	 * @return Selected parents and children as an array of references.
	 */
	private static Chromosome[] selection(List<Chromosome> population) {
		Chromosome familiy[] = new Chromosome[3];
		while (true) {
			familiy[0] = population.get(PRNG.nextInt(population.size()));
			familiy[1] = population.get(PRNG.nextInt(population.size()));
			familiy[2] = population.get(PRNG.nextInt(population.size()));

			/* Parent should be different from the child. */
			if (familiy[0] == familiy[2]) {
				continue;
			}

			/* Parent should be different from the child. */
			if (familiy[1] == familiy[2]) {
				continue;
			}

			/* Parents should be different. */
			if (familiy[0] == familiy[1]) {
				continue;
			}

			/*
			 * Appointed for a child individual in the genetic algorithm
			 * population will replace the previous one that is why the weakest
			 * should be selected.
			 */
			if (familiy[2].fitness() > familiy[0].fitness()) {
				continue;
			}
			if (familiy[2].fitness() > familiy[1].fitness()) {
				continue;
			}

			/*
			 * If parents are different from the child, different from each
			 * other, and the weakest is chosen for population removal go on.
			 */
			break;
		}

		return familiy;
	}

	/**
	 * A simple form of genetic algorithm.
	 * 
	 * @param reel
	 *            Single reel as an array of numbers.
	 */
	private static void simpleGeneticAlgorithm(int[] reel) {
		System.err.println("=== OPTIMIZATION START ===");
		Chromosome original = Chromosome.initializeOriginal(reel, CHUNKS_SIZE,
				HISTOGRAM_THRESHOLD);
		// System.err.println(original);

		List<Chromosome> population = initializeRandomPopulation(reel, original,
				POPULATION_SIZE);
		// System.err.println(population);

		/*
		 * Do an evolutionary optimization. Each loop only a single genetic
		 * algorithm child is created that is why population size should be
		 * multiplied by the number of required generations.
		 */
		for (long g = EVOLUTION_EPOCHS * population.size(); g > 0; g--) {
			/* Select parents and a child slot. */
			Chromosome familiy[] = selection(population);

			/* Stronger parent is the first one. */
			Chromosome parent1 = (familiy[0].fitness() > familiy[1].fitness())
					? familiy[0]
					: familiy[1];
			Chromosome parent2 = (familiy[0].fitness() < familiy[1].fitness())
					? familiy[0]
					: familiy[1];

			/* Crossover between parents. */
			Chromosome child = parent1.crossover(parent2);

			/*
			 * Mutation done according to original chunks available values.
			 */
			child.mutate(original, MUTATION_RATE);

			/*
			 * Evaluate fitness value of the newly created child.
			 * 
			 * Distance is taken with a negative sign because if the candidate
			 * solution is farther away from the original the solution is worse.
			 * 
			 * With such an evaluation of the fitness, all values will be
			 * negative, but the smallest distance gives the best-found
			 * candidate solution.
			 */
			child.sampling(original);
			child.fitness(-child.distance(original));

			/* The new generation replaces the old generation. */
			if (child.fitness() > familiy[2].fitness()) {
				population.remove(familiy[2]);
				population.add(child);
			}

			/* Report optimization progress at each generation. */
			if (g % population.size() == 0) {
				System.err.print(g);
				System.err.print("\t");
				System.err.println(bestFound(population).fitness());
			}
		}

		/* Print the original. */
		System.out.println("=== ORIGIANL ===");
		System.out.println(original);
		System.out.println();

		/* Print the best-found solution. */
		System.out.println("=== BEST FOUND ===");
		System.out.println(bestFound(population));
		System.out.println();

		System.err.println("=== OPTIMIZATION END ===");
	}

	/**
	 * A recursive descent form of genetic algorithm.
	 * 
	 * @param depth
	 *            Level of recursive descent (zero is the bottom).
	 * 
	 * @param original
	 *            The chromosome of the original sequence.
	 * 
	 * @return The best-found solution.
	 */
	private static Chromosome recursiveOptimalSolution(int depth,
			Chromosome original) {
		/*
		 * Recursive depth is identical to the population size. If the recursive
		 * level is below or equal to zero, there is an best-found solution.
		 */
		if (depth <= 0) {
			return null;
		}

		/*
		 * If the recursive level is one there will be only one random
		 * individual and it will be returned.
		 */
		if (depth == 1) {
			/* Create a random solution. */
			Chromosome child = Chromosome.initializeRandom(original);
			
			/* Evaluate the random solution.. */
			child.sampling(original);
			child.fitness(-child.distance(original));
			
			/* Return newly created random solution. */
			return child;
		}

		/*
		 * Build the local population on the specified recursive level according
		 * to best-found individuals from the sub-levels.
		 */
		List<Chromosome> population = new ArrayList<Chromosome>(depth);
		for (int i = 0; i < depth; i++) {
			population.add(recursiveOptimalSolution(depth - 1, original));
		}

		/*
		 * Apply local search until better solutions are found in the local
		 * recursive level population.
		 */
		boolean stop = false;
		Chromosome result = bestFound(population);
		while (stop == false) {
			stop = true;

			/* Crossover and mutation with each other. */
			for (Chromosome first : population) {
				for (Chromosome second : population) {
					/* Crossover. */
					Chromosome child = first.crossover(second);

					/* Mutation. */
					child.mutate(original, MUTATION_RATE);

					/* Evaluation. */
					child.sampling(original);
					child.fitness(-child.distance(original));

					/* Selection. */
					if (child.fitness() > result.fitness()) {
						result = child;
						stop = false;
					}
				}
			}
		}

		/*
		 * Return the best-found solution from the local search on the current
		 * recursive node.
		 */
		return result;
	}

	/**
	 * A hierarchical form of genetic algorithm.
	 * 
	 * @param reel
	 *            Single reel as an array of numbers.
	 */
	private static void hierarchicalGeneticAlgorithm(int[] reel) {
		System.err.println("=== OPTIMIZATION START ===");

		/* Creation of the chromosome with chunks from the original reel. */
		Chromosome original = Chromosome.initializeOriginal(reel, CHUNKS_SIZE,
				HISTOGRAM_THRESHOLD);

		/* Get a recursive optimal solution. */
		Chromosome best = recursiveOptimalSolution(RECURSION_DEPTH, original);

		/* Print the original. */
		System.out.println("=== ORIGIANL ===");
		System.out.println(original);
		System.out.println();

		/* Print the best-found solution. */
		System.out.println("=== BEST FOUND ===");
		System.out.println(best);
		System.out.println();

		System.err.println("=== OPTIMIZATION END ===");
	}

	/**
	 * Application single entry point method.
	 * 
	 * @param args
	 *            Command line arguments.
	 */
	public static void main(String[] args) {
		/* Handle each virtual reel separate. */
		for (int reels[][] : ORIGINAL_STRIPS) {
			System.err.println("=== RELLS ===");
			System.out.println("=== RELLS ===");
			System.out.println();

			for (int reel[] : reels) {
				// simpleGeneticAlgorithm(reel);
				hierarchicalGeneticAlgorithm(reel);
			}
		}
	}

}
