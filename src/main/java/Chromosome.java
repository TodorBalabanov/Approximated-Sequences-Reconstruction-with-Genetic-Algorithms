import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Genetic algorithm chromosome representation.
 * 
 * @author Todor Balabanov
 */
class Chromosome {

	/**
	 * Sequence chunks comparator.
	 */
	private static final Comparator<List<Integer>> CHUNKS_COMAPARATOR = new Comparator<List<Integer>>() {
		@Override
		public int compare(List<Integer> first, List<Integer> second) {
			if (first == null) {
				throw new RuntimeException(
						"First chunks should not be null pointer!");
			}

			if (second == null) {
				throw new RuntimeException(
						"Second chunks should not be null pointer!");
			}

			if (first.size() != second.size()) {
				throw new RuntimeException("Chunks should be with equal size!");
			}

			int size = Math.min(first.size(), second.size());
			for (int i = 0; i < size; i++) {
				if (first.get(i) - second.get(i) == 0) {
					continue;
				}

				return first.get(i) - second.get(i);
			}

			return 0;
		}
	};

	/**
	 * Pseudo-random number generator.
	 */
	private static final Random PRNG = new Random();

	/**
	 * Chunks histogram is used to estimate how often chunks are met in the
	 * original sequence. This estimation is very useful for estimating how many
	 * chunks to be presented in the original chromosome.
	 */
	private static Map<List<Integer>, Integer> histogram = new HashMap<List<Integer>, Integer>();

	/**
	 * Chromosome minimal length.
	 */
	private static int minLength = 0;

	/**
	 * Chromosome maximal length.
	 */
	private static int maxLength = 0;

	/**
	 * Optimal sequence candidate.
	 */
	private int[] sequence = {};

	/**
	 * List of randomly generated chunks.
	 */
	private List<List<Integer>> chunks = new ArrayList<List<Integer>>();

	/**
	 * The chromosome fitness value is used during the parents selection
	 * process.
	 */
	private double fitness = 0.0;

	/**
	 * Creates a chromosome from an original sequence pattern.
	 * 
	 * @param reel
	 *            Pattern as numbers.
	 * 
	 * @param chunkSize
	 *            The size of the chunks represented in the chromosome.
	 * 
	 * @param histogramThreshold
	 *            Minimum count of the least probable chunk from the sample.
	 * 
	 * @return An original reel chromosome representation.
	 */
	public static Chromosome initializeOriginal(int[] reel, int chunkSize,
			int histogramThreshold) {
		/* Build a chunks histogram. */
		int leastProbable = 0;
		histogram = new HashMap<List<Integer>, Integer>();
		while (leastProbable < histogramThreshold) {
			/* Form a single chunk. */
			List<Integer> chunk = new ArrayList<Integer>();
			int position = PRNG.nextInt(reel.length);
			for (int i = 0; i < chunkSize; i++) {
				chunk.add(reel[(position + i) % reel.length]);
			}

			/*
			 * If the chunk is not presented it appears for the first time in
			 * the current sample.
			 */
			if (histogram.containsKey(chunk) == false) {
				/*
				 * If the chunk is not presented it appears for the first time
				 * in the current sample.
				 */
				leastProbable = 1;
				histogram.put(chunk, 1);
			} else {
				/*
				 * If the chunk is presented its counter should be increased.
				 */
				int count = histogram.get(chunk) + 1;
				histogram.put(chunk, count);

				/* Check what is the current least probable chunk. */
				leastProbable = Integer.MAX_VALUE;
				for (List<Integer> key : histogram.keySet()) {
					if (histogram.get(key) < leastProbable) {
						leastProbable = histogram.get(key);
					}
				}
			}
		}
		// System.err.println(histogram);

		/* Form list of chunks according the amount of their appearance. */
		List<List<Integer>> chunks = new ArrayList<List<Integer>>();
		for (List<Integer> chunk : histogram.keySet()) {
			/*
			 * The count can be reduced with least probable value plus one in
			 * order sequences for chromosomes to be shorter.
			 */
			for (int count = histogram.get(chunk); count > 0; count--) {
				chunks.add(chunk);
			}
		}

		/* Create and initialize original. */
		Chromosome result = new Chromosome();
		result.sequence(reel);
		result.chunks(chunks);

		/* Estimation of the unique chunks and unique values amount. */
		int chunksTotalLength = 0;
		Set<Integer> uniqueValues = new HashSet<Integer>();
		Set<List<Integer>> uniqueChunks = new HashSet<List<Integer>>();
		for (List<Integer> chunk : result.chunks()) {
			chunksTotalLength += chunk.size();

			/* Update set of unique chunks. */
			uniqueChunks.add(chunk);

			for (Integer value : chunk) {
				/* Update set of unique values. */
				uniqueValues.add(value);
			}
		}
		minLength = uniqueValues.size();
		maxLength = chunksTotalLength / uniqueChunks.size();
		// System.err.println(uniqueValues);

		return result;
	}

	/**
	 * Construct randomly generated chromosomes according to a given sample.
	 * 
	 * @param sample
	 *            A sample chromosome which is used during chromosome creation.
	 * 
	 * @return Randomly initialized chromosome.
	 */
	public static Chromosome initializeRandom(Chromosome sample) {
		/*
		 * When sequence size is not known in advance it is difficult to guess
		 * the real size. Random size between the number of unique values and
		 * total length of the chunks is used.
		 */
		int sequence[] = new int[minLength
				+ PRNG.nextInt(maxLength - minLength + 1)];

		/*
		 * Fill the candidate sequence with values from the original chunks.
		 */
		for (int j = 0; j < sequence.length; j++) {
			sequence[j] = sample.randomValue();
		}
		// System.err.println(Arrays.toString(sequence));

		/* Form chromosome. */
		Chromosome result = new Chromosome();
		result.sequence(sequence);
		result.sampling(sample);

		return result;
	}

	/**
	 * Constructor without parameters.
	 */
	public Chromosome() {
		super();
	}

	/**
	 * An integer array representation of the sequence getter.
	 * 
	 * @return Sequence as a string reference.
	 */
	public int[] sequence() {
		return sequence;
	}

	/**
	 * An integer array representation of the sequence setter.
	 * 
	 * @param sequence
	 *            Sequence as a string reference.
	 */
	public void sequence(int[] sequence) {
		this.sequence = sequence;
	}

	/**
	 * A list of chunks of the sequence getter.
	 * 
	 * @return Randomly observed chunks of the sequence.
	 */
	public List<List<Integer>> chunks() {
		return chunks;
	}

	/**
	 * A list of chunks of the sequence getter.
	 * 
	 * @param chunks
	 *            Randomly observed chunks of the sequence.
	 */
	public void chunks(List<List<Integer>> chunks) {
		this.chunks = chunks;

		/* Chunks should be sorted when fitness value is estimated. */
		Collections.sort(chunks, CHUNKS_COMAPARATOR);
	}

	/**
	 * Chromosome fitness value getter.
	 * 
	 * @return Chromosome fitness value.
	 */
	public double fitness() {
		return fitness;
	}

	/**
	 * Chromosome fitness value setter.
	 * 
	 * @param fitness
	 *            Chromosome fitness value.
	 */
	public void fitness(double fitness) {
		this.fitness = fitness;
	}

	/**
	 * Provides random value from the chunks.
	 * 
	 * @return Randomly selected value from a randomly selected chunk.
	 */
	private int randomValue() {
		List<Integer> chunk = chunks().get(PRNG.nextInt(chunks().size()));
		return chunk.get(PRNG.nextInt(chunk.size()));
	}

	/**
	 * Do sampling from the available sequence with parameters for sampling
	 * taken from the template chromosome.
	 * 
	 * @param sample
	 *            A sample chromosome which is used during chromosome creation.
	 */
	void sampling(Chromosome sample) {
		/* Generate chunks for the candidate sequence. */
		List<List<Integer>> chunks = new ArrayList<List<Integer>>();
		for (int j = 0; j < sample.chunks().size(); j++) {
			/* Form a single chunk. */
			List<Integer> chunk = new ArrayList<Integer>();
			int position = PRNG.nextInt(sequence.length);
			for (int k = 0; k < sample.chunks().get(j).size(); k++) {
				chunk.add(sequence[(position + k) % sequence.length]);
			}

			/* Add new chunk to the chunks list. */
			chunks.add(chunk);
		}

		chunks(chunks);
	}

	/**
	 * Calculate Euclidean distance between chunks.
	 * 
	 * @param first
	 *            First chunk.
	 * @param second
	 *            Second chunk.
	 * 
	 * @return Distance calculated.
	 */
	private double euclidean(List<Integer> first, List<Integer> second) {
		int chunkSize = Math.min(first.size(), second.size());

		double distance = 0;
		for (int j = 0; j < chunkSize; j++) {
			distance += (first.get(j) - second.get(j))
					* (first.get(j) - second.get(j));
		}

		/* Square root as it is in the Euclidean norm. */
		return Math.sqrt(distance);
	}

	/**
	 * Calculate Levenshtein distance between chunks.
	 * 
	 * @param first
	 *            First chunk.
	 * @param second
	 *            Second chunk.
	 * 
	 * @return Distance calculated.
	 */
	private double levenshtein(List<Integer> first, List<Integer> second) {
		int[][] matrix = new int[first.size() + 1][second.size() + 1];

		for (int i = 0; i <= first.size(); i++) {
			for (int j = 0; j <= second.size(); j++) {
				if (i == 0) {
					matrix[i][j] = j;
				} else if (j == 0) {
					matrix[i][j] = i;
				} else {
					int cost = (first.get(i - 1) == second.get(j - 1)) ? 0 : 1;

					int a = matrix[i - 1][j - 1] + cost;
					int b = matrix[i - 1][j] + 1;
					int c = matrix[i][j - 1] + 1;

					matrix[i][j] = (a < b)
							? ((a < c) ? a : c)
							: ((b < c) ? b : c);
				}
			}
		}

		return matrix[first.size()][second.size()];
	}

	/**
	 * Calculates the distance between two chromosomes.
	 * 
	 * Weighted Euclidean distance between lists of chunks is calculated, but
	 * other distances may be more proper.
	 * 
	 * @param sample
	 *            Sample chromosome to compare with.
	 * 
	 * @return Distance calculated between the two chromosomes.
	 */
	public double distance(Chromosome sample) {
		/* Chunks lists should be with equal length. */
		if (chunks.size() != sample.chunks.size()) {
			throw new RuntimeException(
					"The distance can be calculated only between lists of chunks with equal length!");
		}

		double result = 0;
		int listSize = Math.min(chunks.size(), sample.chunks.size());
		for (int i = 0; i < listSize; i++) {
			List<Integer> first = chunks.get(i);
			List<Integer> second = sample.chunks.get(i);

			/* Chunks lists should be with equal length. */
			if (first.size() != second.size()) {
				throw new RuntimeException(
						"Chunks should be with equal sizes!");
			}

			// result += euclidean(first, second);
			result += levenshtein(first, second);
		}

		/* Do an average. */
		result /= listSize;

		return result;
	}

	/**
	 * Chromosome mutation with a certain probability and source of mutation
	 * information.
	 * 
	 * @param sample
	 *            Source of mutation information.
	 * 
	 * @param rate
	 *            Mutation rate between 0 and 1.
	 */
	public void mutate(Chromosome sample, double rate) {
		for (List<Integer> chunk : chunks()) {
			for (int i = 0; i < chunk.size(); i++) {
				if (PRNG.nextDouble() >= rate) {
					continue;
				}

				/* Mutate only with the proper rate. */
				chunk.set(i, sample.randomValue());
			}
		}
	}

	/**
	 * Crossover with a mate.
	 * 
	 * @param mate
	 *            Mating chromosome.
	 * 
	 * @return Child chromosome after mating.
	 */
	public Chromosome crossover(Chromosome mate) {
		/* Mating threshold is around half of the genes. */
		double threshold = 0.5 + PRNG.nextGaussian() * 0.2;

		/* Child has variable length. */
		int sequence[] = new int[minLength
				+ PRNG.nextInt(maxLength - minLength + 1)];

		/*
		 * Fill the candidate sequence with values from the original chunks.
		 */
		int[] first = sequence();
		int[] second = mate.sequence();
		for (int i = 0; i < sequence.length; i++) {
			if (PRNG.nextDouble() < threshold) {
				sequence[i] = first[i % first.length];
			} else {
				sequence[i] = second[i % second.length];
			}
		}
		// System.err.println(Arrays.toString(sequence));

		/* Form chromosome. */
		Chromosome result = new Chromosome();
		result.sequence(sequence);

		return result;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "Chromosome [sequence=" + Arrays.toString(sequence) + ", chunks="
				+ chunks + ", fitness=" + fitness + "]";
	}

}