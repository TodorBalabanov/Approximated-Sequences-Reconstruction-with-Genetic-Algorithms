import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
	 * Chunks comparator instance.
	 */
	private static final ChunksComparator CHUNKS_COMAPARATOR = new ChunksComparator();

	/**
	 * Original sequences which should be reconstructed.
	 */
	private static final int[][] ORIGINAL_STRIPS = {
			{3, 7, 4, 7, 11, 12, 6, 11, 7, 7, 9, 4, 10, 12, 7, 7, 5, 9, 8, 1, 9,
					8, 10, 9, 7, 5, 5, 7, 9, 10, 10, 12, 6, 6, 10, 10, 8, 11, 8,
					7, 4, 3, 5, 12, 6, 9, 8, 1, 9, 8, 9, 8, 4, 3, 11, 11, 7, 9,
					8, 11, 7, 3, 11},
			{3, 9, 4, 11, 7, 4, 11, 9, 10, 5, 5, 11, 3, 9, 8, 1, 9, 8, 11, 6, 6,
					6, 10, 10, 12, 12, 7, 6, 11, 10, 10, 10, 6, 11, 7, 4, 3, 10,
					5, 6, 9, 4, 6, 7, 7, 7, 12, 12, 11, 3, 10, 11, 11, 8, 12, 9,
					9, 9, 11, 7, 4, 11, 10},
			{8, 11, 10, 7, 9, 11, 3, 8, 8, 7, 11, 3, 8, 8, 3, 3, 11, 6, 7, 4,
					11, 6, 16, 12, 5, 12, 7, 4, 10, 10, 6, 8, 8, 6, 6, 6, 7, 3,
					6, 10, 10, 6, 3, 10, 9, 7, 7, 4, 10, 3, 10, 8, 10, 16, 10,
					5, 9, 5, 4, 9, 10, 3, 10},
			{3, 8, 7, 5, 16, 6, 3, 5, 3, 9, 3, 7, 7, 12, 3, 7, 9, 4, 10, 6, 5,
					8, 4, 7, 9, 11, 7, 6, 6, 5, 16, 8, 7, 7, 4, 8, 10, 9, 16,
					10, 3, 5, 3, 8, 4, 9, 3, 4, 9, 10, 10, 7, 5, 5, 6, 10, 3,
					12, 16, 10, 10, 5, 3},
			{10, 3, 5, 3, 10, 8, 10, 7, 5, 6, 5, 5, 7, 3, 10, 3, 4, 8, 3, 12, 5,
					4, 3, 8, 10, 4, 6, 3, 9, 8, 10, 10, 5, 10, 6, 3, 8, 5, 8, 8,
					4, 5, 5, 7, 3, 4, 6, 6, 5, 3, 3, 5, 6, 3, 3, 12, 5, 4, 4, 3,
					3, 5, 5},};

	/**
	 * Chunk size is related with the size of the visible part of the sequence.
	 */
	private static final int CHUNKS_SIZE = 3;

	/**
	 * The histogram threshold is the minimum number of the appearance of the
	 * less probable unique chunk. When the size of the sequence is unknown by
	 * estimating the less probable unique chunk appearance the amount of chunks
	 * sample can be estimated.This estimation should be done with confidence
	 * level of about 95% and according to the rules of the normal probability
	 * distribution.
	 */
	private static final int HISTOGRAM_THRESHOLD = 10;

	/** How many genetic algorithm generations to be evolved. */
	private static final int EVOLUTION_EPOCHS = 1000;

	/**
	 * Sequence chunks comparator.
	 * 
	 * @author Todor Balabanov
	 */
	private static class ChunksComparator implements Comparator<List<Integer>> {
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
	}

	/**
	 * Genetic algorithm chromosome representation.
	 * 
	 * @author Todor Balabanov
	 */
	private static class Chromosome {
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
		 * A constructor used to generate chromosome by getting only the
		 * sequence. Chunks and fitness are initialized after that.
		 * 
		 * @param sequence
		 *            Sequence represented as a single string.
		 */
		public Chromosome(int[] sequence) {
			super();
			this.sequence = sequence;
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

		@Override
		public String toString() {
			return "Chromosome [sequence=" + Arrays.toString(sequence)
					+ ", chunks=" + chunks + ", fitness=" + fitness + "]";
		}

	}

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
	private static Chromosome initializeOriginal(int[] reel, int chunkSize,
			int histogramThreshold) {
		/* Build a chunks histogram. */
		int leastProbable = 0;
		Map<List<Integer>, Integer> histogram = new HashMap<List<Integer>, Integer>();
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
				/* If the chunk is presented its counter should be increased. */
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

		/* Form list of chunks according the amount of their appearance. */
		List<List<Integer>> chunks = new ArrayList<List<Integer>>();
		for (List<Integer> chunk : histogram.keySet()) {
			for (int count = histogram.get(chunk); count > 0; count--) {
				chunks.add(chunk);
			}
		}

		Chromosome result = new Chromosome(reel);
		result.chunks(chunks);

		return result;
	}

	/**
	 * Application single entry point method.
	 * 
	 * @param args
	 *            Command line arguments.
	 */
	public static void main(String[] args) {
		/* Handle each virtual separate. */
		for (int reel[] : ORIGINAL_STRIPS) {
			Chromosome original = initializeOriginal(reel, CHUNKS_SIZE,
					HISTOGRAM_THRESHOLD);
			System.err.println(original);
		}
	}

}
