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
					3, 5, 5},},

			{{7, 7, 4, 16, 3, 12, 6, 9, 7, 3, 9, 10, 10, 3, 7, 3, 7, 5, 8, 8, 9,
					8, 3, 9, 7, 3, 11, 16, 9, 10, 10, 3, 6, 6, 10, 16, 8, 8, 8,
					7, 10, 12, 4, 4, 6, 8, 4, 8, 7, 8, 9, 8, 10, 7, 11, 9, 7, 9,
					8, 9, 7, 4, 4},
					{12, 9, 4, 11, 9, 8, 1, 9, 8, 5, 3, 11, 3, 7, 7, 8, 11, 3,
							3, 6, 6, 6, 8, 10, 4, 12, 7, 6, 3, 4, 10, 4, 5, 3,
							8, 8, 4, 12, 12, 9, 8, 1, 9, 8, 7, 7, 12, 4, 11, 5,
							10, 12, 11, 8, 8, 9, 9, 9, 10, 7, 11, 11, 10},
					{8, 11, 10, 7, 9, 11, 5, 5, 3, 7, 11, 3, 8, 8, 5, 11, 11, 6,
							7, 11, 5, 6, 11, 12, 5, 12, 7, 10, 10, 10, 8, 8, 8,
							6, 6, 6, 7, 10, 9, 10, 10, 6, 10, 10, 9, 7, 7, 4,
							10, 9, 10, 8, 10, 10, 10, 5, 9, 5, 10, 9, 10, 10,
							10},
					{8, 8, 7, 9, 10, 6, 5, 11, 5, 9, 3, 12, 7, 12, 7, 7, 9, 12,
							10, 6, 8, 8, 12, 7, 9, 11, 12, 6, 6, 8, 8, 8, 7, 16,
							4, 4, 4, 9, 8, 10, 3, 5, 12, 8, 4, 9, 3, 11, 9, 10,
							16, 7, 11, 11, 6, 10, 12, 3, 3, 10, 5, 11, 5},
					{10, 3, 8, 8, 10, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 10, 8,
							8, 12, 5, 4, 10, 4, 10, 8, 12, 10, 3, 3, 12, 10, 5,
							10, 6, 11, 8, 11, 11, 8, 8, 10, 10, 7, 11, 4, 3, 11,
							11, 3, 3, 12, 6, 12, 3, 12, 10, 4, 12, 11, 11, 11,
							12},},

			{{16, 7, 7, 16, 10, 12, 16, 9, 7, 16, 11, 4, 16, 12, 7, 16, 7, 12,
					16, 7, 10, 16, 10, 9, 16, 11, 12, 16, 7, 10, 16, 12, 6, 16,
					10, 10, 16, 8, 7, 16, 4, 12, 16, 12, 7, 16, 10, 10, 16, 10,
					12, 16, 12, 8, 16, 8, 7, 16, 12, 12, 16, 11, 11},
					{9, 8, 16, 3, 8, 16, 5, 8, 16, 5, 9, 16, 8, 9, 16, 8, 9, 16,
							8, 9, 16, 8, 8, 16, 3, 9, 16, 6, 4, 16, 4, 4, 16, 3,
							3, 16, 8, 3, 16, 6, 9, 16, 4, 7, 16, 3, 4, 16, 4, 5,
							16, 4, 3, 16, 9, 9, 16, 9, 3, 16, 4, 3, 16},
					{10, 10, 12, 11, 8, 11, 8, 12, 12, 6, 9, 12, 11, 7, 12, 7,
							12, 3, 10, 9, 10, 8, 8, 12, 5, 4, 10, 4, 10, 8, 12,
							10, 9, 12, 11, 10, 5, 10, 6, 11, 8, 8, 3, 7, 8, 8,
							11, 11, 11, 6, 7, 6, 11, 16, 12, 12, 7, 6, 11, 7, 4,
							10, 7},
					{9, 12, 16, 11, 7, 16, 5, 5, 16, 5, 11, 16, 3, 7, 16, 7, 11,
							16, 11, 6, 16, 6, 12, 16, 12, 12, 16, 6, 11, 16, 10,
							8, 16, 11, 8, 16, 8, 12, 16, 6, 9, 16, 12, 7, 16, 7,
							12, 16, 11, 5, 16, 12, 11, 16, 11, 8, 16, 9, 10, 16,
							6, 12, 16},
					{10, 3, 8, 8, 8, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 10, 8,
							8, 12, 5, 4, 10, 4, 10, 8, 12, 10, 9, 12, 12, 10, 5,
							10, 6, 11, 8, 11, 11, 8, 8, 10, 10, 7, 11, 4, 9, 11,
							11, 8, 8, 12, 6, 12, 3, 12, 10, 4, 16, 11, 11, 11,
							12},},

			{{7, 7, 4, 12, 3, 12, 6, 9, 7, 7, 9, 4, 10, 12, 7, 7, 7, 5, 8, 9, 9,
					8, 10, 9, 7, 11, 11, 12, 9, 10, 10, 12, 6, 6, 10, 8, 8, 8,
					8, 7, 8, 12, 12, 12, 8, 5, 8, 7, 9, 8, 9, 8, 8, 7, 11, 9, 7,
					9, 8, 9, 7, 11, 11},
					{12, 9, 4, 11, 7, 11, 5, 5, 10, 5, 11, 11, 3, 7, 6, 11, 1,
							6, 11, 6, 6, 6, 16, 10, 12, 12, 7, 6, 11, 1, 6, 11,
							5, 11, 8, 8, 8, 16, 12, 6, 6, 11, 1, 6, 11, 7, 12,
							12, 11, 5, 10, 12, 11, 8, 16, 9, 9, 9, 10, 7, 11,
							11, 10},
					{8, 11, 10, 7, 9, 11, 12, 8, 8, 7, 11, 3, 8, 8, 11, 11, 11,
							6, 7, 11, 11, 6, 16, 12, 5, 12, 7, 10, 11, 10, 8, 8,
							8, 6, 6, 6, 7, 10, 9, 10, 10, 6, 10, 10, 9, 7, 7, 4,
							10, 9, 10, 8, 10, 16, 10, 5, 12, 12, 10, 9, 10, 10,
							10},
					{8, 11, 7, 9, 10, 6, 11, 11, 11, 9, 3, 12, 6, 11, 1, 6, 11,
							12, 10, 6, 8, 12, 12, 6, 11, 1, 6, 11, 6, 8, 8, 8,
							7, 7, 4, 4, 4, 9, 8, 10, 3, 5, 12, 8, 4, 9, 11, 11,
							9, 10, 11, 7, 11, 6, 11, 1, 6, 11, 3, 10, 5, 11,
							12},
					{10, 3, 8, 8, 9, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 10, 8,
							8, 12, 5, 4, 10, 4, 10, 8, 12, 10, 9, 12, 11, 10, 5,
							10, 6, 11, 8, 11, 12, 8, 11, 11, 10, 7, 12, 4, 9,
							11, 11, 12, 8, 12, 6, 12, 3, 12, 10, 4, 10, 11, 11,
							11, 12},},

			{{3, 7, 10, 10, 7, 10, 3, 3, 7, 10, 10, 7, 10, 3, 3, 7, 10, 10, 7,
					10, 3, 3, 7, 10, 10, 7, 10, 3, 3, 7, 10, 10, 7, 10, 3, 3, 7,
					10, 10, 7, 10, 3, 3, 7, 10, 10, 7, 10, 3, 3, 7, 10, 10, 7,
					10, 3, 3, 7, 10, 10, 7, 10, 3},
					{4, 4, 8, 16, 4, 8, 8, 4, 8, 16, 4, 8, 16, 8, 4, 4, 8, 16,
							4, 8, 8, 4, 8, 16, 4, 8, 16, 8, 4, 4, 8, 16, 4, 8,
							8, 4, 8, 16, 4, 8, 16, 8, 4, 4, 8, 16, 4, 8, 8, 4,
							8, 16, 4, 8, 16, 8, 4, 8, 16, 4, 8, 16, 8},
					{5, 5, 7, 10, 7, 10, 7, 10, 10, 5, 7, 10, 5, 7, 5, 5, 7, 10,
							7, 10, 7, 10, 10, 5, 7, 10, 5, 7, 5, 5, 7, 10, 7,
							10, 7, 10, 10, 5, 7, 10, 5, 7, 5, 5, 7, 10, 7, 10,
							7, 10, 10, 5, 7, 10, 5, 7, 10, 10, 5, 7, 10, 5, 7},
					{9, 6, 9, 16, 6, 6, 6, 9, 16, 9, 9, 6, 9, 16, 9, 6, 9, 16,
							6, 6, 6, 9, 16, 9, 9, 6, 9, 16, 9, 6, 9, 16, 6, 6,
							6, 9, 16, 9, 9, 6, 9, 16, 9, 6, 9, 16, 6, 6, 6, 9,
							16, 9, 9, 6, 16, 9, 9, 16, 9, 9, 6, 9, 16},
					{12, 9, 8, 12, 9, 8, 9, 8, 12, 8, 8, 9, 8, 8, 12, 9, 8, 12,
							9, 8, 9, 8, 12, 8, 8, 9, 8, 8, 12, 9, 8, 12, 9, 8,
							9, 8, 12, 8, 8, 9, 8, 8, 12, 9, 8, 12, 9, 8, 9, 8,
							12, 8, 8, 9, 8, 8, 8, 12, 8, 8, 9, 8, 8},},

			{{3, 7, 4, 11, 11, 12, 6, 11, 7, 3, 9, 4, 10, 12, 16, 7, 5, 9, 8, 1,
					9, 8, 10, 9, 7, 16, 11, 7, 9, 10, 10, 12, 6, 6, 10, 16, 8,
					8, 8, 7, 4, 12, 12, 12, 6, 9, 8, 1, 9, 8, 9, 8, 4, 7, 11, 4,
					7, 9, 8, 11, 7, 16, 11},
					{12, 9, 4, 11, 7, 4, 11, 9, 10, 5, 11, 11, 3, 7, 7, 8, 3,
							11, 3, 6, 6, 6, 9, 8, 1, 9, 8, 6, 6, 7, 4, 11, 6,
							11, 8, 8, 8, 8, 12, 6, 9, 4, 12, 7, 16, 7, 12, 12,
							11, 9, 10, 12, 11, 8, 4, 9, 9, 9, 10, 7, 4, 11, 10},
					{8, 11, 10, 7, 9, 11, 12, 8, 8, 7, 11, 3, 8, 8, 11, 16, 11,
							6, 7, 4, 11, 6, 16, 12, 5, 12, 7, 4, 10, 10, 8, 8,
							8, 6, 6, 6, 7, 10, 9, 10, 10, 6, 10, 10, 9, 7, 16,
							4, 10, 9, 10, 8, 10, 11, 10, 5, 9, 5, 4, 9, 10, 10,
							10},
					{8, 8, 7, 9, 10, 6, 11, 11, 11, 9, 3, 12, 7, 12, 7, 16, 9,
							4, 10, 6, 8, 8, 4, 7, 9, 16, 12, 6, 6, 8, 8, 8, 7,
							16, 4, 8, 10, 9, 8, 10, 3, 5, 12, 8, 4, 9, 11, 4, 9,
							10, 8, 7, 11, 11, 6, 10, 12, 16, 3, 10, 10, 11, 12},
					{10, 3, 16, 8, 11, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 4, 8,
							8, 12, 5, 4, 10, 8, 10, 4, 12, 16, 9, 12, 8, 10, 5,
							10, 6, 11, 8, 16, 8, 8, 8, 5, 5, 7, 11, 4, 9, 11,
							11, 16, 8, 12, 6, 12, 3, 7, 5, 4, 8, 11, 11, 11,
							12},},

			{{3, 7, 4, 16, 11, 12, 6, 11, 7, 7, 9, 4, 10, 12, 7, 7, 5, 9, 8, 1,
					9, 8, 10, 9, 7, 11, 11, 7, 9, 10, 10, 12, 6, 6, 10, 10, 8,
					8, 8, 7, 4, 12, 12, 12, 6, 9, 8, 1, 9, 8, 9, 8, 4, 7, 11, 9,
					7, 9, 8, 11, 7, 11, 11},
					{12, 9, 4, 11, 7, 4, 16, 9, 10, 5, 11, 11, 3, 7, 7, 8, 16,
							11, 11, 6, 6, 6, 8, 10, 3, 12, 7, 6, 11, 10, 10, 10,
							9, 11, 8, 8, 8, 9, 12, 6, 9, 4, 12, 7, 7, 7, 12, 16,
							11, 9, 10, 12, 11, 8, 9, 9, 9, 9, 10, 7, 4, 11, 10},
					{8, 11, 10, 7, 9, 11, 12, 8, 8, 7, 11, 3, 8, 8, 11, 11, 11,
							6, 7, 4, 16, 6, 8, 12, 5, 12, 7, 4, 10, 16, 8, 8, 8,
							6, 6, 6, 3, 10, 9, 10, 16, 6, 10, 10, 9, 7, 7, 4,
							10, 9, 10, 8, 9, 9, 10, 5, 9, 5, 4, 9, 10, 10, 10},
					{8, 8, 7, 9, 10, 6, 11, 11, 11, 9, 3, 12, 7, 12, 16, 7, 9,
							4, 10, 6, 8, 8, 4, 7, 9, 11, 12, 6, 6, 8, 8, 8, 7,
							16, 4, 8, 10, 9, 8, 10, 3, 5, 12, 8, 4, 9, 16, 4, 9,
							10, 9, 7, 9, 9, 6, 10, 12, 11, 3, 10, 16, 11, 12},
					{10, 3, 8, 8, 16, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 4, 8,
							8, 16, 5, 4, 10, 8, 10, 8, 3, 10, 9, 12, 9, 10, 16,
							10, 6, 11, 8, 11, 5, 8, 8, 5, 5, 7, 11, 8, 9, 11,
							11, 8, 8, 12, 6, 12, 3, 12, 5, 4, 9, 11, 11, 11,
							12},},

			{{3, 7, 4, 6, 11, 12, 6, 11, 7, 7, 9, 4, 10, 12, 7, 7, 5, 9, 8, 1,
					9, 8, 10, 9, 7, 11, 11, 7, 9, 10, 10, 12, 6, 6, 10, 6, 8, 8,
					8, 7, 4, 12, 12, 12, 6, 9, 8, 1, 9, 8, 9, 8, 4, 7, 11, 6, 7,
					9, 8, 11, 7, 11, 11},
					{12, 9, 4, 11, 7, 4, 11, 9, 10, 5, 11, 11, 3, 7, 7, 8, 11,
							11, 11, 6, 6, 6, 6, 10, 12, 12, 7, 6, 11, 10, 10,
							10, 9, 11, 8, 8, 8, 6, 12, 6, 9, 4, 12, 7, 7, 7, 12,
							12, 11, 9, 10, 12, 11, 8, 6, 9, 9, 9, 10, 7, 4, 11,
							10},
					{8, 11, 10, 7, 9, 11, 12, 8, 8, 7, 11, 3, 8, 8, 11, 11, 11,
							6, 7, 4, 11, 6, 6, 12, 5, 12, 7, 4, 10, 10, 8, 8, 8,
							6, 6, 6, 7, 10, 9, 10, 10, 6, 10, 10, 9, 7, 7, 4,
							10, 9, 10, 8, 10, 6, 10, 5, 9, 5, 4, 9, 10, 10, 10},
					{8, 8, 7, 9, 10, 6, 11, 11, 11, 9, 3, 12, 7, 12, 7, 7, 9, 4,
							10, 6, 16, 8, 4, 7, 9, 11, 12, 6, 6, 8, 16, 8, 7, 6,
							4, 8, 10, 9, 8, 10, 3, 5, 12, 8, 4, 9, 11, 4, 9, 10,
							6, 7, 11, 11, 6, 10, 16, 11, 3, 10, 10, 11, 12},
					{10, 3, 8, 8, 6, 8, 10, 7, 9, 6, 8, 8, 7, 3, 10, 9, 4, 8, 8,
							12, 5, 4, 10, 16, 10, 4, 12, 10, 9, 12, 6, 10, 5,
							10, 6, 11, 8, 11, 6, 8, 8, 5, 5, 7, 11, 4, 9, 11,
							11, 16, 8, 12, 6, 12, 3, 12, 5, 4, 6, 11, 11, 11,
							12},},};

	/**
	 * Genetic algorithm chromosome representation.
	 * 
	 * @author Todor Balabanov
	 */
	private static class Chromosome {

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
					throw new RuntimeException(
							"Chunks should be with equal size!");
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
		 * original sequence. This estimation is very useful for estimating how
		 * many chunks to be presented in the original chromosome.
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
				 * If the chunk is not presented it appears for the first time
				 * in the current sample.
				 */
				if (histogram.containsKey(chunk) == false) {
					/*
					 * If the chunk is not presented it appears for the first
					 * time in the current sample.
					 */
					leastProbable = 1;
					histogram.put(chunk, 1);
				} else {
					/*
					 * If the chunk is presented its counter should be
					 * increased.
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
				 * The count can be reduced with least probable value plus one
				 * in order sequences for chromosomes to be shorter.
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
		 *            A sample chromosome which is used during chromosome
		 *            creation.
		 * 
		 * @return Randomly initialized chromosome.
		 */
		public static Chromosome initializeRandom(Chromosome sample) {
			/*
			 * When sequence size is not known in advance it is difficult to
			 * guess the real size. Random size between the number of unique
			 * values and total length of the chunks is used.
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
		 *            A sample chromosome which is used during chromosome
		 *            creation.
		 */
		private void sampling(Chromosome sample) {
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

		private double levenshtein(List<Integer> first, List<Integer> second) {
			int[][] matrix = new int[first.size() + 1][second.size() + 1];
 
			for (int i = 0; i <= first.size(); i++) {			
				for (int j = 0; j <= second.size(); j++) {
					if (i == 0) {
						matrix[i][j] = j;
					}
					else if (j == 0) {
						matrix[i][j] = i;
					} else {
						int cost = (first.get(i - 1) == second.get(j - 1)) ? 0 : 1;
						int a = matrix[i - 1][j - 1] + cost; 
						int b = matrix[i - 1][j] + 1; 
						int c = matrix[i][j - 1] + 1;
						matrix[i][j] = (a < b) ? ((a<c)?a:c) : ((b<c)?b:c);
					}
				}
			}
 
			return matrix[first.size()][second.size()];
		}

		/**
		 * Calculates the distance between two chromosomes.
		 * 
		 * Weighted Euclidean distance between lists of chunks is calculated,
		 * but other distances may be more proper.
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
			return "Chromosome [sequence=" + Arrays.toString(sequence)
					+ ", chunks=" + chunks + ", fitness=" + fitness + "]";
		}

	}

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
	 * Finds the best-found solution.
	 * 
	 * @param population
	 *            Current generation as population of individuals.
	 * 
	 * @return A reference to the best-found solution into the population.
	 */
	private static Chromosome bestFound(List<Chromosome> population) {
		Chromosome result = population.get(0);

		/* The best-found solution is the one with the highest fitness value. */
		for (Chromosome candidate : population) {
			if (candidate.fitness() > result.fitness()) {
				result = candidate;
			}
		}

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
		for (int reels[][] : ORIGINAL_STRIPS) {
			System.err.println(reels);
			System.out.println("=== NEW RELLS ===");
			System.out.println();


			for (int reel[] : reels) {
				System.err.println(reel);
				Chromosome original = Chromosome.initializeOriginal(reel,
						CHUNKS_SIZE, HISTOGRAM_THRESHOLD);
				// System.err.println(original);

				List<Chromosome> population = initializeRandomPopulation(reel,
						original, POPULATION_SIZE);
				// System.err.println(population);

				/*
				 * Do an evolutionary optimization. Each loop only a single
				 * genetic algorithm child is created that is why population
				 * size should be multiplied by the number of required
				 * generations.
				 */
				for (long g = EVOLUTION_EPOCHS
						* population.size(); g > 0; g--) {
					/* Select parents and a child slot. */
					Chromosome familiy[] = selection(population);

					/* Stronger parent is the first one. */
					Chromosome parent1 = (familiy[0].fitness() > familiy[1]
							.fitness()) ? familiy[0] : familiy[1];
					Chromosome parent2 = (familiy[0].fitness() < familiy[1]
							.fitness()) ? familiy[0] : familiy[1];

					/* Crossover between parents. */
					Chromosome child = parent1.crossover(parent2);

					/*
					 * Mutation done according to original chunks available
					 * values.
					 */
					child.mutate(original, MUTATION_RATE);

					/*
					 * Evaluate fitness value of the newly created child.
					 * 
					 * Distance is taken with a negative sign because if the
					 * candidate solution is farther away from the original the
					 * solution is worse.
					 * 
					 * With such an evaluation of the fitness, all values will
					 * be negative, but the smallest distance gives the
					 * best-found candidate solution.
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
			}
		}
	}

}
