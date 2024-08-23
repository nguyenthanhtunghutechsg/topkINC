package BaseLine;

import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;

public class MainBaseLine {

	public static void main(String[] arg) throws IOException {

		String input = "DB.txt";
		String output = ".//output.txt";
		int k=15;

		// the number of updates to be performed
		int numberOfUpdates = 1;

		// scan the database to count the number of lines
		// for our test purpose
		int linecount = countLines(input);

		double addedratio = 1d / ((double) numberOfUpdates);
		int linesForeEachUpdate = (int) (addedratio * linecount);

		// Apply the algorithm several times
		AlgoTopKINC algo = new AlgoTopKINC();
		int firstLine = 0;
		for (int i = 0; i < numberOfUpdates; i++) {
			int lastLine = firstLine + linesForeEachUpdate;
			//

			// Applying the algorithm
			// If this is the last update, we make sure to run until the last line
			if (i == numberOfUpdates - 1) {
				System.out.println("" + i + ") Run the algorithm using line " + firstLine + " to before line "
						+ linecount + " of the input database.");
				algo.runAlgorithm(input,output, k, firstLine, linecount);
			} else {
				// If this is not the last update
				System.out.println("" + i + ") Run the algorithm using line " + firstLine + " to before line "
						+ lastLine + " of the input database.");
				algo.runAlgorithm(input,output, k, firstLine, lastLine);
			}
			algo.printStats();

			firstLine = lastLine;
		}

	}

	/**
	 * This methods counts the number of lines in a text file.
	 * 
	 * @param filepath the path to the file
	 * @return the number of lines as an int
	 * @throws IOException Exception if error reading/writting file
	 */
	public static int countLines(String filepath) throws IOException {
		LineNumberReader reader = new LineNumberReader(new FileReader(filepath));
		while (reader.readLine() != null) {
		}
		int count = reader.getLineNumber();
		reader.close();
		return count;
	}

}
