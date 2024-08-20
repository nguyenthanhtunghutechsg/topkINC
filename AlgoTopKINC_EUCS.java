package TKHUI_INC_EUCS;

/* This file is copyright (c) 2008-2015 Philippe Fournier-Viger
* 
* This file is part of the SPMF DATA MINING SOFTWARE
* (http://www.philippe-fournier-viger.com/spmf).
* 
* SPMF is free software: you can redistribute it and/or modify it under the
* terms of the GNU General Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your option) any later
* version.
* 
* SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
* A PARTICULAR PURPOSE. See the GNU General Public License for more details.
* You should have received a copy of the GNU General Public License along with
* SPMF. If not, see <http://www.gnu.org/licenses/>.
* 
*/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import java.util.PriorityQueue;

/**
 * This is an implementation of the "FHM" algorithm for High-Utility Itemsets
 * Mining as described in the conference paper : <br/>
 * <br/>
 * 
 * Fournier-Viger, P., Wu, C.-W., Zida, S., Tseng, V. (2014) FHM: A Faster
 * High-Utility Itemset Mining Algorithm using Estimated Utility Co-occurrence
 * Pruning. Proc. 21st International Symposium on Methodologies for Intelligent
 * Systems (ISMIS 2014), Springer, LNAI, 12 pages (to appear).
 *
 * @see UtilityList
 * @see Element
 * @author Philippe Fournier-Viger
 */
public class AlgoTopKINC_EUCS {

	/** the time at which the algorithm started */
	public long startTimestamp = 0;

	/** the time at which the algorithm ended */
	public long endTimestamp = 0;

	int k = 0;
	/** the number of candidate high-utility itemsets */
	public long candidateCount = 0;

	/** Map to remember the TWU of each item */
	Map<Integer, Long> mapItemToTWU;

	/** writer to write the output file */
	BufferedWriter writer = null;

	/** The eucs structure: key: item key: another item value: twu */
	Map<Integer, Map<Integer, PairItem>> mapEUCSLeaf;

	/** enable LA-prune strategy */
	boolean ENABLE_LA_PRUNE = true;

	/** variable for debug mode */
	boolean DEBUG = false;

	/** this class represent an item and its utility in a transaction */
	class Item {
		int item;
		int utility;

		Item(int item, int utility) {
			this.item = item;
			this.utility = utility;
		}
	}

	class PairItem {
		long twu = 0;
		int utility = 0;
	}

	/**
	 * Default constructor
	 */
	public AlgoTopKINC_EUCS() {

	}

	/**
	 * Run the algorithm
	 * 
	 * @param input      the input file path
	 * @param output     the output file path
	 * @param minUtility the minimum utility threshold
	 * @throws IOException exception if error while writing the file
	 */

	int firstLine;
	long totalDBUtility = 0;
	long min_utility = 0;
	long first_min_utility = 0;
	Map<Integer, UtilityList> mapItemToUtilityList;
	Map<Integer, Integer> mapItemToUtility;
	PriorityQueue<Itemset> kItemsets;
	boolean EUCS_PRUNING;
	boolean LEAF_PRUNING;

	public void runAlgorithm(String input, String output, int k, int firstLine, int lastLine, boolean EUCS,
			boolean LEAF) throws IOException {
		// reset maximum
		MemoryLogger.getInstance().reset();

		this.firstLine = firstLine;
		writer = new BufferedWriter(new FileWriter(output));
		EUCS_PRUNING = EUCS;
		LEAF_PRUNING = LEAF;
		boolean firstTime = (mapEUCSLeaf == null);
		startTimestamp = System.currentTimeMillis();
		BufferedReader myInput = null;
		String thisLine;
		kItemsets = new PriorityQueue<Itemset>();
		if (firstTime) {
			mapEUCSLeaf = new HashMap<Integer, Map<Integer, PairItem>>();
			mapItemToUtilityList = new HashMap<Integer, UtilityList>();
			totalDBUtility = 0;
			mapItemToTWU = new HashMap<Integer, Long>();
			mapItemToUtility = new HashMap<Integer, Integer>();
			this.k = k;
		}
		int tid = 0;
		try {
			// prepare the object for reading the file
			myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(input))));
			// for each line (transaction) until the end of file
			while ((thisLine = myInput.readLine()) != null && tid < lastLine) {
				// if the line is a comment, is empty or is a
				// kind of metadata
				if (tid >= firstLine) {
					if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%'
							|| thisLine.charAt(0) == '@') {
						continue;
					}

					// split the transaction according to the : separator
					String split[] = thisLine.split(":");
					// the first part is the list of items
					String items[] = split[0].split(" ");
					int[] arrayItems = Arrays.stream(items).mapToInt(Integer::parseInt).toArray();
					String utilities[] = split[2].split(" ");
					int[] arrayUtilities = Arrays.stream(utilities).mapToInt(Integer::parseInt).toArray();
					// the second part is the transaction utility
					Integer transactionUtility = Integer.parseInt(split[1]);
					// for each item, we add the transaction utility to its TWU
					List<Item> transactionItems = new ArrayList<>();
					for (int i = 0; i < arrayItems.length; i++) {
						Integer itemInTran = arrayItems[i];
						Integer utilityInTran = arrayUtilities[i];
						Item itemObj = new Item(itemInTran,utilityInTran);
						transactionItems.add(itemObj);
					}
					transactionItems.sort(new Comparator<Item>() {
						@Override
						public int compare(Item o1, Item o2) {
							// TODO Auto-generated method stub
							return o1.item-o2.item;
						}
					});
//						// convert item to integer
//						Integer item = arrayItems[i];
//						Integer utility = arrayUtilities[i];
//						// get the current TWU of that item
//						Long twu = mapItemToTWU.get(item);
//						Integer getedUtility = mapItemToUtility.get(item);
//						// add the utility of the item in the current transaction to its twu
//						if (twu == null) {
//							twu = (long)transactionUtility;
//							getedUtility = utility;
//							UtilityList newUL = new UtilityList(item);
//							Element element = new Element(tid, utility, 0);
//							newUL.addElement(element);
//							mapItemToUtilityList.put(item, newUL);
//
//						} else {
//							Element element = new Element(tid, utility, 0);
//							UtilityList uLItem = mapItemToUtilityList.get(item);
//							uLItem.addElement(element);
//							mapItemToUtilityList.put(item, uLItem);
//							twu = twu + transactionUtility;
//							getedUtility = getedUtility + utility;
//						}
//						mapItemToTWU.put(item, twu);
//						mapItemToUtility.put(item, getedUtility);
//						// EUCS
//						Map<Integer, Long> mapEUCSItem = mapEUCS.get(item);
//						if (mapEUCSItem == null) {
//							mapEUCSItem = new HashMap<Integer, Long>();
//							mapEUCS.put(item, mapEUCSItem);
//						}
//						for (int j = i + 1; j < arrayUtilities.length; j++) {
//							int ItemAfter = arrayItems[j];
//							Long twuSum = mapEUCSItem.get(ItemAfter);
//							if (twuSum == null) {
//								mapEUCSItem.put(ItemAfter, (long)transactionUtility);
//							} else {
//								mapEUCSItem.put(ItemAfter, twuSum + transactionUtility);
//							}
//						}
//					}
					totalDBUtility += transactionUtility;
				}
				tid++;
			}
		} catch (Exception e) {
			// catches exception if error while reading the input file
			e.printStackTrace();
		} finally {
			if (myInput != null) {
				myInput.close();
			}
		}

		List<Integer> listAllItemIntegers = new ArrayList<Integer>();
		for (Entry<Integer, Integer> entry : mapItemToUtility.entrySet()) {
			listAllItemIntegers.add(entry.getKey());
		}
		Collections.sort(listAllItemIntegers, new Comparator<Integer>() {
			@Override
			public int compare(Integer item1, Integer item2) {
				// TODO Auto-generated method stub
				return mapItemToUtility.get(item2) - mapItemToUtility.get(item1);
			}
		});

		if (firstTime) {
			if (k > listAllItemIntegers.size()) {
				min_utility = 1;
			} else {
				int itemk = listAllItemIntegers.get(k - 1);
				min_utility = mapItemToUtility.get(itemk);
			}
		} else {
			if (k <= listAllItemIntegers.size()) {
				int itemk = listAllItemIntegers.get(k - 1);
				int newMin_utility = mapItemToUtility.get(itemk);
				if (newMin_utility > min_utility) {
					min_utility = newMin_utility;
				}
			}
		}
//		if (firstTime) {
//			if (k > listAllItemIntegers.size()) {
//				min_utility = 1;
//			} else {
//				int itemk = listAllItemIntegers.get(k - 1);
//				min_utility = mapItemToUtility.get(itemk);
//			}
//		} else {
//			if (k > listAllItemIntegers.size()) {
//				min_utility = 1;
//			} else {
//				int itemk = listAllItemIntegers.get(k - 1);
//				min_utility = mapItemToUtility.get(itemk);
//			}
//		}
		first_min_utility = min_utility;
		List<UtilityList> listOfUtilityLists = new ArrayList<UtilityList>();
		for (Entry<Integer, UtilityList> entry : mapItemToUtilityList.entrySet()) {
			if (mapItemToTWU.get(entry.getKey()) >= min_utility) {
				listOfUtilityLists.add(entry.getValue());
			}
		}

		Collections.sort(listOfUtilityLists, new Comparator<UtilityList>() {
			public int compare(UtilityList o1, UtilityList o2) {
				return compareItems(o1.item, o2.item);
			}
		});

		int arrayRu[] = new int[tid + 1];
		for (int i = listOfUtilityLists.size() - 1; i >= 0; i--) {
			UtilityList ul = listOfUtilityLists.get(i);
			int newRemain = 0;
			for (int j = 0; j < ul.elements.size(); j++) {
				Element element = ul.elements.get(j);
				element.rutils = arrayRu[element.tid];
				arrayRu[element.tid] += element.iutils;
				newRemain += element.rutils;
			}
			ul.sumRutils = newRemain;
		}
//		for (UtilityList ul : listOfUtilityLists) {
//			System.out.println(ul.toString());
//		}

		// check the memory usage
		MemoryLogger.getInstance().checkMemory();
		System.out.println("mining... " + min_utility);
		// Mine the database recursively
		fhm(new int[0], 0, null, listOfUtilityLists);

		// check the memory usage again and close the file.
		MemoryLogger.getInstance().checkMemory();
		// close output file
		writer.close();
		// record end time
		endTimestamp = System.currentTimeMillis();
	}

	/**
	 * Method to compare items by their TWU
	 * 
	 * @param item1 an item
	 * @param item2 another item
	 * @return 0 if the same item, >0 if item1 is larger than item2, <0 otherwise
	 */
	private int compareItems(int item1, int item2) {
		int compare = (int) (mapItemToTWU.get(item1) - mapItemToTWU.get(item2));
		// if the same, use the lexical order otherwise use the TWU
		return (compare == 0) ? item1 - item2 : compare;
	}

	/**
	 * This is the recursive method to find all high utility itemsets. It writes the
	 * itemsets to the output file.
	 * 
	 * @param prefix       This is the current prefix. Initially, it is empty.
	 * @param pUL          This is the Utility List of the prefix. Initially, it is
	 *                     empty.
	 * @param ULs          The utility lists corresponding to each extension of the
	 *                     prefix.
	 * @param minUtility   The minUtility threshold.
	 * @param prefixLength The current prefix length
	 * @throws IOException
	 */
	private void fhm(int[] prefix, int prefixLength, UtilityList pUL, List<UtilityList> ULs) throws IOException {

		// For each extension X of prefix P
		for (int i = 0; i < ULs.size(); i++) {
			UtilityList X = ULs.get(i);
			candidateCount++;
			// If pX is a high utility itemset.
			// we save the itemset: pX
			if (X.sumIutils >= min_utility) {
				// save to file
				writeOut(prefix, prefixLength, X.item, X.sumIutils);
			}

			// If the sum of the remaining utilities for pX
			// is higher than minUtility, we explore extensions of pX.
			// (this is the pruning condition)
			if (X.sumIutils + X.sumRutils >= min_utility) {
				// This list will contain the utility lists of pX extensions.
				List<UtilityList> exULs = new ArrayList<UtilityList>();
				// For each extension of p appearing
				// after X according to the ascending order
				for (int j = i + 1; j < ULs.size(); j++) {
					UtilityList Y = ULs.get(j);

					// ======================== NEW OPTIMIZATION USED IN FHM
					Map<Integer, Long> mapTWUF = mapEUCS.get(X.item);
					Map<Integer, Long> mapTWUFRev = mapEUCS.get(Y.item);
					long SumEUCS = 0l;
					if (mapTWUF != null) {
						Long twuF = mapTWUF.get(Y.item);
						if (twuF != null) {
							SumEUCS += twuF;
						}
					}
					if (mapTWUFRev != null) {
						Long twuF = mapTWUFRev.get(X.item);
						if (twuF != null) {
							SumEUCS += twuF;
						}
					}
					if (SumEUCS < min_utility) {
						continue;
					}

					// =========================== END OF NEW OPTIMIZATION

					// we construct the extension pXY
					// and add it to the list of extensions of pX
					UtilityList temp = construct(pUL, X, Y);
					if (temp != null) {
						exULs.add(temp);
					}
				}
				// We create new prefix pX
				int[] newPrefix = new int[prefix.length + 1];
				System.arraycopy(prefix, 0, newPrefix, 0, prefix.length);
				newPrefix[prefixLength] = X.item;
				// We make a recursive call to discover all itemsets with the prefix pXY
				fhm(newPrefix, prefixLength + 1, X, exULs);
			}
		}
		MemoryLogger.getInstance().checkMemory();
	}

	/**
	 * This method constructs the utility list of pXY
	 * 
	 * @param P  : the utility list of prefix P.
	 * @param px : the utility list of pX
	 * @param py : the utility list of pY
	 * @return the utility list of pXY
	 */
	private UtilityList construct(UtilityList P, UtilityList px, UtilityList py) {
		// create an empy utility list for pXY
		UtilityList pxyUL = new UtilityList(py.item);

		// == new optimization - LA-prune == /
		// Initialize the sum of total utility
		long totalUtility = px.sumIutils + px.sumRutils;
		// ================================================

		// for each element in the utility list of pX
		for (Element ex : px.elements) {
			// do a binary search to find element ey in py with tid = ex.tid
			Element ey = findElementWithTID(py, ex.tid);
			if (ey == null) {
				// == new optimization - LA-prune == /
				if (ENABLE_LA_PRUNE) {
					totalUtility -= (ex.iutils + ex.rutils);
					if (totalUtility < min_utility) {
						return null;
					}
				}
				// =============================================== /
				continue;
			}
			// if the prefix p is null
			if (P == null) {
				// Create the new element
				Element eXY = new Element(ex.tid, ex.iutils + ey.iutils, ey.rutils);
				// add the new element to the utility list of pXY
				pxyUL.addElement(eXY);

			} else {
				// find the element in the utility list of p wih the same tid
				Element e = findElementWithTID(P, ex.tid);
				if (e != null) {
					// Create new element
					Element eXY = new Element(ex.tid, ex.iutils + ey.iutils - e.iutils, ey.rutils);
					// add the new element to the utility list of pXY
					pxyUL.addElement(eXY);
				}
			}
		}
		// return the utility list of pXY.
		return pxyUL;
	}

	/**
	 * Do a binary search to find the element with a given tid in a utility list
	 * 
	 * @param ulist the utility list
	 * @param tid   the tid
	 * @return the element or null if none has the tid.
	 */
	private Element findElementWithTID(UtilityList ulist, int tid) {
		List<Element> list = ulist.elements;

		// perform a binary search to check if the subset appears in level k-1.
		int first = 0;
		int last = list.size() - 1;

		// the binary search
		while (first <= last) {
			int middle = (first + last) >>> 1; // divide by 2

			if (list.get(middle).tid < tid) {
				first = middle + 1; // the itemset compared is larger than the subset according to the lexical order
			} else if (list.get(middle).tid > tid) {
				last = middle - 1; // the itemset compared is smaller than the subset is smaller according to the
									// lexical order
			} else {
				return list.get(middle);
			}
		}
		return null;
	}

	/**
	 * Method to write a high utility itemset to the output file.
	 * 
	 * @param the          prefix to be writent o the output file
	 * @param an           item to be appended to the prefix
	 * @param utility      the utility of the prefix concatenated with the item
	 * @param prefixLength the prefix length
	 */
	private void writeOut(int[] prefix, int prefixLength, int item, long utility) throws IOException {
		Itemset itemset = new Itemset(prefix, item, utility);
		kItemsets.add(itemset);
		if (kItemsets.size() > k) {
			if (utility > this.min_utility) {
				Itemset lower;
				do {
					lower = kItemsets.peek();
					if (lower == null) {
						break;
					}
					kItemsets.remove(lower);
				} while (kItemsets.size() > k);
				this.min_utility = kItemsets.peek().utility;
			}
		}

//		StringBuilder buffer = new StringBuilder();
//		// append the prefix
//		for (int i = 0; i < prefixLength; i++) {
//			buffer.append(prefix[i]);
//			buffer.append(' ');
//		}
//		// append the last item
//		buffer.append(item);
//		// append the utility value
//		buffer.append(" #UTIL: ");
//		buffer.append(utility);
//
//		// write to file
//		System.out.println(buffer.toString());
////		writer.newLine();
	}

	/**
	 * Print statistics about the latest execution to System.out.
	 * 
	 * @throws IOException
	 */
	public void printStats() throws IOException {
		System.out.println("=============  TOP-K INCRE - EUCS- SPMF 0.97e - STATS =============");
		System.out.println(" Total time ~ " + (endTimestamp - startTimestamp) + " ms");
		System.out.println(" Memory ~ " + MemoryLogger.getInstance().getMaxMemory() + " MB");
		System.out.println(" minU : " + min_utility);
		System.out.println(" FirstMinU : " + first_min_utility);
		System.out.println(" totalDBUtility : " + totalDBUtility);
		System.out.println(" Candidate count : " + candidateCount);

		if (DEBUG) {
			int pairCount = 0;
			double maxMemory = getObjectSize(mapEUCS);
			for (Entry<Integer, Map<Integer, Long>> entry : mapEUCS.entrySet()) {
				maxMemory += getObjectSize(entry.getKey());
				for (Entry<Integer, Long> entry2 : entry.getValue().entrySet()) {
					pairCount++;
					maxMemory += getObjectSize(entry2.getKey()) + getObjectSize(entry2.getValue());
				}
			}
			System.out.println("CMAP size " + maxMemory + " MB");
			System.out.println("PAIR COUNT " + pairCount);
		}
		System.out.println("===================================================");

//		Iterator<Itemset> iter = kItemsets.iterator();
//		while (iter.hasNext()) {
//			StringBuffer buffer = new StringBuffer();
//			Itemset itemset = (Itemset) iter.next();
//
//			// append the prefix
//			for (int i = 0; i < itemset.getItemset().length; i++) {
//				buffer.append(itemset.getItemset()[i]);
//				buffer.append(' ');
//			}
//			buffer.append(itemset.item);
//
//			// append the utility value
//			buffer.append(" #UTIL: ");
//			buffer.append(itemset.utility);
//
//			// write to file
//			System.out.println(buffer.toString());
//			if (iter.hasNext()) {
//
//			}
//		}
	}

	/**
	 * Get the size of a Java object (for debugging purposes)
	 * 
	 * @param object the object
	 * @return the size in MB
	 * @throws IOException
	 */
	private double getObjectSize(Object object) throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		ObjectOutputStream oos = new ObjectOutputStream(baos);
		oos.writeObject(object);
		oos.close();
		double maxMemory = baos.size() / 1024d / 1024d;
		return maxMemory;
	}
}