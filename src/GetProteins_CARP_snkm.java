/***
 * 2026-09-22 Update for CARP snakemake file rule get_protein_info_from_consensus
 * Delete
 * private static String dDir = "../";
	private static String data = "output/notKnown.fa";
	private static String ipFile = dDir + "output/notKnown.fa.spwb.gff";
	private static String opFile = dDir + "output/protein.txt";
	private static String nknpFile = dDir + "output/notKnownNotProtein.fa";
 * 
 */
import java.io.*;
import java.util.*;

public class GetProteins_CARP_snkm {

	private static int minLength = 21;

	
	public static void main(String[] args) {
      if (args.length != 4) {
         System.err.println(
            "Usage: java GetProteins_CARP_snkm " +
            "<not_known.fa> <protein_hits.gff> " +
            "<protein.txt> <not_protein.fa>"
         );
         System.exit(1);
      }

      String data = args[0];
      String ipFile = args[1];
      String opFile = args[2];
      String nknpFile = args[3];

		Hashtable<String, String> proteins = getProteins(ipFile);
		writeProteinFamilies(opFile, proteins);
		writeNonProteinSequences(data, nknpFile, proteins);
	}



	private static void addProtein (Hashtable<String, String> proteins, String line) {
		StringTokenizer st = new StringTokenizer(line);
		String family = st.nextToken();
		for (int i=1; i<3; i++) st.nextToken();
		try {
			int start = Integer.parseInt(st.nextToken());
			int end = Integer.parseInt(st.nextToken());
			if ((end - start + 1) >= minLength) {
				for (int i=1; i<5; i++) st.nextToken();
				String protein = st.nextToken();
				if (proteins.containsKey(family)) {
					String prots = proteins.get(family);
					prots += ":" + protein;
					proteins.put(family, prots);
					}
				else proteins.put(family, protein);
				}
			}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);}
		}

	private static Hashtable<String, String> getProteins(String inFile) {
		Hashtable<String, String> proteins = new Hashtable<String, String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			String line = null;
			while ((line = in.readLine()) != null) addProtein (proteins, line);
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		System.out.println("" + proteins.size() + " consensus sequences have been identified as proteins");
		return proteins;
		}

	private static void writeNonProteinSequences (String inFile, String outFile, 
		Hashtable<String, String> proteins) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			boolean wanted = false;
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.charAt(0) == '>') 
					wanted = !proteins.containsKey(line.substring(1).split(" ")[0]);
				if (wanted) out.write(line + "\n");
				}
			in.close();
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
		
	private static void writeProteinFamilies (String outFile, Hashtable<String, String> proteins) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			out.write("Sequence MappedTo\n");
			for (Enumeration<String> e=proteins.keys(); e.hasMoreElements();) {
				String family = e.nextElement();
				out.write(family + " " + proteins.get(family) + "\n");
				}
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}		
		
	}
