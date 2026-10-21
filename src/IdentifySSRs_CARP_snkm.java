/***
 * 2026-09-22 Update for CARP snakemake identify_ssr_from_phobos_output block
 * Delete: 
 * private static String dDir = "../";
	private static String inPrefix = dDir + "output/notKnownNotProtein";
	private static String pFile = inPrefix + ".phobos";
 */
import java.io.*;
import java.util.*;


public class IdentifySSRs_CARP_snkm {
	
	public static void main(String[] args) {
		if (args.length != 4) {
			System.err.println(
            "Usage: java IdentifySSRs_CARP_snkm " +
            "<phobos> <not_protein.fa> <ssr.txt> <not_ssr.fa>"
			);
			System.exit(1);
		}

		String pFile = args[0];
		String notProtein = args[1];
		String ssrOutput = args[2];
		String notSsrOutput = args[3];

		Hashtable<String, String> ssrs = getSSRs(pFile);

		outputSSRs(ssrs, ssrOutput);
		outputUnknownFasta(ssrs, notProtein, notSsrOutput);
	}
		

	private static double addFamilySSR (String name, BufferedReader in, Hashtable<String, String> ssrs) {
		String line = null;
		double maxCovered = 0.0;
		try {
			line = in.readLine();
			int length = Integer.parseInt(line.substring(line.indexOf(":")+2));
			double lowerLimit = length * .1;
			double upperLimit = length * .9;
			in.readLine();
			StringTokenizer st = null;
			while ((line = in.readLine()) != null && line.charAt(0) != '#') {
				st = new StringTokenizer(line);
				st.nextToken();
				int start = Integer.parseInt(st.nextToken());
				if (start <= lowerLimit) { 
					st.nextToken();
					int end = Integer.parseInt(st.nextToken());
					double covered = (end - start + 1) * 100./length;
					if (covered > maxCovered) maxCovered = covered;
					if (end >= upperLimit) {
						for (int i=1; i<15; i++) st.nextToken();
						ssrs.put(name, st.nextToken());
						}
					}
				}
			}
		catch (IOException ie) {ie.printStackTrace();}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);}
		return maxCovered;
		}
		 		
	private static Hashtable<String, String> getSSRs (String inFile) {
		Hashtable<String, String> ssrs = new Hashtable<String, String> ();
		double maxCovered = 0.0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			String line = null;
			double covered = 0.0;
			while ((line = in.readLine()) != null) 
				if (line.charAt(0) == '>') {
					covered = addFamilySSR(line.substring(1), in, ssrs); 
					if (covered > maxCovered) maxCovered = covered;
					}
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		System.out.println("There are " + ssrs.size() + " families that are SSRs, maximum coverage was " + 
			maxCovered + " percent");
		return ssrs;
		}

	private static void outputSSRs (Hashtable<String, String> ssrs, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			out.write("Sequence MappedTo\n");
			for (Enumeration<String> e=ssrs.keys(); e.hasMoreElements();) {
				String key = e.nextElement();
				out.write(key + " " + ssrs.get(key) + "\n");
				}
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
		
	private static void outputUnknownFasta (Hashtable<String, String> ssrs, String inFile, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			boolean wanted = false;
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.charAt(0) == '>') 
					wanted = !ssrs.containsKey(line.substring(1).split(" ")[0]);
				if (wanted) out.write(line + "\n");
				}
			in.close();
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		} 
}
