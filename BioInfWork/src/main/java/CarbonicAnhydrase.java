import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;
import javafx.util.Pair;

import java.io.*;
import java.util.*;

/**
 * Builds the alignment of the carbonic anhydrase.
 */
public class CarbonicAnhydrase {
    /**
     * Extracts the name of protein.
     * @param base given set of lines from database
     * @param line number of line
     * @return pair of protein name and sequence
     */
    private static Pair<String, String> extractProtein(ArrayList<String> base, Integer line) {
        StringBuilder proteinName = new StringBuilder();
        StringBuilder proteinSeq = new StringBuilder();

        Integer index = 0;

        while (base.get(line).charAt(index) != '>') {
            index++;
        }
        while (base.get(line).charAt(index) == '>') {
            index++;
        }
        while (base.get(line).charAt(index) != ' ') {
            proteinName.append(base.get(line).charAt(index));
            index++;
        }
        line++;

        while (!base.get(line).equals("")) {
            proteinSeq.append(base.get(line));
            line++;
        }

        return new Pair<>(proteinName.toString(), proteinSeq.toString());
    }

    /**
     * Gets pre-amino acid.
     * @param protein given protein
     * @return pre-amino acid of the peptide
     */
    static Character getPre(String protein) {
        return protein.charAt(protein.indexOf("pre=") + 4);
    }

    /**
     * Gets post-amino acid.
     * @param protein given protein
     * @return pre-amino acid of the peptide
     */
    private static Character getPost(String protein) {
        return protein.charAt(protein.indexOf("post=") + 5);
    }

    /**
     * Makes set of spaces to make an alignment.
     * @param n number of spaces
     * @return set of spaces to make an alignment
     */
    private static String getSpaces(int n) {
        StringBuilder spaces = new StringBuilder();
        for (int i = 0; i < n; i++) {
            spaces.append(" ");
        }
        return spaces.toString();
    }

    /**
     * Removes modification of peptide to get only peptide sequence.
     * @param peptide given peptide
     * @return peptide without modifications
     */
    static String removeModifications(String peptide) {
        StringBuilder peptideNew = new StringBuilder();
        for (int j = 0; j < peptide.length(); j++) {
            if (Character.isLetter(peptide.charAt(j))) {
                peptideNew.append(peptide.charAt(j));
            }
        }

        return peptideNew.toString();
    }

    /**
     * Makes an alignment of the protein.
     */
    public static void main(String[] args) throws IOException {
        TsvParserSettings settings = new TsvParserSettings();
        settings.getFormat().setLineSeparator("\n");
        TsvParser parser = new TsvParser(settings);

        // parses all rows in one go.
        List<String[]> allRows = parser.parseAll(new FileReader("../ChymotrypsinNTT1.tsv"));
        ArrayList<Pair<String, String>> peptides = new ArrayList<>();

        for (String[] allRow : allRows) {
            String peptideDraft = allRow[8];
            String peptide = removeModifications(peptideDraft);

            String protein = allRow[9];
            String evalue = allRow[13];
            Integer posMinus = evalue.indexOf('-');
            if (posMinus == -1) {
                continue;
            }
            Integer magnitude = Integer.parseInt(evalue.substring(posMinus + 1, evalue.length()));

            if (magnitude > 9) {
                peptides.add(new Pair<>(peptide, protein));
            }
        }

        String everything = null;

        ArrayList<String> lines = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader("../cah.fasta"))) {
            StringBuilder sb = new StringBuilder();
            String line = br.readLine();

            while (line != null) {
                lines.add(line);
                sb.append(line);
                sb.append(System.lineSeparator());
                line = br.readLine();
            }
            everything = sb.toString();
        }

        Integer line = 0;

        Pair<String, String> protein = extractProtein(lines, line);

        PrintWriter out = new PrintWriter("alignment_ChymotrypsinNTT1Maxcharge3_Len30.txt");
        out.println(protein.getKey());
        out.println(protein.getValue());

        ArrayList<TreeSet<String>> positions = new ArrayList<>();
        HashMap<String, Integer> peptideNumber = new HashMap<>();

        for (int i = 0; i < protein.getValue().length(); i++) {
            positions.add(new TreeSet<>((o1, o2) -> o2.length() - o1.length()));
        }

        for (Pair<String, String> pair : peptides) {
            String curPeptide = pair.getKey();
            String curProtein = pair.getValue();

            Integer number = peptideNumber.get(curPeptide);
            if (number == null) {
                number = 1;
            } else {
                number++;
            }
            peptideNumber.put(curPeptide, number);

            if (!curProtein.contains(protein.getKey())) {
                continue;
            }

            String pre = Character.toString(getPre(curProtein));
            String post = Character.toString(getPost(curProtein));

            if (pre.equals("-")) {
                positions.get(0).add(curPeptide);
            } else if (post.equals("-")) {
                positions.get(protein.getValue().length() - curPeptide.length()).add(curPeptide);
            } else {
                Integer ind = protein.getValue().indexOf(pre + curPeptide + post);
                positions.get(ind + 1).add(curPeptide);
            }
        }

        for (int i = 0; i < protein.getValue().length(); i++) {
            for (String peptide : positions.get(i)) {
                out.println(getSpaces(i) + peptide + " (" + peptideNumber.get(peptide) + ")");
            }
        }

        out.close();
    }
}
