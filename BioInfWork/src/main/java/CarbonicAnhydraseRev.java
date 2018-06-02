import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;
import javafx.util.Pair;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * Builds the alignment of the carbonic anhydrase to the right side.
 */
public class CarbonicAnhydraseRev {
    public static Pair<String, String> extractProtein(ArrayList<String> base, Integer line) {
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

    public static Character getPre(String protein) {
        return protein.charAt(protein.indexOf("pre=") + 4);
    }
    public static Character getPost(String protein) {
        return protein.charAt(protein.indexOf("post=") + 5);
    }

    public static String getSpaces(int n) {
        StringBuilder spaces = new StringBuilder();
        for (int i = 0; i < n; i++) {
            spaces.append(" ");
        }
        return spaces.toString();
    }

    public static void main(String[] args) throws IOException {
        TsvParserSettings settings = new TsvParserSettings();
        settings.getFormat().setLineSeparator("\n");
        TsvParser parser = new TsvParser(settings);

        // parses all rows in one go.
        List<String[]> allRows = parser.parseAll(new FileReader("../ChymotrypsinNTT1.tsv"));
        ArrayList<Pair<String, String>> peptides = new ArrayList<>();

        for (String[] allRow : allRows) {
            String peptideDraft = allRow[8];
            StringBuilder peptide = new StringBuilder();
            for (int j = 0; j < peptideDraft.length(); j++) {
                if (Character.isLetter(peptideDraft.charAt(j))) {
                    peptide.append(peptideDraft.charAt(j));
                }
            }
            String protein = allRow[9];
            String evalue = allRow[13];
            Integer posMinus = evalue.indexOf('-');
            if (posMinus == -1) {
                continue;
            }
            Integer magnitude = Integer.parseInt(evalue.substring(posMinus + 1, evalue.length()));

            if (magnitude > 9) {
                peptides.add(new Pair<>(peptide.toString(), protein));
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

        PrintWriter out = new PrintWriter("carboanhydrase_rev.txt");
        out.println(protein.getKey());
        out.println(protein.getValue());

        ArrayList<TreeSet<String>> positions = new ArrayList<TreeSet<String>>();
        for (int i = 0; i < protein.getValue().length(); i++) {
            positions.add(new TreeSet<>((o1, o2) -> o2.length() - o1.length()));
        }

        for (Pair<String, String> pair : peptides) {
            String curPeptide = pair.getKey();
            String curProtein = pair.getValue();

            if (!curProtein.contains(protein.getKey())) {
                continue;
            }

            String pre = Character.toString(getPre(curProtein));
            String post = Character.toString(getPost(curProtein));

            if (pre.equals("-")) {
                positions.get(curProtein.length() - 1).add(curPeptide);
            } else if (post.equals("-")) {
                positions.get(protein.getValue().length() - 1).add(curPeptide);
            } else {
                Integer ind = protein.getValue().indexOf(pre + curPeptide + post);
                positions.get(ind + curPeptide.length()).add(curPeptide);
            }
        }

        for (int i = 0; i < protein.getValue().length(); i++) {
            for (String peptide : positions.get(i)) {
                out.println(getSpaces(i - peptide.length() + 1) + peptide);
            }
        }

        out.close();
    }
}
