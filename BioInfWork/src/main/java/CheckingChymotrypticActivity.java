import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;

import java.io.*;
import java.util.HashMap;
import java.util.List;

public class CheckingChymotrypticActivity {
    private static final int MAGNITUDE_CONST = 10;

    public static void main(String[] args) throws IOException {
        TsvParserSettings settings = new TsvParserSettings();
        settings.getFormat().setLineSeparator("\n");
        TsvParser parser = new TsvParser(settings);

        // parses all rows in one go.
        List<String[]> allRows = parser.parseAll(new FileReader("../TrypsinNTT2Charge2-5Maxlen50.tsv"));
        List<String[]> allRows2 = parser.parseAll(new FileReader("../ChymotrypsinNTT1Charge2-5_Len50.tsv"));

        PrintStream st = new PrintStream(new FileOutputStream("difference_charge2-5_maxlen50_all.txt"));
        System.setOut(st);

        HashMap<String, String> peptidesFrom1 = new HashMap<>();
        HashMap<String, String> scanEvalueFrom1 = new HashMap<>();
        HashMap<String, String> scanProteinFrom1 = new HashMap<>();

        for (String[] allRow : allRows) {
            Character pre = CarbonicAnhydrase.getPre(allRow[9]);

            peptidesFrom1.put(allRow[1], pre + "|" + CarbonicAnhydrase.removeModifications(allRow[8]));

            scanEvalueFrom1.put(allRow[1], allRow[13]);

            String colProtein = allRow[9];
            String proteinName = colProtein.split("\\(")[0];

            scanProteinFrom1.put(allRow[1], proteinName);
        }

        for (String[] anAllRows2 : allRows2) {
            String curPeptide = CarbonicAnhydrase.removeModifications(anAllRows2[8]);
            Character pre = CarbonicAnhydrase.getPre(anAllRows2[9]);
            String curProtein = anAllRows2[9];
            String curProteinName = curProtein.split("\\(")[0];

            String combination = pre + "|" + curPeptide;

            String evalue = anAllRows2[13];
            Integer posMinus = evalue.indexOf('-');
            if (posMinus == -1) {
                continue;
            }
            Integer magnitude = Integer.parseInt(evalue.substring(posMinus + 1, evalue.length()));

            if (peptidesFrom1.containsKey(anAllRows2[1])) {
                String valuePeptideAndPre = peptidesFrom1.get(anAllRows2[1]);
                if (!valuePeptideAndPre.equals(combination) && magnitude > MAGNITUDE_CONST) {
                    System.out.println(anAllRows2[1] + " " + valuePeptideAndPre +
                            " (" + scanProteinFrom1.get(anAllRows2[1]) + ")" +
                            " (" + scanEvalueFrom1.get(anAllRows2[1]) + ") " + " -> " + combination +
                            " (" + curProteinName + ")" +
                            " (" + evalue + ")");
                } else if (valuePeptideAndPre.equals(combination) && magnitude > MAGNITUDE_CONST) {
                    String evalueFrom1 = scanEvalueFrom1.get(anAllRows2[1]);
                    Integer posMinusFrom1 = evalueFrom1.indexOf('-');
                    if (posMinusFrom1 == -1) {
                        System.out.println("* " + anAllRows2[1] + " " + valuePeptideAndPre +
                                " (" + scanProteinFrom1.get(anAllRows2[1]) + " -> " + curProteinName + ")");
                        continue;
                    }
                    Integer magnitudeFrom1 = Integer.parseInt(evalueFrom1.substring(posMinusFrom1 + 1, evalueFrom1.length()));
                    if (magnitudeFrom1 <= MAGNITUDE_CONST) {
                        System.out.println("* " + anAllRows2[1] + " " + valuePeptideAndPre +
                                " (" + scanProteinFrom1.get(anAllRows2[1]) + " -> " + curProteinName + ")");
                    }
                }
            }
        }

        System.out.println();
        System.out.println("New:");

        for (String[] anAllRows2 : allRows2) {
            String curPeptide = CarbonicAnhydrase.removeModifications(anAllRows2[8]);
            Character pre = CarbonicAnhydrase.getPre(anAllRows2[9]);
            String curProtein = anAllRows2[9];
            String curProteinName = curProtein.split("\\(")[0];

            String combination = pre + "|" + curPeptide;

            String evalue = anAllRows2[13];
            Integer posMinus = evalue.indexOf('-');
            if (posMinus == -1) {
                continue;
            }
            Integer magnitude = Integer.parseInt(evalue.substring(posMinus + 1, evalue.length()));

            if (!peptidesFrom1.containsKey(anAllRows2[1]) && magnitude > MAGNITUDE_CONST) {
                System.out.println(anAllRows2[1] + " " + combination + " (" + curProteinName + ")" +
                        " (" + evalue + ")");
            }
        }

    }
}
