import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;
import com.univocity.parsers.tsv.TsvWriter;
import com.univocity.parsers.tsv.TsvWriterSettings;

import java.io.*;
import java.util.*;

/**
 * Shows peptides which have requested properties.
 */
public class BioInf {
    private static List<Character> trypticAmino = Arrays.asList('R', 'K');
    private static List<Character> chymotrypticAmino = Arrays.asList('L', 'F', 'Y', 'W');

    @SuppressWarnings("unused")
    public Reader getReader(String relativePath) throws UnsupportedEncodingException {

        return new InputStreamReader(this.getClass().getResourceAsStream(relativePath), "UTF-8");
    }

    /**
     * Checks if the peptide has corresponding ends.
     * @param peptide given peptide
     * @param pre given pre-acid
     * @return true if the peptide has corresponding ends
     */
    public static boolean checkRightTrypticLeftChymotryptic(String peptide, Character pre) {
        Set<Character> tryptic = new HashSet<>(trypticAmino);
        Set<Character> chymotryptic = new HashSet<>(chymotrypticAmino);

        if (tryptic.contains(peptide.charAt(peptide.length() - 1)) && chymotryptic.contains(pre)) {
            return true;
        }

        return false;
    }

    /**
     * Checks if the peptide has corresponding ends.
     * @param peptide given peptide
     * @param pre given pre-acid
     * @return true if the peptide has corresponding ends
     */
    public static boolean checkLeftTrypticRightChymotryptic(String peptide, Character pre) {
        Set<Character> tryptic = new HashSet<>(trypticAmino);
        Set<Character> chymotryptic = new HashSet<>(chymotrypticAmino);

        if (chymotryptic.contains(peptide.charAt(peptide.length() - 1)) && tryptic.contains(pre)) {
            return true;
        }

        return false;
    }

    /**
     * Gets pre-amino acid.
     * @param protein given protein
     * @return pre-amino acid of the peptide
     */
    public static Character getPre(String protein) {
        return protein.charAt(protein.indexOf("pre=") + 4);
    }

    public static void main(String[] args) throws FileNotFoundException {
        if (args.length < 2) {
            System.out.println("Not enough arguments.");
            return;
        }

        String fileNameIn = args[0];
        String fileNameOut = args[1];

        TsvParserSettings settings = new TsvParserSettings();
        settings.getFormat().setLineSeparator("\n");
        TsvParser parser = new TsvParser(settings);

        // parses all rows in one go.
        List<String[]> allRows = parser.parseAll(new FileReader(fileNameIn));

        // Writing to an in-memory byte array. This will be printed out to the standard output so you can easily see the result.
        ByteArrayOutputStream tsvResult = new ByteArrayOutputStream();

        // TsvWriter (and all other file writers) work with an instance of java.io.Writer
        Writer outputWriter = new OutputStreamWriter(tsvResult);

        // As with the CsvWriter, all you need is to create an instance of TsvWriter with the default TsvWriterSettings.
        TsvWriter writer = new TsvWriter(outputWriter, new TsvWriterSettings());

        // Write the record headers of this file
        writer.writeHeaders(allRows.get(0));

        List<Object[]> resRows = new ArrayList<>();

        for (int i = 1; i < allRows.size(); i++) {
            if (checkRightTrypticLeftChymotryptic(allRows.get(i)[8], getPre(allRows.get(i)[9]))) {
                resRows.add(allRows.get(i));
            }
        }

        // Here we just tell the writer to write everything and close the given output Writer instance.
        writer.writeRowsAndClose(resRows);

        try(OutputStream outputStream = new FileOutputStream(fileNameOut)) {
            tsvResult.writeTo(outputStream);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
