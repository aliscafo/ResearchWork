import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Counts the number of MS/MS spectrum.
 */
public class CountMSMSScans {
    public static void main(String[] args) throws IOException {
        String experiment;

        if (args.length == 1) {
            experiment = args[0];
        } else {
            experiment = "../140509QXc1_car_anh_tryp_001.mzXML";
        }

        String content = readFile(experiment);

        int scanNum = 0;
        int lastIndex = 0;
        String findStr = "msLevel=\"2\"";

        while(lastIndex != -1){

            lastIndex = content.indexOf(findStr, lastIndex);

            if(lastIndex != -1){
                scanNum++;
                lastIndex += findStr.length();
            }
        }

        System.out.println(scanNum);
    }

    private static String readFile(String path) throws IOException {
        byte[] encoded = Files.readAllBytes(Paths.get(path));
        return new String(encoded);
    }
}
