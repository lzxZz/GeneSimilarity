


import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;

public class Main {
    public static void main(String ... args){

        try {
//             System.out.println("Start");
// //            ArrayList<TermNode> nodes =  ReadFile.readOboFile("onto.obo");
// //            HashMap<String,TermNode> nodesdict = ReadFile.getNodesDict();
// //            ArrayList<AnnoNode> list =  ReadFile.getGeneNodes("gene.gaf");
//             ReadFile.getCoFunctionNet("net.txt");
           
             Pjj.getSimilarity("GO0006606", "GO0071035");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
