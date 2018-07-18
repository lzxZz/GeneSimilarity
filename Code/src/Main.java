


import com.sun.org.apache.bcel.internal.generic.Select;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Main {
    public static void main(String ... args){

        try {
//             System.out.println("Start");
           ArrayList<TermNode> nodes =  ReadFile.readOboFile("onto.obo");
// //            HashMap<String,TermNode> nodesdict = ReadFile.getNodesDict();
            ArrayList<AnnoNode> list =  ReadFile.getGeneNodes("gene.gaf");
             ReadFile.getCoFunctionNet("net.txt");

            HashSet<String> bioprocess = new HashSet<>();
            HashSet<String> cellcomponet = new HashSet<>();
            HashSet<String> ff = new HashSet<>();
            for (AnnoNode node:list)
            {
                switch (node.NameSpace){
                    case "P":
                        bioprocess.add(node.ID);
                        break;
                    case "F":
                        ff.add(node.ID);
                        break;
                    case "C":
                        cellcomponet.add(node.ID);
                        break;
                }
            }


             Pjj.init();

            CalcThread th = new CalcThread();
            th.start(cellcomponet);

//             for (TermNode i: nodes){
//                 for (TermNode j: nodes){
////                     System.out.println(i.ID + "__" +j.ID + "\t---\t"+ Pjj.getSimilarity(i.ID,j.ID));
//                     System.out.println(i.ID + "\t" +j.ID);
//                     Files.write(Paths.get("D:/pjjlog.txt"),(i.ID + "\t" +j.ID + "\t" + Pjj.getSimilarity(i.ID,j.ID)+"\n").getBytes(),StandardOpenOption.APPEND);
//                 }
//             }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static class CalcThread extends Thread{
        HashSet<String> ids;
        public void run(){
            int i =0 ;
            StringBuffer sb ;
            for (String id1:ids){
                sb= new StringBuffer();
                try {


                    for (String id2 : ids) {

                        double value = Pjj.getSimilarity(id1, id2);


                        sb.append(id1 + "\t" + id2 + "\t" + value + "\n");

                    }
                    Files.write(Paths.get("D:/pjjcc.txt"), (sb.toString()).getBytes(), StandardOpenOption.APPEND);
                    System.out.println(++i);
                }catch (IOException e)
                {

                }
            }
        }

        public void start(HashSet set){
            ids  = set;
            run();
        }
    }
}

