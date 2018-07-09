import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

 public  class ReadFile {


   static private HashMap<String,TermNode> NodesDict = new HashMap<String,TermNode>();

   public static ArrayList<TermNode> readOboFile(String oboFile) throws IOException {
       ArrayList<TermNode> AllNode  =   new ArrayList<TermNode>();
       String content               =   new String(Files .readAllBytes(Paths.get(oboFile)));

        //正则匹配每一个术语及其属性
        Pattern termPattern         =   Pattern.compile("\\[Term\\][\\s\\S.+?]+?\\n\\n");
        Pattern idPattern           =   Pattern.compile("(?<=id: )GO:\\d{7}");
        Pattern namePattern         =   Pattern.compile("(?<=name: ).+?\\n");
        Pattern namespacePattern    =   Pattern.compile("(?<=namespace: ).+?\\n");
        Pattern defPattern          =   Pattern.compile("(?<=def: ).+?\\n");
        Pattern isaPattern          =   Pattern.compile("(?<=is_a: )GO:\\d{7}");
        Pattern partofPattern       =   Pattern.compile("(?<=relationship: part_of )GO:\\d{7}");
        Pattern obsoletePattern     =   Pattern.compile("is_obsolete: true");

        Matcher m = termPattern.matcher(content);
        while (m.find()){

            TermNode    node        = new TermNode();
            String      term        = m.group();
            String      id          = "";
            String      name        = "";
            String      namespace   = "";
            String      def         = "";
            String      isa         = "";
            String      partof      = "";
            boolean     isObsolete  = false;
            System.out.println(term);

            Matcher     idm         = idPattern.matcher(term);
            if(idm.find()){
                System.out.println("ID : " + idm.group().replace(":",""));
                id = idm.group().replace(":","");
                node.ID = id;
            }

            Matcher namem = namePattern.matcher(term);
            if (namem.find()){
                System.out.println("name : " + namem.group());
                name = namem.group();
                node.name = name;
            }

            Matcher nsm = namespacePattern.matcher(term);
            if (nsm.find()){
                System.out.println("NameSpace : " + nsm.group());
                namespace  = nsm.group();
                node.namespace =  namespace;
            }

            Matcher defm  = defPattern.matcher(term);
            if (defm.find()){
                System.out.println("Def : " + defm.group());
                def = defm.group();
                node.def = def;
            }

            Matcher isam = isaPattern.matcher(term);
            while (isam.find()) {
                System.out.println("Is_a Parent ID : " + isam.group().replace(":",""));
                isa = isam.group().replace(":","");
                node.IParentIDs.add(isa);
            }


            Matcher partofm = partofPattern.matcher(term);
            while (partofm.find()){
                System.out.println("Part_of Parent ID : " + partofm.group());
                partof = partofm.group();
                node.PParentIDs.add(partof);
            }

            if (obsoletePattern.matcher(term).find()){
                isObsolete = true;
                node.isObsolete = isObsolete;
            }

            AllNode.add(node);
            NodesDict.put(node.ID,node);
        }

        return AllNode;
    }
    static HashMap<String,TermNode> getNodesDict(){
        return NodesDict;
     }


   static ArrayList<AnnoNode> getGeneNodes(String gafFile) throws IOException {
       ArrayList<AnnoNode> genes = new ArrayList<AnnoNode>();
       for (String line : Files.readAllLines(Paths.get(gafFile))){
           String[] ss = line.split("\t");
//           System.out.println(ss.length);
           if (ss.length > 11) {
//                System.out.println("Gene ID : " + ss[2]  + "\t\t Go ID : " + ss[4]);
//               System.out.println(line.replace("\t","---"));

               String geneSymbol = "";//col3
               String goID = "";//col5
               ArrayList<String>  refID = new ArrayList<String>();//col6
               String evidCode = "";//col7
               String nameSpace = "";//col9  PFC
               ArrayList<String>  synonym = new ArrayList<String>();//col11
               String type = "";//col12 默认为gene

               geneSymbol   = ss[2];
               goID         = ss[4].replace(":", "");
               evidCode     = ss[6];
               nameSpace    = ss[8];
               type         = ss[11];
               refID.addAll(Arrays.asList(ss[5].split("\\|")));
               synonym.addAll(Arrays.asList(ss[10].split("\\|")));
               synonym.remove(synonym.size()-1);


                   AnnoNode node  = new AnnoNode();
                   node.GeneName = geneSymbol;
                   node.ID  = goID;
                   node.Synonym.addAll(synonym);
                   node.EvidenceCode = evidCode;
                   node.NameSpace = nameSpace;
                   node.RefIDS.addAll(refID);
                   node.Type = type;

                    genes.add(node);

           }
       }

       return  genes;
   }

    static  ArrayList<NetEdge> getCoFunctionNet(String netFile) throws IOException {
       ArrayList<NetEdge> edges = new ArrayList<NetEdge>();

       for (String line: Files.readAllLines(Paths.get(netFile))){
           String[] ss = line.split("\t");
           NetEdge edge = new NetEdge();
           edge.Gene1 = ss[0];
           edge.Gene2 = ss[1];
           edge.Weight = Double.parseDouble(ss[2]);
           edges.add(edge);
       }


       return edges;
    }

}
