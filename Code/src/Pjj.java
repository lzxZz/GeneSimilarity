import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


//彭佳杰博士论文实现，目前只实现到计算功能距离
public class Pjj{

    private static ArrayList<AnnoNode> annoNodes;
    private static ArrayList<TermNode> ontoDict;
    private static ArrayList<NetEdge> netEdges;
    private static HashMap<String,Double> netHashMap;

    //数据初始化操作
    public static void init()throws IOException{
        annoNodes  = ReadFile.getGeneNodes("gene.gaf");
        ontoDict = ReadFile.readOboFile("onto.obo");
        netEdges = ReadFile.getCoFunctionNet("net.txt");
        netHashMap = new HashMap<String,Double>();
        
        //由于边有方向，为了简化操作，使用key为name1-name2，value为weight/10这样的键值对存储基因功能网络，并且提前做好归一化处理
        for (NetEdge edge : netEdges) {
            String name = edge.Gene1 + "-" + edge.Gene2;
            if (! netHashMap.containsKey(name)){
                netHashMap.put(name, edge.Weight/10);//做归一化处理
            }
        }
    }

    //return the D(go1,go2) 公式2-1
    public static double getSimilarity(String go1,String go2)throws IOException{
        init();
        //根据给定的基因本体，获取对应的基因集合
        ArrayList<String> geneSet1 = new ArrayList<String>();
        ArrayList<String> geneSet2 = new ArrayList<String>();
        //遍历所有注释节点，将指定的GO术语对应的基因添加到集合中去
        for(AnnoNode node: annoNodes){
            if (node.ID.equals(go1)){
                geneSet1.add(node.GeneName);   
                geneSet1.addAll(node.Synonym);
            }
            if (node.ID.equals(go2)){
                geneSet2.add(node.GeneName);   
                geneSet2.addAll(node.Synonym);
            }
        }
        //去除重复元素
        HashSet set = new HashSet(geneSet1);
        geneSet1.clear();
        geneSet1.addAll(set);
        geneSet1.trimToSize();//去除空值

        set = new HashSet(geneSet2);
        geneSet2.clear();
        geneSet2.addAll(set);
        geneSet2.trimToSize();//去除空值
        
        //数据输出到控制台，用于观察
        // for (String var : geneSet1) {
        //     System.out.println("Go Term " + go1 + ":" + var);
        // }

        //计算公式 L(set1,set2) + L(set2,set1)/(2(num of set1 + num of set2) - L(set1,set2) - L(set2,set1))
        double value = 0.0;
        double l12 = doubleSetAB(geneSet1, geneSet2);
        double l21 = doubleSetAB(geneSet2, geneSet1);

        value = (l12 + l21) / (2 * (geneSet1.size() + geneSet2.size()) - l12 - l21);

        return value;
    }

    public static double doubleSetAB(ArrayList<String> set1,ArrayList<String> set2){
        double sum = 0.0;
        for (String gene1:set1){
            for (String gene2:set2){
                String key = gene1 + "-" + gene2;
                if (gene1.equals(gene2)){
                    sum += 1;
                }else if (netHashMap.containsKey(key)){
                    sum += (1 -  netHashMap.get(key));
                }
            }
        }

        return sum;
    }

}