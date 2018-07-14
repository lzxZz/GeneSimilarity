import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

//彭佳杰博士论文实现，目前只实现到计算功能距离
public class Pjj {

    private static ArrayList<AnnoNode> annoNodes;
    private static ArrayList<TermNode> ontoNodes;
    private static ArrayList<NetEdge> netEdges;
    private static HashMap<String,TermNode> ontoDict;
    private static HashMap<String, Double> netHashMap;
    private static HashMap<String,AnnoNode> annoDict;
    private static HashMap<String,PTermNode> pontoDict;
    // 数据初始化操作
    public static void init() throws IOException {
        annoNodes = ReadFile.getGeneNodes("gene.gaf");
        ontoNodes = ReadFile.readOboFile("onto.obo");
        ontoDict = ReadFile.getNodesDict();
        netEdges = ReadFile.getCoFunctionNet("net.txt");
        
        netHashMap = new HashMap<String, Double>();
        annoDict = new HashMap<>();
        for (AnnoNode node:annoNodes){
            annoDict.put(node.ID, node);
        }

        pontoDict = new HashMap<>();
        for (TermNode node:ontoNodes){
            PTermNode pnode = new PTermNode();
            pnode.ID = node.ID;
            pnode.def =  node.def;
            pnode.IParentIDs =  node.IParentIDs;
            pnode.PParentIDs = node.PParentIDs;
            pnode.IParentNodes = node.IParentNodes;
            pnode.PParentNodes = node.PParentNodes;
            pnode.isObsolete =  node.isObsolete;;
            pnode.name = node.name;
            pnode.namespace = node.namespace;
            
            
            
            pontoDict.put(node.ID, pnode);
        }
        for (TermNode node:ontoNodes){
            for (String id:node.PParentIDs){

                pontoDict.get(node.ID).childNodes.add(pontoDict.get(id));
            }

            for (String id:node.IParentIDs){
                pontoDict.get(node.ID).childNodes.add(pontoDict.get(id));
            }
        }

        // 由于边有方向，为了简化操作，使用key为name1-name2，value为weight/10这样的键值对存储基因功能网络，并且提前做好归一化处理
        for (NetEdge edge : netEdges) {
            String name = edge.Gene1 + "-" + edge.Gene2;
            if (!netHashMap.containsKey(name)) {
                netHashMap.put(name, edge.Weight / 10);// 做归一化处理
            }
        }

    }

    //返回-200表示两个术语非同分支
    public static double getSimilarity(String ta, String tb) throws IOException {
        //获取G
        init();
        int countG = 0;
        if (annoDict.containsKey(ta) && annoDict.containsKey(tb)){
            if ((! annoDict.get(ta).NameSpace.equals(  annoDict.get(tb).NameSpace))){
                return -200;
            }else{
                String ns = annoDict.get(ta).NameSpace;
                for (AnnoNode node:annoNodes){
                    if (ns.equals(node.NameSpace)){
                        countG++;
                    }

                }
            }
        }

        // 获取Ga，Gb，P
        ArrayList<String> Ga = getGeneSetByGoId(ta);
        ArrayList<String> Gb = getGeneSetByGoId(tb);

        double d = getSimilarity2_1(ta, tb);

        ArrayList<String> P = getPublicParent(ta, tb);

        double similarityValue = 0.0;
        for (String p : P){

        }

        return similarityValue;
    }

    public static ArrayList<TermNode> getParentNode(String id) {
        ArrayList<TermNode> parents = new ArrayList<>();
        if (ontoDict.containsKey(id)) {
            TermNode node = ontoDict.get(id);
            for (String pid:node.PParentIDs){
                parents.add(ontoDict.get(pid));
            }
            
            for(String innerid: node.PParentIDs){
                parents.addAll(getParentNode(innerid));
            }

            for (String pid:node.IParentIDs){
                parents.add(ontoDict.get(pid));
            }

            for(String innerid: node.IParentIDs){
                parents.addAll(getParentNode(innerid));
            }
        }
        
        return parents;
    }

    public static ArrayList<String> getPublicParent(String ta, String tb) {
        ArrayList<TermNode> set1 = getParentNode(ta);
        ArrayList<TermNode> set2 = getParentNode(tb);

        ArrayList<String> publicNode = new ArrayList<>();

        for (TermNode node: set1){
            if (set2.contains(node)){
                publicNode.add(node.ID);
            }
        }

        return publicNode;
    }

    // return the D(go1,go2) 公式2-1
    public static double getSimilarity2_1(String go1, String go2) throws IOException {
        init();
        // 根据给定的基因本体，获取对应的基因集合
        ArrayList<String> geneSet1 = getGeneSetByGoId(go1);
        ArrayList<String> geneSet2 = getGeneSetByGoId(go2);

        // 数据输出到控制台，用于观察
        // for (String var : geneSet1) {
        // System.out.println("Go Term " + go1 + ":" + var);
        // }

        // 计算公式 L(set1,set2) + L(set2,set1)/(2(num of set1 + num of set2) - L(set1,set2)
        // - L(set2,set1))
        double value = 0.0;
        double l12 = doubleSetAB(geneSet1, geneSet2);
        double l21 = doubleSetAB(geneSet2, geneSet1);

        value = (l12 + l21) / (2 * (geneSet1.size() + geneSet2.size()) - l12 - l21);

        return value;
    }

    public static ArrayList<String> getGeneSet2_3(String termA, String termB, String termP) {
        ArrayList<String> result = new ArrayList<>();

        result.addAll(getGeneSetByGoId(termA));
        result.addAll(getGeneSetByGoId(termB));
        result.addAll(getGeneSetByGoId(termP));

        // ta到p上路径所有术语的注释
        // tb到p上路径所有术语的注释

        return result;
    }

    // 计算AB两集合的功能距离
    public static double doubleSetAB(ArrayList<String> set1, ArrayList<String> set2) {
        double sum = 0.0;
        for (String gene1 : set1) {
            for (String gene2 : set2) {
                String key = gene1 + "-" + gene2;
                if (gene1.equals(gene2)) {
                    sum += 1;
                } else if (netHashMap.containsKey(key)) {
                    sum += (1 - netHashMap.get(key));
                }
            }
        }

        return sum;
    }

    private static ArrayList<String> getGeneSetByGoId(String id) {
        ArrayList<String> genes = new ArrayList<>();
        for (AnnoNode node : annoNodes) {
            if (id.equals(node.ID)) {
                genes.add(node.GeneName);
                genes.addAll(node.Synonym);
            }
        }
        // 去除重复元素
        HashSet set = new HashSet(genes);
        genes.clear();
        genes.addAll(set);
        genes.trimToSize();// 去除空值

        return genes;

    }
}