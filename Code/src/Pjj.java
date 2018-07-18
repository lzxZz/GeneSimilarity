import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

//彭佳杰博士论文实现，目前只实现到计算功能距离
public class Pjj {

    private static ArrayList<AnnoNode> annoNodes;
    private static ArrayList<TermNode> ontoNodes;
    private static ArrayList<NetEdge> netEdges;
    private static HashMap<String, TermNode> ontoDict;
    private static HashMap<String, Double> netHashMap;
    private static HashMap<String, AnnoNode> annoDict;
    private static HashMap<String, PTermNode> pontoDict;
    private static final String Log = "D:/pjjlog.txt";

    // 数据初始化操作
    public static void init() throws IOException {
        annoNodes = ReadFile.getGeneNodes("gene.gaf");
        System.out.println("read gene finished");
        ontoNodes = ReadFile.readOboFile("onto.obo");
        System.out.println("read onto finished");
        ontoDict = ReadFile.getNodesDict();
        netEdges = ReadFile.getCoFunctionNet("net.txt");
        System.out.println("read net finished");

        netHashMap = new HashMap<String, Double>();
        annoDict = new HashMap<>();
        for (AnnoNode node : annoNodes) {
            annoDict.put(node.ID, node);
        }

        pontoDict = new HashMap<>();
        for (TermNode node : ontoNodes) {
            PTermNode pnode = new PTermNode();
            pnode.ID = node.ID;
            pnode.def = node.def;
            pnode.IParentIDs = node.IParentIDs;
            pnode.PParentIDs = node.PParentIDs;
            pnode.IParentNodes = node.IParentNodes;
            pnode.PParentNodes = node.PParentNodes;
            pnode.isObsolete = node.isObsolete;
            ;
            pnode.name = node.name;
            pnode.namespace = node.namespace;

            pontoDict.put(node.ID, pnode);
        }
        for (TermNode node : ontoNodes) {
            for (String id : node.PParentIDs) {

                pontoDict.get(node.ID).childNodes.add(pontoDict.get(id));
            }

            for (String id : node.IParentIDs) {
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


    // 返回-200表示两个术语非同分支
    public static double getSimilarity(String ta, String tb) throws IOException {
        // 获取G
        //init();
        int countG = 0;
        if (annoDict.containsKey(ta) && annoDict.containsKey(tb)) {
            if ((!annoDict.get(ta).NameSpace.equals(annoDict.get(tb).NameSpace))) {
                return -200;
            } else {
                String ns = annoDict.get(ta).NameSpace;
                for (AnnoNode node : annoNodes) {
                    if (ns.equals(node.NameSpace)) {
                        countG++;
                    }

                }
            }
        }

        // 获取Ga，Gb，P
        ArrayList<String> Ga = getGeneSetByGoId(ta);
        ArrayList<String> Gb = getGeneSetByGoId(tb);

        if (Ga.size()==0 || Gb.size()==0){
            System.out.println("该术语对没有注释基因");
            return -100;
        }


        ArrayList<String> PNodeIds = getPublicParent(ta, tb);

        HashSet set = new HashSet(PNodeIds);
        PNodeIds.clear();
        PNodeIds.addAll(set);


        double similarityValue = 0.0;

        if (PNodeIds.size() == 0) {

            System.out.println("没有公共祖先节点");
            return -300;
        }
        double d = getSimilarity2_1(ta, tb);

        ArrayList<String> uabpGenes = new ArrayList<>();
        for (String pid : PNodeIds) {
            // u
            uabpGenes.clear();
            uabpGenes = getUabp(pontoDict.get(ta), pontoDict.get(tb), pontoDict.get(pid));
            // f
            double f = d * d * uabpGenes.size() + (1 - d * d) * Math.sqrt(Ga.size() * Gb.size());
            // h
            double h = d * d * countG + (1 - d * d) * Math.max(Ga.size(), Gb.size());
            // s


            double sim = 0.0;
            sim = (2 * Math.log(countG) - 2 * Math.log(f))
                    / (2 * Math.log(countG) - Math.log(Ga.size()) - Math.log(Gb.size()))
                    * (1 - h / countG * PNodeIds.size() / countG);
           // System.out.println(similarityValue);
            similarityValue = sim > similarityValue ? sim : similarityValue;
        }

        return similarityValue;
    }

    // 计算U(ta,tb,p)的注释基因集合
    public static ArrayList<String> getUabp(PTermNode ta, PTermNode tb, PTermNode tp) {
        ArrayList<PTermNode> nodes = new ArrayList<>();
        nodes.addAll(getPathWayNode(tp, ta));
        nodes.addAll(getPathWayNode(tp, tb));
        nodes.add(tp);
        nodes.addAll(getAllAnnoNodeByNode(ta));
        nodes.addAll(getAllAnnoNodeByNode(tb));

        ArrayList<String> genes = new ArrayList<>();
        for (PTermNode node : nodes) {
            if (node != null) {
                AnnoNode anno = annoDict.get(node.ID);
                if (anno != null) {
                    genes.add(anno.GeneName);
                    genes.addAll(anno.Synonym);
                }
            }
        }

        // 去重复
        HashSet set = new HashSet<>(genes);
        genes.clear();
        genes.addAll(set);
        genes.trimToSize();

        return genes;
    }

    private static ArrayList<PTermNode> getAllAnnoNodeByNode(PTermNode node) {
        ArrayList<PTermNode> nodes = new ArrayList<>();
        if (node.childNodes != null) {
            nodes.addAll(node.childNodes);
            for (PTermNode tNode : node.childNodes) {
                if (tNode != null)
                    nodes.addAll(getAllAnnoNodeByNode(tNode));
            }
        } else {
            return null;
        }
        return nodes;
    }

    public static ArrayList<TermNode> getParentNode(String id) {
        ArrayList<TermNode> parents = new ArrayList<>();
        if (ontoDict.containsKey(id)) {
            TermNode node = ontoDict.get(id);
            for (String pid : node.PParentIDs) {
                parents.add(ontoDict.get(pid));
            }

            for (String innerid : node.PParentIDs) {
                parents.addAll(getParentNode(innerid));
            }

            for (String pid : node.IParentIDs) {
                parents.add(ontoDict.get(pid));
            }

            for (String innerid : node.IParentIDs) {
                parents.addAll(getParentNode(innerid));
            }
        }

        return parents;
    }

    public static ArrayList<String> getPublicParent(String ta, String tb) {
        ArrayList<TermNode> set1 = getParentNode(ta);
        ArrayList<TermNode> set2 = getParentNode(tb);

        ArrayList<String> publicNode = new ArrayList<>();
        try {
            for (TermNode node : set1) {
                if (node != null) {


                    if (set2.contains(node)) {
                        publicNode.add(node.ID);
                    }
                }
            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }


        return publicNode;
    }

    // 获取节点之间所有路径的节点
    private static ArrayList<PTermNode> getPathWayNode(PTermNode pnode, PTermNode cnode) {
        ArrayList<PTermNode> nodes = new ArrayList<>();
        for (PTermNode node : pnode.childNodes) {
            if (node !=null) {
                if (hasChildNode(node, cnode)) {
                    nodes.add(node);
                    nodes.addAll(getPathWayNode(node, cnode));
                }
            }
        }
        return nodes;
    }

    // 递归算法，递归次数会很多，严重影响效率
    private static boolean hasChildNode(PTermNode pnode, PTermNode cnode) {

        if (pnode.childNodes.contains(cnode)) {
            return true;
        }
        if (pnode.childNodes.size() == 0) {
            return false;
        } else {
            for (PTermNode node : pnode.childNodes) {
                if (node != null) {
                    return hasChildNode(node, cnode);
                }else{
                    return false;
                }
            }
        }

        return false;
    }

    // return the D(go1,go2) 公式2-1
    public static double getSimilarity2_1(String go1, String go2) throws IOException {

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

        if (geneSet1.size() == 0 || geneSet2.size() == 0) {
            System.out.println("该术语对没有注释基因");
        }
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
            double value = 1;
            for (String gene2 : set2) {
                String key = gene1 + "-" + gene2;
//                try {
                    if (gene1.equals(gene2)) {
                        value *= 0;
//                        Files.write(Paths.get(Log), (key + "\t---\t0\r\n").getBytes(), StandardOpenOption.APPEND);
                    } else if (netHashMap.containsKey(key)) {
                         value *= (1 - netHashMap.get(key));


//                        Files.write(Paths.get(Log), (key + "\t---\t" + (1 - netHashMap.get(key)) + "\r\n").getBytes(), StandardOpenOption.APPEND);

                    } else {

//                        Files.write(Paths.get(Log), (key + "\t---\t1\r\n").getBytes(), StandardOpenOption.APPEND);
                    }
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
            }
            sum += value;
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