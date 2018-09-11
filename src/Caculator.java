import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

import javax.xml.crypto.dom.DOMCryptoContext;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.lang.Math.*;

import static java.lang.Math.log;

public class Caculator {


    //数据结构
    static HashMap<String, Term> idTermHashMap;       //id到节点的hash
    static HashMap<String, ArrayList<String>> idParentsHashMap;    //id到所有祖先节点ID
    static HashMap<String, ArrayList<String>> idChildesHashMap;    //id到所有的子孙节点ID
    static HashMap<String, ArrayList<String>> idSonHashMap;        //id到儿子节点的ID
    static ArrayList<Anno> AnnoList;           //基因ID到注释节点的
    static HashMap<String, ArrayList<Anno>> idAnnoHashmap;      //GoID到注释条目
    static HashMap<String, ArrayList<String>> geneTermsHashMap;   //基因到术语集合
    static HashMap<String, ArrayList<String>> idGenesHashMap;     //术语到基因集合
    static HashMap<String, Double> netHashMap;         //gene1:gene2到边权重，（已经归一化）
    static HashMap<String,HashSet<String>> ecnumGeneSetHashMap;      //EC编号到基因集合

    static boolean isInit;

    //计算两个基因的相似度
    public static double getGeneSimilarity(String gene1, String gene2) throws Exception {
        try {

                data_initializer();



        } catch (IOException e) {
            System.out.println("数据初始化错误");
            System.out.println(e.getMessage());
            e.printStackTrace();
            return -1000;
        }
        //获取两个基因对应的术语集合
        ArrayList<String> terms1 = geneTermsHashMap.get(gene1);
        ArrayList<String> terms2 = geneTermsHashMap.get(gene2);

        if (terms1 == null || terms2 == null) {
            return 0;
        }

        double geneSimilarity = 0.0;
        double termSimilaritySum1 = 0.0;
        for (String id : terms1) {
            double maxvalue = getMaxSimilarityBetweenGoAndGoSet(id, terms2, gene1, gene2);
            termSimilaritySum1 += maxvalue;
            //System.out.println("MaxSimilarity:" + maxvalue);
        }

        double termSimilaritySum2 = 0.0;
        for (String id : terms2) {
            termSimilaritySum2 += getMaxSimilarityBetweenGoAndGoSet(id, terms1, gene1, gene2);
        }

        geneSimilarity = (termSimilaritySum1 + termSimilaritySum2) / (terms1.size() + terms2.size());


        //获取注释基因的术语集合
        return geneSimilarity;
    }

    //计算一个术语到另一个术语集合最大的相似度,忽略gene1和gene2
    private static double getMaxSimilarityBetweenGoAndGoSet(String goId, ArrayList<String> goset, String gene1, String gene2)throws Exception {
        double maxSimilarity = 0.0;
        for (String id : goset) {
            double tmpSimilarity = getTermSimilarutyIgnorGene(goId, id, gene1, gene2);
            maxSimilarity = tmpSimilarity > maxSimilarity ? tmpSimilarity : maxSimilarity;
        }

        return maxSimilarity;
    }

    //根据给定的基因名称，获取注释术语集合
    private ArrayList<String> getTermAnnoListFromGene(String gene) {
        return null;
    }

    //计算两个术语间的相似度，忽略指定的基因
    //**其中所有的变量名都来自于论文中伪代码**
    private static double getTermSimilarutyIgnorGene(String term1, String term2, String gene1, String gene2)  throws Exception{
        if (idTermHashMap.get(term1) == null || idTermHashMap.get(term2) == null) {
            return 0.0;
        }
        if (!idTermHashMap.get(term1).namespace.equals(idTermHashMap.get(term2).namespace)) {
            return 0.0;
        }

        ArrayList<String> P = new ArrayList<>();
        for (String id : idParentsHashMap.get(term1)) {
            if (idParentsHashMap.get(term2).contains(id)) {
                P.add(id);
            }
        }

        ArrayList<String> Ga = idGenesHashMap.get(term1);
        ArrayList<String> Gb = idGenesHashMap.get(term2);

        int countG = 1; //论文中只使用到了该分支注释的所有基因的数量
        // for (Anno node : AnnoList) {
        //     if (node.nameSpace.equals(idAnnoHashmap.get(term1).get(0).nameSpace)) {
        //         countG++;
        //     }
        // }
        // double Dab = get2_1(term1, term2, gene1, gene2);
        double Dab = 1.1;
        //System.out.println(Dab);
        double Dab2 = Dab * Dab; //论文中大量使用Dab的平方，这里做预处理
        double maxSimilarity = 0.0;
        P =  arrayDistinct(P);
        for (String id : P) {
            int U = getUabp(term1, term2, id, gene1, gene2);

            // double f = Dab2 * U + (1 - Dab2) * Math.sqrt(Ga.size() * Gb.size());

            // int maxGab = 0;
            // maxGab = Ga.size() > Gb.size() ? Ga.size() : Gb.size();//公式2-5计算h时使用到了Ga与Gb数量的最大值

            // double h = Dab2 * countG + (1 - Dab2) * maxGab;

            // int Gp = 0;
            // ArrayList<String> pgenes = new ArrayList<>();
            // if (idChildesHashMap.get(id) != null) {
            //     for (String childid : idChildesHashMap.get(id)) {
            //         ArrayList<String> tmppgenes = idGenesHashMap.get(childid);

            //         if (tmppgenes != null) {
            //             pgenes.addAll(tmppgenes);
            //         }
            //     }
            //     if (pgenes != null) {
            //         arrayDistinct(pgenes);
            //         Gp = pgenes.size();
            //     }
            // }


            // double similarity = ((2 * log(countG) - 2 * log(f)) / (2 * log(countG) - (log(Ga.size()) + log(Gb.size())))) * (1 - ((h * Gp) / (countG * countG)));

            //maxSimilarity = maxSimilarity > similarity ? maxSimilarity : similarity;
        }

        maxSimilarity = 1;
        return maxSimilarity;
    }

    static StringBuilder sb = new StringBuilder();
    static int cccc = 0;
    //获取UABP注释基因的数量
    //**论文中只使用了数量**
    private static int getUabp(String term1, String term2, String parent, String gene1, String gene2) throws Exception {
        HashSet<String> genes = new HashSet<>();
        //分别获取，a,b,p对应的注释基因。
        ArrayList<String> agenes = idGenesHashMap.get(term1);
        if (agenes != null) {
            genes.addAll(agenes);
        }
        ArrayList<String> bgenes = idGenesHashMap.get(term2);
        if (bgenes != null) {
            genes.addAll(bgenes);
        }


        ArrayList<String> pgenes = idGenesHashMap.get(parent);
        if (pgenes != null) {
            genes.addAll(pgenes);
        }


        // System.out.println("要计算的路径为：" + term1 + "---" + parent);
        // System.out.println("要计算的路径为：" + term2 + "---" + parent);
        
        sb.append(term1 + ":" + parent + "\n");
        sb.append(term2 + ":" + parent + "\n");
        cccc++;
        System.out.println("计算出第" + cccc + "次") ;
        if (cccc > 1000)
        {
            Files.write(Paths.get("path.buf"), sb.toString().getBytes(), StandardOpenOption.APPEND);

            
            sb.delete(0, sb.length());
        }



        // //获取a,p路径上的术语注释基因
        // //ArrayList<String> pathAP = getPathWayTerm(term1, parent);
        // ArrayList<String> pathAP = new ArrayList<>();
        // for (String node : pathAP) {
        //     try {
        //         ArrayList<String> tmpgenes = idGenesHashMap.get(node);
        //         if (tmpgenes != null) {
        //             genes.addAll(tmpgenes);
        //         }
        //     } catch (NullPointerException e) {
        //         System.out.println(e.getMessage());
        //     }
        // }


        // //获取b，p路径上术语的注释基因
        // ArrayList<String> pathBP = getPathWayTerm(term1, parent);
        // for (String node : pathBP) {
        //     ArrayList<String> tmpgenes = idGenesHashMap.get(node);
        //     if (tmpgenes != null) {
        //         genes.addAll(tmpgenes);
        //     }
        // }
        // genes.remove(gene1);
        // genes.remove(gene2);
        return 10;
    }

    private static ArrayList<String> getPathWayTerm(String child, String parent) {
        ArrayList<String> nodes = new ArrayList<>();
        if (child.equals(parent)) {
            nodes.add(child);
            return nodes;
        }
        if (idSonHashMap.get(parent) != null) {
            for (String node : idSonHashMap.get(parent)) {
                ArrayList<String> childes = idChildesHashMap.get(node);
                if (childes != null) {
                    if (childes.contains(child)) {
                        nodes.add(node);
                        nodes.addAll(getPathWayTerm(child, node));
                    }
                }
            }
        }

        return nodes;
    }


    private static double get2_1(String term1, String term2, String gene1, String gene2) {
        //获取注释基因集合
        ArrayList<String> set1 = new ArrayList<>(idGenesHashMap.get(term1));
        ArrayList<String> set2 = new ArrayList<>(idGenesHashMap.get(term2));
        //忽略要计算的基因对
        set1.remove(gene1);
        set1.remove(gene2);
        set2.remove(gene1);
        set2.remove(gene2);
        //去重复
        set1 = arrayDistinct(set1);
        set2 = arrayDistinct(set2);

        double l12 = geneSetDistance(set1, set2);
        double l21 = geneSetDistance(set2, set1);

       // System.out.println("L12 : " +  l12 + "\tL21:" + l21);
        double value = 0.0;
        value = (l12 + l21) / (2 * (set1.size() + set2.size()) - l12 - l21);
        return value;
    }

    //公式2中累乘求和
    private static double geneSetDistance(ArrayList<String> set1, ArrayList<String> set2) {
        double sum = 0.0;
        for (String gene1 : set1) {
            double value = 1;
            for (String gene2 : set2) {
                //System.out.println("distance"  + value);
                String key = gene1 + ":" + gene2;
                if (gene1.equals(gene2)) {
                    value *= 0;
                } else {
                    if (netHashMap.containsKey(key)) {
                        value *= (1 - netHashMap.get(key));
                    } else {
                        value *= 1;
                    }
                }
            }
            sum += value;
        }

        return sum;
    }

    private static ArrayList<String> arrayDistinct(ArrayList<String> arr) {
        HashSet<String> set = new HashSet<>(arr);
        arr.clear();
        arr.addAll(set);
        return new ArrayList<>(set);
    }

    public static void data_initializer() throws IOException {
        if (isInit){
            return;
        }

        //读取术语列表，转化为hashmap,便于操作
        idTermHashMap = new HashMap<>();
        for (Term term : Reader.readOboFile("onto.obo")) {
            //忽略过时术语
            if (!term.isObsolote) {
                idTermHashMap.put(term.id, term);
            }

        }


        //id到所有祖先节点
        idParentsHashMap = new HashMap<>();
        for (String id : idTermHashMap.keySet()) {
            ArrayList<String> parents = new ArrayList<>();
            idParentsHashMap.put(id, parents);

            //已经添加的术语
            HashSet<String> added = new HashSet<>();
            //等待添加的术语
            ArrayList<String> adding = new ArrayList<>();
            adding.addAll(idTermHashMap.get(id).partId);
            adding.addAll(idTermHashMap.get(id).isId);
            while (adding.size() > 0) {

                try {
                    //判断是否添加过这一节点
                    if (!added.contains(adding.get(0))) {
                        parents.add(adding.get(0));
                        Term t = idTermHashMap.get(adding.get(0));
                        if (t.partId != null && t.partId.size() > 0) {
                            adding.addAll(t.partId);
                        }
                        if (t.isId != null && t.isId.size() > 0) {
                            adding.addAll(t.isId);
                        }
                    }
                    adding.remove(0);
                } catch (NullPointerException e) {
                    System.out.println(adding.get(0));
                }
            }


        }

//        for (String id : idTermHashMap.keySet()) {
//            if (idParentsHashMap.get(id).size() == 0)
//                System.out.println(idTermHashMap.get(id).name);
//        }

        //id到所有的子孙节点
        idChildesHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("child.txt"))) {
            ArrayList<String> childs = new ArrayList<>();
            if (line.split(":").length > 1) {
                for (String c : line.split(":")[1].split("\t")) {
                    childs.add(c);
                }
            }
            idChildesHashMap.put(line.split(":")[0], childs);
        }

//        for (String id:idTermHashMap.keySet()){
//            idChildesHashMap.put(id,new ArrayList<>());
//        }//初始化数据结构
//
//        for (String id:idTermHashMap.keySet()){
//            for (String pid:idTermHashMap.keySet()){
//                if (idParentsHashMap.get(id).contains(pid)){
//                    idChildesHashMap.get(pid).add(id);
//                }
//            }
//        }

        //子孙节点存盘
//        int count = 0;
//        StringBuilder sb = new StringBuilder();
//        for (String id:idChildesHashMap.keySet()){
//
//            sb.append(id + ":");
//            for (String child:idChildesHashMap.get(id)){
//                sb.append(child + "\t");
//            }
//            sb.append("\n");
//            count++;
//            if (count > 4000){
//                Path childfile = Paths.get("child.txt");
//                if (!Files.exists(childfile)){
//                    Files.createFile(childfile);
//                }
//                Files.write(childfile,sb.toString().getBytes(),StandardOpenOption.APPEND);
//                count = 0;
//                sb = new StringBuilder();
//            }
//        }

        //id到所有的子节点
        idSonHashMap = new HashMap<>();

        for (String line : Files.readAllLines(Paths.get("son.txt"))) {
            ArrayList<String> sons = new ArrayList<>();
            if (line.split(":").length > 1) {
                for (String c : line.split(":")[1].split("\t")) {
                    sons.add(c);
                }
            }
            idSonHashMap.put(line.split(":")[0], sons);
        }

//        for (String id :idTermHashMap.keySet()){
//            idSonHashMap.put(id,new ArrayList<>());
//        }
//        for (String id:idTermHashMap.keySet()){
//            for (String pid:idTermHashMap.keySet()){
//                if (idTermHashMap.get(id).partId.contains(pid)){
//                    idSonHashMap.get(pid).add(id);
//                }
//                if (idTermHashMap.get(id).isId.contains(pid)){
//                    idSonHashMap.get(pid).add(id);
//                }
//            }
//        }
//
//
//        int count = 0;
//        StringBuilder sb = new StringBuilder();
//        for (String id:idSonHashMap.keySet()){
//
//            sb.append(id + ":");
//            for (String child:idSonHashMap.get(id)){
//                sb.append(child + "\t");
//            }
//            sb.append("\n");
//            count++;
//            if (count > 4000){
//                Path sonfile = Paths.get("son.txt");
//                if (!Files.exists(sonfile)){
//                    Files.createFile(sonfile);
//                }
//                Files.write(sonfile,sb.toString().getBytes(),StandardOpenOption.APPEND);
//                count = 0;
//                sb = new StringBuilder();
//            }
//        }


        AnnoList = Reader.readGafFile("gene.gaf");


        //读取基因注释列表，转化为hashmap，便于操作
        idAnnoHashmap = new HashMap<>();
        for (Anno anno : AnnoList) {
            if (!idAnnoHashmap.containsKey(anno.goID)) {
                idAnnoHashmap.put(anno.goID, new ArrayList<>());
            }
            idAnnoHashmap.get(anno.goID).add(anno);
        }


        geneTermsHashMap = new HashMap<>();
        idGenesHashMap = new HashMap<>();


        for (Anno anno : AnnoList) {
            if (!idGenesHashMap.containsKey(anno.goID)) {
                idGenesHashMap.put(anno.goID, new ArrayList<>());
            }
            idGenesHashMap.get(anno.goID).add(anno.geneID);
            if (anno.sym != null) {
                idGenesHashMap.get(anno.goID).addAll(anno.sym);
            }
        }

        for (Anno anno : AnnoList) {
            if (!geneTermsHashMap.containsKey(anno.geneID)) {
                geneTermsHashMap.put(anno.geneID, new ArrayList<>());
            }
            geneTermsHashMap.get(anno.geneID).add(anno.goID);
            if (anno.sym != null) {
                for (String gene : anno.sym) {
                    if (!geneTermsHashMap.containsKey(gene)) {
                        geneTermsHashMap.put(gene, new ArrayList<>());
                    }
                    geneTermsHashMap.get(gene).add(anno.goID);
                }
            }
        }
        for (String gene : geneTermsHashMap.keySet()) {
            arrayDistinct(geneTermsHashMap.get(gene));
        }
        for (String id : idGenesHashMap.keySet()) {
            arrayDistinct(idGenesHashMap.get(id));
        }

        //初始化基因功能网路数据
        netHashMap = new HashMap<>();
        ArrayList<String> genenames = new ArrayList<>(4172);
        for (String name: Files.readAllLines(Paths.get("goldlist.txt"))){
            genenames.add(name);
        }
        int i =0, j = 0;
        for(String line:Files.readAllLines(Paths.get("result.mat"))){


            for (String value:line.split("\t")){
                if (value.equals("0.0")){
                    continue;
                }

                String key = genenames.get(i) + ":" + genenames.get(j);
                netHashMap.put(key,Double.parseDouble(value));
                j++;
            }
            i++;
            j = 0;
        }

//        for (FunctionNet edge : Reader.readNetFle("net.txt")) {
//            //做归一化处理
//            netHashMap.put(edge.gene1 + ":" + edge.gene2, edge.value / 10);
//        }

        ecnumGeneSetHashMap = new HashMap<>();
        for (String line:Files.readAllLines(Paths.get("ec.tab"))){
            String[] values = line.split("\t");
            if (values.length < 4){
                continue;
            }
            String ec = values[2];
            String gene = values[3];
            if(!ecnumGeneSetHashMap.containsKey(ec)){
                ecnumGeneSetHashMap.put(ec,new HashSet<>());
            }

            ecnumGeneSetHashMap.get(ec).add(gene);
        }
        HashMap<String,HashSet<String>> tmp = new HashMap<>();
        for(String ec:ecnumGeneSetHashMap.keySet()){
            if (ecnumGeneSetHashMap.get(ec).size() > 1){
                tmp.put(ec,ecnumGeneSetHashMap.get(ec));
            }
        }

        ecnumGeneSetHashMap.clear();
         ecnumGeneSetHashMap.putAll(tmp);
        ecnumGeneSetHashMap.remove("");


        isInit = true;
        System.out.println("init finished");
    }
    private static ArrayList<String> goldList;
    private static ArrayList<ArrayList<Double>> netMatrix;

    public static void storageMatrix() throws Exception {
        try {

                data_initializer();


        } catch (IOException e) {
            System.out.println("数据初始化错误");
            System.out.println(e.getMessage());
            e.printStackTrace();
            return;
        }

//        HashSet<String> goldGene = new HashSet<>();
//        for (String line : Files.readAllLines(Paths.get("gold.txt"))) {
//            goldGene.add(line.split("\t")[0]);
//            goldGene.add(line.split("\t")[1]);
//        }

        goldList = new ArrayList<>();
        for (String line: Files.readAllLines(Paths.get("goldlist.txt"))){
            goldList.add(line);
        }
        StringBuffer sb = new StringBuffer();

        for(String g1:goldList){
            for(String g2:goldList) {
                if (!netHashMap.containsKey(g1+":"+g2)){
                    sb.append(0.0 + "\t");
                }else {
                    sb.append(netHashMap.get(g1 + ":" + g2) + "\t");

                }
            }
            sb.append("\n");
            if (!Files.exists(Paths.get("goldmat.txt"))) {
                Files.createFile(Paths.get("goldmat.txt"));
            }
            Files.write(Paths.get("goldmat.txt"),sb.toString().getBytes(),StandardOpenOption.APPEND);
            sb = new StringBuffer();
        }



        netMatrix = new ArrayList<>();
        netMatrix.ensureCapacity(4172);
        for (int i = 0; i < goldList.size(); i++) {
            netMatrix.add(new ArrayList<>());
            netMatrix.get(i).ensureCapacity(4172);
            for (int j = 0; j <goldList.size() ; j++) {
                String key = goldList.get(i) + ":" +goldList.get(j);
                if (netHashMap.containsKey(key  )) {
                    netMatrix.get(i).add(netHashMap.get(key));
                }else{
                    netMatrix.get(i).add(0.0);
                }
            }
        }

        ArrayList<ArrayList<Double>> pm =  getPMatrix();
        //StringBuilder sb = new StringBuilder();
        for (ArrayList<Double> line:pm){

            for (double value:line){
                sb.append(value + "\t");
            }
            sb.append("\n");
        }

        if (!Files.exists(Paths.get("pMatrix.txt"))){
            Files.createFile(Paths.get("pMatrix.txt"));
        }
        Files.write(Paths.get("pMatrix.txt"),sb.toString().getBytes());
    }


    private static ArrayList<ArrayList<Double>> getWMatrix() throws Exception {
        double o2 = geto2();
        System.out.println("开始计算W");
        long start  = System.currentTimeMillis();
        ArrayList<ArrayList<Double>> wM = new ArrayList<>();
        wM.ensureCapacity(4172);
        for (int i = 0; i < netMatrix.size(); i++) {
            wM.add(new ArrayList<>());
            wM.get(i).ensureCapacity(4172);
            for (int j = 0; j < netMatrix.size(); j++) {
                wM.get(i).add(Math.exp(-(oDistan(netMatrix.get(i),netMatrix.get(j)) / o2)));
            }
        }
        long end = System.currentTimeMillis();
        System.out.println("计算W矩阵耗时" + (end-start) + "毫秒");
        return wM;
    }


    private static double oDistan(ArrayList<Double> l1,ArrayList<Double> l2) throws Exception {
        double value = 0.0;
        if (l1.size() != l2.size()){
            throw  new Exception("向量不等长");
        }
        for (int i = 0; i < l1.size(); i++) {
            value += (l1.get(i) - l2.get(i)) * (l1.get(i) - l2.get(i));
        }
        value = Math.sqrt(value);
        return value;
    }



    private static double geto2() throws Exception {
        System.out.println("开始计算O2");
        long start  = System.currentTimeMillis();
        double sum = 0.0;
        for (int i = 0; i <netMatrix.size() ; i++) {
            for (int j = 0; j < netMatrix.size(); j++) {
                if(i==j){
                    continue;
                }
                sum += oDistan(netMatrix.get(i),netMatrix.get(j));

            }
        }
        sum /= netMatrix.size();
        sum /= (netMatrix.size() -1);

        sum = sum * sum;

        long end = System.currentTimeMillis();
        System.out.println("计算O2耗时" + (end-start) + "毫秒");
        return sum;
    }

    private static ArrayList<ArrayList<Double>> getPMatrix() throws Exception {
        ArrayList<ArrayList<Double>> pMatrix = new ArrayList<>();
        ArrayList<ArrayList<Double>> wMatrix = getWMatrix();
        pMatrix.ensureCapacity(4172);
        System.out.println("开始计算P矩阵");
        long start  = System.currentTimeMillis();
        for (int i = 0; i < netMatrix.size() ; i++) {

            pMatrix.add(new ArrayList<>());
            pMatrix.get(i).ensureCapacity(4172);
            for (int j = 0; j < netMatrix.size() ; j++) {
                //String key = goldList.get(i) + ":" + goldList.get(j);
                if (netMatrix.get(i).get(j) !=  0.0){
                    double value = wMatrix.get(i).get(j);

                    double sum = 0.0;
                    for (int k = 0; k < netMatrix.size(); k++) {
                        sum += wMatrix.get(i).get(k);
                    }
                    value = value / sum;
                    pMatrix.get(i).add(value);
                }else{
                    pMatrix.get(i).add(0.0);
                }
            }
        }

        long end = System.currentTimeMillis();
        System.out.println("计算P矩阵耗时" + (end-start) + "毫秒");
        return pMatrix;
    }

//
//    static class Evalutor{
//        //计算差异倍数
//        public static void LFC(){
//            try {
//                data_initializer();
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//
//            for (String ec:ecnumGeneSetHashMap.keySet()){
//                System.out.println(getLFC(ec));
//            }
//        }
//
//        private static ArrayList<HashSet<String>> getNotInterSet(HashSet<String> eiset){
//            ArrayList<HashSet<String>> result = new ArrayList<>();
//
//            for (String ec : ecnumGeneSetHashMap.keySet()){
//                HashSet<String> ej =  ecnumGeneSetHashMap.get(ec);
//                boolean flag = true;
//                for (String gene:ej){
//                    if (eiset.contains(gene)){
//                        flag = false;
//                        break;
//                    }
//                }
//                if (flag){
//                    result.add(ej);
//                }
//            }
//
//            return result;
//        }
//
//
//        //根据指定的EC号，计算其LFC得分
//        private static double getLFC(String ecNumber){
//            HashSet<String> eiSet = ecnumGeneSetHashMap.get(ecNumber);
//            //System.out.println("开始计算EC编号为：" + ecNumber + "的LFC得分");
//            long start = System.currentTimeMillis();
//            //获取不相交集
//            ArrayList<HashSet<String>> ecSet = getNotInterSet(eiSet);
//
//
//            double value = 0.0;
//
//            int t = 0;
//            for (HashSet<String> ejSet:ecSet){
//                t++;
//                //System.out.println("正在计算" + t + "/" + ecSet.size() );
//
//                //如果ejSet与eiSet相交则跳过
//
//                double sumValue = 0.0;
//                for (String gene:eiSet){
//                    sumValue += getDiff(gene,eiSet,ejSet);
//                    //System.out.println("已使用时间："  + ((System.currentTimeMillis() - start)/1000) + "秒");
//
//                }
//                value += sumValue/eiSet.size();
//
//            }
//            double  EC_NUMBER = ecnumGeneSetHashMap.size();
//
//            //System.out.println("计算EC编号为：" + ecNumber + "的LFC得分总耗时为：" + ((System.currentTimeMillis() - start)/1000) + "秒");
//            return value / EC_NUMBER;
//        }
//        //获取对数差异
//        private static double getDiff(String gene,HashSet<String> ei, HashSet<String > ej){
//            //分子
//            double topValue = 0.0;
//            //分母
//            double bottomValue = 0.0;
//            //拉普拉斯平滑参数
//            double c = 1.0E-10;
//            StringBuilder sb = new StringBuilder();
//
//
//
//            for(String gi:ei){
//                //topValue += getGeneSimilarity(gene,gi) + c;
//                sb.append(gene + "\t" + gi + "\n");
//                topValue+=1;
//            }
//            topValue *= ej.size();
//
//            for (String gj:ej){
//                //bottomValue += getGeneSimilarity(gene,gj) + c;
//                sb.append(gene + "\t" + gj + "\n");
//                bottomValue += 1;
//            }
//            bottomValue *=ei.size();
//
//            try {
//                Files.write(Paths.get("lfc.txt"),sb.toString().getBytes(),StandardOpenOption.APPEND);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//            return Math.log(topValue/bottomValue);
//        }
//    }

}

