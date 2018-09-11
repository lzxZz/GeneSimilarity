import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Evaluator {
    static HashMap<String,HashSet<String>> ecnumGeneSetHashMap;
    static HashMap<String,Double> geneSimHashMap;
    private static void data_initializer(String gene3Wsimfile) throws IOException {
        ecnumGeneSetHashMap = new HashMap<>();
        for (String line:Files.readAllLines(Paths.get("ec.tab"))) {
            String[] values = line.split("\t");
            if (values.length < 4) {
                continue;
            }
            String ec = values[2];
            String gene = values[3];
            if (!ecnumGeneSetHashMap.containsKey(ec)) {
                ecnumGeneSetHashMap.put(ec, new HashSet<>());
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


        geneSimHashMap=  new HashMap<>();
        for (String line: Files.readAllLines(Paths.get(gene3Wsimfile))) {
            String g1, g2;
            if (line.split("\t").length == 3) {


                double value;
                g1 = line.split("\t")[0];
                g2 = line.split("\t")[1];
                value = Double.parseDouble(line.split("\t")[2]);
                geneSimHashMap.put(g1 + ":" + g2, value);
            }
        }


    }

    //计算差异倍数
    public static void LFC(String genesim){
        try {
            data_initializer(genesim);
        } catch (IOException e) {
            e.printStackTrace();
        }

        for (String ec:ecnumGeneSetHashMap.keySet()){
            System.out.println(ec + "\t" + getLFC(ec));
        }
    }

    private static ArrayList<HashSet<String>> getNotInterSet(HashSet<String> eiset){
        ArrayList<HashSet<String>> result = new ArrayList<>();

        for (String ec : ecnumGeneSetHashMap.keySet()){
            HashSet<String> ej =  ecnumGeneSetHashMap.get(ec);
            boolean flag = true;
            for (String gene:ej){
                if (eiset.contains(gene)){
                    flag = false;
                    break;
                }
            }
            if (flag){
                result.add(ej);
            }
        }

        return result;
    }


    //根据指定的EC号，计算其LFC得分
    private static double getLFC(String ecNumber){
        HashSet<String> eiSet = ecnumGeneSetHashMap.get(ecNumber);
        //System.out.println("开始计算EC编号为：" + ecNumber + "的LFC得分");
        long start = System.currentTimeMillis();
        //获取不相交集
        ArrayList<HashSet<String>> ecSet = getNotInterSet(eiSet);


        double value = 0.0;

        int t = 0;
        for (HashSet<String> ejSet:ecSet){
            t++;
            //System.out.println("正在计算" + t + "/" + ecSet.size() );

            //如果ejSet与eiSet相交则跳过

            double sumValue = 0.0;
            for (String gene:eiSet){
                sumValue += getDiff(gene,eiSet,ejSet);
                //System.out.println("已使用时间："  + ((System.currentTimeMillis() - start)/1000) + "秒");

            }
            value += sumValue/eiSet.size();

        }
        double  EC_NUMBER = ecnumGeneSetHashMap.size();
        //EC_NUMBER = 20;

        //System.out.println("计算EC编号为：" + ecNumber + "的LFC得分总耗时为：" + ((System.currentTimeMillis() - start)/1000) + "秒");
        return value / EC_NUMBER;
    }
    //获取对数差异
    private static double getDiff(String gene,HashSet<String> ei, HashSet<String > ej){
        //分子
        double topValue = 0.0;
        //分母
        double bottomValue = 0.0;
        //拉普拉斯平滑参数
        double c = 1.0E-10;




        for(String gi:ei){
            //topValue += getGeneSimilarity(gene,gi) + c;
            String key1,key2;
            key1 = gene +":"+ gi;
            key2 = gi +":"+ gene;
            if (geneSimHashMap.containsKey(key1)){
                topValue+=geneSimHashMap.get(key1);
            }else if(geneSimHashMap.containsKey(key2)){
               topValue += geneSimHashMap.get(key2);
            }
            else{
                topValue +=0;
            }
//            sb.append(gene + "\t" + gi + "\n");
//            topValue+=1;
            topValue +=c;
        }
        topValue *= ej.size();

        for (String gj:ej){
            //bottomValue += getGeneSimilarity(gene,gj) + c;
            String key1,key2;
            key1 = gene +":"+ gj;
            key2 = gj +":"+ gene;
            if (geneSimHashMap.containsKey(key1)){
                bottomValue+=geneSimHashMap.get(key1);
            }else if(geneSimHashMap.containsKey(key2)){
                bottomValue += geneSimHashMap.get(key2);
            }
            else{
                bottomValue +=0;
            }
            bottomValue +=c;
//            sb.append(gene + "\t" + gj + "\n");
//            bottomValue += 1;
        }
        bottomValue *=ei.size();


        return Math.log(topValue/bottomValue);
    }
}
