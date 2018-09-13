import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Main {
//     public static void main2(String... args){
//         if(args.length != 1){
//             System.out.printf("请只输入3W对基因的文件名");
//         }
//         if (!Files.exists(Paths.get(args[0]))){
//             System.out.println("请输入正确的3W对基因的文件名");
//         }else{


//             Evaluator.LFC(args[0]);
//         }
//     }

//     public static void main(String... args)  throws Exception {
//         System.out.println("start");
//         long start = System.currentTimeMillis();
//         int time = 1;
        
//         for (String line: Files.readAllLines(Paths.get("lfcpair.txt"))){
//             System.out.println(line);
//             String gene1,gene2;
//             gene1 = line.split(":")[0];
//             gene2 = line.split(":")[1];
//             double sim  = Caculator.getGeneSimilarity(gene1,gene2);
//             System.out.println("第" + time++ + "组计算\t共用时：" +((System.currentTimeMillis() - start)/1000) + "秒");
//             System.out.println( gene1 + "\t" + gene2 + "\t" + sim);

//             //if (!Files.exists(Paths.get(output))){
//                 //Files.createFile(Paths.get(output));
//             //}
//             //Files.write(Paths.get(output),( gene1 + "\t" + gene2 + "\t" + sim + "\n").getBytes(),StandardOpenOption.APPEND);
//         }
//     }

//     public  static void   main1(String... args) throws Exception{
//         String input;
// //
//         long start = System.currentTimeMillis();
//         input = "lfcpair.txt";
//         String output = "lfcresult.txt";
//         Caculator.data_initializer();

//         String[] lines = Files.readAllLines(Paths.get(input)).toArray(new String[0]);
//         int thCnt = 4;
//         if (args.length == 1){
//             thCnt = Integer.parseInt(args[0]);
//         }

//         for (int i = 0; i < thCnt; i++) {
//             final int startIndex = i;
//             int finalThCnt = thCnt;
//             Thread th = new Thread(new Runnable() {
//                 @Override
//                 public void run() {
//                     long sss = System.currentTimeMillis();
//                     for (int j = startIndex; j < lines.length; j+= finalThCnt) {
//                         String g1,g2;
//                         g1 = lines[j].split(":")[0];
//                         g2 = lines[j].split(":")[1];
//                         //double sim = Caculator.getGeneSimilarity(g1,g2);
//                         double sim = 0;
//                         System.out.print("Thread:" + (startIndex+1) + "\t" + g1 + "\t" + g2 + "\t" + sim);
//                         //System.out.println("Current Time Use:" + (System.currentTimeMillis()-sss));
//                         System.out.println("\tAll Time Use:" + (System.currentTimeMillis()-start));
//                         if (!Files.exists(Paths.get(output))){
//                             try {
//                                 Files.createFile(Paths.get(output));
//                             } catch (IOException e) {
//                                 e.printStackTrace();
//                             }
//                         }
//                         try {
//                             Files.write(Paths.get(output),( g1 + "\t" + g2 + "\t" + sim + "\n").getBytes(),StandardOpenOption.APPEND);
//                         } catch (IOException e) {
//                             e.printStackTrace();
//                         }

//                     }
//                 }
//             });
//             th.start();
//         }


// //        for (String line: Files.readAllLines(Paths.get(input))){
// //            String gene1,gene2;
// //            gene1 = line.split(":")[0];
// //            gene2 = line.split(":")[1];
// //            double sim  = Caculator.getGeneSimilarity(gene1,gene2);
// //            System.out.println("第" + time++ + "组计算\t共用时：" +((System.currentTimeMillis() - start)/1000) + "秒");
// //            System.out.println( gene1 + "\t" + gene2 + "\t" + sim);
// //
// //            if (!Files.exists(Paths.get(output))){
// //                //Files.createFile(Paths.get(output));
// //            }
// //            //Files.write(Paths.get(output),( gene1 + "\t" + gene2 + "\t" + sim + "\n").getBytes(),StandardOpenOption.APPEND);
// //        }
// ////        long start = System.currentTimeMillis();
// ////        for (int i = 0; i < 1; i++) {
// ////            Matrix.times(Matrix.getIMatrix(4000),Matrix.getIMatrix(4000));
// ////        }
// ////        long end = System.currentTimeMillis();
// ////
// ////        System.out.println(end-start);
// //
// //        //Caculator.storageMatrix();
// //

// //            double arg = Double.parseDouble(args[0]);
// //            if (arg >=1 || arg < 0){
// //                System.out.println("请输入迭代的阿尔法参数，范围为0-1之间");
// //            }
// //            RandWalk.walk(arg);
// //
// //        //Caculator.Evalutor.LFC();
// //
// //



//         }
//     }

// //    public static void main(String... args){
// //        printHelp();
// //    }
// //
// //    private static void printHelp(){
// //
// //        String info1 = "支持的计算一共有三种，\n\t一是进行指定参数阿尔法的随机游走迭代过程，\n\t二是指定随机游走矩阵计算LFC所需要的3W基因对的相似度，\n\t三是指定3W基因对计算LFC的值，";
// //        String info2 = "支持的参数一共有三种，\n\t一是线程数-th value来指定,默认为2，\n\t二是类别-type value，必须指定，\n\t三是参数-arg 对应于不同的类型有不同的要求";
// //        System.out.println(info1);
// //        System.out.println(info2);
// //
// //    }



//main函数负责组合各个方法，来完成lfc的计算，其中各个步骤分开运行。
// 其中会变化的数据有
    // 矩阵/网络数据
    // 基因集合之间的功能距离dab
    // 术语相似度
    // 基因相似度
    // lfc
// 计算步骤分为
//     网络数据的获取（矩阵完备化）
//     术语相似度计算（避免计算机因相似度时的重复计算）
//     基因相似度计算（避免计算lfc时候的重复计算）
//     lfc计算（获取最终的结果）
// 项目的伸缩性只需要在矩阵完备化处能够支持接入不同的方法即可


    public static void main(String... args) throws IOException
    {
        // 前三步的计算都能够通过多线程来进行加速
        // 读取参数，选择要进行的计算，输出文件如果已经存在，则不会进行计算，避免覆盖掉已有数据


        //矩阵完备化 MatrixComple
        //输入参数,  类型选择，输出文件，参数列表(字符串，逗号分隔多个参数)，参数解析由完备化的具体类进行解析， 线程数，默认为2
        // 只依赖于net数据

        // 术语相似度计算 TermSimCalc
        // 输入参数， 网络数据文件， 输出文件，线程数，默认为2
        // 依赖于网络数据，本体结构，注释信息
        TermSim.Calculator("net_file", "out_file", 2);


        // 基因相似度计算   GeneSim
        // 输入参数， 术语相似度文件，输出文件，线程数，默认为2
        // 依赖于术语相似度
        GeneSim.calculator("term_sim_file", "out_file", 2);

        // lfc计算  Evaluator
        // 输入参数， 基因相似度文件，输出文件，是否输出到控制台
        // 依赖于基因相似度
        
        LFC.Calculator("./result/gene.result", "./result/lfc.result", true);
    }
}




