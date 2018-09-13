import java.io.File;
import java.nio.file.Path;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.spec.ECGenParameterSpec;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.ArrayList;

class LFC {
    static HashSet<String> ec_numbers = new HashSet<String>();
    static HashMap<String, HashSet<String>> ec_genes = new HashMap<>();

    private static void init_data(String gene_file) {
        // 获取ec号，并移除掉只含有单个基因的
        try {
            // HashMap<String, String> ec_gene = new HashMap<>();
            for (String line : Files.readAllLines(Paths.get("./data/ec.tab"))) {
                String infos[] = line.split("\t");
                if (infos.length != 5) {
                    continue;
                }
                
                //忽略掉无编号的路径，和最后一位不是数字的路径
                if (infos[2].equals("") || infos[2].endsWith("-") ) {
                    continue;
                }

                String ec = infos[2];
                ec_numbers.add(ec);
                if (ec_genes.containsKey(ec)) {
                    ec_genes.get(ec).add(infos[3]);
                } else {
                    HashSet<String> tmp_set = new HashSet<>();
                    tmp_set.add(infos[3]);
                    ec_genes.put(ec, tmp_set);
                }

            }
            System.out.println(ec_numbers.size());
            System.out.println(ec_genes.size());

            for (String key : ec_genes.keySet()) {
                if (ec_genes.get(key).size() < 2) {
                    ec_numbers.remove(key);
                }
            }
            HashMap<String, HashSet<String>> tmp_map = new HashMap<>();
            for (String key : ec_numbers) {
                if (ec_genes.containsKey(key)) {
                    tmp_map.put(key, ec_genes.get(key));
                }
            }
            ec_genes.clear();
            for (String key : tmp_map.keySet()) {
                ec_genes.put(key, tmp_map.get(key));
            }

            System.out.println(ec_numbers.size());
            System.out.println(ec_genes.size());
            // for (String key : ec_genes.keySet())
            // {
            // System.out.println(key);
            // }




            for (String line : Files.readAllLines(Paths.get(gene_file)))
            {
                if (line.split("\t").length != 3)
                {
                    continue;
                }
                String key;
                double value;
                key = line.split("\t")[0] + ":" + line.split("\t")[1];
                value = Double.parseDouble(line.split("\t")[2]);
                gene_value.put(key, value);
                
            }


        } catch (Exception e) {
            System.err.println("读取文件失败");
            e.printStackTrace();
        }

    }

    public static void Calculator(String gene_file, String out_file, boolean out_to_console) throws IOException {
        init_data(gene_file);

        if (Files.exists(Paths.get(out_file)))
        {
            System.err.println("输出文件已经存在，请重新指定文件名！");
            return;
        }

        Path out_path = Paths.get(out_file);
        Files.createFile(out_path);

        StringBuilder sb =  new StringBuilder();

        for (String ei : ec_numbers) {
            double lfc = 0;

            for (String ej : ec_numbers) {
                if (!is_interact_by_keys(ei, ej)) {
                    double diff_value = 0;
                    for (String gene : ec_genes.get(ei)) {
                        diff_value += get_diff_value(gene, ei, ej);
                    }

                    diff_value /= ec_genes.get(ei).size();

                    lfc += diff_value;
                }
            }

            lfc /= ec_numbers.size();

            sb.append(ei);
            sb.append("\t\t");
            sb.append(lfc);
            sb.append("\n");

            if (out_to_console)
            {
            System.out.println(lfc);
            }
        }

        Files.write(out_path, sb.toString().getBytes());
        

    }

    static HashMap<String, Double> gene_value = new HashMap<>();

    private static boolean is_interact_by_keys(String key1, String key2) {
        try {
            for (String gene1 : ec_genes.get(key1)) {
                for (String gene2 : ec_genes.get(key2)) {
                    if (gene1.equals(gene2)) {
                        return true;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return false;
    }

    private static double get_diff_value(String gene, String ei, String ej) {
        double c = 1E-10;

        double top_value = 0, bottom_value = 0;

        for (String gi : ec_genes.get(ei)) {
            if (gi.equals(""))
            {
                continue;
            }

            String key = gene + ":" + gi;
            bottom_value += (1 - get_gene_value_by_key(key) + c);
        }
        bottom_value *= ec_genes.get(ej).size();

        for (String gj : ec_genes.get(ej)) {

            if (gj.equals(""))
            {
                continue;
            
            }

            String key = gene + ":" + gj;
            top_value += (1 - get_gene_value_by_key(key) + c);
        }
        top_value *= ec_genes.get(ei).size();

        return Math.log(top_value / bottom_value);
    }

    
    

    private static double get_gene_value_by_key(String key) {
        if (key.split(":").length != 2)
        {
            return 0;
        }
        
        String g1,g2;
        g1 = key.split(":")[0];
        g2 = key.split(":")[1];

        if (g1.compareTo(g2) > 0)
        {
            key = g2 + ":" + g1;
        }


        // gene_keys.add(key);  //读取基因对时使用，目前不使用
        
        if (gene_value.containsKey(key)) {
            return gene_value.get(key);
        }
        return 0;
    }

    // static HashSet<String> gene_keys = new HashSet<>();
    // private static void get_gene_pair()
    // {
    //     try {
    //         StringBuilder sb = new StringBuilder();
    //         File file = new File("./buf/genepair.buf");
    //         if (!Files.exists(file.toPath())) {
    //             Files.createFile(file.toPath());
    //             int buf_count = 1;
    //             for (String pair : gene_keys) {
    //                 buf_count++;
    //                 sb.append(pair);
    //                 sb.append("\n");
    //                 if (buf_count > 5000) {
    //                     Files.write(file.toPath(), sb.toString().getBytes(), StandardOpenOption.APPEND);
    //                     buf_count = 0;
    //                     sb.delete(0, sb.length());
    //                 }
    //             }
    //         }
    // 
    //     } catch (Exception e) {
    // 
    //     }
    // }
}