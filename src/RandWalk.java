import com.sun.jmx.snmp.internal.SnmpSubSystem;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

public class RandWalk {


    static private Matrix readPMatrix() throws IOException {
        ArrayList<ArrayList<Double>> plist = new ArrayList<>();
        String[] lines = Files.readAllLines(Paths.get("pMatrix.txt")).toArray(new String[0]);
        for (int i = 0; i < lines.length; i++) {
            plist.add(new ArrayList<>(4172));
            String[] values = lines[i].split("\t");

            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);
                plist.get(i).add(value);
            }

        }
        for (int i = 0; i < 4172 ; i++) {
            double sum = 0.0;

            for (int j = 0; j <4172 ; j++) {
                sum += plist.get(j).get(i) * plist.get(j).get(i);
            }
            for (int j = 0; j <4172 ; j++) {
                if (plist.get(j).get(i) == 0.0){
                    continue;
                }
                plist.get(j).set(i,plist.get(j).get(i)/sum);
            }

        }

        return new Matrix(plist);
    }

    static private Matrix readSimMatrix() throws IOException {
        ArrayList<ArrayList<Double>> plist = new ArrayList<>();
        String[] lines = Files.readAllLines(Paths.get("goldmat.txt")).toArray(new String[0]);

        for (int i = 0; i < lines.length; i++) {
            plist.add(new ArrayList<>(4172));

            String[] values = lines[i].split("\t");


            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);

                plist.get(i).add(value);
            }

        }
        for (int i = 0; i < 4172 ; i++) {
            double sum = 0.0;

            for (int j = 0; j <4172 ; j++) {
                sum += plist.get(j).get(i) * plist.get(j).get(i);
            }
            for (int j = 0; j <4172 ; j++) {
                if (plist.get(j).get(i) == 0.0){
                    continue;
                }
                plist.get(j).set(i,plist.get(j).get(i)/sum);
            }

        }


        return new Matrix(plist);
    }

    static public void walk(double arg) throws Exception {
        Matrix simm = readSimMatrix();
        Matrix pm = readPMatrix();
        pm = readSimMatrix();

        double c = arg;
        int times = 1;//迭代轮次，存盘

        long s = System.currentTimeMillis();
        for (int k = 0; k < simm.getRowcount(); k++) {
            long start = System.currentTimeMillis();
            for (int l = 0; l < 15; l++) {


                Matrix newRow = Matrix.times(simm.getRow(k), pm);
                newRow = Matrix.dot(newRow, c);

                newRow = Matrix.plus(newRow, Matrix.dot(Matrix.getEVector(pm.getColcount(), k), (1 - c)));

                Matrix.setRow(k, newRow, simm);

                //newRow.print();
            }
            long end = System.currentTimeMillis();
            System.out.println("第" + k + "次，耗时：" +  (end-start) + "总耗时:" + (end-s));
        }
        if (!Files.exists(Paths.get("18result"+ arg +".mat"))){
            Files.createFile(Paths.get("18result"+ arg +".mat"));
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i <4172 ; i++) {

            for (int j = 0; j <4172 ; j++) {
                //Files.write(Paths.get("result"+ arg +".mat"),(simm.getMatriValue(i,j) + "\t").getBytes(),StandardOpenOption.APPEND);
                sb.append(simm.getMatriValue(i,j) + "\t");
            }
            sb.append("\n");
            //Files.write(Paths.get("result"+ arg +".mat"),("\n").getBytes(),StandardOpenOption.APPEND);
            Files.write(Paths.get("18result"+ arg +".mat"),sb.toString().getBytes(),StandardOpenOption.APPEND);
            sb = new StringBuilder();
        }


    }
}

