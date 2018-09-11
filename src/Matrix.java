import java.util.ArrayList;

public class Matrix {
    private int rowcount;
    private int colcount;

    private double matrixValue[][];

    public int getColcount() {
        return colcount;
    }

    public int getRowcount() {
        return rowcount;
    }
    public Matrix(int r, int c){
        rowcount = r;
        colcount = c;
        matrixValue = new double[r][c];
    }
    public Matrix getRow(int r){
        double[][] row = new double[1][];
        row[0] = matrixValue[r];


        return  new Matrix(row);
    }
    public static void setRow(int r,Matrix row,Matrix matrix) throws Exception {
        if (row.getRowcount() != 1){
            throw new Exception("不是行向量");
        }
        if (row.getColcount() != matrix.getColcount()){
            throw new Exception("列数不匹配");
        }
        if (r >= matrix.rowcount){
            throw new Exception("没有那么多行");
        }
        for (int i = 0; i < matrix.colcount; i++) {
            matrix.setMatrixValue(r,i, row.getMatriValue(0,i));
        }
    }

    public void print(){
        System.out.println("Row:" + getRowcount() + "\tCol:" + getColcount());
        for (int i = 0; i < getRowcount(); i++) {
            for (int j = 0; j < getColcount() ; j++) {
                try {
                    System.out.print(getMatriValue(i,j)+"\t");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            System.out.println("");
        }
    }

    public Matrix(int r,int c,double value){
        rowcount = r;
        colcount = c;
        matrixValue = new double[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                matrixValue[i][j] = value;
            }
        }
    }
    public Matrix(double[][] values){
        rowcount = values.length;
        colcount = values[0].length;
        matrixValue = new double[rowcount][colcount];
        for (int i = 0; i < rowcount; i++) {
            for (int j = 0; j < colcount; j++) {
                matrixValue[i][j] = values[i][j];
            }
        }
    }

    public Matrix(ArrayList<ArrayList<Double>> value){
        matrixValue = new double[value.size()][value.get(0).size()];
        rowcount = value.size();
        colcount = value.get(0).size();
        for (int i = 0; i < rowcount; i++) {
            for (int j = 0; j < colcount; j++) {
                matrixValue[i][j] = value.get(i).get(j);
            }
        }
    }

    public double getMatriValue(int r, int c) throws Exception {
        if (r >= rowcount || c >= colcount){
            throw new Exception("Out of bound");
        }
        return matrixValue[r][c];
    }

    public void setMatrixValue(int r, int c, double value) throws Exception {
        if (r > rowcount || c > colcount){
            throw new Exception("Out of bound");
        }
        matrixValue[r][c] = value;
    }

    public static Matrix getIMatrix(int n){
        Matrix result = new Matrix(n,n);
        for (int i = 0 ;i < n ;i++){
            try {
                result.setMatrixValue(i,i,1);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return result;
    }

    public static Matrix getEVector(int dimension,int index){
        Matrix result = new Matrix(1,dimension);
        try {
            result.setMatrixValue(0,index,1);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return  result;
    }

    public static Matrix times(Matrix m1,Matrix m2) throws Exception {
        if (m1.getColcount() != m2.getRowcount()){
            throw new Exception("维度不匹配，无法进行矩阵乘法");
        }
        Matrix result = new Matrix(m1.getRowcount(),m2.getColcount());
        for (int i = 0 ; i < result.getRowcount();i++){
            for (int j = 0; j < result.getColcount() ; j++) {
                for (int k = 0 ; k < m1.getColcount();k++)
                    result.setMatrixValue(i,j,result.getMatriValue(i,j)  + m1.getMatriValue(i,k) * m2.getMatriValue(k,j));
            }
        }

        return result;
    }

    public static Matrix times(double[][] m1,double[][] m2) throws Exception {
        return times(conver(m1),conver(m2));
    }

    public static Matrix transpose(Matrix m){
        Matrix newMatrix = new Matrix(m.getColcount(),m.getRowcount());
        for (int i = 0;i < m.getRowcount();i++){
            for (int j = 0; j < m.getColcount(); j++) {
                try {
                    newMatrix.setMatrixValue(j,i,m.getMatriValue(i,j));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        return  newMatrix;
    }

    public static Matrix conver(double m[][]) throws Exception {
        int r = m.length;
        if (r < 1){
            throw new Exception("NULL ARRAY");
        }
        int c = m[0].length;
        if (c < 1){
            throw new Exception("NULL ARRAY");
        }

        Matrix result = new Matrix(r,c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                result.setMatrixValue(i,j,m[i][j]);
            }
        }


        return result;
    }

    public static Matrix plus(Matrix m1,Matrix m2) throws Exception {
        if (m1.getRowcount() != m2.getRowcount() || m1.getColcount() != m2.getColcount()){
            throw new Exception("Can't plus,Matrix's shape are diffrent!");
        }
        Matrix result = new Matrix(m1.getRowcount(),m1.getColcount());
        for (int i = 0; i < result.getRowcount(); i++) {
            for (int j = 0; j < result.getColcount(); j++) {
                result.setMatrixValue(i,j,m1.getMatriValue(i,j) + m2.getMatriValue(i,j));
            }
        }
        return result;
    }

    public static Matrix plus(double[][] m1, double[][] m2) throws Exception {
        return  plus(conver(m1),conver(m2));
    }

    public static Matrix dot(Matrix m,double value){
        for (int i = 0; i < m.getRowcount(); i++) {
            for (int j = 0; j <m.getColcount() ; j++) {
                try {
                    m.setMatrixValue(i,j,value * m.getMatriValue(i,j));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        return m;
    }
}
