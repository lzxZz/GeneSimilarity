import java.util.ArrayList;

class Term {
    String name;
    String namespace;
    String id;
    boolean isObsolote;
    ArrayList<String> partId;
    ArrayList<String> isId;
    public Term(){
        partId = new ArrayList<>();
        isId = new ArrayList<>();
    }
}
